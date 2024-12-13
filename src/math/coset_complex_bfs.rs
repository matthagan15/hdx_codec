use core::panic;
use std::{
    collections::{HashMap, VecDeque},
    fs::{self},
    io::{Read, Write},
    path::{Path, PathBuf},
    sync::Arc,
    thread,
    time::Instant,
};

use fxhash::{FxHashMap, FxHashSet};
use log::trace;
use mhgl::{ConGraph, HGraph, HyperGraph};
use rand::{seq::SliceRandom, thread_rng, Rng};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use serde::ser::SerializeStruct;
use serde::{Deserialize, Serialize};
use serde_json::Value;

use super::{
    coset_complex_subgroups::{CosetRep, KTypeSubgroup},
    finite_field::{FFRep, FiniteField},
    galois_field::GaloisField,
    polynomial::{FFPolynomial, PolyDegree},
};
use crate::matrices::{galois_matrix::GaloisMatrix, polymatrix::PolyMatrix};

#[derive(Debug, Clone, Serialize, Deserialize)]
struct GroupBFSConfig {
    quotient: String,
    dim: usize,
}

/// Helper struct for BFS
#[derive(Debug, Clone, Serialize, Deserialize)]
struct GroupBFSNode {
    distance: u32,
    hg_node: Vec<u32>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct GroupBFSCache {
    quotient: FFPolynomial,
    frontier: VecDeque<(GaloisMatrix, u32)>,
    visited: FxHashMap<GaloisMatrix, GroupBFSNode>,
    coset_to_node: HashMap<CosetRep, u32>,
    current_distance: u32,
    last_flushed_distance: u32,
    num_matrices_completed: u32,
    last_cached_matrices_done: u64,
    cumulative_time_ms: u128,
    directory: PathBuf,
    filename: String,
}

// #[derive(Debug, Clone, Serialize, Deserialize)]
// pub struct NodeData {
//     type_ix: u16,
// }
pub type NodeData = u16;

#[derive(Debug, Clone)]
pub struct BFSState {
    frontier: VecDeque<(GaloisMatrix, u32)>,
    visited: HashMap<GaloisMatrix, GroupBFSNode>,
    current_bfs_distance: u32,
    last_flushed_distance: u32,
    num_matrices_completed: usize,
    last_cached_matrices_done: u64,
}

impl Serialize for BFSState {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        let mut visited_string = String::new();
        visited_string.push_str("[");
        for (k, v) in self.visited.iter() {
            let key_string = serde_json::to_string(k).expect("Could not serialize matrix");
            let val_string = serde_json::to_string(v).expect("Could not serialize GroupBFSNode");
            visited_string.push_str(&key_string[..]);
            visited_string.push_str("+");
            visited_string.push_str(&val_string[..]);
            visited_string.push_str(";");
        }
        visited_string.pop();
        visited_string.push_str("]");

        let mut s = serializer.serialize_struct("BFSState", 6)?;
        s.serialize_field("frontier", &self.frontier)?;
        s.serialize_field("visited", &visited_string[..])?;
        s.serialize_field("current_bfs_distance", &self.current_bfs_distance)?;
        s.serialize_field("last_flushed_distance", &self.last_flushed_distance)?;
        s.serialize_field("num_matrices_completed", &self.num_matrices_completed)?;
        s.serialize_field("last_cached_matrices_done", &self.last_cached_matrices_done)?;
        s.end()
    }
}

impl<'de> Deserialize<'de> for BFSState {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let mut data = <serde_json::Value>::deserialize(deserializer)?;
        let frontier: Result<VecDeque<(GaloisMatrix, u32)>, _> =
            serde_json::from_value(data["frontier"].take());

        let current_bfs_distance: Result<u32, _> =
            serde_json::from_value(data["current_bfs_distance"].take());
        let last_flushed_distance: Result<u32, _> =
            serde_json::from_value(data["last_flushed_distance"].take());
        let num_matrices_completed: Result<usize, _> =
            serde_json::from_value(data["num_matrices_completed"].take());
        let last_cached_matrices_done: Result<u64, _> =
            serde_json::from_value(data["last_cached_matrices_done"].take());

        let mut visited_string =
            String::from(data["visited"].as_str().expect("Cannot get string."));
        if visited_string.starts_with('[') {
            visited_string.remove(0);
        }
        if visited_string.ends_with(']') {
            visited_string.pop();
        }
        let mut visited: HashMap<GaloisMatrix, GroupBFSNode> = HashMap::new();
        for key_val_string in visited_string.split(';') {
            let key_val: Vec<&str> = key_val_string.split('+').collect();
            if key_val.len() != 2 {
                panic!("Cannot deserialize visited map, does not have key_value pair in string.")
            }
            let key: GaloisMatrix =
                serde_json::from_str(key_val[0]).expect("Could not deserialize visited key");
            let val: GroupBFSNode =
                serde_json::from_str(key_val[1]).expect("Could not deserialize visited val");
            visited.insert(key, val);
        }
        Ok(BFSState {
            frontier: frontier.expect("Could not deserialize frontier"),
            visited,
            current_bfs_distance: current_bfs_distance.unwrap(),
            last_flushed_distance: last_flushed_distance.unwrap(),
            num_matrices_completed: num_matrices_completed.unwrap(),
            last_cached_matrices_done: last_cached_matrices_done.unwrap(),
        })
    }
}

impl BFSState {
    pub fn new(cache_file: Option<&Path>) -> Self {
        match cache_file {
            Some(path) => {
                if path.is_file() == false {
                    log::trace!("Provided cache file does not exist yet. Creating new BFSState.");
                    return BFSState::new(None);
                }
                let file_data = fs::read_to_string(path).expect("Could not read from file.");
                match serde_json::from_str(&file_data) {
                    Ok(bfs_state) => {
                        log::trace!("BFSState retrieved from cache.");
                        bfs_state
                    }
                    Err(e) => {
                        log::error!("Could not deserialize cache from file although a file was present. Cache is invalid, fix it!");
                        log::error!("serde_json error = {:}", e);
                        panic!()
                    }
                }
            }
            None => {
                let mut frontier = VecDeque::new();
                frontier.push_back((GaloisMatrix::id(3), 0));
                Self {
                    frontier,
                    visited: HashMap::new(),
                    current_bfs_distance: 0,
                    last_flushed_distance: 0,
                    num_matrices_completed: 0,
                    last_cached_matrices_done: 0,
                }
            }
        }
    }

    pub fn cache(&self, cache_file: &Path) {
        match serde_json::to_string(&self) {
            Ok(data) => {
                let out = std::fs::write(cache_file, data);
                if out.is_err() {
                    log::error!("Error on writing cache: {:?}", out.unwrap());
                    panic!()
                }
            }
            Err(e) => {
                log::error!("Error on serializing data: {:?}", e);
                // dbg!(self);
                panic!()
            }
        }
    }

    /// Returns newly discovered edge_id's.
    pub fn step(
        &mut self,
        subgroup_generators: &KTypeSubgroup,
        hg: &mut HGraph<u16, ()>,
    ) -> Vec<u64> {
        if self.frontier.is_empty() {
            return Vec::new();
        }
        let mut new_edges = Vec::new();
        let x = self.frontier.pop_front().unwrap();
        // process this matrix first, compute the cosets and triangles it can
        // be a part of.
        let neighbors = subgroup_generators.generate_right_mul(&x.0);
        let (c0, c1, c2) = subgroup_generators.coset_reps(&neighbors[..]);

        // flush visited and coset to node
        self.current_bfs_distance = x.1;

        for neighbor in neighbors {
            let neighbor_bfs = GroupBFSNode {
                distance: x.1 + 1,
                hg_node: Vec::new(),
            };
            if self.visited.contains_key(&neighbor) == false {
                self.visited.insert(neighbor.clone(), neighbor_bfs);
                self.frontier.push_back((neighbor, x.1 + 1));
            }
        }
        self.num_matrices_completed += 1;
        let n0 = if self.visited.contains_key(&c0.rep) {
            let v: Vec<u32> = self
                .visited
                .get(&c0.rep)
                .unwrap()
                .hg_node
                .iter()
                .filter(|n| {
                    if let Some(x) = hg.get_node(n) {
                        *x == 0
                    } else {
                        false
                    }
                })
                .cloned()
                .collect();
            if v.len() == 0 {
                let new_node = hg.add_node(0);
                self.visited
                    .get_mut(&c0.rep)
                    .unwrap()
                    .hg_node
                    .push(new_node);
                new_node
            } else if v.len() == 1 {
                v[0]
            } else {
                panic!("Why do I have two nodes from the same matrix with the same type?")
            }
        } else {
            panic!("Why have I computed a coset but not added the matrix to visited yet?")
        };
        let n1 = if self.visited.contains_key(&c1.rep) {
            let v: Vec<u32> = self
                .visited
                .get(&c1.rep)
                .unwrap()
                .hg_node
                .iter()
                .filter(|n| {
                    if let Some(x) = hg.get_node(n) {
                        *x == 1
                    } else {
                        false
                    }
                })
                .cloned()
                .collect();
            if v.len() == 0 {
                let new_node = hg.add_node(1);
                self.visited
                    .get_mut(&c1.rep)
                    .unwrap()
                    .hg_node
                    .push(new_node);
                new_node
            } else if v.len() == 1 {
                v[0]
            } else {
                panic!("Why do I have two nodes from the same matrix with the same type?")
            }
        } else {
            panic!("Why have I computed a coset but not added the matrix to visited yet?")
        };
        let n2 = if self.visited.contains_key(&c2.rep) {
            let v: Vec<u32> = self
                .visited
                .get(&c2.rep)
                .unwrap()
                .hg_node
                .iter()
                .filter(|n| {
                    if let Some(x) = hg.get_node(n) {
                        *x == 2
                    } else {
                        false
                    }
                })
                .cloned()
                .collect();
            if v.len() == 0 {
                let new_node = hg.add_node(2);
                self.visited
                    .get_mut(&c2.rep)
                    .unwrap()
                    .hg_node
                    .push(new_node);
                new_node
            } else if v.len() == 1 {
                v[0]
            } else {
                panic!("Why do I have two nodes from the same matrix with the same type?")
            }
        } else {
            panic!("Why have I computed a coset but not added the matrix to visited yet?")
        };
        let check_1 = if let Some(id) = hg.find_id(&[n0, n1]) {
            id
        } else {
            hg.add_edge(&[n0, n1], ())
        };
        let check_2 = if let Some(id) = hg.find_id(&[n0, n2]) {
            id
        } else {
            hg.add_edge(&[n0, n2], ())
        };
        let check_3 = if let Some(id) = hg.find_id(&[n1, n2]) {
            id
        } else {
            hg.add_edge(&[n1, n2], ())
        };
        let message_id = if let Some(id) = hg.find_id(&[n0, n1, n2]) {id} else {
            hg.add_edge(&[n0, n1, n2], ())
        };
        new_edges
    }
}

/// Returns the number of maximal faces that will be in the coset complex, aka
/// computes the size of SL(dim, q = p^n).
pub fn size_of_coset_complex(quotient: &FFPolynomial, dim: usize) -> usize {
    let p = quotient.field_mod;
    let n = quotient.degree();
    let q = p.pow(n) as usize;
    let mut prod: usize = 1;
    let q_dim = q
        .checked_pow(dim as u32)
        .expect("Coset Complex too big to even compute size of with fixed width unsigned int.");
    for k in 0..dim {
        prod = prod.checked_mul(q_dim
            - q.checked_pow(k as u32).expect(
                "Coset Complex too big to even compute size of with fixed width unsigned int.",
            )).expect("Coset Complex too big to compute size of with fixed width ints.");
    }
    prod / (q - 1)
}

/// Returns all newly filled out parity checks for iterative rank estimation.
pub fn bfs(
    quotient: FFPolynomial,
    matrix_dim: usize,
    truncation: Option<usize>,
    cache_file: Option<PathBuf>,
    num_cache_checkpoints: Option<usize>,
) -> (HGraph<u16, ()>, Vec<u64>) {
    let _maximum_number_matrices = size_of_coset_complex(&quotient, matrix_dim);
    let mut new_edges = Vec::new();
    let truncation = truncation
        .unwrap_or(usize::MAX);
    if matrix_dim != 3 {
        panic!("Only dimension 3 matrices are currently supported.")
    }
    let mut bfs_state = BFSState::new(cache_file.as_ref().map(|x| x.as_path()));
    log::trace!("-------------- BFS -------------");
    log::trace!(
        "Frontier: {:}, visited: {:}, current_distance: {:}, num_matrices_completed: {:}, truncation: {:}",
        bfs_state.frontier.len(),
        bfs_state.visited.len(),
        bfs_state.current_bfs_distance,
        bfs_state.num_matrices_completed, 
        truncation,
    );
    let mut hg = if bfs_state.num_matrices_completed > 0 {
        // The only way this will be positive is if a cache was successfully retrieved, so
        // cache_file must be Some()
        let hg_path = cache_file.as_ref().unwrap().with_extension("hg");
        match HGraph::from_file(hg_path.as_path()) {
            Some(hg) => {
                log::trace!("HGraph retrieved from cache file.");
                hg
            }
            None => {
                panic!("Could not retrieve HGraph from cache file but bfs cache was present. Fix!")
            }
        }
    } else {
        log::trace!("Initial BFS Step, creating new HGraph.");
        HGraph::new()
    };
    let mut cache_points = Vec::new();
    if let Some(num_cache_checkpoints) = num_cache_checkpoints {
        let cache_rate = truncation / (num_cache_checkpoints);
        log::trace!("truncation: {:}, cache_rate: {:}", truncation, cache_rate);
        let mut cur_cache_point = 0;
        for _ in 0..num_cache_checkpoints {
            cur_cache_point += cache_rate;
            if cur_cache_point > bfs_state.num_matrices_completed {
                cache_points.push(cur_cache_point);
            }
        }
        log::trace!("Caching checkpoints: {:?}", cache_points);
    }
    let mut num_steps = 0;
    let mut time_in_step = 0.0;
    let logging_rate = 100_000;
    let lookup = Arc::new(GaloisField::new(quotient.clone()));
    let subgroup_generators = KTypeSubgroup::new(&lookup);
    while bfs_state.num_matrices_completed < truncation && bfs_state.frontier.len() > 0 {
        let start = Instant::now();
        let new_step_edges = bfs_state.step(&subgroup_generators, &mut hg);
        time_in_step += start.elapsed().as_secs_f64();
        num_steps += 1;
        for new_edge in new_step_edges {
            new_edges.push(new_edge);
        }
        if let Some(cache_file) = cache_file.clone() {
            if cache_points.len() > 0 && cache_points[0] == bfs_state.num_matrices_completed {
                log::trace!("Reached Caching checkpoint. Here are some stats. num_matrices_completed: {:}, frontier: {:}, bfs_distance: {:}", bfs_state.num_matrices_completed, bfs_state.frontier.len(),bfs_state.current_bfs_distance );
                bfs_state.cache(cache_file.as_path());
                hg.to_disk(cache_file.with_extension("hg").as_path());
                cache_points.remove(0);
            }
        }
        if num_steps % logging_rate == 0 {
            log::trace!("Reached Logging checkpoint. Here are some stats. num_matrices_completed: {:}, frontier: {:}, bfs_distance: {:}", bfs_state.num_matrices_completed, bfs_state.frontier.len(),bfs_state.current_bfs_distance );
            let time_per_step = time_in_step / (num_steps as f64);
            let num_steps_remaining = truncation - bfs_state.num_matrices_completed;
            let time_remaining = time_per_step * num_steps_remaining as f64;
            log::trace!("Estimated {:} seconds remaining. {:} seconds per step", time_remaining, time_per_step);
        }
    }
    log::trace!("BFS complete!");
    if let Some(cache_file) = cache_file {
        log::trace!("Caching bfs state and hgraph.");
        bfs_state.cache(cache_file.as_path());
        hg.to_disk(cache_file.with_extension("hg").as_path());
    }
    (hg, new_edges)
}

mod tests {
    use std::{
        collections::{HashMap, HashSet},
        io::Write,
        path::PathBuf,
        rc::Rc,
        str::FromStr,
        sync::Arc,
    };

    use simple_logger::SimpleLogger;

    use crate::math::{
        coset_complex_bfs::BFSState, galois_field::GaloisField, polynomial::FFPolynomial,
    };
    use crate::matrices::galois_matrix::GaloisMatrix;

    use super::{bfs};

    fn simple_quotient_and_field() -> (u32, FFPolynomial) {
        let p = 3_u32;
        let primitive_coeffs = [(2, (1, p).into()), (1, (2, p).into()), (0, (2, p).into())];
        let q = FFPolynomial::from(&primitive_coeffs[..]);
        (p, q)
    }




    #[test]
    fn as_function() {
        let logger = SimpleLogger::new().init().unwrap();
        let p = 3_u32;
        let primitive_coeffs = [(2, (1, p).into()), (1, (2, p).into()), (0, (2, p).into())];
        let q = FFPolynomial::from(&primitive_coeffs[..]);
        let cache_file = PathBuf::from("/Users/matt/repos/qec/tmp/new_cache");
        if cache_file.is_file() {
            std::fs::remove_file(cache_file.as_path()).unwrap();
        }
        let first_step_edges = vec![0, 1, 2, 3];
        let second_step_edges = vec![4, 5, 6];
        let (hg, first_step_computed) = bfs(q.clone(), 3, Some(1), Some(cache_file.clone()), None);
        assert_eq!(first_step_computed, first_step_edges);
        let (hg2, second_step_computed) = bfs(q.clone(), 3, Some(3), Some(cache_file.clone()), None);
        assert_eq!(second_step_computed, second_step_edges);
        let (hg3, third_step) = bfs(q.clone(), 3, Some(10_001), Some(cache_file.clone()), Some(10));
        let (hg4, fourth_step) = bfs(q, 3, Some(11_001), Some(cache_file), Some(10));
    }
}
