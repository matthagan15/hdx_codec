use core::panic;
use std::{
    collections::{HashMap, VecDeque},
    fs::{self},
    path::{Path, PathBuf},
    str::FromStr,
    sync::{Arc, RwLock},
    time::Instant,
};

use mhgl::{HGraph, HyperGraph};
use serde::{ser::SerializeStruct, Deserialize, Serialize};

use super::{
    coset_complex_subgroups::KTypeSubgroup, galois_field::GaloisField, polynomial::FFPolynomial,
};
use crate::matrices::galois_matrix::GaloisMatrix;

pub type NodeData = u16;

pub fn get_first_node_complete_star(quotient: FFPolynomial, dim: usize) -> HGraph<NodeData, ()> {
    let mut bfs_manager = BFSState::new(quotient.clone(), dim, vec![usize::MAX]);
    let lookup = Arc::new(RwLock::new(GaloisField::new(quotient.clone())));
    let mut hg = HGraph::new();
    while bfs_manager.current_bfs_distance < 3 {
        bfs_manager.step(&KTypeSubgroup::new(dim, lookup.clone()), &mut hg);
    }
    dbg!(hg.containing_edges_of_nodes([0]));
    hg.star([0])
}

pub fn bfs_benchmark() {
    let q = FFPolynomial::from_str("1*x^2 + 2 * x^1 + 2 * x^0 % 3").unwrap();
    let dim = 3;
    let lookup = Arc::new(RwLock::new(GaloisField::new(q.clone())));
    let subgroups = KTypeSubgroup::new(dim, lookup.clone());
    let mut bfs = BFSState::new(q.clone(), dim, vec![1000]);
    let mut hg = HGraph::new();
    let start = Instant::now();
    for _ in 0..1000 {
        bfs.step(&subgroups, &mut hg);
    }
    println!("time elapsed: {:}", start.elapsed().as_secs_f64());
}

/// Helper struct for BFS
#[derive(Debug, Clone, Serialize, Deserialize)]
struct GroupBFSNode {
    distance: u32,
    hg_node: Vec<u32>,
}

#[derive(Debug, Clone)]
pub struct BFSState {
    frontier: VecDeque<(GaloisMatrix, u32)>,
    visited: HashMap<GaloisMatrix, GroupBFSNode>,
    current_bfs_distance: u32,
    last_flushed_distance: u32,
    num_matrices_completed: usize,
    last_cached_matrices_done: u64,
    quotient: FFPolynomial,
    dim: usize,
    checkpoints: Vec<usize>,
    next_checkpoint: usize,
    prev_completed_nodes: Vec<Vec<u32>>,
    hgraph: HGraph<NodeData, ()>,
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
        s.serialize_field("quotient", &self.quotient)?;
        s.serialize_field("dim", &self.dim)?;
        s.serialize_field("checkpoints", &self.checkpoints[..])?;
        s.serialize_field("next_checkpoint", &self.next_checkpoint)?;
        // let prev_completed_nodes: Vec<Vec<u32>> = self.prev_completed_nodes.iter().cloned().collect();
        s.serialize_field("prev_completed_nodes", &self.prev_completed_nodes[..])?;
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
        let quotient: Result<FFPolynomial, _> = serde_json::from_value(data["quotient"].take());
        let dim: Result<usize, _> = serde_json::from_value(data["dim"].take());
        let checkpoints: Result<Vec<usize>, _> = serde_json::from_value(data["checkpoints"].take());
        let next_checkpoint: Result<usize, _> =
            serde_json::from_value(data["next_checkpoint"].take());
        let prev_completed_nodes: Result<Vec<Vec<u32>>, _> =
            serde_json::from_value(data["prev_completed_nodes"].take());
        let hg: Result<HGraph<NodeData, ()>, _> = serde_json::from_value(data["hgraph"].take());
        Ok(BFSState {
            frontier: frontier.expect("Could not deserialize frontier"),
            visited,
            current_bfs_distance: current_bfs_distance.unwrap(),
            last_flushed_distance: last_flushed_distance.unwrap(),
            num_matrices_completed: num_matrices_completed.unwrap(),
            last_cached_matrices_done: last_cached_matrices_done.unwrap(),
            quotient: quotient.unwrap(),
            dim: dim.unwrap(),
            checkpoints: checkpoints.unwrap(),
            next_checkpoint: next_checkpoint.unwrap(),
            prev_completed_nodes: prev_completed_nodes.unwrap(),
            hgraph: hg.unwrap(),
        })
    }
}

impl BFSState {
    pub fn new(quotient: FFPolynomial, dim: usize, checkpoints: Vec<usize>) -> Self {
        let mut frontier = VecDeque::new();
        frontier.push_back((GaloisMatrix::id(3), 0));
        let first_checkpoint = *checkpoints
            .first()
            .expect("No checkpoints provided for BFS.");
        Self {
            frontier,
            visited: HashMap::new(),
            current_bfs_distance: 0,
            last_flushed_distance: 0,
            num_matrices_completed: 0,
            last_cached_matrices_done: 0,
            quotient,
            dim,
            checkpoints,
            next_checkpoint: first_checkpoint,
            prev_completed_nodes: Vec::new(),
            hgraph: HGraph::new(),
        }
    }

    pub fn from_cache(cache_file: Option<&Path>) -> Option<Self> {
        cache_file.map(|path| {
            if path.is_file() == false {
                None
            } else {
                let file_data = fs::read_to_string(path).expect("Could not read from file.");
                match serde_json::from_str::<BFSState>(&file_data) {
                    Ok(bfs_state) => {
                        log::trace!("BFSState retrieved from cache.");
                        Some(bfs_state)
                    }
                    Err(e) => {
                        log::error!("Could not deserialize cache from file although a file was present. Cache is invalid, fix it!");
                        log::error!("serde_json error = {:}", e);
                        panic!()
                    }
                }
            }
        }).flatten()
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
        let new_edges = Vec::new();
        let x = self.frontier.pop_front().unwrap();
        // process this matrix first, compute the cosets and triangles it can
        // be a part of.
        let neighbors = subgroup_generators.generate_right_mul(&x.0);
        let coset_reps = subgroup_generators.coset_reps(&neighbors[..]);
        // let (c0, c1, c2) = (coset_reps[0], coset_reps[1], coset_reps[2]);

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
        let n0 = if self.visited.contains_key(&coset_reps[0].rep) {
            let v: Vec<u32> = self
                .visited
                .get(&coset_reps[0].rep)
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
                    .get_mut(&coset_reps[0].rep)
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
        let n1 = if self.visited.contains_key(&coset_reps[1].rep) {
            let v: Vec<u32> = self
                .visited
                .get(&coset_reps[1].rep)
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
                    .get_mut(&coset_reps[1].rep)
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
        let n2 = if self.visited.contains_key(&coset_reps[2].rep) {
            let v: Vec<u32> = self
                .visited
                .get(&coset_reps[2].rep)
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
                    .get_mut(&coset_reps[2].rep)
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
        let _check_1 = if let Some(id) = hg.find_id(&[n0, n1]) {
            id
        } else {
            hg.add_edge(&[n0, n1], ())
        };
        let _check_2 = if let Some(id) = hg.find_id(&[n0, n2]) {
            id
        } else {
            hg.add_edge(&[n0, n2], ())
        };
        let _check_3 = if let Some(id) = hg.find_id(&[n1, n2]) {
            id
        } else {
            hg.add_edge(&[n1, n2], ())
        };
        let _message_id = if let Some(id) = hg.find_id(&[n0, n1, n2]) {
            id
        } else {
            hg.add_edge(&[n0, n1, n2], ())
        };
        new_edges
    }

    fn get_next_checkpoint(&self) -> Option<usize> {
        self.checkpoints
            .get(self.prev_completed_nodes.len())
            .cloned()
    }

    pub fn run(&self) {}
}

/// Returns the number of maximal faces that will be in the coset complex, aka
/// computes the size of SL(dim, q = p^n).
pub fn size_of_coset_complex(quotient: &FFPolynomial, dim: usize) -> usize {
    let p = quotient.field_mod;
    let n = quotient.degree();
    let q = p.pow(n.try_into().unwrap()) as usize;
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
    // let _maximum_number_matrices = size_of_coset_complex(&quotient, matrix_dim);
    let mut new_edges = Vec::new();
    let truncation = truncation.unwrap_or(usize::MAX);
    if matrix_dim != 3 {
        panic!("Only dimension 3 matrices are currently supported.")
    }
    let mut bfs_state = BFSState::from_cache(cache_file.as_ref().map(|x| x.as_path())).unwrap();
    // log::trace!("-------------- BFS -------------");
    // log::trace!(
    //     "Frontier: {:}, visited: {:}, current_distance: {:}, num_matrices_completed: {:}, truncation: {:}",
    //     bfs_state.frontier.len(),
    //     bfs_state.visited.len(),
    //     bfs_state.current_bfs_distance,
    //     bfs_state.num_matrices_completed,
    //     truncation,
    // );
    let mut hg = if bfs_state.num_matrices_completed > 0 {
        // The only way this will be positive is if a cache was successfully retrieved, so
        // cache_file must be Some()
        let hg_path = cache_file.as_ref().unwrap().with_extension("hg");
        match HGraph::from_file(hg_path.as_path()) {
            Some(hg) => {
                // log::trace!("HGraph retrieved from cache file.");
                hg
            }
            None => {
                panic!("Could not retrieve HGraph from cache file but bfs cache was present. Fix!")
            }
        }
    } else {
        // log::trace!("Initial BFS Step, creating new HGraph.");
        HGraph::new()
    };
    let mut cache_points = Vec::new();
    if let Some(num_cache_checkpoints) = num_cache_checkpoints {
        let cache_rate = truncation / (num_cache_checkpoints);
        // log::trace!("truncation: {:}, cache_rate: {:}", truncation, cache_rate);
        let mut cur_cache_point = 0;
        for _ in 0..num_cache_checkpoints {
            cur_cache_point += cache_rate;
            if cur_cache_point > bfs_state.num_matrices_completed {
                cache_points.push(cur_cache_point);
            }
        }
        // log::trace!("Caching checkpoints: {:?}", cache_points);
    }
    let mut num_steps = 0;
    let mut time_in_step = 0.0;
    let logging_rate = 100_000;
    let lookup_start = Instant::now();
    let lookup = Arc::new(RwLock::new(GaloisField::new(quotient.clone())));
    println!(
        "Galois lookup table creation time: {:}",
        lookup_start.elapsed().as_secs_f64()
    );
    let subgroup_generators = KTypeSubgroup::new(3, lookup);
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
                // log::trace!("Reached Caching checkpoint. Here are some stats. num_matrices_completed: {:}, frontier: {:}, bfs_distance: {:}", bfs_state.num_matrices_completed, bfs_state.frontier.len(),bfs_state.current_bfs_distance );
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
            log::trace!(
                "Estimated {:} seconds remaining. {:} seconds per step",
                time_remaining,
                time_per_step
            );
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

#[cfg(test)]
mod tests {
    use super::bfs;
    use crate::{
        math::{
            coset_complex_bfs::BFSState, coset_complex_subgroups::KTypeSubgroup,
            finite_field::FFRep, galois_field::GaloisField, polynomial::FFPolynomial,
        },
        matrices::galois_matrix::GaloisMatrix,
        quantum,
    };
    use mhgl::HGraph;
    use simple_logger::SimpleLogger;
    use std::{
        path::PathBuf,
        sync::{Arc, RwLock},
    };

    #[test]
    fn small_checkpoints() {
        let p: FFRep = 3;
        let dim = 3;
        let quotient = FFPolynomial::from(
            &vec![(0, (2, p).into()), (1, (2, p).into()), (2, (1, p).into())][..],
        );
        let mut bfs_manager = BFSState::new(quotient.clone(), dim, vec![usize::MAX]);
        let lookup = Arc::new(RwLock::new(GaloisField::new(quotient.clone())));
        let subgroups = KTypeSubgroup::new(dim, lookup.clone());

        println!("monomial: {:}", FFPolynomial::monomial((0, p).into(), 2));

        for s in subgroups.generate_left_mul(&GaloisMatrix::id(dim)) {
            if s == GaloisMatrix::id(dim) {
                println!("{:}", s.pretty_print(lookup.clone()));
            }
        }
        let mut hg = HGraph::new();
        for _ in 0..3 {
            bfs_manager.step(&subgroups, &mut hg);
        }

        println!("hg: {:}", hg);
    }

    #[test]
    fn as_function() {
        let _logger = SimpleLogger::new().init().unwrap();
        let p = 3_u32;
        let primitive_coeffs = [(2, (1, p).into()), (1, (2, p).into()), (0, (2, p).into())];
        let q = FFPolynomial::from(&primitive_coeffs[..]);
        let cache_file = PathBuf::from("/Users/matt/repos/qec/tmp/new_cache");
        if cache_file.is_file() {
            std::fs::remove_file(cache_file.as_path()).unwrap();
        }
        let first_step_edges = vec![0, 1, 2, 3];
        let second_step_edges = vec![4, 5, 6];
        let (_hg, first_step_computed) = bfs(q.clone(), 3, Some(1), Some(cache_file.clone()), None);
        assert_eq!(first_step_computed, first_step_edges);
        let (_hg2, second_step_computed) =
            bfs(q.clone(), 3, Some(3), Some(cache_file.clone()), None);
        assert_eq!(second_step_computed, second_step_edges);
        let (_hg3, _third_step) = bfs(
            q.clone(),
            3,
            Some(10_001),
            Some(cache_file.clone()),
            Some(10),
        );
        let (_hg4, _fourth_step) = bfs(q, 3, Some(11_001), Some(cache_file), Some(10));
    }
}
