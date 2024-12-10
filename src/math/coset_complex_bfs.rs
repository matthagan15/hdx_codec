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
    let lookup = Arc::new(GaloisField::new(quotient.clone()));
    let subgroup_generators = KTypeSubgroup::new(&lookup);
    while bfs_state.num_matrices_completed < truncation && bfs_state.frontier.len() > 0 {
        let new_step_edges = bfs_state.step(&subgroup_generators, &mut hg);
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
    }
    log::trace!("BFS complete!");
    if let Some(cache_file) = cache_file {
        log::trace!("Caching bfs state and hgraph.");
        bfs_state.cache(cache_file.as_path());
        hg.to_disk(cache_file.with_extension("hg").as_path());
    }
    (hg, new_edges)
}

#[derive(Debug)]
/// Generates the group GL_n(F_q) via a cached BFS.
/// The primay usage is
/// ```
///    let mut bfs = Bfs::new()
/// ```
pub struct GroupBFS {
    subgroup_generators: KTypeSubgroup,
    quotient: FFPolynomial,
    frontier: VecDeque<(GaloisMatrix, u32)>,
    // visited: FxHashSet<GroupBFSNode>,
    visited: FxHashMap<GaloisMatrix, GroupBFSNode>,
    hg: HGraph<NodeData, ()>,
    // coset_to_node: HashMap<CosetRep, u32>,
    current_bfs_distance: u32,
    last_flushed_distance: u32,
    num_matrices_completed: u32,
    last_cached_matrices_done: u64,
    cumulative_time_ms: u128,
    directory: PathBuf,
    filename: String,
    cache: bool,
    lookup: Arc<GaloisField>,
}

impl Serialize for GroupBFS {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        let mut s = serializer.serialize_struct("GroupBFS", 11)?;
        s.serialize_field("quotient", &self.quotient)?;
        s.serialize_field("frontier", &self.frontier)?;
        s.serialize_field("visited", &self.visited)?;
        // s.serialize_field("coset_to_node", &self.coset_to_node)?;
        s.serialize_field("current_distance", &self.current_bfs_distance)?;
        s.serialize_field("last_flushed_distance", &self.last_flushed_distance)?;
        s.serialize_field("num_matrices_completed", &self.num_matrices_completed)?;
        s.serialize_field("last_cached_matrices_done", &self.last_cached_matrices_done)?;
        s.serialize_field("cumulative_time_ms", &self.cumulative_time_ms)?;
        s.serialize_field("directory", &self.directory)?;
        s.serialize_field("filename", &self.filename)?;
        s.end()
    }
}

// impl<'de> Deserialize<'de> for Mat {
//     /// Note: will default to Edge::Undirected
//     fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
//     where
//         D: serde::Deserializer<'de>,
//     {
//         let mut data = <serde_json::Value>::deserialize(deserializer)?;
//         println!("data: {:?}", data);
//         let mut vec_data = String::from(data["data"].take().as_str().unwrap());
//         let dim = data["dim"].take().as_u64().unwrap() as usize;
//         dbg!(&dim);
//         dbg!(&vec_data);
//         if vec_data.starts_with("[") {
//             vec_data.remove(0);
//         }
//         if vec_data.ends_with("]") {
//             vec_data.remove(vec_data.len() - 1);
//         }
//         if vec_data.contains(",") {
//             let mut v: Vec<u32> = vec_data
//                 .split(',')
//                 .filter_map(|x| -> Option<u32> {
//                     if let Ok(number) = x.parse() {
//                         Some(number)
//                     } else {
//                         None
//                     }
//                 })
//                 .collect();
//             v.sort();
//             Ok(Mat { data: v, dim })
//         } else {
//             if let Ok(n) = vec_data.parse::<u32>() {
//                 Ok(Mat { data: vec![n], dim })
//             } else {
//                 if vec_data.len() == 0 {
//                     Ok(Mat { data: vec![], dim })
//                 } else {
//                     println!("vec_data: {:?}", vec_data);
//                     panic!("Could not parse single input.");
//                 }
//             }
//         }
//     }
// }

// impl Serialize for Mat {
//     fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
//     where
//         S: serde::Serializer,
//     {
//         let mut s = String::new();
//         if self.data.len() == 0 {
//             s.push_str("[]");
//         } else {
//             s.push_str("[");
//             for node in self.data.iter() {
//                 s.push_str(&format!("{:?},", node));
//             }
//             if s.ends_with(',') {
//                 s.remove(s.len() - 1);
//             }
//             s.push_str("]");
//         }

//         let mut ser = serializer.serialize_struct("Mat", 2)?;
//         ser.serialize_field("data", &s[..])?;
//         ser.serialize_field("dim", &self.dim)?;
//         ser.end()
//     }
// }

impl GroupBFS {
    pub fn new(directory: &Path, filename: String, quotient: &FFPolynomial, cache: bool) -> Self {
        println!("Checking for existing cache.");
        if let Some(cache) = GroupBFS::from_cache(directory, filename.clone()) {
            println!(
                "Successfully loaded GroupBFS from cache in directory: {:}",
                directory.to_str().unwrap()
            );
            if cache.quotient != *quotient {
                println!(
                    "Inconsistent quotients found. Provided: {:}. Cached: {:}.",
                    quotient, cache.quotient
                );
                panic!()
            }
            println!("size of frontier found: {:}", cache.frontier.len());
            println!("Size of visited found: {:}", cache.visited.len());
            cache
        } else {
            println!("No cache found, creating new GroupBFS.");
            let mut pbuf = PathBuf::new();
            pbuf.push(directory);
            let lookup = Arc::new(GaloisField::new(quotient.clone()));

            let mut ret = Self {
                subgroup_generators: KTypeSubgroup::new(&lookup),
                quotient: quotient.clone(),
                frontier: VecDeque::new(),
                visited: FxHashMap::default(),
                // visited: HashSet::new(),
                hg: HGraph::new(),
                // coset_to_node: HashMap::new(),
                current_bfs_distance: 0,
                last_flushed_distance: 0,
                num_matrices_completed: 0,
                last_cached_matrices_done: 0,
                cumulative_time_ms: 0,
                directory: pbuf,
                filename,
                cache,
                lookup,
            };
            ret.frontier.push_front((GaloisMatrix::id(3), 0));
            // ret.visited.insert(bfs_start);
            ret
        }
    }

    pub fn hgraph(&self) -> &HGraph<NodeData, ()> {
        &self.hg
    }

    pub fn field_mod(&self) -> FFRep {
        self.lookup.field_mod
    }

    fn from_cache(directory: &Path, filename: String) -> Option<Self> {
        // check if directory is a directory
        println!("retrieving cache from: {:}", directory.to_str().unwrap());
        println!("is_dir: {:}", directory.is_dir());

        if directory.is_dir() {
            let mut cache_path = PathBuf::new();
            cache_path.push(directory);
            cache_path.push("hdx.cache");
            println!("cache_path: {:}", cache_path.to_str().unwrap());
            if cache_path.is_file() {
                // todo: serialize this to JsonObject
                let file_data = fs::read_to_string(&cache_path).expect("Could not read file.");
                println!(
                    "Attempting to read GroupBFSCache from {:}",
                    cache_path.to_str().unwrap()
                );
                let cache: GroupBFSCache =
                    serde_json::from_str(&file_data).expect("could not parse serde");
                println!("GroupBFSCache parse successful.");
                let mut hg_path = PathBuf::new();
                hg_path.push(directory);
                hg_path.push(&filename[..]);
                let hg = HGraph::from_file(&hg_path).expect("Could not get hgraph.");
                let lookup = Arc::new(GaloisField::new(cache.quotient.clone()));
                let h_gens = KTypeSubgroup::new(&lookup);
                let bfs = GroupBFS {
                    subgroup_generators: h_gens,
                    quotient: cache.quotient,
                    frontier: cache.frontier,
                    visited: cache.visited,
                    hg,
                    // coset_to_node: cache.coset_to_node,
                    current_bfs_distance: cache.current_distance,
                    last_flushed_distance: cache.last_flushed_distance,
                    num_matrices_completed: cache.num_matrices_completed,
                    last_cached_matrices_done: cache.last_cached_matrices_done,
                    cumulative_time_ms: cache.cumulative_time_ms,
                    directory: cache.directory,
                    filename: cache.filename,
                    cache: true,
                    lookup,
                };
                return Some(bfs);
            }
        }
        None
    }

    pub fn cache(&self) {
        // serialiaze self first, then leave lookup table out and store hgraph
        // to disk.
        let mut tmp_path = self.directory.clone();
        tmp_path.push("tmp.cache");
        println!("Caching.");
        match serde_json::to_string(self) {
            Ok(serialized_self) => {
                println!("Serialized data, writing to disk.");
                let mut old_cache_path = self.directory.clone();
                old_cache_path.push("hdx.cache");
                thread::spawn(move || match fs::write(&tmp_path, serialized_self) {
                    Ok(_) => {
                        fs::rename(tmp_path, old_cache_path).expect("Could not rename file");
                        println!("Succesfully cached.");
                    }
                    Err(_) => {
                        println!("Failed to write cache.");
                    }
                });
            }
            Err(e) => {
                println!("Could not serialize GroupBFS. Serde error: {:}", e);
            }
        }
        println!("Serializing Graph");
        let mut hg_path = self.directory.clone();
        hg_path.push(&self.filename[..]);
        self.hg.to_disk(&hg_path);
        println!("Graph saved to disk.");
        println!("Done caching.");
        // todo: if I can't serialize the graph I should probably delete the cache?
    }

    fn clear_cache(&self) {
        let mut cache_path = self.directory.clone();
        cache_path.push("hdx.cache");
        let remove_result = fs::remove_file(cache_path);
        if remove_result.is_ok() {
            println!("Removed cache file.");
        } else {
            println!(
                "Cache file could not be removed. Error: {:?}",
                remove_result
            );
        }
    }

    fn flush(&mut self) {
        // let mut to_flush = Vec::new();
        let mut to_save = Vec::new();
        for bfs_node in self.visited.drain() {
            if bfs_node.1.distance >= self.current_bfs_distance - 1 {
                to_save.push(bfs_node);
            }
        }
        for t in to_save.into_iter() {
            self.visited.insert(t.0, t.1);
        }

        // for t in to_flush {
        //     let (c0, c1, c2) = self.subgroup_generators.get_coset_reps(&t.mat);
        //     self.coset_to_node.remove(&c0);
        //     self.coset_to_node.remove(&c1);
        //     self.coset_to_node.remove(&c2);
        // }
        self.last_flushed_distance = self.current_bfs_distance;
    }
    /// return the full path of the file used to store the
    /// resulting hgraph
    pub fn get_hgraph_file_path(&self) -> PathBuf {
        let mut filename = String::new();
        filename.push_str("p_");
        filename.push_str(&self.quotient.field_mod.to_string());
        filename.push_str("_dim_");
        filename.push_str("3");
        filename.push_str("_deg_");
        filename.push_str(&self.quotient.degree().to_string());
        filename.push_str(".hg");
        let mut file_path = self.directory.clone();
        file_path.push(&filename);
        file_path
    }

    pub fn bfs(&mut self, num_steps: usize) {
        let p = self.quotient.field_mod;
        let dim = 3;
        let deg = self.quotient.degree();
        let estimated_num_matrices = (p as u64).pow(deg * (dim * dim - 1));

        let cache_step_size = estimated_num_matrices / 32;
        let trace_step_size = 40_000;
        let start_time = Instant::now();

        let mut counter = 0;
        log::trace!("Starting BFS.");
        log::trace!("Estimated number matrices: {:}", estimated_num_matrices);
        log::trace!("cache step_size: {:}", cache_step_size);
        while self.frontier.is_empty() == false && counter < num_steps {
            self.step();
            counter += 1;
            if self.current_bfs_distance > self.last_flushed_distance + 1 {
                log::trace!("Flushing, current distance: {:}", self.current_bfs_distance);
                log::trace!("Completed {:} matrices", self.num_matrices_completed);
                self.flush();
            }
            if self.last_cached_matrices_done + cache_step_size < self.num_matrices_completed.into()
                && self.cache
            {
                log::trace!("Caching.");
                self.cache();
                self.last_cached_matrices_done = self.num_matrices_completed as u64;
            }
            if self.num_matrices_completed % trace_step_size == 0 {
                let time_since_start = start_time.elapsed().as_secs_f64();
                let time_per_matrix = time_since_start / self.num_matrices_completed as f64;
                let estimated_time_remaining =
                    (estimated_num_matrices - self.num_matrices_completed as u64) as f64
                        * time_per_matrix;
                log::info!("Time elapsed: {:} seconds", time_since_start);
                log::info!("Matrices processed: {:}", self.num_matrices_completed);
                log::info!("Time per matrix: {:}", time_per_matrix);
                log::info!("Estimated time remaining: {:}", estimated_time_remaining);
                log::info!("{:}", "$".repeat(65));
            }
        }
        let time_taken = start_time.elapsed().as_secs_f64();
        log::trace!("{:}", "@".repeat(65));
        log::trace!("Succesfully completed BFS!");
        log::trace!("Time taken (secs): {:}", time_taken);
        log::trace!("Matrices processed: {:}", self.num_matrices_completed);
        log::trace!(
            "Seconds per matrix: {:}",
            time_taken / self.num_matrices_completed as f64
        );

        // TODO: Need to get rid of saving this to disk during the BFS call.
        // should be a separate call so that way bfs can be used without side-effects.
        log::trace!("Saving complex to disk");
        let mut hg_path = self.directory.clone();
        hg_path.push(&self.filename[..]);
        self.hg.to_disk(&hg_path);
        // crate::coset_complex_to_disk(&self.hg, hg_path);
        if self.cache {
            log::trace!("Clearing cache.");
            self.clear_cache();
        }
        log::trace!("All done.");
    }

    /// Returns a partially completed, or 'trimmed', breadth first search
    /// of the coset complex for GL_n(F_q).
    pub fn trimmed_bfs(&mut self, num_steps: usize) {
        self.bfs(num_steps);
        let num_triangles_complete_border = self.field_mod() as usize;
        let mut rng = thread_rng();
        let mut nodes = self.hg.nodes();
        nodes.shuffle(&mut rng);
        let get_gluers = |hg: &HGraph<u16, ()>| {
            for ix in 0..nodes.len() {
                let ix_type = hg.get_node(&nodes[ix]);
                if ix_type.is_none() {
                    continue;
                }
                if *ix_type.unwrap() != 0 {
                    continue;
                }
                let node_ix_link = hg.link_of_nodes(&[nodes[ix]]);
                for (triangle_id, border_nodes_ix) in node_ix_link {
                    let border_id_ix = hg.find_id(&border_nodes_ix[..]);
                    if border_id_ix.is_none() {
                        continue;
                    }
                    let border_id_ix = border_id_ix.unwrap();
                    let link_border_ix = hg.link_of_nodes(&border_nodes_ix[..]);
                    if link_border_ix.len() < num_triangles_complete_border {
                        'jx_loop: for jx in ix..nodes.len() {
                            let jx_type = hg.get_node(&nodes[jx]);
                            if jx_type.is_none() {
                                continue 'jx_loop;
                            }
                            if *jx_type.unwrap() != 0 {
                                continue;
                            }
                            for (_, border_link_nodes) in link_border_ix.iter() {
                                if border_link_nodes.contains(&nodes[jx]) {
                                    continue 'jx_loop;
                                }
                            }
                            let node_jx_link = hg.link_of_nodes([nodes[jx]]);
                            for (_, jx_border_nodes) in node_jx_link {
                                let link_border_jx = hg.link_of_nodes(&jx_border_nodes[..]);
                                if link_border_ix.len() + link_border_jx.len()
                                    <= num_triangles_complete_border
                                {
                                    let border_id_jx = hg.find_id(&jx_border_nodes[..]);
                                    if border_id_jx.is_none() {
                                        continue;
                                    }
                                    let border_id_jx = border_id_jx.unwrap();
                                    let l1 = self.hg.link_of_nodes([border_nodes_ix[0]]);
                                    let l2 = self.hg.link_of_nodes([border_nodes_ix[1]]);
                                    let l3 = self.hg.link_of_nodes([jx_border_nodes[0]]);
                                    let l4 = self.hg.link_of_nodes([jx_border_nodes[1]]);
                                    return Some((border_id_ix, border_id_jx));
                                }
                            }
                        }
                    }
                }
            }
            None
        };

        loop {
            let e = get_gluers(&self.hg);
            if e.is_none() {
                break;
            }
            let (edge1, edge2) = e.unwrap();
            let nodes1 = self.hg.query_edge(&edge1);
            if nodes1.is_none() {
                log::trace!("Current nodes1 for {:} is busted?", edge1);
            }
            let nodes1 = nodes1.unwrap();
            assert!(nodes1.len() == 2);
            let (f1, f2) = (nodes1[0], nodes1[1]);
            let t1 = self.hg.get_node(&f1).unwrap();
            let t2 = self.hg.get_node(&f2).unwrap();
            let (f1, f2) = match (t1, t2) {
                (1, 2) => (f1, f2),
                (2, 1) => (f2, f1),
                _ => {
                    panic!("Incorrect types discovered on border check nodes.");
                }
            };

            let nodes2 = self.hg.query_edge(&edge2);
            if nodes2.is_none() {
                log::trace!("current_glue_nodes for {:} is busted?", edge2);
                continue;
            }
            let nodes2 = nodes2.unwrap();
            assert!(nodes2.len() == 2);
            let (g1, g2) = (nodes2[0], nodes2[1]);
            let t1 = *self.hg.get_node(&g1).unwrap();
            let t2 = *self.hg.get_node(&g2).unwrap();
            let (g1, g2) = match (t1, t2) {
                (1, 2) => (g1, g2),
                (2, 1) => (g2, g1),
                _ => {
                    panic!("Incorrect types discovered on border check nodes.");
                }
            };
            log::trace!("Gluing together edges {:} and {:}", edge1, edge2);
            println!("{:?}", self.hg.link_of_nodes([f1]));
            println!("{:?}", self.hg.link_of_nodes([f2]));
            // self.hg.concatenate_nodes(&g1, &f1);
            // self.hg.concatenate_nodes(&g2, &f2);
            println!("after{:}", "&".repeat(40));
            println!("{:?}", self.hg.link_of_nodes([f1]));
            println!("{:?}", self.hg.link_of_nodes([f2]));
            let num_triangles = self.hg.link_of_nodes([f1, f2]).len();
            if num_triangles == num_triangles_complete_border {
                log::trace!("Completed a border check between nodes {:} - {:}", f1, f2);
            }
        }

        log::trace!("Saving complex to disk");
        let mut hg_path = self.directory.clone();
        hg_path.push(&self.filename[..]);
        self.hg.to_disk(&hg_path);
    }

    /// Computes the border checks for all fully discovered local checks.
    pub fn get_border_checks(&self) -> Vec<u64> {
        let mut ret = FxHashSet::default();
        for node in self.hg.nodes() {
            if *self.hg.get_node(&node).unwrap() == 0 {
                println!("Found node: {:}", node);
                let link = self.hg.link_of_nodes([node]);
                let num_borders = link.iter().filter(|(_, nodes)| nodes.len() == 2).count();
                if num_borders == self.field_mod().pow(3) as usize {
                    link.into_iter()
                        .filter_map(|(id, nodes)| {
                            if nodes.len() == 2 {
                                self.hg.find_id(&nodes[..])
                            } else {
                                None
                            }
                        })
                        .for_each(|id| {
                            ret.insert(id);
                        });
                }
            }
        }
        ret.into_iter().collect()
    }
    /// Returns the nodes associated with the triangle found during the step. Will
    /// return an empty vec if the frontier is empty
    pub fn step(&mut self) -> Vec<u32> {
        if self.frontier.is_empty() {
            return Vec::new();
        }
        let x = self.frontier.pop_front().unwrap();
        // process this matrix first, compute the cosets and triangles it can
        // be a part of.
        let neighbors = self.subgroup_generators.generate_right_mul(&x.0);
        let (c0, c1, c2) = self.subgroup_generators.coset_reps(&neighbors[..]);

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
                    if let Some(x) = self.hg.get_node(n) {
                        *x == 0
                    } else {
                        false
                    }
                })
                .cloned()
                .collect();
            if v.len() == 0 {
                let new_node = self.hg.add_node(0);
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
                    if let Some(x) = self.hg.get_node(n) {
                        *x == 1
                    } else {
                        false
                    }
                })
                .cloned()
                .collect();
            if v.len() == 0 {
                let new_node = self.hg.add_node(1);
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
                    if let Some(x) = self.hg.get_node(n) {
                        *x == 2
                    } else {
                        false
                    }
                })
                .cloned()
                .collect();
            if v.len() == 0 {
                let new_node = self.hg.add_node(2);
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

        self.hg.add_edge(&[n0, n1], ());
        self.hg.add_edge(&[n0, n2], ());
        self.hg.add_edge(&[n1, n2], ());
        self.hg.add_edge(&[n0, n1, n2], ());
        vec![n0, n1, n2]
    }
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

    use super::{bfs, GroupBFS, GroupBFSNode};

    fn simple_quotient_and_field() -> (u32, FFPolynomial) {
        let p = 3_u32;
        let primitive_coeffs = [(2, (1, p).into()), (1, (2, p).into()), (0, (2, p).into())];
        let q = FFPolynomial::from(&primitive_coeffs[..]);
        (p, q)
    }

    fn simple_group_bfs() -> GroupBFS {
        let p = 3_u32;
        let primitive_coeffs = [(2, (1, p).into()), (1, (2, p).into()), (0, (2, p).into())];
        let q = FFPolynomial::from(&primitive_coeffs[..]);
        let directory = PathBuf::from_str("/Users/matt/repos/qec/tmp").unwrap();
        GroupBFS::new(&directory, String::from("tester"), &q, true)
    }

    #[test]
    fn test_group_bfs_manager_new() {
        let p = 3_u32;
        let primitive_coeffs = [(2, (1, p).into()), (1, (2, p).into()), (0, (2, p).into())];
        let q = FFPolynomial::from(&primitive_coeffs[..]);
        let directory = PathBuf::from_str("/Users/matt/repos/qec/tmp/").unwrap();
        let mut bfs_manager = GroupBFS::new(&directory, String::from("tester"), &q, true);
        bfs_manager.bfs(usize::MAX);
    }

    #[test]
    fn trimmer_test() {
        let logger = SimpleLogger::new().init().unwrap();
        let p = 3_u32;
        let primitive_coeffs = [(2, (1, p).into()), (1, (2, p).into()), (0, (2, p).into())];
        let q = FFPolynomial::from(&primitive_coeffs[..]);
        let directory = PathBuf::from_str("/Users/matt/repos/qec/tmp/").unwrap();
        let mut bfs_manager = GroupBFS::new(&directory, String::from("tester"), &q, true);
        bfs_manager.bfs(50);
        // bfs_manager.trimmed_bfs(250);
    }

    #[test]
    fn test_short_walk() {
        let dir = PathBuf::from("/Users/matt/repos/qec/tmp");
        let p = 3_u32;
        let primitive_coeffs = [(2, (1, p).into()), (1, (2, p).into()), (0, (2, p).into())];
        let q = FFPolynomial::from(&primitive_coeffs[..]);
        let mut bfs = GroupBFS::new(&dir, String::from("tester"), &q, false);
        bfs.bfs((2 as usize).pow(10));
        println!("graph: {:}", bfs.hg);
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
