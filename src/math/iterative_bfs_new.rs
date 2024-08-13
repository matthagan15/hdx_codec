use core::panic;
use std::{
    collections::{HashMap, HashSet, VecDeque},
    fs::{self, File},
    hash::Hash,
    io::{Read, Write},
    ops::Mul,
    path::{Path, PathBuf},
    rc::Rc,
    sync::Arc,
    thread,
    time::{Duration, Instant},
};

use log::trace;
use mhgl::{ConGraph, HGraph, HyperGraph};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use serde::ser::SerializeStruct;
use serde::{Deserialize, Serialize};
use serde_json::Value;

use crate::{
    matrices::sparse_ffmatrix::SparseVector, reed_solomon::ReedSolomon, tanner_code::TannerCode,
};

use super::{
    coset_complex_subgroups::{CosetRep, KTypeSubgroup},
    finite_field::{FFRep, FiniteField},
    galois_field::GaloisField,
    polynomial::{FFPolynomial, PolyDegree},
};
use crate::matrices::{galois_matrix::GaloisMatrix, polymatrix::PolyMatrix};

const BFS_FILENAME: &str = "hdx_bfs.cache";

/// Helper struct for BFS
#[derive(Debug, Clone, Serialize, Deserialize)]
struct GroupBFSNode {
    mat: GaloisMatrix,
    distance: u32,
}
impl PartialEq for GroupBFSNode {
    fn eq(&self, other: &Self) -> bool {
        self.mat == other.mat
    }
}
impl Eq for GroupBFSNode {}
impl Hash for GroupBFSNode {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.mat.hash(state);
    }
}
#[derive(Debug, Serialize, Deserialize)]
pub struct GroupBFSCache {
    quotient: FFPolynomial,
    frontier: VecDeque<GroupBFSNode>,
    visited: HashSet<GroupBFSNode>,
    coset_to_node: HashMap<CosetRep, u32>,
    current_distance: u32,
    last_flushed_distance: u32,
    num_matrices_completed: u32,
    last_cached_matrices_done: u64,
    cumulative_time_ms: u128,
    directory: PathBuf,
    filename: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NodeData {
    type_ix: u16,
}

#[derive(Debug)]
pub struct GroupBFS {
    subgroup_generators: KTypeSubgroup,
    quotient: FFPolynomial,
    frontier: VecDeque<GroupBFSNode>,
    visited: HashSet<GroupBFSNode>,
    hg: HGraph<NodeData, ()>,
    coset_to_node: HashMap<CosetRep, u32>,
    current_bfs_distance: u32,
    last_flushed_distance: u32,
    num_matrices_completed: u32,
    last_cached_matrices_done: u64,
    cumulative_time_ms: u128,
    directory: PathBuf,
    filename: String,
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
        s.serialize_field("coset_to_node", &self.coset_to_node)?;
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

impl GroupBFS {
    pub fn new(directory: &Path, filename: String, quotient: &FFPolynomial) -> Self {
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
                visited: HashSet::new(),
                hg: HGraph::new(),
                coset_to_node: HashMap::new(),
                current_bfs_distance: 0,
                last_flushed_distance: 0,
                num_matrices_completed: 0,
                last_cached_matrices_done: 0,
                cumulative_time_ms: 0,
                directory: pbuf,
                filename,
                lookup,
            };
            let bfs_start = GroupBFSNode {
                mat: GaloisMatrix::id(3),
                distance: 0,
            };
            ret.frontier.push_front(bfs_start.clone());
            ret.visited.insert(bfs_start);
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
                    coset_to_node: cache.coset_to_node,
                    current_bfs_distance: cache.current_distance,
                    last_flushed_distance: cache.last_flushed_distance,
                    num_matrices_completed: cache.num_matrices_completed,
                    last_cached_matrices_done: cache.last_cached_matrices_done,
                    cumulative_time_ms: cache.cumulative_time_ms,
                    directory: cache.directory,
                    filename: cache.filename,
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
        let mut to_flush = Vec::new();
        let mut to_save = HashSet::new();
        for bfs_node in self.visited.drain() {
            if bfs_node.distance < self.current_bfs_distance - 1 {
                to_flush.push(bfs_node);
            } else {
                to_save.insert(bfs_node);
            }
        }
        for t in to_save.into_iter() {
            self.visited.insert(t);
        }

        for t in to_flush {
            let (c0, c1, c2) = self.subgroup_generators.get_coset_reps(&t.mat);
            self.coset_to_node.remove(&c0);
            self.coset_to_node.remove(&c1);
            self.coset_to_node.remove(&c2);
        }
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

    pub fn bfs(&mut self, max_num_steps: usize) {
        let p = self.quotient.field_mod;
        let dim = 3;
        let deg = self.quotient.degree();
        let estimated_num_matrices = (p as u64).pow(deg * (dim * dim - 1));

        let cache_step_size = estimated_num_matrices / 32;
        let trace_step_size = 40_000;
        let start_time = Instant::now();
        let mut counter = 0;
        println!("Starting BFS.");
        println!("Estimated number matrices: {:}", estimated_num_matrices);
        println!("cache step_size: {:}", cache_step_size);
        let mut thread_handles = Vec::new();
        while self.frontier.is_empty() == false && counter < max_num_steps {
            self.step();
            counter += 1;
            if self.current_bfs_distance > self.last_flushed_distance + 1 {
                println!("Flushing, current distance: {:}", self.current_bfs_distance);
                println!("Completed {:} matrices", self.num_matrices_completed);
                self.flush();
            }
            if self.last_cached_matrices_done + cache_step_size < self.num_matrices_completed.into()
            {
                println!("Caching.");
                self.cache();
                self.last_cached_matrices_done = self.num_matrices_completed as u64;
            }
            if self.num_matrices_completed % trace_step_size == 0 {
                let time_since_start = start_time.elapsed().as_secs_f64();
                let time_per_matrix = time_since_start / self.num_matrices_completed as f64;
                let estimated_time_remaining =
                    (estimated_num_matrices - self.num_matrices_completed as u64) as f64
                        * time_per_matrix;
                println!("Time elapsed: {:} seconds", time_since_start);
                println!("Matrices processed: {:}", self.num_matrices_completed);
                println!("Time per matrix: {:}", time_per_matrix);
                println!("Estimated time remaining: {:}", estimated_time_remaining);
                println!("{:}", "$".repeat(65));
            }
            if self.num_matrices_completed % 100000 == 0 {
                println!("{:}", "*".repeat(70));
                println!("Generating TannerCode.");
                let new_hg = self.hg.clone();
                let mut filename = self.directory.clone();
                let mut name = String::from("bfs_");
                name.push_str(&self.num_matrices_completed.to_string());
                name.push_str(".hg");
                filename.push(name);
                let th = thread::spawn(move || {
                    let num_nodes = new_hg.num_nodes();
                    let num_edges = new_hg.edges_of_size(2).len();
                    let num_triangles = new_hg.edges_of_size(3).len();
                    new_hg.to_disk(&filename);

                    println!(
                        "ConGraph saved to disk. Number of nodes, edges, triangles: {:}, {:}, {:}.",
                        num_nodes, num_edges, num_triangles
                    );
                });
                thread_handles.push(th);
                println!("Matrix thread returned.");
                println!("{:}", "*".repeat(70));
            }
        }
        let time_taken = start_time.elapsed().as_secs_f64();
        println!("{:}", "@".repeat(65));
        println!("Succesfully completed BFS!");
        println!("Time taken (secs): {:}", time_taken);
        println!("Matrices processed: {:}", self.num_matrices_completed);
        println!(
            "Seconds per matrix: {:}",
            time_taken / self.num_matrices_completed as f64
        );
        println!("Waiting for writer threads:");
        for th in thread_handles {
            th.join().expect("Thread did not join?");
        }
        println!("Saving ConGraph to disk");
        let mut hg_path = self.directory.clone();
        hg_path.push(&self.filename[..]);
        self.hg.to_disk(&hg_path);

        trace!("Clearing cache.");
        self.clear_cache();
        println!("All done.");
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
        // TODO: Can I get by with using just one subgroup generator pass?
        // Currently we do two multiplications, one to get the coset reps
        // and another to generate the neighbors. The question is then can
        // we get by with generating neighbors using right multiplication?
        let neighbors = self.subgroup_generators.generate_right_mul(&x.mat);
        let (c0, c1, c2) = self.subgroup_generators.coset_reps(&neighbors[..]);
        let n0 = if self.coset_to_node.contains_key(&c0) {
            *self.coset_to_node.get(&c0).unwrap()
        } else {
            let new_node = self.hg.add_node(NodeData { type_ix: 0 });
            self.coset_to_node.insert(c0, new_node);
            new_node
        };
        let n1 = if self.coset_to_node.contains_key(&c1) {
            *self.coset_to_node.get(&c1).unwrap()
        } else {
            let new_node = self.hg.add_node(NodeData { type_ix: 1 });
            self.coset_to_node.insert(c1, new_node);
            new_node
        };
        let n2 = if self.coset_to_node.contains_key(&c2) {
            *self.coset_to_node.get(&c2).unwrap()
        } else {
            let new_node = self.hg.add_node(NodeData { type_ix: 2 });
            self.coset_to_node.insert(c2, new_node);
            new_node
        };

        self.hg.add_edge(&[n0, n1], ());
        self.hg.add_edge(&[n0, n2], ());
        self.hg.add_edge(&[n1, n2], ());
        self.hg.add_edge(&[n0, n1, n2], ());

        // flush visited and coset to node
        self.current_bfs_distance = x.distance;
        // let neighbors = self.subgroup_generators.generate_left_mul(&x.mat);

        for neighbor in neighbors {
            let neighbor_bfs = GroupBFSNode {
                mat: neighbor,
                distance: x.distance + 1,
            };
            if self.visited.contains(&neighbor_bfs) == false {
                self.visited.insert(neighbor_bfs.clone());
                self.frontier.push_back(neighbor_bfs);
            }
        }
        self.num_matrices_completed += 1;
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

    use crate::math::{galois_field::GaloisField, polynomial::FFPolynomial};
    use crate::matrices::galois_matrix::GaloisMatrix;

    use super::{GroupBFS, GroupBFSNode};

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
        GroupBFS::new(&directory, String::from("tester"), &q)
    }

    #[test]
    fn test_group_bfs_manager_new() {
        let p = 3_u32;
        let primitive_coeffs = [(2, (1, p).into()), (1, (2, p).into()), (0, (2, p).into())];
        let q = FFPolynomial::from(&primitive_coeffs[..]);
        let directory = PathBuf::from_str("/Users/matt/repos/qec/tmp/").unwrap();
        let mut bfs_manager = GroupBFS::new(&directory, String::from("tester"), &q);
        bfs_manager.bfs(usize::MAX);
    }

    #[test]
    fn test_short_walk() {
        let dir = PathBuf::from("/Users/matt/repos/qec/tmp");
        let p = 3_u32;
        let primitive_coeffs = [(2, (1, p).into()), (1, (2, p).into()), (0, (2, p).into())];
        let q = FFPolynomial::from(&primitive_coeffs[..]);
        let mut bfs = GroupBFS::new(&dir, String::from("tester"), &q);
        bfs.bfs((2 as usize).pow(10));
        println!("graph: {:}", bfs.hg);
    }
}
