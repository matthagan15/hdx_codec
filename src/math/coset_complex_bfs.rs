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

use fxhash::{FxHashMap, FxHashSet};
use log::trace;
use mhgl::{ConGraph, HGraph, HyperGraph};
use rand::{seq::SliceRandom, thread_rng};
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

const BFS_FILENAME: &str = "hdx_bfs.cache";

/// Helper struct for BFS
#[derive(Debug, Clone, Serialize, Deserialize)]
struct GroupBFSNode {
    distance: u32,
    hg_node: Vec<u32>,
}
// impl PartialEq for GroupBFSNode {
//     fn eq(&self, other: &Self) -> bool {
//         self.mat == other.mat
//     }
// }
// impl Eq for GroupBFSNode {}
// impl Hash for GroupBFSNode {
//     fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
//         self.mat.hash(state);
//     }
// }
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

        // To get the incomplete border checks I first want to filter
        // only the border checks that have completely filled out local checks.
        // need to take all existing border checks and randomly glue them together.
        let border_checks = self.get_border_checks();
        // let num_triangles_complete_border = (2 * self.field_mod().pow(2)) as usize;
        let num_triangles_complete_border = self.field_mod() as usize;
        let mut incomplete_border_checks: Vec<(u64, usize)> = border_checks
            .into_iter()
            .filter_map(|check| {
                let num_triangles = self.hg.maximal_edges(&check).len();
                if num_triangles < num_triangles_complete_border {
                    Some((check, num_triangles))
                } else {
                    None
                }
            })
            .collect();
        let mut rng = thread_rng();
        incomplete_border_checks.shuffle(&mut rng);
        let incomplete_copy = incomplete_border_checks.clone();
        // used to keep track of which border checks have been glued or fixed already.
        let mut already_used_borders: HashSet<u64> = HashSet::new();
        let mut newly_completed_borders = Vec::new();
        while incomplete_border_checks.is_empty() == false {
            let (current_fixer, mut fixer_num_triangles) = incomplete_border_checks.pop().unwrap();
            if already_used_borders.contains(&current_fixer) {
                continue;
            }
            already_used_borders.insert(current_fixer);
            for ix in 0..incomplete_copy.len() {
                if already_used_borders.contains(&incomplete_copy[ix].0) {
                    continue;
                }
                let (to_glue, glue_num_triangles) = &incomplete_copy[ix];
                // need to check to make sure that to_glue and current_fixer do not share a
                // green vertex in common
                let fix_link = self.hgraph().link(&current_fixer);
                let glue_link = self.hgraph().link(to_glue);
                let glue_link_set = glue_link
                    .into_iter()
                    .map(|(id, nodes)| {
                        if nodes.len() != 1 {
                            panic!("link of a border check has more than one node.")
                        }
                        nodes[0]
                    })
                    .fold(HashSet::new(), |mut acc, x| {
                        acc.insert(x);
                        acc
                    });
                let mut skip_glue = false;
                for (_, nodes) in fix_link {
                    if nodes.len() != 1 {
                        panic!("link of a border check has more than one node.")
                    }
                    if glue_link_set.contains(&nodes[0]) {
                        skip_glue = true;
                        break;
                    }
                }
                if skip_glue {
                    continue;
                }
                if glue_num_triangles + fixer_num_triangles <= num_triangles_complete_border {
                    // glue
                    log::trace!(
                        "current_fixer {:} has {:} / {num_triangles_complete_border} triangles",
                        current_fixer,
                        fixer_num_triangles
                    );
                    let current_fixer_nodes = self.hg.query_edge(&current_fixer).unwrap();
                    assert!(current_fixer_nodes.len() == 2);
                    let (f1, f2) = (current_fixer_nodes[0], current_fixer_nodes[1]);
                    let t1 = self.hg.get_node(&f1).unwrap();
                    let t2 = self.hg.get_node(&f2).unwrap();
                    let (f1, f2) = match (t1, t2) {
                        (1, 2) => (f1, f2),
                        (2, 1) => (f2, f1),
                        _ => {
                            panic!("Incorrect types discovered on border check nodes.");
                        }
                    };

                    let current_glue_nodes = self.hg.query_edge(&to_glue);
                    if current_glue_nodes.is_none() {
                        log::trace!("current_glue_nodes for {:} is busted?", to_glue);
                        continue;
                    }
                    let current_glue_nodes = current_glue_nodes.unwrap();
                    assert!(current_glue_nodes.len() == 2);
                    let (g1, g2) = (current_glue_nodes[0], current_glue_nodes[1]);
                    let t1 = self.hg.get_node(&g1).unwrap();
                    let t2 = self.hg.get_node(&g2).unwrap();
                    let (g1, g2) = match (t1, t2) {
                        (1, 2) => (g1, g2),
                        (2, 1) => (g2, g1),
                        _ => {
                            panic!("Incorrect types discovered on border check nodes.");
                        }
                    };
                    log::trace!("Gluing together edges {:} and {:}", current_fixer, to_glue);
                    self.hg.concatenate_nodes(&g1, &f1);
                    self.hg.concatenate_nodes(&g2, &f2);
                    fixer_num_triangles += glue_num_triangles;
                    already_used_borders.insert(*to_glue);
                    log::trace!(
                        "num triangles of fixer: {:} / {num_triangles_complete_border}",
                        fixer_num_triangles
                    );
                    if fixer_num_triangles == num_triangles_complete_border {
                        newly_completed_borders.push(current_fixer);
                        log::trace!("Completed fixer {:}!", current_fixer);
                        break;
                    }
                } else {
                    continue;
                }
            }
        }
        log::trace!(
            "was able to fill out {:} incomplete borders",
            newly_completed_borders.len()
        );
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
        bfs_manager.trimmed_bfs(500);
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
}
