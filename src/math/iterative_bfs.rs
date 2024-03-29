use core::panic;
use std::{
    collections::{HashMap, HashSet, VecDeque},
    fs,
    hash::Hash,
    io::Write,
    path::{Path, PathBuf},
    rc::Rc,
    thread,
    time::Instant,
};

use mhgl::HGraph;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use serde::{Deserialize, Serialize};

use super::{
    finite_field::FiniteField,
    galois_matrix::GaloisMatrix,
    polymatrix::PolyMatrix,
    polynomial::{FiniteFieldPolynomial, PolyDegree},
};

const BFS_FILENAME: &str = "hdx_bfs.cache";

/// Helper function for creating subgroups for BFS generation.
pub fn compute_deg(dim: usize, type_ix: usize, row_ix: usize, col_ix: usize) -> PolyDegree {
    let mut how_many_over: HashMap<usize, usize> = HashMap::new();
    let mut row = type_ix;
    let mut hoppable_cols_left = (dim - 1) as i32;
    for _ in 0..dim {
        let e = how_many_over.entry(row).or_default();
        *e = hoppable_cols_left as usize;
        hoppable_cols_left -= 1;
        row += 1;
        row %= dim;
    }
    let hoppable_cols = how_many_over[&row_ix];
    let mut dist_from_diag = col_ix as i32 - row_ix as i32;
    while dist_from_diag < 0 {
        dist_from_diag += dim as i32;
    }
    dist_from_diag %= dim as i32;
    if dist_from_diag > (hoppable_cols as i32) {
        0
    } else {
        dist_from_diag.try_into().unwrap()
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct CosetRepNew {
    rep: GaloisMatrix,
    type_ix: u16,
}

impl Serialize for CosetRepNew {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        let mut s = String::new();
        let rep_serde = serde_json::to_string(&self.rep);
        if let Ok(rep_s) = rep_serde {
            s.push_str(&rep_s);
        } else {
            panic!("cannot serialize matrix rep")
        }
        s.push('@');
        s.push_str(&self.type_ix.to_string());
        serializer.serialize_str(&s)
    }
}

impl<'de> Deserialize<'de> for CosetRepNew {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let data = <String>::deserialize(deserializer)?;
        let splitted: Vec<&str> = data.split('@').collect();
        if splitted.len() != 2 {
            panic!("improper CosetRep encountered.")
        }
        let matrix = splitted[0];
        let type_ix = splitted[1].parse::<u16>();
        let mat_out = serde_json::from_str::<GaloisMatrix>(matrix);
        if mat_out.is_err() {
            panic!("Could not deserialize matrix");
        }
        let mat = mat_out.unwrap();
        Ok(CosetRepNew {
            rep: mat,
            type_ix: type_ix.expect("Could not parse CosetRep type"),
        })
    }
}

/// The canonical representative of a coset, meaning we can
/// hash it directly without reference to subgroups.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct CosetRep {
    rep: PolyMatrix,
    type_ix: usize,
}

impl Serialize for CosetRep {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        let mut s = String::new();
        let rep_serde = serde_json::to_string(&self.rep);
        if let Ok(rep_s) = rep_serde {
            s.push_str(&rep_s);
        } else {
            panic!("cannot serialize matrix rep")
        }
        s.push('@');
        s.push_str(&self.type_ix.to_string());
        serializer.serialize_str(&s)
    }
}

impl<'de> Deserialize<'de> for CosetRep {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let data = <String>::deserialize(deserializer)?;
        let splitted: Vec<&str> = data.split('@').collect();
        if splitted.len() != 2 {
            panic!("improper CosetRep encountered.")
        }
        let matrix = splitted[0];
        let type_ix = splitted[1].parse::<usize>();
        let mat_out = serde_json::from_str::<PolyMatrix>(matrix);
        if mat_out.is_err() {
            panic!("Could not deserialize matrix");
        }
        let mat = mat_out.unwrap();
        Ok(CosetRep {
            rep: mat,
            type_ix: type_ix.expect("Could not parse CosetRep type"),
        })
    }
}

#[derive(Debug, Clone, Deserialize, Serialize)]
struct HTypeSubgroup {
    gens: Vec<PolyMatrix>,
    type_zero_end: usize,
    type_one_end: usize,
}

impl HTypeSubgroup {
    pub fn new(quotient: &FiniteFieldPolynomial) -> Self {
        let mut type_zero = h_type_subgroup(0, quotient);
        let mut type_one = h_type_subgroup(1, quotient);
        let mut type_two = h_type_subgroup(2, quotient);
        let mut gens = Vec::with_capacity(type_zero.len() + type_one.len() + type_two.len());
        gens.append(&mut type_zero);
        let type_zero_end = gens.len() - 1;
        gens.append(&mut type_one);
        let type_one_end = gens.len() - 1;
        gens.append(&mut type_two);
        Self {
            gens,
            type_zero_end,
            type_one_end,
        }
    }

    pub fn generate_left_mul(&self, mat: &PolyMatrix) -> Vec<PolyMatrix> {
        self.gens.par_iter().map(|h| h * mat).collect()
    }

    pub fn generate_right_mul(&self, mat: &PolyMatrix) -> Vec<PolyMatrix> {
        self.gens.par_iter().map(|h| mat * h).collect()
    }
}

fn h_type_subgroup(type_ix: usize, quotient: &FiniteFieldPolynomial) -> Vec<PolyMatrix> {
    let mut ret = Vec::new();
    let dim = 3;
    let p = quotient.field_mod;
    let id = PolyMatrix::id(dim, quotient.clone());
    let mut row_ix = type_ix as i32 - 1;
    while row_ix <= 0 {
        row_ix += dim as i32;
    }
    row_ix %= dim as i32;
    for a in 0..p {
        let mut tmp = id.clone();
        let e = tmp.get_mut(row_ix as usize, type_ix);
        *e = FiniteFieldPolynomial::monomial(FiniteField::new(a, p), 1);
        ret.push(tmp);
    }
    ret
}
/// Helper struct for BFS
#[derive(Debug, Clone, Serialize, Deserialize)]
struct GroupBFSNode {
    mat: PolyMatrix,
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
pub struct GroupBFS {
    subgroups: ParallelSubgroups,
    h_gens: HTypeSubgroup,
    quotient: FiniteFieldPolynomial,
    frontier: VecDeque<GroupBFSNode>,
    visited: HashSet<GroupBFSNode>,
    hg: HGraph,
    coset_to_node: HashMap<CosetRep, u32>,
    current_distance: u32,
    last_flushed_distance: u32,
    num_matrices_completed: u32,
    last_cached_matrices_done: u32,
    directory: PathBuf,
}

impl GroupBFS {
    pub fn new(directory: &Path, quotient: &FiniteFieldPolynomial) -> Self {
        let previously_cached = GroupBFS::load_from_disk(directory);
        if previously_cached.is_some() {
            println!(
                "Successfully loaded GroupBFS from cache in directory: {:}",
                directory.to_str().unwrap()
            );
            let ret = previously_cached.unwrap();
            println!("size of frontier found: {:}", ret.frontier.len());
            println!("Size of visited found: {:}", ret.visited.len());
            ret
        } else {
            println!("No cache found, creating new GroupBFS.");
            let mut pbuf = PathBuf::new();
            pbuf.push(directory);
            let mut ret = Self {
                subgroups: ParallelSubgroups::new(3, quotient),
                h_gens: HTypeSubgroup::new(quotient),
                quotient: quotient.clone(),
                frontier: VecDeque::new(),
                visited: HashSet::new(),
                hg: HGraph::new(),
                coset_to_node: HashMap::new(),
                current_distance: 0,
                last_flushed_distance: 0,
                num_matrices_completed: 0,
                last_cached_matrices_done: 0,
                directory: pbuf,
            };
            let e = PolyMatrix::id(3, quotient.clone());
            let bfs_start = GroupBFSNode {
                mat: e,
                distance: 0,
            };
            ret.frontier.push_front(bfs_start.clone());
            ret.visited.insert(bfs_start);
            ret.cache();
            ret
        }
    }

    fn load_from_disk(directory: &Path) -> Option<Self> {
        // check if directory is a directory
        if directory.is_dir() {
            let mut cache_path = PathBuf::new();
            cache_path.push(directory);
            cache_path.push(BFS_FILENAME);
            if cache_path.is_file() {
                let file_data = fs::read_to_string(cache_path).expect("Could not read file.");
                let serde_out = serde_json::from_str::<GroupBFS>(&file_data);
                if serde_out.is_ok() {
                    return Some(serde_out.unwrap());
                }
            }
        }
        None
    }

    fn cache(&self) {
        let mut tmp_path = self.directory.clone();
        tmp_path.push("tmp_bfs.cache");
        println!("Attempting to write to: {:}", tmp_path.display());
        let serde_out = serde_json::to_string(self);
        if serde_out.is_ok() {
            let mut old_cache = self.directory.clone();
            thread::spawn(move || {
                let file_out = fs::write(&tmp_path, serde_out.unwrap());
                if file_out.is_ok() {
                    // successfully wrote to disk
                    // rely on rename to delete the old cache and update name
                    old_cache.push(BFS_FILENAME);
                    fs::rename(tmp_path, old_cache).expect("Could not rename file");
                    println!("Succesfully cached.");
                } else {
                    println!("Could not write to temporary cache file.")
                }
            });
            println!("main thread returning to computation.");
        } else {
            println!("Could not serialize self. Serde error: {:?}", serde_out);
        }
    }

    fn clear_cache(&self) {
        let mut cache_path = self.directory.clone();
        cache_path.push(BFS_FILENAME);
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
            if bfs_node.distance < self.current_distance - 1 {
                to_flush.push(bfs_node);
            } else {
                to_save.insert(bfs_node);
            }
        }
        for t in to_save.into_iter() {
            self.visited.insert(t);
        }

        for t in to_flush {
            let (c0, c1, c2) = self.subgroups.get_coset_reps(&t.mat);
            self.coset_to_node.remove(&c0);
            self.coset_to_node.remove(&c1);
            self.coset_to_node.remove(&c2);
        }
        self.last_flushed_distance = self.current_distance;
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

    fn save_hgraph(&self) {
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

        if let Ok(s) = serde_json::to_string(&self.hg) {
            let mut file = fs::File::create(&file_path).expect("Could not create hgraph file.");
            file.write_all(s.as_bytes())
                .expect("Could not write hgraph to file.");
            println!(
                "HGraph successfully saved at: {:}",
                file_path.to_str().unwrap()
            );
        }
    }

    pub fn bfs(&mut self, max_num_steps: usize) {
        let p = self.quotient.field_mod;
        let dim = 3;
        let deg = self.quotient.degree();
        // let estimated_num_matrices = (p.pow(deg as u32)).pow(dim * dim - 1);
        let estimated_num_matrices = 1000;
        let cache_step_size = estimated_num_matrices / 32;
        let start_time = Instant::now();
        let mut counter = 0;
        println!("Starting BFS.");
        println!("Estimated number matrices: {:}", estimated_num_matrices);
        println!("cache step_size: {:}", cache_step_size);
        while self.frontier.is_empty() == false && counter < max_num_steps {
            self.step();
            counter += 1;
            if self.current_distance > self.last_flushed_distance + 1 {
                println!("Flushing, current distance: {:}", self.current_distance);
                println!("Completed {:} matrices", self.num_matrices_completed);
                self.flush();
            }
            if self.last_cached_matrices_done + cache_step_size < self.num_matrices_completed {
                // println!("Caching.");
                // self.cache();
                self.last_cached_matrices_done = self.num_matrices_completed;
            }
            if self.num_matrices_completed % 40_000 == 0 {
                let time_since_start = start_time.elapsed().as_secs_f64();
                let time_per_matrix = time_since_start / self.num_matrices_completed as f64;
                let estimated_time_remaining =
                    (estimated_num_matrices - self.num_matrices_completed) as f64 * time_per_matrix;
                println!("Time elapsed: {:} seconds", time_since_start);
                println!("Matrices processed: {:}", self.num_matrices_completed);
                println!("Time per matrix: {:}", time_per_matrix);
                println!("Estimated time remaining: {:}", estimated_time_remaining);
                println!("{:}", "$".repeat(65));
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
        println!("Saving HGraph to disk");
        self.save_hgraph();
        println!("Clearing Cache.");
        self.clear_cache();
        println!("All done.");
    }

    fn step(&mut self) {
        let x = self.frontier.pop_front().unwrap();
        // process this matrix first, compute the cosets and triangles it can
        // be a part of.
        let (c0, c1, c2) = self.subgroups.get_coset_reps(&x.mat);
        let n0 = if self.coset_to_node.contains_key(&c0) {
            *self.coset_to_node.get(&c0).unwrap()
        } else {
            let new_node = self.hg.add_node();
            self.coset_to_node.insert(c0, new_node);
            new_node
        };
        let n1 = if self.coset_to_node.contains_key(&c1) {
            *self.coset_to_node.get(&c1).unwrap()
        } else {
            let new_node = self.hg.add_node();
            self.coset_to_node.insert(c1, new_node);
            new_node
        };
        let n2 = if self.coset_to_node.contains_key(&c2) {
            *self.coset_to_node.get(&c2).unwrap()
        } else {
            let new_node = self.hg.add_node();
            self.coset_to_node.insert(c2, new_node);
            new_node
        };

        self.hg.create_edge_no_dups(&[n0, n1]);
        self.hg.create_edge_no_dups(&[n0, n2]);
        self.hg.create_edge_no_dups(&[n2, n1]);
        self.hg.create_edge_no_dups(&[n0, n1, n2]);

        // flush visited and coset to node
        self.current_distance = x.distance;

        let neighbors = self.h_gens.generate_left_mul(&x.mat);

        // how do I answer the question: have I seen this matrix before?
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
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct ParallelSubgroups {
    generators: Vec<PolyMatrix>,
    type_zero_end: usize,
    type_one_end: usize,
}

impl ParallelSubgroups {
    pub fn new(dim: usize, quotient: &FiniteFieldPolynomial) -> Self {
        let p = quotient.field_mod;
        let id = PolyMatrix::id(dim, quotient.clone());
        let mut type_to_gens = HashMap::new();

        for type_ix in 0..dim {
            let ret_type_i: &mut HashSet<PolyMatrix> = type_to_gens.entry(type_ix).or_default();
            ret_type_i.insert(id.clone());
            for row_ix in 0..dim {
                for col_ix in 0..dim {
                    if row_ix == col_ix {
                        continue;
                    }
                    let deg = compute_deg(dim, type_ix, row_ix, col_ix);
                    if deg > 0 {
                        let matrices = type_to_gens
                            .get_mut(&type_ix)
                            .expect("couldn't get matrices");
                        let mut new_matrices = HashSet::new();
                        for matrix in matrices.iter() {
                            for c in 1..p {
                                let mut new_matrix = matrix.clone();
                                let monomial = FiniteFieldPolynomial::monomial((c, p).into(), deg);
                                new_matrix.set_entry(row_ix, col_ix, &monomial);
                                new_matrices.insert(new_matrix);
                            }
                        }
                        for matrix in new_matrices {
                            matrices.insert(matrix);
                        }
                    }
                }
            }
        }
        let mut gens = Vec::new();
        let zero_gens = type_to_gens.remove(&0).expect("No type zero gens");
        let one_gens = type_to_gens.remove(&1).expect("No type zero gens");
        let two_gens = type_to_gens.remove(&2).expect("No type zero gens");
        for m in zero_gens.into_iter() {
            gens.push(m);
        }
        let type_zero_end = gens.len() - 1;
        for m in one_gens.into_iter() {
            gens.push(m);
        }
        let type_one_end = gens.len() - 1;
        for m in two_gens.into_iter() {
            gens.push(m);
        }
        Self {
            generators: gens,
            type_zero_end,
            type_one_end,
        }
    }

    /// Returns the sorted coset of the given representative (member of the coset) and the type of the coset.
    pub fn get_coset(&self, coset_rep: &PolyMatrix, type_ix: usize) -> Vec<PolyMatrix> {
        if type_ix == 0 {
            let subs = &self.generators[0..=self.type_zero_end];
            let mut coset: Vec<PolyMatrix> = subs.par_iter().map(|k| coset_rep * k).collect();
            coset.sort();
            coset
        } else if type_ix == 1 {
            let subs = &self.generators[self.type_zero_end + 1..=self.type_one_end];
            let mut coset: Vec<PolyMatrix> = subs.par_iter().map(|k| coset_rep * k).collect();
            coset.sort();
            coset
        } else if type_ix == 2 {
            let subs = &self.generators[self.type_one_end + 1..self.generators.len()];
            let mut coset: Vec<PolyMatrix> = subs.par_iter().map(|k| coset_rep * k).collect();
            coset.sort();
            coset
        } else {
            panic!("This type is not implemented yet.")
        }
    }

    /// Returns the sorted cosets of a given matrix of all 3 types.
    /// Returns all 3 at once as this can be more easily parallelized than
    /// getting each coset one at a time.
    pub fn get_all_cosets(
        &self,
        rep: &PolyMatrix,
    ) -> (Vec<PolyMatrix>, Vec<PolyMatrix>, Vec<PolyMatrix>) {
        let v: Vec<PolyMatrix> = self.generators.par_iter().map(|k| rep * k).collect();
        let (mut v0, mut v1, mut v2) = (
            v[..=self.type_zero_end].to_vec(),
            v[self.type_zero_end + 1..=self.type_one_end].to_vec(),
            v[self.type_one_end + 1..].to_vec(),
        );
        v0.sort();
        v1.sort();
        v2.sort();
        (v0, v1, v2)
    }

    pub fn get_canonical_rep(&self, rep: &PolyMatrix, type_ix: usize) -> PolyMatrix {
        let coset = self.get_coset(rep, type_ix);
        coset[0].clone()
    }

    pub fn get_coset_reps(&self, mat: &PolyMatrix) -> (CosetRep, CosetRep, CosetRep) {
        let (coset0, coset1, coset2) = self.get_all_cosets(mat);
        let c0 = CosetRep {
            rep: coset0[0].clone(),
            type_ix: 0,
        };
        let c1 = CosetRep {
            rep: coset1[0].clone(),
            type_ix: 1,
        };
        let c2 = CosetRep {
            rep: coset2[0].clone(),
            type_ix: 2,
        };
        (c0, c1, c2)
    }
}

mod tests {
    use std::{
        collections::{HashMap, HashSet},
        io::Write,
        path::PathBuf,
        rc::Rc,
        str::FromStr,
    };

    use crate::math::{
        ffmatrix::FFMatrix, iterative_bfs::HTypeSubgroup, polymatrix::PolyMatrix,
        polynomial::FiniteFieldPolynomial,
    };

    use super::{GroupBFS, GroupBFSNode, ParallelSubgroups};

    fn simple_quotient_and_field() -> (u32, FiniteFieldPolynomial) {
        let p = 3_u32;
        let primitive_coeffs = [(2, (1, p).into()), (1, (2, p).into()), (0, (2, p).into())];
        let q = FiniteFieldPolynomial::from(&primitive_coeffs[..]);
        (p, q)
    }

    #[test]
    fn test_group_bfs_manager() {
        let p = 3_u32;
        let primitive_coeffs = [(2, (1, p).into()), (1, (2, p).into()), (0, (2, p).into())];
        let q = FiniteFieldPolynomial::from(&primitive_coeffs[..]);
        let directory = PathBuf::from_str("/Users/matt/repos/qec/tmp2/").unwrap();
        let mut bfs_manager = GroupBFS::new(&directory, &q);
        bfs_manager.bfs((2 as usize).pow(12));
    }

    #[test]
    fn test_bfs_nodes() {
        let p = 3_u32;
        let primitive_coeffs = [(2, (1, p).into()), (1, (2, p).into()), (0, (2, p).into())];
        let q = FiniteFieldPolynomial::from(&primitive_coeffs[..]);
        let e = PolyMatrix::id(3, q.clone());
        let sg = ParallelSubgroups::new(3, &q);
        let coset = sg.get_coset(&e, 0);
        let m1 = coset[3].clone();
        let m2 = coset[4].clone();
        let bfs_node1 = GroupBFSNode {
            mat: m1.clone(),
            distance: 0,
        };
        let bfs_node2 = GroupBFSNode {
            mat: m2.clone(),
            distance: 1,
        };
        let bfs_node3 = GroupBFSNode {
            mat: m1,
            distance: 2,
        };
        let set = HashSet::from([bfs_node1, bfs_node2]);
        println!("constains 3? {:}", set.contains(&bfs_node3));
    }

    #[test]
    fn test_cereal() {
        let p = 3_u32;
        let primitive_coeffs = [(2, (1, p).into()), (1, (2, p).into()), (0, (2, p).into())];
        let q = FiniteFieldPolynomial::from(&primitive_coeffs[..]);
        let e = PolyMatrix::id(3, q.clone());
        let sg = ParallelSubgroups::new(3, &q);
        let coset = sg.get_coset(&e, 0);
        let m1 = coset[3].clone();
        let m2 = coset[4].clone();
        let (c0, c1, c2) = sg.get_coset_reps(&m2);
        let mapper = HashMap::from([(c0, 12_u32), (c1, 10), (c2, 99)]);
        let serde_out = serde_json::to_string(&mapper);
        println!("{:}", serde_out.unwrap());
    }

    #[test]
    fn test_short_walk() {
        let dir = PathBuf::from("/Users/matt/repos/qec/tmp");
        let p = 3_u32;
        let primitive_coeffs = [(2, (1, p).into()), (1, (2, p).into()), (0, (2, p).into())];
        let q = FiniteFieldPolynomial::from(&primitive_coeffs[..]);
        let mut bfs = GroupBFS::new(&dir, &q);
        bfs.bfs((2 as usize).pow(10));
        println!("graph: {:}", bfs.hg);
    }

    #[test]
    fn test_old_subgroups() {
        let p = 3_u32;
        let primitive_coeffs = [(2, (1, p).into()), (1, (2, p).into()), (0, (2, p).into())];
        let q = FiniteFieldPolynomial::from(&primitive_coeffs[..]);
        let subs = ParallelSubgroups::new(3, &q);
        let h_type_subs = HTypeSubgroup::new(&q);
        let id = PolyMatrix::id(3, q.clone());
        let out = h_type_subs.generate_left_mul(&id);
        println!("matrix mul output");
        for o in out.iter() {
            println!("{:}", o);
        }
        let out2 = h_type_subs.generate_left_mul(&out[0]);
        println!("out two!");
        for o in out2 {
            println!("{:}", o);
        }
    }
}
