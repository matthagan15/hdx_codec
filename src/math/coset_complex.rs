use core::num;
use std::{
    borrow::{Borrow, BorrowMut}, collections::{HashMap, HashSet, VecDeque}, io::{Read, Write}, os::unix::process::parent_id, path::PathBuf, str::FromStr, sync::{Arc, Mutex, RwLock}, time
};

use bitvec::ptr::read;
use mhgl::HGraph;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use serde::{Deserialize, Serialize};
use uuid::Uuid;

use super::{
    finite_field::FiniteField,
    group_ring_field::Group,
    polymatrix::{dim_three_det, PolyMatrix},
    polynomial::FiniteFieldPolynomial,
};

const SUBGROUP_FILE_EXTENSION: &str = ".subs";
const GROUP_FILE_EXTENSION: &str = ".group";
const HGRAPH_FILE_EXTENSION: &str = ".hgraph";
const NODE_TO_COSET_FILE_EXTENSION: &str = ".nodes";

pub const SUBGROUP_FILENAME: &str = "subgroups.json";
pub const GROUP_FILENAME: &str = "group.json";
pub const HGRAPH_FILENAME: &str = "hgraph.json";
const NODE_TO_COSET_FILENAME: &str = "node_to_coset.json";

fn generate_all_polys(
    field_mod: u32,
    max_degree: usize,
) -> HashMap<usize, HashSet<FiniteFieldPolynomial>> {
    if max_degree == 0 {
        let mut ret = HashSet::new();
        for a in 0..field_mod {
            ret.insert(FiniteFieldPolynomial::constant(a, field_mod));
        }
        HashMap::from([(0, ret)])
    } else {
        let smalls = generate_all_polys(field_mod, max_degree - 1);
        let mut ret = HashMap::new();
        let monomials: Vec<FiniteFieldPolynomial> = (1..field_mod)
            .into_iter()
            .map(|x| FiniteFieldPolynomial::monomial(FiniteField::new(x, field_mod), max_degree))
            .collect();
        for (deg, polys) in smalls.into_iter() {
            for poly in polys.iter() {
                for mono in monomials.iter() {
                    let with_mono: &mut HashSet<FiniteFieldPolynomial> =
                        ret.entry(max_degree).or_default();
                    with_mono.insert(mono + poly);

                    let without_mono = ret.entry(poly.degree()).or_default();
                    without_mono.insert(poly.clone());
                }
            }
        }
        ret
    }
}

fn compute_deg(dim: usize, type_ix: usize, row_ix: usize, col_ix: usize) -> usize {
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
        dist_from_diag as usize
    }
}

fn compute_subgroups(dim: usize, quotient: FiniteFieldPolynomial) -> CosetGenerators {
    let p = quotient.field_mod;
    let id = PolyMatrix::id(dim, quotient.clone());
    // let polys = generate_all_polys(p, dim - 1);
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
                // TODO: This control logic is obfuscated. I don't think
                // there is any need to drain the matrices, I think
                // you can just iterate over them.
                if deg > 0 {
                    let matrices = type_to_gens
                        .get_mut(&type_ix)
                        .expect("couldn't get matrices");
                    let mut new_matrices = HashSet::new();
                    for mut matrix in matrices.drain() {
                        new_matrices.insert(matrix.clone());

                        // TODO: This is the way claimed in the paper that
                        // should generate SL_3(R_n), but my test seems to
                        // indicate that it doesn't...
                        for c in 1..p {
                            let entry = matrix.get_mut(row_ix, col_ix);
                            *entry = FiniteFieldPolynomial::monomial((c, p).into(), deg);
                            new_matrices.insert(matrix.clone());
                        }
                        // for d in 0..=deg {

                        //     // TODO: This is wrong, the polynomial should be a monomial over the degree, not contain all polynomials of that degree..
                        //     for poly in polys.get(&d).expect("no poly?") {
                        //         let entry = matrix.get_mut(row_ix, col_ix);
                        //         *entry = poly.clone();
                        //         new_matrices.insert(matrix.clone());
                        //     }
                        // }
                    }
                    type_to_gens.insert(type_ix, new_matrices);
                }
            }
        }
    }
    CosetGenerators {
        type_to_generators: type_to_gens,
        dim,
        quotient,
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct CosetGenerators {
    pub type_to_generators: HashMap<usize, HashSet<PolyMatrix>>,
    dim: usize,
    quotient: FiniteFieldPolynomial,
}

/// Currently comptes the entire group using Breadth-First-Search
/// starting at the identity matrix over the generators provided.
fn compute_group(generators: &CosetGenerators) -> HashSet<PolyMatrix> {
    println!("Computing group with the following parameters:");
    println!("quotient: {:}", generators.quotient);
    println!("dim: {:}", generators.dim);
    for (t, set) in generators.type_to_generators.iter() {
        println!("{:}", "*".repeat(50));
        println!("Type: {:}", t);
        for s in set.iter() {
            println!("{:}", s);
        }
    }
    let e = PolyMatrix::id(generators.dim, generators.quotient.clone());


    // TODO: Currently a matrix is being stored twice while it is in
    // the frontier as we also put it in visited. Instead just keep track
    // of an extra bit if the matrix is visited or not.
    // Also, do not store completed in RAM. come up with some way of storing
    // them on disk in the meanwhile.
    let mut completed = HashSet::new();
    let mut frontier = VecDeque::from([e.clone()]);
    let mut visited = HashSet::from([e.clone()]);
    let gens = generators.type_to_generators.clone();
    
    if gens.is_empty() {
        completed.insert(e);
        return completed;
    }
   

    let mut counter = 0;
    while frontier.len() > 0 {
        counter += 1;
        if (counter % 100) == 0 {
            println!("{:}", ".".repeat(50));
            println!(
                "frontier length: {:} , {:.4}% of visitied.",
                frontier.len(),
                frontier.len() as f64 / visited.len() as f64
            );
            println!("completed length: {:}", completed.len());
            println!("visited length: {:}", visited.len());
        }
        let x = frontier.pop_front().expect("no frontier?");
        for (j, gen_list) in gens.iter() {
            for g in gen_list {
                let new = g * &x;
                if visited.contains(&new) == false {
                    visited.insert(new.clone());
                    frontier.push_back(new);
                }
            }
        }
        completed.insert(x);
    }
    completed
}

/// Computes the coset of a provided `start` group element and
fn compute_coset(start: &PolyMatrix, coset_type: usize, coset_gens: &CosetGenerators) -> Coset {
    let subgroup = coset_gens
        .type_to_generators
        .get(&coset_type)
        .expect("Could not find coset of given type.");
    let mut completed = HashSet::new();
    for elem in subgroup.iter() {
        completed.insert(start * elem);
    }
    let mut ret: Vec<PolyMatrix> = completed.into_iter().collect();
    // for m in ret.iter_mut() {
    //     m.clean()
    // }
    ret.sort();
    Coset {
        type_ix: coset_type,
        set: ret,
    }
}

#[derive(Hash, Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
struct Coset {
    type_ix: usize,
    set: Vec<PolyMatrix>,
}

fn compute_vertices(
    group: &HashSet<PolyMatrix>,
    subgroups: &CosetGenerators,
    hgraph: &mut HGraph,
) -> HashMap<u32, Coset> {
    let mut node_to_coset = HashMap::new();
    let mut cosets = HashSet::new();
    let mut coset_size = None;
    let mut counter = 0;
    let tot = group.len();
    for g in group.iter() {
        let percent_done = counter as f64 / tot as f64;
        counter += 1;
        if counter % 50 == 0 {
            println!("{:.4}% done", percent_done);
        }
        for coset_type in subgroups.type_to_generators.keys() {
            let coset = compute_coset(g, *coset_type, subgroups);
            if cosets.contains(&coset) == false {
                let nodes = hgraph.add_nodes(1);
                let node = nodes[0];

                if coset_size.is_none() {
                    coset_size = Some(coset.set.len());
                }

                cosets.insert(coset.clone());
                node_to_coset.insert(node, coset);
            }
        }
    }
    println!("{:}", "*".repeat(75));
    println!("size of group: {:}", group.len());
    println!("number of cosets: {:}", cosets.len());
    println!("coset size: {:?}", coset_size);
    node_to_coset
}

fn compute_edges(nodes_to_coset: &HashMap<u32, Coset>, hgraph: &mut HGraph) {
    let edge_set: HashSet<(u32, u32)> = HashSet::new();
    let mut locker = Arc::new(RwLock::new(edge_set));
    for (n1, coset1) in nodes_to_coset.iter() {
        nodes_to_coset.par_iter().for_each(|(n2, coset2)| {
            if n1 == n2 {
                return;
            }
            let reader = locker.read().expect("cannot read lock.");
            let contains_edge = reader.contains(&(*n1.min(n2), *n1.max(n2)));
            drop(reader);
            if contains_edge == false {
                let set1: HashSet<PolyMatrix> = coset1.set.clone().into_iter().collect();
                for m2 in coset2.set.iter() {
                    if set1.contains(m2) {
                        let mut writer = locker
                            .write()
                            .expect("could not get write access to edge creator.");
                        writer.insert((*n1.min(n2), *n1.max(n2)));
                        break;
                    }
                }
            }
        });
    }
    let reader = locker.read().expect("Could not read locker.");
    for (n1, n2) in reader.iter() {
        hgraph.create_edge(&[*n1, *n2]);
    }
}

fn compute_triangles(nodes_to_coset: &HashMap<u32, Coset>, hgraph: &mut HGraph) {
    let edge_set: HashSet<Vec<u32>> = HashSet::new();
    let locker = Arc::new(Mutex::new(edge_set));
    // nodes_to_coset.par_iter().for_each(|(n1, coset1)| {
    let mut counter = 0;
    for (n1, coset1) in nodes_to_coset.iter() {
        if counter % 100 == 0 {
            let percent_done = counter as f64 / nodes_to_coset.len() as f64;
            println!("{:.4} % done.", percent_done);
        }
        counter += 1;
        for (n2, coset2) in nodes_to_coset.iter() {
            if n1 == n2 {
                continue;
            }
            if hgraph.query_edge(&[*n1, *n2]) == false {
                continue;
            }
            for (n3, coset3) in nodes_to_coset.iter() {
                if n1 == n3 || n2 == n3 {
                    continue;
                }
                if hgraph.query_edge(&[*n1, *n3]) == false
                    || hgraph.query_edge(&[*n2, *n3]) == false
                {
                    continue;
                }
                let set1: HashSet<PolyMatrix> = coset1.set.clone().into_iter().collect();
                let set2: HashSet<PolyMatrix> = coset2.set.clone().into_iter().collect();
                let intersection: HashSet<PolyMatrix> = set1.intersection(&set2).cloned().collect();
                for m3 in coset3.set.iter() {
                    if intersection.contains(m3) {
                        let mut node_vec = vec![*n1, *n2, *n3];
                        node_vec.sort();
                        let mut write_lock = locker.lock().expect("couldn't lock mutex.");
                        write_lock.insert(node_vec);
                    }
                }
            }
        }
        // });
    }
    let triangle_set = locker.lock().expect("could not lock triangle set.");
    for triangle in triangle_set.iter() {
        hgraph.create_edge(&triangle[..]);
    }
}


#[derive(Debug, Clone)]
pub struct CosetComplex {
    file_base: String,
    dim: usize,
    quotient: FiniteFieldPolynomial,
    group: Option<HashSet<PolyMatrix>>,
    subgroups: Option<CosetGenerators>,
    hgraph: HGraph,
    node_to_coset: Option<HashMap<u32, Coset>>,
}

impl CosetComplex {
    pub fn new(file_base: String, dim: usize, quotient: &FiniteFieldPolynomial) -> Self {
        CosetComplex {
            file_base: file_base.trim().to_string(),
            dim,
            quotient: quotient.clone(),
            group: None,
            subgroups: None,
            hgraph: HGraph::new(),
            node_to_coset: None,
        }
    }

    pub fn print_hgraph(&self) {
        if self.node_to_coset.is_none() {
            println!("Need to compute vertices first, graph is empty.")
        } else {
            println!("{:}", self.hgraph);
        }
    }

    pub fn hgraph_ref(&self) -> &HGraph {
        &self.hgraph
    }

    pub fn print_subgroups(&self) {
        if let Some(coset_gens) = &self.subgroups {
            for (type_ix, gens) in coset_gens.type_to_generators.iter() {
                println!("{:}", "-".repeat(50));
                println!("Subgroup of type {:}", type_ix);
                for gen in gens.iter() {
                    println!("{:}", gen);
                }
            }
            println!("{:}", "-".repeat(50));
        }
    }

    pub fn get_triangles_containing_edge(&self, edge: (u32, u32)) -> Vec<Uuid> {
        self.hgraph.get_containing_edges(&[edge.0, edge.1])
    }

    pub fn load_subgroups_from_disk(&mut self) {
        let mut subgroup_file_path = PathBuf::from_str(&self.file_base).expect("No pathbuf?");
        subgroup_file_path.push(SUBGROUP_FILENAME);
        println!("Attempting to load subgroups from: {:}", subgroup_file_path.display());
        if let Ok(mut subgroup_file) =
            std::fs::File::open(&subgroup_file_path) {
                let mut subgroup_file_string = String::new();
            let read_ok = subgroup_file
                .read_to_string(&mut subgroup_file_string);
            if read_ok.is_err() {
                println!("Could not read subgroup file");
                return;
            }
            let subgroup_serialization = serde_json::from_str(&subgroup_file_string);
            if subgroup_serialization.is_ok() {
                self.subgroups = Some(subgroup_serialization.unwrap());
            } else {
                println!("Could not deserialize subgroups.");
            }
        }
    }

    pub fn load_groups_from_disk(&mut self) {
        let mut group_file_path = PathBuf::from_str(&self.file_base).expect("No pathbuf?");
        group_file_path.push(GROUP_FILENAME);
        if let Ok(mut group_file) = std::fs::File::open(&group_file_path) {
            let mut group_file_string = String::new();
            let read_ok = group_file
                .read_to_string(&mut group_file_string);
            if read_ok.is_err() {
                println!("could not read group file.");
                return;
            }
            let groups = serde_json::from_str(&group_file_string);
            if groups.is_ok() {
                self.group = Some(groups.unwrap());
            } else {
                println!("Could not deserialize groups.");
            }
        }
    }

    pub fn load_hgraph_from_disk(&mut self) {
        let mut hgraph_file_path = PathBuf::from_str(&self.file_base).expect("No pathbuf?");
        hgraph_file_path.push(HGRAPH_FILENAME);
        if let Ok(mut hgraph_file) =
            std::fs::File::open(&hgraph_file_path) {
                let mut hgraph_file_string = String::new();
            let read_ok = hgraph_file
                .read_to_string(&mut hgraph_file_string);
            if read_ok.is_err() {
                println!("Could not read hgraph file.");
                return;
            }
            let hgraph = serde_json::from_str(&hgraph_file_string);
            if hgraph.is_ok() {
                self.hgraph = hgraph.unwrap();
            } else {
                println!("Could not deserialie hgraph.");
            }
            }
    }

    pub fn load_node_to_coset_from_disk(&mut self) {
        let mut n_to_c_file_path = PathBuf::from_str(&self.file_base).expect("No pathbuf?");
        n_to_c_file_path.push(NODE_TO_COSET_FILENAME);
        if let Ok(mut n_to_c_file) = std::fs::File::open(&n_to_c_file_path) {
            let mut n_to_c_string = String::new();
            let read_ok = n_to_c_file
                .read_to_string(&mut n_to_c_string);
            if read_ok.is_err() {
                println!("could not read node_to_coset file.");
                return;
            }
            let node_to_coset = serde_json::from_str(&n_to_c_string);
            if node_to_coset.is_ok() {
                self.node_to_coset = Some(node_to_coset.unwrap());
            } else {
                println!("Could not deserialize node_to_coset map.");
            }
        }
    }

    pub fn load_from_disk(&mut self) {
        println!("Loading subgroups from disk.");
        self.load_subgroups_from_disk();
        self.load_groups_from_disk();
        self.load_hgraph_from_disk();
        self.load_node_to_coset_from_disk();   
    }

    pub fn generate_subgroups(&mut self) {
        if self.subgroups.is_none() {
            let gens = compute_subgroups(self.dim, self.quotient.clone());
            self.subgroups = Some(gens);
        }
    }

    pub fn generate_group(&mut self) {
        if self.subgroups.is_none() {
            self.generate_subgroups();
        }
        if self.group.is_none() {
            if let Some(subs) = &self.subgroups {
                let g = compute_group(subs);
                self.group = Some(g);
            }
        }
    }

    pub fn compute_vertices(&mut self) {
        if self.subgroups.is_none() {
            self.generate_group()
        }
        if self.group.is_none() {
            self.generate_group();
        }
        let n_to_c = compute_vertices(
            self.group.as_ref().unwrap(),
            self.subgroups.as_ref().unwrap(),
            &mut self.hgraph,
        );
        self.node_to_coset = Some(n_to_c);
    }

    pub fn compute_edges(&mut self) {
        if self.subgroups.is_none() {
            self.generate_subgroups();
        }
        if self.group.is_none() {
            self.generate_group();
        }
        if self.node_to_coset.is_none() {
            self.compute_vertices();
        }
        compute_edges(self.node_to_coset.as_ref().unwrap(), &mut self.hgraph);
    }

    pub fn compute_triangles(&mut self) {
        if self.subgroups.is_none() {
            println!("Did not find subgroups, computing now.");
            self.generate_subgroups();
            println!("Saving to disk.");
            self.save_subgroups_to_disk();
        }
        if self.group.is_none() {
            self.generate_group();
            self.save_group_to_disk();
        }
        if self.node_to_coset.is_none() {
            self.compute_vertices();
            self.save_node_to_coset_to_disk();
        }
        if self.hgraph.edges_of_size(2).len() == 0 {
            self.compute_edges();
        }
        if self.hgraph.edges_of_size(3).len() == 0 {
            compute_triangles(self.node_to_coset.as_ref().unwrap(), &mut self.hgraph);
        }
    }

    /// Computes the link degrees of each vertex and edge.
    pub fn check_degrees(&self) {
        if self.node_to_coset.is_none() {
            println!("No nodes.");
            return;
        }
        let nodes: Vec<&u32> = self.node_to_coset.as_ref().unwrap().keys().collect();
        let mut node_to_node_stats = HashMap::new();
        let mut node_to_edges_stats = HashMap::new();
        let edges = self.hgraph.edges_of_size(2);
        let mut edge_stats = HashMap::new();
        for node in nodes {
            let link = self.hgraph.link_as_vec(&[*node]);
            let mut num_nodes = 0;
            let mut num_edges = 0;
            for (set, weight) in link.into_iter() {
                if set.len() == 1 {
                    num_nodes += 1;
                    continue;
                }
                if set.len() == 2 {
                    num_edges += 1;
                }
            }
            let e: &mut usize = node_to_node_stats.entry(num_nodes).or_default();
            *e += 1;
            let ee: &mut usize = node_to_edges_stats.entry(num_edges).or_default();
            *ee += 1;
        }
        for edge in edges {
            if let Some(edge_set) = self.hgraph.query_edge_id(&edge) {
                let mut num_nodes = 0;
                for (set, weight) in self.hgraph.link_as_vec(&edge_set[..]) {
                    if set.len() == 1 {
                        num_nodes += 1;
                    }
                }
                let e: &mut usize = edge_stats.entry(num_nodes).or_default();
                *e += 1;
            }
        }
        println!("Node statistics. The node -> node degree stats are:");
        println!("{:?}", node_to_node_stats);
        println!("node -> edge degree stats are:");
        println!("{:?}", node_to_edges_stats);
        println!("And the edge -> node degree stats are:");
        println!("{:?}", edge_stats);
    }

    pub fn save_subgroups_to_disk(&self) {
        let mut subgroup_file_path = PathBuf::from_str(&self.file_base).expect("No pathbuf?");
        subgroup_file_path.push(SUBGROUP_FILENAME);
        let mut subgroup_file =
            std::fs::File::create(subgroup_file_path).expect("could not open file for writing.");
        if let Some(subs) = &self.subgroups {
            let subs_string = serde_json::to_string(subs).expect("Could not serialize subgroups.");
            subgroup_file
                .write(subs_string.as_bytes())
                .expect("Failed to write subgroups");
        } else {
            println!("Did not have subgroup to save to disk.");
        }
    }

    pub fn save_group_to_disk(&self) {
        let mut group_file_path = PathBuf::from_str(&self.file_base).expect("No pathbuf?");
        group_file_path.push(GROUP_FILENAME);
        let mut group_file =
            std::fs::File::create(group_file_path).expect("could not create group file.");
        if let Some(group) = &self.group {
            let group_string = serde_json::to_string(group).expect("could not serialize group");
            group_file
                .write(group_string.as_bytes())
                .expect("Failed to write groups.");
        } else {
            println!("Tried to write non-existent group to disk.");
        }
    }

    pub fn save_hgraph_to_disk(&self) {
        let mut hgraph_file_path = PathBuf::from_str(&self.file_base).expect("No pathbuf?");
        hgraph_file_path.push(HGRAPH_FILENAME);
        let mut hgraph_file =
            std::fs::File::create(hgraph_file_path).expect("Could not create hgraph file.");
        let hgraph_string =
            serde_json::to_string(&self.hgraph).expect("Could not serialize hgraph.");
        hgraph_file
            .write(hgraph_string.as_bytes())
            .expect("Could not write hgraph to file.");
    }

    pub fn save_node_to_coset_to_disk(&self) {
        if let Some(n_to_c) = &self.node_to_coset {
            let mut n_to_c_file_path = PathBuf::from_str(&self.file_base).expect("No pathbuf?");
            n_to_c_file_path.push(NODE_TO_COSET_FILENAME);
            let mut n_to_c_file = std::fs::File::create(n_to_c_file_path)
                .expect("Could not create node_to_coset file.");
            let n_to_c_string =
                serde_json::to_string(n_to_c).expect("Could not serialize node_to_cosest.");
            n_to_c_file
                .write(n_to_c_string.as_bytes())
                .expect("Could not write node_to_coset to file.");
        } else {
            println!("Did not have node_to_coset complex to save to disk.");
        }
    }

    pub fn save_to_disk(&self) {
        self.save_subgroups_to_disk();
        self.save_group_to_disk();
        self.save_hgraph_to_disk();
        self.save_node_to_coset_to_disk();
    }

    pub fn compute_hgraph_at_all_costs(&mut self) -> HGraph {
        // first attempt to load previous session from disk
        self.load_from_disk();
        self.compute_triangles();
        self.save_hgraph_to_disk();
        self.hgraph.clone()
    }
}

struct GroupIterator {
    indices: Vec<usize>,
    dimension: usize,
    polys: Vec<FiniteFieldPolynomial>,
    quotient: FiniteFieldPolynomial,
}

impl GroupIterator {
    /// returns true if indices can be incremented again, false if it cannot be incremented.
    fn increment_indices(&mut self) -> bool {
        let mut carry_increment = true;
        let mut ix = 0;
        while carry_increment {
            let mut tmp = self.indices[ix];
            tmp += 1;
            tmp %= self.polys.len();
            self.indices[ix] = tmp;
            if tmp != 0 {
                carry_increment = false;
            } else {
                if ix == self.indices.len() - 1 {
                    return false;
                }
                carry_increment = true;
                ix += 1;
            }
        }
        true
    }

    fn construct_matrix(&self) -> PolyMatrix {
        let entries: Vec<FiniteFieldPolynomial> = self
            .indices
            .clone()
            .into_iter()
            .map(|ix| self.polys[ix].clone())
            .collect();
        let p = entries[0].field_mod;
        PolyMatrix {
            entries,
            n_rows: self.dimension,
            n_cols: self.dimension,
            field_mod: p,
            quotient: self.quotient.clone(),
        }
    }

    fn generate_sl3(&mut self) -> HashSet<PolyMatrix> {
        let upper_bound = self.polys.len().pow(self.indices.len() as u32);
        let mut counter = 0;
        let mut can_increment = true;
        let mut sl3 = HashSet::new();
        while can_increment {
            if counter % 10000 == 0 {
                let percent_done = counter as f64 / upper_bound as f64;
                println!("{:.4} % done", percent_done);
            }
            counter += 1;
            let m = self.construct_matrix();
            if dim_three_det(&m).is_one() {
                sl3.insert(m);
            }
            can_increment = self.increment_indices();
        }
        sl3
    }
}

mod tests {
    use std::collections::HashSet;

    use crate::math::{
        coset_complex::{compute_coset, compute_triangles},
        finite_field::FiniteField,
        polymatrix::PolyMatrix,
        polynomial::FiniteFieldPolynomial,
    };

    use super::{
        compute_edges, compute_group, compute_subgroups, compute_vertices, generate_all_polys,
        CosetComplex, CosetGenerators, GroupIterator,
    };

    use mhgl::HGraph;

    #[test]
    fn test_sl3_generation() {
        let (gens, bfs_sol) = simplest_group();
        for (k, v) in gens.type_to_generators.iter() {
            println!("Type: {:}", k);
            for m in v.iter() {
                println!("{:}", m);
            }
        }
        let p = 3_u32;
        let primitive_coeffs = [(2, (1, p).into()), (1, (2, p).into()), (0, (2, p).into())];
        let primitive_poly = FiniteFieldPolynomial::from(&primitive_coeffs[..]);
        let deg_to_polys = generate_all_polys(p, 1);
        let mut all_polys: Vec<FiniteFieldPolynomial> = Vec::new();
        for (_, set) in deg_to_polys.into_iter() {
            let mut vecd = set.into_iter().collect();
            all_polys.append(&mut vecd);
        }

        let mut gi = GroupIterator {
            indices: vec![0, 0, 0, 0, 0, 0, 0, 0, 0],
            dimension: 3,
            polys: all_polys,
            quotient: primitive_poly,
        };
        println!("Generating sl3");
        let exhaustive_sol = gi.generate_sl3();
        println!("Exhaustive # sols = {:}", exhaustive_sol.len());
        println!("Bfs # sols = {:}", bfs_sol.len());
        let mut exhaustive_in_bfs = true;
        let mut bfs_in_exhaustive = true;
        for m in exhaustive_sol.iter() {
            if bfs_sol.contains(m) == false {
                exhaustive_in_bfs = false;
                break;
            }
        }
        for m in bfs_sol.iter() {
            if exhaustive_sol.contains(m) == false {
                bfs_in_exhaustive = false;
                break;
            }
        }
        println!(
            "two sets equal: {:}",
            exhaustive_in_bfs && bfs_in_exhaustive
        );
    }

    fn simplest_group() -> (CosetGenerators, HashSet<PolyMatrix>) {
        let p = 3_u32;
        let primitive_coeffs = [(2, (1, p).into()), (1, (2, p).into()), (0, (2, p).into())];
        let primitive_poly = FiniteFieldPolynomial::from(&primitive_coeffs[..]);
        let dim = 3;
        let gens = compute_subgroups(dim, primitive_poly);
        println!("subgroups computed.");
        let g = compute_group(&gens);
        (gens, g)
    }

    #[test]
    fn test_serialize_groups() {
        let p = 2_u32;
        let primitive_coeffs = [(2, (1, p).into()), (1, (2, p).into()), (0, (2, p).into())];
        let q = FiniteFieldPolynomial::from(&primitive_coeffs[..]);
        let mut gm = CosetComplex::new(
            String::from("/Users/matt/repos/qec/data/groups/p_2_dim_3_deg_2.group"),
            3,
            &q,
        );
        gm.generate_group();
        gm.save_to_disk();
    }

    #[test]
    fn test_deserialize_groups() {
        let p = 2_u32;
        let primitive_coeffs = [(2, (1, p).into()), (1, (2, p).into()), (0, (2, p).into())];
        let q = FiniteFieldPolynomial::from(&primitive_coeffs[..]);
        let mut gm = CosetComplex::new(
            String::from("/Users/matt/repos/qec/data/groups/p_2_dim_3_deg_2.group"),
            3,
            &q,
        );
        gm.load_from_disk();
    }

    fn get_nontrivial_group_manager() -> CosetComplex {
        let p = 2_u32;
        let primitive_coeffs = [(2, (1, p).into()), (1, (2, p).into()), (0, (2, p).into())];
        let q = FiniteFieldPolynomial::from(&primitive_coeffs[..]);
        let mut gm = CosetComplex::new(
            String::from("/Users/matt/repos/qec/data/groups/p_2_dim_3_deg_2.group"),
            3,
            &q,
        );
        gm.load_from_disk();
        gm
    }

    #[test]
    fn test_compute_group() {
        let (gens, g) = simplest_group();
        println!(
            "size of subgroups: {:}",
            gens.type_to_generators.get(&0).unwrap().len()
        );
    }

    #[test]
    fn test_vertex_creation() {
        let gm = get_nontrivial_group_manager();
        let mut hg = HGraph::new();
        let group = gm.group.unwrap();
        let subgroups = gm.subgroups.unwrap();
        let node_to_coset = compute_vertices(&group, &subgroups, &mut hg);
        println!("hg:\n{:}", hg);
    }

    #[test]
    fn test_compute_edges() {
        let gm = get_nontrivial_group_manager();
        let mut hg = HGraph::new();
        let group = gm.group.unwrap();
        let subgroups = gm.subgroups.unwrap();
        let nodes_to_coset = compute_vertices(&group, &subgroups, &mut hg);
        compute_edges(&nodes_to_coset, &mut hg);
        println!("hg:\n{:}", hg);
    }

    #[test]
    fn test_compute_triangles() {
        let gm = get_nontrivial_group_manager();
        let mut hg = HGraph::new();
        let group = gm.group.unwrap();
        let subgroups = gm.subgroups.unwrap();
        let nodes_to_coset = compute_vertices(&group, &subgroups, &mut hg);
        compute_triangles(&nodes_to_coset, &mut hg);
        let edges = hg.edges_of_size(2);
        let triangles = hg.edges_of_size(3);
        println!("number edges: {:}", edges.len());
        println!("number triangles: {:}", triangles.len());
    }

    #[test]
    fn test_coset() {
        let (gens, g) = simplest_group();
        println!("len of g: {:}", g.len());
        for x in g.iter() {
            println!("{:}", "#".repeat(75));
            println!("x: {:}", x);
            let gen_type_0 = gens.type_to_generators.get(&0).unwrap().clone();
            let coset = compute_coset(&x, 0, &gens);
            println!("{:}", ".".repeat(75));
            println!("coset.");
            for h in coset.set.iter() {
                println!("{:}", h);
            }
        }
    }

    #[test]
    fn test_alternative_generators() {
        let p = 3_u32;
        let primitive_coeffs = [(2, (1, p).into()), (1, (2, p).into()), (0, (2, p).into())];
        let primitive_poly = FiniteFieldPolynomial::from(&primitive_coeffs[..]);
        let dim = 3;
        let out = compute_subgroups(dim, primitive_poly);
        for (t, gens) in out.type_to_generators {
            println!("{:}", "*".repeat(75));
            println!("type {:}", t);
            println!("len of gens: {:}", gens.len());
            println!("{:}", "-".repeat(75));
        }
    }

    #[test]
    fn test_all_polys() {
        let out = generate_all_polys(3, 3);
        for (deg, polys) in out {
            println!("degree {:} polynomials", deg);
            println!("{:}", "-".repeat(75));
            for poly in polys {
                println!("{:}", poly);
            }
            println!("{:}", "#".repeat(75));
        }
    }

    #[test]
    fn test_basic_group() {
        let (gens, g) = simplest_group();
        let q = gens.quotient.clone();
        let all_polys = generate_all_polys(q.field_mod, q.degree());
    }
}
