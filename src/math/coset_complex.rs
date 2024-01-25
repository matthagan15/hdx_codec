
use core::num;
use std::{borrow::{Borrow, BorrowMut}, collections::{HashSet, HashMap, VecDeque}, io::{Read, Write}, os::unix::process::parent_id, sync::{Arc, RwLock}, time};

use mhgl::HGraph;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use serde::{Deserialize, Serialize};

use super::{polynomial::FiniteFieldPolynomial, matrix::PolyMatrix, finite_field::FiniteField};

const SUBGROUP_FILE_EXTENSION: &str = ".subs";
const GROUP_FILE_EXTENSION: &str = ".group";
const HGRAPH_FILE_EXTENSION: &str = ".hgraph";
const NODE_TO_COSET_FILE_EXTENSION: &str = ".nodes";

/// Computes the set of all matrices reachable from start by multiplying by 
/// any matrices from the set generators.
fn matrix_bfs(start: &PolyMatrix, generators: &Vec<PolyMatrix>, act_on_left: bool) -> HashSet<PolyMatrix> {
    let mut completed = HashSet::new();
    let mut frontier = VecDeque::from([start.clone()]);
    let mut visited = HashSet::from([start.clone()]);
    while frontier.len() > 0 {
        let x = frontier.pop_front().expect("no frontier?");
        for g in generators {
            let new = if act_on_left {
                g * &x
            } else {
                &x * g
            };
            if visited.contains(&new) == false {
                visited.insert(new.clone());
                frontier.push_back(new);
            }
        }
        completed.insert(x);
    }
    completed.insert(start.clone());
    completed
}

fn generate_all_polys(field_mod: u32, max_degree: usize) -> HashMap<usize, HashSet<FiniteFieldPolynomial>> {
    if max_degree == 0 {
        let mut ret = HashSet::new();
        for a in 0..field_mod {
            ret.insert(FiniteFieldPolynomial::constant(a, field_mod));
        }
        HashMap::from([(0, ret)])
    } else {
        let smalls = generate_all_polys(field_mod, max_degree - 1);
        let mut ret = HashMap::new();
        let monomials: Vec<FiniteFieldPolynomial> = (1..field_mod).into_iter().map(|x| FiniteFieldPolynomial::monomial(FiniteField::new(x, field_mod), max_degree)).collect();
        for (deg, polys) in smalls.into_iter() {
            for poly in polys.iter() {
                for mono in monomials.iter() {
                    let with_mono: &mut HashSet<FiniteFieldPolynomial> = ret.entry(max_degree).or_default();
                    with_mono.insert(mono + poly);

                    let without_mono = ret.entry(poly.degree()).or_default();
                    without_mono.insert(poly.clone());
                }
            }
        }
        ret
    }
}



fn compute_subgroups(dim: usize, quotient: FiniteFieldPolynomial) -> CosetGenerators {
    let p = quotient.field_mod;
    let compute_deg = |ix: usize, j: usize, k: usize| {
        let mut how_many_over: HashMap<usize, usize> = HashMap::new();
        let mut row = ix;
        let mut hoppable_cols_left = (dim - 1) as i32;
        for _ in 0..dim {
            let e = how_many_over.entry(row).or_default();
            *e = hoppable_cols_left as usize;
            hoppable_cols_left -= 1;
            row += 1;
            row %= dim;
        }
        let hoppable_cols = how_many_over[&j];
        let mut dist_from_diag = k as i32 - j as i32;
        while dist_from_diag < 0 {
            dist_from_diag += dim as i32;
        }
        dist_from_diag %= dim as i32;
        if dist_from_diag > (hoppable_cols as i32) {
            0
        } else {
            dist_from_diag as usize
        }
    };
    let id = PolyMatrix::id(dim, quotient.clone());
    let polys = generate_all_polys(p, dim - 1);
    let mut ret = HashMap::new();
    for i in 0..dim {
        let ret_type_i: &mut HashSet<PolyMatrix> = ret.entry(i).or_default();
        ret_type_i.insert(id.clone());
        for j in 0..dim {
            for k in 0..dim {
                if j == k {
                    continue;
                }
                let deg = compute_deg(i, j, k);
                if deg > 0 {
                    let matrices = ret.get_mut(&i).expect("couldn't get matrices");
                    let mut new_matrices = HashSet::new();
                    for mut matrix in matrices.drain() {
                        new_matrices.insert(matrix.clone());
                        for d in 0..=deg {
                            for poly in polys.get(&d).expect("no poly?") {
                                let entry = matrix.get_mut(j, k);
                                *entry = poly.clone();
                                new_matrices.insert(matrix.clone());
                            }
                        }
                    }
                    ret.insert(i, new_matrices);
                }
            }
        }
    }
    CosetGenerators {
        type_to_generators: ret,
        dim,
        quotient,
    }
    
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct CosetGenerators {
    type_to_generators: HashMap<usize, HashSet<PolyMatrix>>,
    dim: usize,
    quotient: FiniteFieldPolynomial,
}


/// Currently comptes the entire group using Breadth-First-Search 
/// starting at the identity matrix over the generators provided.
fn compute_group(generators: &CosetGenerators) -> HashSet<PolyMatrix> {
    let mut completed = HashSet::new();
    let gens = generators.type_to_generators.clone();
    let e = PolyMatrix::id(generators.dim, generators.quotient.clone());
    if gens.is_empty() {
        completed.insert(e);
        return completed;
    }
    let mut frontier = VecDeque::from([e.clone()]);
    let mut visited = HashSet::from([e.clone()]);

    let mut counter = 0;
    while frontier.len() > 0 {
        counter += 1;
        if (counter % 100) == 0 {
            println!("{:}", ".".repeat(50));
            println!("frontier length: {:}", frontier.len());
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
    let subgroup = coset_gens.type_to_generators.get(&coset_type).expect("Could not find coset of given type.");
    let mut completed = HashSet::new();
    for elem in subgroup.iter() {
        completed.insert(start * elem);
    }
    let mut ret: Vec<PolyMatrix> = completed.into_iter().collect();
    // for m in ret.iter_mut() {
    //     m.clean()
    // }
    ret.sort();
    Coset { type_ix: coset_type, set: ret }
}

#[derive(Hash, Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
struct Coset {
    type_ix: usize,
    set: Vec<PolyMatrix>,
}

fn compute_vertices(group: &HashSet<PolyMatrix>, subgroups: &CosetGenerators, hgraph: &mut HGraph) -> HashMap<u32, Coset> {
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
    // println!("printing cosets.");
    // for c in cosets.iter() {
    //     println!("type: {:}", c.type_ix);
    //     for m in c.set.iter() {
    //         println!("{:}", m);
    //     }
    // }
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
                        let mut writer = locker.write().expect("could not get write access to edge creator.");
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
    let mut num_edges = 0;
    let mut num_triangles = 0;
    let mut counter1 = 0;
    let mut counter2 = 0;
    for (n1, coset1) in nodes_to_coset.iter() {
        let percent_done = counter1 as f64 / nodes_to_coset.len() as f64;
        counter1 += 1;
        println!("{:.4}% done loop 1.", percent_done);

        for (n2, coset2) in nodes_to_coset.iter() {
            if n1 == n2  {
                continue;
            }
            let percent_done = counter2 as f64 / nodes_to_coset.len() as f64;
            counter2 += 1;
            println!("{:.4}% done loop 2.", percent_done);

            if hgraph.query_edge(&[*n1, *n2]) == false {
                let set1: HashSet<PolyMatrix> = coset1.set.clone().into_iter().collect();
                for m2 in coset2.set.iter() {
                    if set1.contains(m2) {
                        hgraph.create_edge(&[*n1, *n2]);
                        num_edges += 1;
                        break;
                    }
                }
            }

            for (n3, coset3) in nodes_to_coset.iter() {
                if n1 == n3 || n2 == n3 {
                    continue;
                }

                if hgraph.query_edge(&[*n1, *n3]) == false {
                    let set1: HashSet<PolyMatrix> = coset1.set.clone().into_iter().collect();
                    for m3 in coset3.set.iter() {
                        if set1.contains(m3) {
                            hgraph.create_edge(&[*n1, *n3]);
                            num_edges += 1;
                            break;
                        }
                    }
                }
                if hgraph.query_edge(&[*n2, *n3]) == false {
                    let set2: HashSet<PolyMatrix> = coset2.set.clone().into_iter().collect();
                    for m3 in coset3.set.iter() {
                        if set2.contains(m3) {
                            hgraph.create_edge(&[*n2, *n3]);
                            num_edges += 1;
                            break;
                        }
                    }
                }
                if hgraph.query_edge(&[*n1, *n2, *n3]) == false {
                    let set1: HashSet<PolyMatrix> = coset1.set.clone().into_iter().collect();
                    let set2: HashSet<PolyMatrix> = coset2.set.clone().into_iter().collect();
                    for m3 in coset3.set.iter() {
                        if set1.contains(m3) && set2.contains(m3) {
                            hgraph.create_edge(&[*n1, *n2, *n3]);
                            num_triangles +=1;
                            break;
                        }
                    }
                }
            }
        }
    }
    dbg!(num_edges);
    dbg!(num_triangles);
}

pub struct DiskManager {
    file_base: String,
    dim: usize,
    quotient: FiniteFieldPolynomial,
    field_mod: u32,
    group: Option<HashSet<PolyMatrix>>,
    subgroups: Option<CosetGenerators>,
    hgraph: HGraph,
    node_to_coset: Option<HashMap<u32, Coset>>,
}

impl DiskManager {
    pub fn new(file_base: String, dim: usize, quotient: &FiniteFieldPolynomial) -> Self {
        DiskManager {
            file_base,
            dim,
            quotient: quotient.clone(),
            field_mod: quotient.field_mod,
            group: None,
            subgroups: None,
            hgraph: HGraph::new(),
            node_to_coset: None,
        }
    }
    pub fn load_from_disk(&mut self) {
        let mut subgroup_file_path = self.file_base.clone();
        subgroup_file_path.push_str(SUBGROUP_FILE_EXTENSION);
        let mut subgroup_file = std::fs::File::open(&subgroup_file_path).expect("Could not load file.");
        let mut subgroup_file_string = String::new();
        subgroup_file.read_to_string(&mut subgroup_file_string).expect("Could not read file");
        let subgroup_serialization = serde_json::from_str(&subgroup_file_string);
        if subgroup_serialization.is_ok() {
            self.subgroups = Some(subgroup_serialization.unwrap());
        } else {
            println!("Could not deserialize subgroups.");
        }

        let mut group_file_path = self.file_base.clone();
        group_file_path.push_str(GROUP_FILE_EXTENSION);
        let mut group_file = std::fs::File::open(&group_file_path).expect("Could not load file.");
        let mut group_file_string = String::new();
        group_file.read_to_string(&mut group_file_string).expect("Could not read file");
        let groups = serde_json::from_str(&group_file_string);
        if groups.is_ok() {
            self.group = Some(groups.unwrap());
        } else {
            println!("Could not deserialize groups.");
        }

        let mut hgraph_file_path = self.file_base.clone();
        hgraph_file_path.push_str(HGRAPH_FILE_EXTENSION);
        let mut hgraph_file = std::fs::File::open(&hgraph_file_path).expect("Could not open hgraph file.");
        let mut hgraph_file_string = String::new();
        hgraph_file.read_to_string(&mut hgraph_file_string).expect("Could not read hgraph file.");
        let hgraph = serde_json::from_str(&hgraph_file_string);
        if hgraph.is_ok() {
            self.hgraph = hgraph.unwrap();
        } else {
            println!("Could not deserialie hgraph.");
        }

        let mut n_to_c_file_path = self.file_base.clone();
        n_to_c_file_path.push_str(NODE_TO_COSET_FILE_EXTENSION);
        let mut n_to_c_file = std::fs::File::open(&n_to_c_file_path).expect("Could not open node_to_coset file for read.");
        let mut n_to_c_string = String::new();
        n_to_c_file.read_to_string(&mut n_to_c_string).expect("Could not read node_to_coset file.");
        let node_to_coset = serde_json::from_str(&n_to_c_string);
        if node_to_coset.is_ok() {
            self.node_to_coset = Some(node_to_coset.unwrap());
        } else {
            println!("Could not deserialize node_to_coset map.");
        }
    }

    pub fn generate_group(&mut self) {
        if self.subgroups.is_none() {
            let gens = compute_subgroups(self.dim, self.quotient.clone());
            self.subgroups = Some(gens);
        }
        if self.group.is_none() {
            if let Some(subs) = &self.subgroups {
                let g = compute_group(subs);
                self.group = Some(g);
            }
        }
    }

    pub fn compute_vertices(&mut self) {}

    pub fn save_to_disk(&self) {
        let mut subgroup_file_path = self.file_base.clone();
        subgroup_file_path.push_str(SUBGROUP_FILE_EXTENSION);
        let mut subgroup_file = std::fs::File::create(&self.file_base).expect("could not open file for writing.");
        if let Some(subs) = &self.subgroups {
            let subs_string = serde_json::to_string(subs).expect("Could not serialize subgroups.");
            subgroup_file.write(subs_string.as_bytes()).expect("Failed to write subgroups");
        } else {
            println!("Did not have subgroup to save to disk.");
        }

        let mut group_file_path = self.file_base.clone();
        group_file_path.push_str(GROUP_FILE_EXTENSION);
        let mut group_file = std::fs::File::create(&group_file_path).expect("could not create group file.");
        if let Some(group) = &self.group {
            let group_string = serde_json::to_string(group).expect("could not serialize group");
            group_file.write(group_string.as_bytes()).expect("Failed to write groups.");
        } else {
            println!("Tried to write non-existent group to disk.");
        }

        let mut hgraph_file_path = self.file_base.clone();
        hgraph_file_path.push_str(HGRAPH_FILE_EXTENSION);
        let mut hgraph_file = std::fs::File::create(hgraph_file_path).expect("Could not create hgraph file.");
        let hgraph_string = serde_json::to_string(&self.hgraph).expect("Could not serialize hgraph.");
        hgraph_file.write(hgraph_string.as_bytes()).expect("Could not write hgraph to file.");

        if let Some(n_to_c) = &self.node_to_coset {
            let mut n_to_c_file_path = self.file_base.clone();
            n_to_c_file_path.push_str(NODE_TO_COSET_FILE_EXTENSION);
            let mut n_to_c_file = std::fs::File::create(n_to_c_file_path).expect("Could not create node_to_coset file.");
            let n_to_c_string = serde_json::to_string(n_to_c).expect("Could not serialize node_to_cosest.");
            n_to_c_file.write(n_to_c_string.as_bytes()).expect("Could not write node_to_cosest to file.");
        } else {
            println!("Did not have node_to_coset complex to save to disk.");
        }
    }
}

mod tests {
    use std::collections::HashSet;

    use crate::math::{coset_complex::compute_coset, finite_field::FiniteField, matrix::PolyMatrix, polynomial::FiniteFieldPolynomial};

    use super::{compute_edges, compute_group, compute_subgroups, compute_triangles, compute_vertices, generate_all_polys, CosetGenerators, DiskManager};

    use deepsize::DeepSizeOf;
    use mhgl::HGraph;

    fn simplest_group() -> (CosetGenerators, HashSet<PolyMatrix>) {
        let p = 2_u32;
        let primitive_coeffs = [
            (2, (1, p).into()),
            (1, (2, p).into()),
            (0, (2, p).into())];
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
        let primitive_coeffs = [
            (2, (1, p).into()),
            (1, (2, p).into()),
            (0, (2, p).into())];
        let q = FiniteFieldPolynomial::from(&primitive_coeffs[..]);
        let mut gm = DiskManager::new(String::from("/Users/matt/repos/qec/data/groups/p_2_dim_3_deg_2.group"), 3, &q);
        gm.generate_group();
        gm.save_to_disk();
    }

    #[test]
    fn test_deserialize_groups() {
        let p = 2_u32;
        let primitive_coeffs = [
            (2, (1, p).into()),
            (1, (2, p).into()),
            (0, (2, p).into())];
        let q = FiniteFieldPolynomial::from(&primitive_coeffs[..]);
        let mut gm = DiskManager::new(String::from("/Users/matt/repos/qec/data/groups/p_2_dim_3_deg_2.group"), 3, &q);
        gm.load_from_disk();
    }

    fn get_nontrivial_group_manager() -> DiskManager {
        let p = 2_u32;
        let primitive_coeffs = [
            (2, (1, p).into()),
            (1, (2, p).into()),
            (0, (2, p).into())];
        let q = FiniteFieldPolynomial::from(&primitive_coeffs[..]);
        let mut gm = DiskManager::new(String::from("/Users/matt/repos/qec/data/groups/p_2_dim_3_deg_2.group"), 3, &q);
        gm.load_from_disk();
        gm
    }


    #[test]
    fn test_compute_group() {
        let (gens, g) = simplest_group();
        println!("size of subgroups: {:}", gens.type_to_generators.get(&0).unwrap().len());

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
        // println!("hg:\n{:}", hg);
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
        let primitive_coeffs = [
            (2, (1, p).into()),
            (1, (2, p).into()),
            (0, (2, p).into())];
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
}