use core::panic;
use std::{
    collections::{HashMap, VecDeque},
    path::PathBuf,
    str::FromStr,
    sync::{Arc, RwLock},
    time::Instant,
    usize,
};

use mhgl::{HGraph, HyperGraph};
use serde::{ser::SerializeStruct, Deserialize, Serialize};

use super::{
    coset_complex_subgroups::KTypeSubgroup, galois_field::GaloisField, polynomial::FFPolynomial,
};
use crate::{matrices::galois_matrix::GaloisMatrix, powerset};

pub type NodeData = u16;

/// Returns the number of maximal faces that will be in the coset complex, aka
/// computes the size of SL(dim, q = p^n). Returns `None` if the size is larger than
/// `usize::MAX`.
pub fn size_of_coset_complex(quotient: &FFPolynomial, dim: usize) -> Option<usize> {
    let p = quotient.field_mod;
    let n = quotient.degree();
    let q = p.pow(n.try_into().unwrap()) as usize;
    let mut prod: usize = 1;
    let q_dim = q.checked_pow(dim as u32);
    if q_dim.is_none() {
        return None;
    }
    let q_dim = q_dim.unwrap();
    for k in 0..dim {
        prod = prod.checked_mul(q_dim
            - q.checked_pow(k as u32).expect(
                "Coset Complex too big to even compute size of with fixed width unsigned int.",
            )).expect("Coset Complex too big to compute size of with fixed width ints.");
    }
    Some(prod / (q - 1))
}

pub fn bfs_benchmark() {
    let q = FFPolynomial::from_str("1*x^2 + 2 * x^1 + 2 * x^0 % 3").unwrap();
    let dim = 3;
    let lookup = Arc::new(RwLock::new(GaloisField::new(q.clone())));
    let subgroups = KTypeSubgroup::new(dim, lookup.clone());
    let mut bfs = BFSState::new(q.clone(), dim, Some(1000));
    let mut hg = HGraph::new();
    let start = Instant::now();
    for _ in 0..1000 {
        bfs.step(&subgroups, &mut hg);
    }
    println!("time elapsed: {:}", start.elapsed().as_secs_f64());
}

/// Filters the nodes in a coset complex to find the type_ix = 0 nodes that have the maximal
/// number of maximally containing faces.
pub fn find_complete_nodes(hg: &HGraph<u16, ()>, quotient: &FFPolynomial, dim: usize) -> Vec<u32> {
    let lookup = Arc::new(RwLock::new(GaloisField::new(quotient.clone())));
    let subgroup_generators = KTypeSubgroup::new(dim, lookup);
    let num_maximal_faces = subgroup_generators.size_of_single_subgroup();
    let filter = |node_id| {
        hg.maximal_edges_of_nodes([node_id]).len() == num_maximal_faces
            && hg.get_node(&node_id) == Some(&0)
    };
    hg.filter_nodes(filter)
}

/// Returns the hgraph associated with the coset complex breadth first search
pub fn bfs(
    quotient: FFPolynomial,
    dim: usize,
    truncation: Option<usize>,
    hgraph_output_file: Option<PathBuf>,
    log_rate: Option<usize>,
) -> HGraph<u16, ()> {
    // if dim != 3 {
    //     panic!("Only dimension 3 matrices are currently supported.")
    // }
    let mut bfs_state = BFSState::new(quotient.clone(), dim, truncation);
    if log_rate.is_some() {
        log::trace!("-------------- BFS -------------");
        log::trace!(
        "Frontier: {:}\nvisited: {:}\ncurrent_distance: {:}\nnum_matrices_completed: {:}\ntruncation: {:}",
        bfs_state.frontier.len(),
        bfs_state.visited.len(),
        bfs_state.current_bfs_distance,
        bfs_state.num_matrices_completed,
        bfs_state.truncation,
    );
    }
    let mut hg = HGraph::new();
    let mut num_steps = 0;
    let mut time_in_step = 0.0;
    let max_number_steps = match (size_of_coset_complex(&quotient, dim), truncation) {
        (None, None) => usize::MAX,
        (None, Some(trunc)) => trunc,
        (Some(coset_complex_size), None) => {
            if coset_complex_size == usize::MAX {
                usize::MAX
            } else {
                coset_complex_size + 1
            }
        }
        (Some(coset_complex_size), Some(trunc)) => coset_complex_size.min(trunc),
    };
    let lookup = Arc::new(RwLock::new(GaloisField::new(quotient.clone())));
    let subgroup_generators = KTypeSubgroup::new(dim, lookup);
    while bfs_state.num_matrices_completed < truncation.unwrap_or(usize::MAX)
        && bfs_state.frontier.len() > 0
    {
        let start = Instant::now();
        let step_nodes_with_types = bfs_state.step(&subgroup_generators, &mut hg);
        time_in_step += start.elapsed().as_secs_f64();
        assert_eq!(step_nodes_with_types.len(), dim);
        num_steps += 1;
        if num_steps % log_rate.unwrap_or(usize::MAX) == 0 {
            log::trace!("Reached Logging checkpoint. Here are some stats.\nnum_matrices_completed: {:}, frontier: {:}, bfs_distance: {:}", bfs_state.num_matrices_completed, bfs_state.frontier.len(),bfs_state.current_bfs_distance );
            let time_per_step = time_in_step / (num_steps as f64);
            let num_steps_remaining = max_number_steps - bfs_state.num_matrices_completed;
            let time_remaining = time_per_step * num_steps_remaining as f64;
            log::trace!(
                "Estimated {:} seconds remaining. {:} seconds per step",
                time_remaining,
                time_per_step
            );
        }
    }
    if log_rate.is_some() {
        log::trace!("BFS complete!");
    }

    if let Some(output_file) = hgraph_output_file {
        log::trace!(
            "Saving bfs state and hgraph to {:}.",
            output_file.with_extension("hg").as_path().display()
        );
        hg.to_disk(output_file.with_extension("hg").as_path());
    }
    hg
}

/// Helper struct for BFS. Each node in a BFS for a coset complex consists of a single
/// matrix over the dimension of the BFS. Each matrix can be associated with a `dim` number
/// of hypergraph nodes, one for each type. `hg_nodes` contains these nodes and `distance`
/// records the number of neighbors away from the identity matrix this BFS node is.
#[derive(Debug, Clone, Serialize, Deserialize)]
struct GroupBFSNode {
    distance: u32,
    // Pairs of nodes and their associated types.
    hg_nodes: Vec<(u32, NodeData)>,
}

#[derive(Debug, Clone)]
pub struct BFSState {
    /// A vector of matrices (which correspond to maximal faces) and their corresponding
    /// distance from the identity matrix.
    frontier: VecDeque<(GaloisMatrix, u32)>,
    visited: HashMap<GaloisMatrix, GroupBFSNode>,
    pub current_bfs_distance: u32,
    last_flushed_distance: u32,
    num_matrices_completed: usize,
    truncation: usize,
    last_cached_matrices_done: u64,
    quotient: FFPolynomial,
    dim: usize,
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
        s.serialize_field("truncation", &self.truncation)?;
        s.serialize_field("last_cached_matrices_done", &self.last_cached_matrices_done)?;
        s.serialize_field("quotient", &self.quotient)?;
        s.serialize_field("dim", &self.dim)?;
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
        let truncation: Result<usize, _> = serde_json::from_value(data["truncation"].take());
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
        Ok(BFSState {
            frontier: frontier.expect("Could not deserialize frontier"),
            visited,
            current_bfs_distance: current_bfs_distance.unwrap(),
            last_flushed_distance: last_flushed_distance.unwrap(),
            num_matrices_completed: num_matrices_completed.unwrap(),
            truncation: truncation.unwrap(),
            last_cached_matrices_done: last_cached_matrices_done.unwrap(),
            quotient: quotient.unwrap(),
            dim: dim.unwrap(),
        })
    }
}

impl BFSState {
    pub fn new(quotient: FFPolynomial, dim: usize, truncation: Option<usize>) -> Self {
        let mut frontier = VecDeque::new();
        frontier.push_back((GaloisMatrix::id(dim), 0));
        Self {
            frontier,
            visited: HashMap::new(),
            current_bfs_distance: 0,
            last_flushed_distance: 0,
            num_matrices_completed: 0,
            truncation: truncation.unwrap_or(usize::MAX),
            last_cached_matrices_done: 0,
            quotient,
            dim,
        }
    }
    /// Returns newly discovered edge_id's.
    pub fn step(
        &mut self,
        subgroup_generators: &KTypeSubgroup,
        hg: &mut HGraph<u16, ()>,
    ) -> Vec<(u32, NodeData)> {
        if self.frontier.is_empty() {
            return Vec::new();
        }
        let x = self.frontier.pop_front().unwrap();
        // process this matrix first, compute the cosets and triangles it can
        // be a part of.
        let neighbors = subgroup_generators.generate_right_mul(&x.0);
        let coset_reps = subgroup_generators.coset_reps(&neighbors[..]);

        // flush visited and coset to node
        self.current_bfs_distance = x.1;

        for neighbor in neighbors {
            let neighbor_bfs = GroupBFSNode {
                distance: x.1 + 1,
                hg_nodes: Vec::new(),
            };
            if self.visited.contains_key(&neighbor) == false {
                self.visited.insert(neighbor.clone(), neighbor_bfs);
                self.frontier.push_back((neighbor, x.1 + 1));
            }
        }
        self.num_matrices_completed += 1;
        // Map each coset rep to the nodes in the hgraph, adding new ones if need be.
        let nodes_and_types = coset_reps
            .iter()
            .map(|coset_rep| {
                let coset_rep_bfs_node = self
                    .visited
                    .get_mut(&coset_rep.rep)
                    .expect("Previously added a BFS node when processing the neighbors.");
                let found_nodes = coset_rep_bfs_node
                    .hg_nodes
                    .iter()
                    .filter_map(|(node, type_ix)| {
                        if *type_ix == coset_rep.type_ix {
                            Some(*node)
                        } else {
                            None
                        }
                    })
                    .collect::<Vec<u32>>();
                if found_nodes.len() == 0 {
                    let new_node = hg.add_node(coset_rep.type_ix);
                    coset_rep_bfs_node
                        .hg_nodes
                        .push((new_node, coset_rep.type_ix));
                    (new_node, coset_rep.type_ix)
                } else if found_nodes.len() == 1 {
                    (found_nodes[0], coset_rep.type_ix)
                } else {
                    panic!("Found more than one node of a given type per coset rep.")
                }
            })
            .collect::<Vec<(u32, NodeData)>>();
        for possible_edge in powerset(
            nodes_and_types
                .iter()
                .map(|(node, _type_ix)| *node)
                .collect(),
        ) {
            if possible_edge.len() < 2 {
                continue;
            }
            if hg.find_id(&possible_edge[..]).is_none() {
                hg.add_edge(possible_edge, ());
            }
        }
        nodes_and_types
    }
}

#[cfg(test)]
mod tests {
    use super::bfs;
    use crate::{
        math::{
            coset_complex_bfs::{find_complete_nodes, BFSState},
            coset_complex_subgroups::KTypeSubgroup,
            finite_field::FFRep,
            galois_field::GaloisField,
            polynomial::FFPolynomial,
        },
        matrices::galois_matrix::GaloisMatrix,
    };
    use mhgl::HGraph;

    use std::{
        str::FromStr,
        sync::{Arc, RwLock},
    };

    #[test]
    fn small_checkpoints() {
        let p: FFRep = 3;
        let dim = 3;
        let quotient = FFPolynomial::from(
            &vec![(0, (2, p).into()), (1, (2, p).into()), (2, (1, p).into())][..],
        );
        let mut bfs_manager = BFSState::new(quotient.clone(), dim, None);
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
    fn small_complete_nodes() {
        let q = FFPolynomial::from_str("x^2 + 2x + 2 % 3").unwrap();
        let dim = 3;
        let hg = bfs(q.clone(), dim, Some(100), None, None);
        dbg!(find_complete_nodes(&hg, &q, dim));
    }
}
