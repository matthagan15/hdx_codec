use core::num;
use std::{clone, collections::HashMap, rc::Rc};

use mhgl::HGraph;
use uuid::Uuid;

use crate::{
    code::Code,
    math::{
        ffmatrix::FFMatrix,
        finite_field::{FFRep, FiniteField as FF},
        polynomial::FiniteFieldPolynomial,
    },
};

#[derive(Debug, Clone, Hash, PartialEq, PartialOrd, Eq, Ord)]
enum Check {
    Node(u32),
    Edge(Uuid),
}

#[derive(Debug, Clone)]
struct TannerCode<C: Code> {
    id_to_message_ix: HashMap<Uuid, usize>,
    check_to_code: HashMap<Check, Rc<C>>,
    graph: HGraph,
    field_mod: FFRep,
}

impl TannerCode<ParityCode> {
    pub fn new(hgraph: HGraph, dim_of_message_symbols: usize, field_mod: FFRep) -> Self {
        if dim_of_message_symbols == 0 {
            panic!("Trying to store messages on nodes, not cool.")
        }
        if dim_of_message_symbols == 1 {
            // Using a graph version
            let nodes = hgraph.nodes();
            let mut deg_to_code: HashMap<usize, Rc<ParityCode>> = HashMap::new();
            let mut node_to_code = HashMap::new();
            for node in nodes.iter() {
                let deg = hgraph.degree(*node);
                if deg_to_code.contains_key(&deg) {
                    node_to_code.insert(Check::Node(*node), deg_to_code.get(&deg).unwrap().clone());
                } else {
                    let new_code = Rc::new(ParityCode::new(field_mod, deg));
                    deg_to_code.insert(deg, new_code.clone());
                    node_to_code.insert(Check::Node(*node), new_code);
                }
            }
            let message_spots = hgraph.edges_of_size(dim_of_message_symbols + 1);
            let id_to_message_ix = (0..message_spots.len())
                .into_iter()
                .map(|ix| (message_spots[ix], ix))
                .collect();
            TannerCode {
                id_to_message_ix,
                check_to_code: node_to_code,
                graph: hgraph,
                field_mod,
            }
        } else {
            // Using a hdx version
            let checks = hgraph.edges_of_size(dim_of_message_symbols);
            let mut deg_to_code: HashMap<usize, Rc<ParityCode>> = HashMap::new();
            let mut id_to_code = HashMap::new();
            for check in checks.iter() {
                let star = hgraph.star_id(check);
                if deg_to_code.contains_key(&star.len()) {
                    id_to_code.insert(
                        Check::Edge(*check),
                        deg_to_code.get(&star.len()).unwrap().clone(),
                    );
                } else {
                    let new_code = Rc::new(ParityCode::new(field_mod, star.len()));
                    deg_to_code.insert(star.len(), new_code.clone());
                    id_to_code.insert(Check::Edge(*check), new_code);
                }
            }
            let message_spots = hgraph.edges_of_size(dim_of_message_symbols + 1);
            let id_to_message_ix = (0..message_spots.len())
                .into_iter()
                .map(|ix| (message_spots[ix], ix))
                .collect();
            TannerCode {
                id_to_message_ix,
                check_to_code: id_to_code,
                graph: hgraph,
                field_mod,
            }
        }
    }

    pub fn get_check_view(&self, check: &Check, message: &Vec<FF>) -> Vec<FF> {
        match check {
            Check::Node(node) => {
                let mut containing_edges = self.graph.get_containing_edges(&[*node]);
                containing_edges.sort();
                // assuming that the user will not specify a sub-maximal dimension for message bits...
                let mut ret = Vec::new();
                for edge in containing_edges {
                    if let Some(ix) = self.id_to_message_ix.get(&edge) {
                        ret.push(
                            message
                                .get(*ix)
                                .expect("Could not find message symbol.")
                                .clone(),
                        );
                    }
                }
                ret
            }
            Check::Edge(id) => {
                let mut star = self.graph.star_id(id);
                star.sort();
                star.into_iter()
                    .map(|id| {
                        let ix = self.id_to_message_ix.get(&id).expect("ID not found");
                        message
                            .get(*ix)
                            .expect("could not find message symbol.")
                            .clone()
                    })
                    .collect()
            }
        }
    }

    fn get_single_parity_check(&self, check: &Check, message: &Vec<FF>) -> Vec<FF> {
        let check_view = self.get_check_view(check, message);
        let code = self
            .check_to_code
            .get(check)
            .expect("Check does not have local code.");
        code.parity_check(&check_view)
    }
}

impl Code for TannerCode<ParityCode> {
    fn encode(&self, message: &Vec<FF>) -> Vec<FF> {
        todo!()
    }

    fn decode(&self, encrypted: &Vec<FF>) -> Vec<FF> {
        todo!()
    }

    fn code_check(&self, encrypted: &Vec<FF>) -> bool {
        let pc = self.parity_check(encrypted);
        let mut are_all_zero = true;
        for symbol in encrypted.iter() {
            if symbol.0 != 0 {
                are_all_zero = false;
                break;
            }
        }
        are_all_zero
    }

    /// Returns the total parity check, local checks are combined into
    /// a final output simply using a sorted order on Uuid's, so it's
    /// essentially a random order.
    fn parity_check(&self, encrypted: &Vec<FF>) -> Vec<FF> {
        // Each local check needs to compute their local parity check
        let mut ret = Vec::new();
        for (check, code) in self.check_to_code.iter() {
            let check_view = self.get_check_view(check, encrypted);
            let local_parity_check = code.parity_check(&check_view);
            ret.push((check, local_parity_check));
        }
        ret.sort_by(|a, b| a.0.cmp(b.0));
        ret.into_iter().fold(Vec::new(), |mut acc, (_, mut v)| {
            acc.append(&mut v);
            acc
        })
    }

    fn parity_check_matrix(&self) -> FFMatrix {
        let message_len = self.id_to_message_ix.len();
        let mut zero: Vec<FF> = (0..message_len)
            .into_iter()
            .map(|_| FF::new(0, self.field_mod))
            .collect();
        let mut pc_cols = Vec::new();
        for ix in 0..message_len {
            zero[ix].0 = 1;
            let pc_col = self.parity_check(&zero);
            pc_cols.push(pc_col);
            zero[ix].0 = 0;
        }
        let mut entries = Vec::with_capacity(message_len * pc_cols[0].len());
        let n_rows = pc_cols[0].len();
        for row_ix in 0..n_rows {
            for col_ix in 0..message_len {
                let col = &pc_cols[col_ix];
                entries.push(col[row_ix]);
            }
        }
        FFMatrix::new(entries, n_rows, message_len)
    }
}

fn cycle_graph(num_nodes: u32) -> HGraph {
    let mut hg = HGraph::new();
    let nodes = hg.add_nodes(num_nodes as usize);
    for ix in 0..nodes.len() {
        hg.create_edge(&[nodes[ix], nodes[(ix + 1) % nodes.len()]]);
    }
    hg
}

/// A variable length parity check, as it's just the sum of all
/// the input messages. This is not really a code but I don't
/// want to make a new trait?
#[derive(Debug, Clone)]
struct ParityCode {
    field_mod: FFRep,
    message_len: usize,
}

impl ParityCode {
    pub fn new(field_mod: FFRep, message_len: usize) -> Self {
        ParityCode {
            field_mod,
            message_len,
        }
    }
}

impl Code for ParityCode {
    /// This doesn't really make sense for variable
    fn encode(&self, message: &Vec<FF>) -> Vec<FF> {
        todo!()
    }

    /// this also doesn't really make sense for parity.
    fn decode(&self, encrypted: &Vec<FF>) -> Vec<FF> {
        todo!()
    }

    fn code_check(&self, encrypted: &Vec<FF>) -> bool {
        let pc = self.parity_check(encrypted);
        pc[0].0 == 0
    }

    fn parity_check(&self, encrypted: &Vec<FF>) -> Vec<FF> {
        vec![FF::new(encrypted.iter().map(|x| x.0).sum(), self.field_mod)]
    }

    fn parity_check_matrix(&self) -> FFMatrix {
        let entries = (0..self.message_len)
            .into_iter()
            .map(|_| FF::new(1, self.field_mod))
            .collect();
        FFMatrix::new(entries, 1, self.message_len)
    }
}

mod tests {
    use std::collections::HashMap;

    use uuid::Uuid;

    use crate::{
        code::{get_generator_from_parity_check, Code},
        lps::compute_lps_graph,
        math::finite_field::FiniteField,
    };

    use super::{cycle_graph, ParityCode, TannerCode};

    #[test]
    fn test_parity_check_code() {
        let parity = ParityCode::new(2, 4);
        let one = FiniteField::new(1, 2);
        let z = FiniteField::new(0, 2);
        dbg!(parity.parity_check(&vec![one.clone(), z.clone(), one.clone(), one.clone()]));
        dbg!(parity.parity_check(&vec![z.clone(), z.clone(), z.clone(), z.clone()]));
    }

    #[test]
    fn test_rank_repitition_code() {
        let mut hg = cycle_graph(5);
        let nodes = hg.nodes();
        dbg!(hg.link_as_vec(&[0]));
        println!("hg: {:}", hg);
        let tc = TannerCode::new(hg, 1, 2);
        let pc = tc.parity_check(&vec![
            FiniteField::new(1, 2),
            FiniteField::new(0, 2),
            FiniteField::new(1, 2),
            FiniteField::new(1, 2),
            FiniteField::new(1, 2),
        ]);
        dbg!(pc);
        let mat = tc.parity_check_matrix();
        let rank = mat.rank();
        println!("Parity check matrix: {:}", mat);
        println!("Parity check rank: {:}", rank);
        println!("Rate: {:}", mat.n_rows - rank);
    }

    #[test]
    fn test_lps_code() {
        let lps = compute_lps_graph(7, 3).unwrap();
        println!("{:}", lps);
        let tc = TannerCode::new(lps, 1, 2);
        let mut mat = tc.parity_check_matrix();
        let rank = mat.rank();
        let rate = mat.n_cols - rank;
        println!("Mat: {:}", mat);
        println!("n_cols, rank, rate: {:}, {:}, {:}", mat.n_cols, rank, rate);
        mat.rref();
        println!("rref: {:}", mat);
        let gen = get_generator_from_parity_check(&mat);
        println!("Gen: {:}", gen);
    }
}
