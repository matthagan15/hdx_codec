use std::{
    collections::{HashMap, HashSet},
    hash::Hash,
    rc::Rc,
};

use log::trace;
use mhgl::{ConGraph, HyperGraph};

use crate::{
    code::Code,
    math::finite_field::{FFRep, FiniteField as FF},
    matrices::{
        ffmatrix::FFMatrix,
        sparse_ffmatrix::{MemoryLayout, SparseFFMatrix},
        sparse_vec::SparseVector,
    },
    reed_solomon::ReedSolomon,
};

#[derive(Debug, Clone, Hash, PartialEq, PartialOrd, Eq, Ord)]
pub enum Check {
    Node(u32),
    Edge(u64),
}

#[derive(Debug, Clone)]
/// Currently only supports parity checks that are all of the same
/// dimension.
pub struct TannerCode<C: Code> {
    id_to_message_ix: HashMap<u64, usize>,
    check_to_code: HashMap<Check, Rc<C>>,
    /// Note the returned ranges are INCLUSIVE, so a check owns outputs (0, 4) means if pc = \[1,2,3,4,5,6,7] then the check produced outputs \[1,2,3,4,5].
    check_to_output_bounds: HashMap<Check, (usize, usize)>,
    #[allow(dead_code)]
    parity_check_len: usize,
    graph: ConGraph,
    field_mod: FFRep,
}

impl TannerCode<ParityCode> {
    pub fn new(hgraph: ConGraph, dim_of_message_symbols: usize, field_mod: FFRep) -> Self {
        if dim_of_message_symbols == 0 {
            panic!("Trying to store messages on nodes, not cool.")
        }
        let mut id_to_message_ix = HashMap::new();
        let mut check_to_code = HashMap::new();
        let mut check_to_output_bounds = HashMap::new();
        if dim_of_message_symbols == 1 {
            // Using a graph version
            let nodes = hgraph.nodes();
            let mut deg_to_code: HashMap<usize, Rc<ParityCode>> = HashMap::new();

            for node in nodes.iter() {
                let link = hgraph.link_of_nodes([*node]);
                let deg = link.len();
                if deg_to_code.contains_key(&deg) {
                    check_to_code
                        .insert(Check::Node(*node), deg_to_code.get(&deg).unwrap().clone());
                } else {
                    let new_code = Rc::new(ParityCode::new(field_mod, deg));
                    deg_to_code.insert(deg, new_code.clone());
                    check_to_code.insert(Check::Node(*node), new_code);
                }
            }
            let message_spots = hgraph.edges_of_size(dim_of_message_symbols + 1);
            for ix in 0..message_spots.len() {
                id_to_message_ix.insert(message_spots[ix], ix);
            }
        } else {
            // Using a hdx version
            let checks = hgraph.edges_of_size(dim_of_message_symbols);
            let mut deg_to_code: HashMap<usize, Rc<ParityCode>> = HashMap::new();
            for check in checks.iter() {
                let star = hgraph.maximal_edges(check);
                if deg_to_code.contains_key(&star.len()) {
                    check_to_code.insert(
                        Check::Edge(*check),
                        deg_to_code.get(&star.len()).unwrap().clone(),
                    );
                } else {
                    let new_code = Rc::new(ParityCode::new(field_mod, star.len()));
                    deg_to_code.insert(star.len(), new_code.clone());
                    check_to_code.insert(Check::Edge(*check), new_code);
                }
            }
            let message_spots = hgraph.edges_of_size(dim_of_message_symbols + 1);
            for ix in 0..message_spots.len() {
                id_to_message_ix.insert(message_spots[ix], ix);
            }
        }
        let mut prev_upper_bound: Option<usize> = Some(0);
        let mut parity_check_len = 0;
        let mut checks: Vec<Check> = check_to_code.keys().cloned().collect();
        checks.sort();
        for check in checks {
            let code = check_to_code
                .get(&check)
                .expect("Just received key from map.");
            parity_check_len += code.parity_check_len();
            let output_bounds = match prev_upper_bound {
                Some(prev_upper) => (prev_upper + 1, prev_upper + 1 + code.parity_check_len()),
                None => (0, code.parity_check_len()),
            };
            prev_upper_bound = Some(output_bounds.1);
            check_to_output_bounds.insert(check, output_bounds);
        }
        TannerCode {
            id_to_message_ix,
            check_to_code,
            check_to_output_bounds,
            parity_check_len,
            graph: hgraph,
            field_mod,
        }
    }
}

impl TannerCode<ReedSolomon> {
    pub fn new(
        hgraph: ConGraph,
        dim_of_message_symbols: usize,
        field_mod: FFRep,
        quotient_degree: usize,
    ) -> Self {
        trace!("Constructing new TannerCode<ReedSolomon>");
        if dim_of_message_symbols == 0 {
            panic!("Trying to store messages on nodes, not cool.")
        }
        let mut id_to_message_ix = HashMap::new();
        let check_to_code: HashMap<Check, Rc<ReedSolomon>>;
        let mut check_to_output_bounds = HashMap::new();

        // currently we restrict checks to those that have
        // maximal degree.
        if dim_of_message_symbols == 1 {
            panic!("I don't want to implement ReedSolomon TannerCode for graphs just yet.")
        } else {
            // Using a hdx version
            // note that edges of size takes a cardinality, so the
            // dim_of_message_symbols is equal to the cardinality of the checks
            // and the cardinality of the messages is dim_of_message_symbols + 1

            // there are a few checks I need to do
            // First we get all the check faces. Then we go through and
            // find the maximal degree of all checks, where the degree
            // is the number of message spots (maximal faces) containing
            // the check. Once we have these, then the message spots is the
            // union of the maximal edges containing all
            let checks = hgraph.edges_of_size(dim_of_message_symbols);
            let mut max_deg = 0;
            let check_to_maximal_edges: HashMap<_, _> = checks
                .into_iter()
                .map(|check| {
                    let max_edges = hgraph.maximal_edges(&check);
                    max_deg = max_deg.max(max_edges.len());
                    (check, max_edges)
                })
                .collect();
            trace!("Maximum check degree: {:}", max_deg);
            let code = Rc::new(ReedSolomon::new(field_mod, quotient_degree));
            let mut message_spots: HashSet<u64> = HashSet::new();
            check_to_code = check_to_maximal_edges
                .into_iter()
                .filter_map(|(check, maximal_edges)| {
                    if maximal_edges.len() == max_deg {
                        for message_spot in maximal_edges.iter() {
                            message_spots.insert(*message_spot);
                        }
                        Some((Check::Edge(check), code.clone()))
                    } else {
                        None
                    }
                })
                .collect();

            let mut count = 0;
            for check_id in message_spots.drain() {
                id_to_message_ix.insert(check_id, count);
                count += 1;
            }
        }
        let mut prev_upper_bound: Option<usize> = Some(0);
        let mut parity_check_len = 0;
        for check in check_to_code.keys() {
            let code = check_to_code
                .get(&check)
                .expect("Just received key from map.");
            parity_check_len += code.parity_check_len();
            let output_bounds = match prev_upper_bound {
                Some(prev_upper) => (prev_upper + 1, prev_upper + 1 + code.parity_check_len()),
                None => (0, code.parity_check_len()),
            };
            prev_upper_bound = Some(output_bounds.1);
            check_to_output_bounds.insert(check.clone(), output_bounds);
        }
        TannerCode {
            id_to_message_ix,
            check_to_code,
            check_to_output_bounds,
            parity_check_len,
            graph: hgraph,
            field_mod,
        }
    }
}

impl<C: Code> TannerCode<C> {
    pub fn get_check_view(&self, check: &Check, message: &Vec<FF>) -> Vec<FF> {
        match check {
            Check::Node(node) => {
                let mut containing_edges = self.graph.containing_edges_of_nodes(&[*node]);
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
                let mut star = self.graph.maximal_edges(id);
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

    /// returns a sorted dense view of a given checks visibility of the sparse
    /// message.
    pub fn get_check_view_sparse(&self, check: &Check, message: &SparseVector) -> Vec<FF> {
        match check {
            Check::Node(node) => {
                let mut containing_edges = self.graph.containing_edges_of_nodes(&[*node]);
                containing_edges.sort();
                // assuming that the user will not specify a sub-maximal dimension for message bits...
                let mut ret = Vec::new();
                for edge in containing_edges {
                    if let Some(ix) = self.id_to_message_ix.get(&edge) {
                        ret.push(FF::new(message.query(ix), self.field_mod));
                    }
                }
                ret
            }
            Check::Edge(id) => {
                let mut star = self.graph.maximal_edges(id);
                star.sort();
                star.into_iter()
                    .map(|id| {
                        let ix = self.id_to_message_ix.get(&id).expect("ID not found");
                        FF::new(message.query(ix), self.field_mod)
                    })
                    .collect()
            }
        }
    }

    pub fn sparse_parity_check_matrix(&self) -> SparseFFMatrix {
        let message_len = self.id_to_message_ix.len();
        let _zero: Vec<FF> = (0..message_len)
            .into_iter()
            .map(|_| FF::new(0, self.field_mod))
            .collect();
        let mut parity_check_cols = Vec::new();
        for ix in 0..message_len {
            let sparse_message = SparseVector::new_with_entries(vec![(ix, 1)]);

            let sparse_check_output = self.parity_check_sparse(&sparse_message);
            parity_check_cols.push(sparse_check_output);
        }
        SparseFFMatrix::new_with_sections(
            self.field_mod,
            MemoryLayout::ColMajor,
            &parity_check_cols,
        )
    }

    /// Returns the total parity check, local checks are combined into
    /// a final output simply using a sorted order on Uuid's, so it's
    /// essentially a random order.
    fn parity_check_sparse(&self, encrypted: &SparseVector) -> SparseVector {
        // Each local check needs to compute their local parity check
        let mut ret = SparseVector::new_empty();
        for (check, code) in self.check_to_code.iter() {
            let check_view = self.get_check_view_sparse(check, encrypted);
            let local_parity_check = code.parity_check(&check_view);
            let ix_bounds = self
                .check_to_output_bounds
                .get(check)
                .expect("Check does not have output upper bounds.");
            let sparse_entries = (ix_bounds.0..=ix_bounds.1)
                .into_iter()
                .zip(local_parity_check.into_iter().map(|ff| ff.0))
                .collect();
            let local_out = SparseVector::new_with_entries(sparse_entries);
            ret.add_to_self(&local_out, self.field_mod);
        }
        ret
    }
}

impl Code for TannerCode<ParityCode> {
    fn encode(&self, _message: &Vec<FF>) -> Vec<FF> {
        todo!()
    }

    fn decode(&self, _encrypted: &Vec<FF>) -> Vec<FF> {
        todo!()
    }

    fn is_codeword(&self, encrypted: &Vec<FF>) -> bool {
        let _pc = self.parity_check(encrypted);
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
                entries.push(col[row_ix].0);
            }
        }
        FFMatrix::new(entries, n_rows, message_len, self.field_mod)
    }

    fn parity_check_len(&self) -> usize {
        todo!()
    }

    fn message_len(&self) -> usize {
        self.id_to_message_ix.len()
    }
}

/// A variable length parity check, as it's just the sum of all
/// the input messages. This is not really a code but I don't
/// want to make a new trait?
#[derive(Debug, Clone)]
pub struct ParityCode {
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
    fn encode(&self, _message: &Vec<FF>) -> Vec<FF> {
        todo!()
    }

    /// this also doesn't really make sense for parity.
    fn decode(&self, _encrypted: &Vec<FF>) -> Vec<FF> {
        todo!()
    }

    fn is_codeword(&self, encrypted: &Vec<FF>) -> bool {
        let pc = self.parity_check(encrypted);
        pc[0].0 == 0
    }

    fn parity_check(&self, encrypted: &Vec<FF>) -> Vec<FF> {
        vec![FF::new(encrypted.iter().map(|x| x.0).sum(), self.field_mod)]
    }

    fn parity_check_matrix(&self) -> FFMatrix {
        let entries = (0..self.message_len).into_iter().map(|_| 1).collect();
        FFMatrix::new(entries, 1, self.message_len, self.field_mod)
    }

    fn parity_check_len(&self) -> usize {
        1
    }

    fn message_len(&self) -> usize {
        self.message_len
    }
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use mhgl::{ConGraph, HyperGraph};

    use crate::{
        code::{get_generator_from_parity_check, Code},
        lps::compute_lps_graph,
        math::{coset_complex_bfs::bfs, finite_field::FiniteField, polynomial::FFPolynomial},
        matrices,
        reed_solomon::ReedSolomon,
    };

    use super::{ParityCode, TannerCode};

    fn cycle_graph(num_nodes: u32) -> ConGraph {
        let mut hg = ConGraph::new();
        let nodes = hg.add_nodes(num_nodes as usize);
        for ix in 0..nodes.len() {
            hg.add_edge(&[nodes[ix], nodes[(ix + 1) % nodes.len()]]);
        }
        hg
    }

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
        let hg = cycle_graph(5);
        let _nodes = hg.nodes();
        dbg!(hg.link_of_nodes(&[0]));
        println!("hg: {:}", hg);
        let tc = TannerCode::<ParityCode>::new(hg, 1, 2);
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
        let tc = TannerCode::<ParityCode>::new(lps, 1, 2);
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

    #[test]
    fn coset_complex_code() {
        let poly = FFPolynomial::from_str("1*x^2 + 2*x^1 + 2*x^0 % 3").unwrap();
        let hg = bfs(poly, 3, Some(1000), None, None);
        // println!("{:}", hg);
        for i in 0..3 {
            let max_edges = hg.maximal_edges_of_nodes([i]);
            let _containing_edges = hg.containing_edges_of_nodes([i]);
            println!("{:} max edges containing {:}", max_edges.len(), i);
            for e in max_edges {
                println!("e: {:?}", hg.query_edge(&e).unwrap());
            }
        }
    }

    #[test]
    fn test_complex_code() {
        let mut hg = ConGraph::new();
        let _nodes = hg.add_nodes(15);
        let _e11 = hg.add_edge(&[0, 1]);
        let _e12 = hg.add_edge(&[0, 2]);
        let _e17 = hg.add_edge(&[0, 3]);
        let _e13 = hg.add_edge(&[0, 4]);
        let _e14 = hg.add_edge(&[1, 2]);
        let _e15 = hg.add_edge(&[1, 3]);
        let _e16 = hg.add_edge(&[1, 4]);
        let _e18 = hg.add_edge(&[2, 3]);
        let _e19 = hg.add_edge(&[3, 4]);
        let _e20 = hg.add_edge(&[0, 5]);
        let _e21 = hg.add_edge(&[0, 6]);
        let _e22 = hg.add_edge(&[0, 7]);
        let _e23 = hg.add_edge(&[0, 8]);
        let _e24 = hg.add_edge(&[0, 9]);
        let _e25 = hg.add_edge(&[5, 6]);
        let _e26 = hg.add_edge(&[5, 7]);
        let _e27 = hg.add_edge(&[8, 9]);

        let _e1 = hg.add_edge(&[0, 1, 2]);
        let _e2 = hg.add_edge(&[0, 1, 3]);
        let _e3 = hg.add_edge(&[0, 1, 4]);
        let _e4 = hg.add_edge(&[0, 5, 6]);
        let _e5 = hg.add_edge(&[0, 5, 7]);
        let _e6 = hg.add_edge(&[0, 8, 9]);
        let _e30 = hg.add_edge(&[0, 5, 11]);
        println!("hg: {:}", hg);
        let tc = TannerCode::<ReedSolomon>::new(hg, 2, 3, 2);

        let _message_raw: Vec<FiniteField> = vec![
            (1, 3).into(),
            (0, 3).into(),
            (0, 3).into(),
            (0, 3).into(),
            (1, 3).into(),
            (0, 3).into(),
        ];
        let _message = matrices::sparse_vec::SparseVector(vec![(0, 1), (3, 1)]);
        let mut mat = tc.sparse_parity_check_matrix();
        mat.swap_layout();
        println!("parity check matrix {:}", mat.clone().to_dense());
        println!("rank per dim: {:}", mat.rank() as f64 / mat.n_cols as f64);
    }
}
