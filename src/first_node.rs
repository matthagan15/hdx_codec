use std::{
    collections::{HashMap, HashSet},
    sync::{Arc, RwLock},
    usize,
};

use mhgl::{HGraph, HyperGraph};

use crate::{
    code::Code,
    math::{
        coset_complex_bfs::BFSState,
        coset_complex_subgroups::KTypeSubgroup,
        finite_field::{FFRep, FiniteField},
        galois_field::GaloisField,
        polynomial::FFPolynomial,
    },
    matrices::sparse_ffmatrix::{MemoryLayout, SparseFFMatrix},
    reed_solomon::ReedSolomon,
};

type NodeData = u16;

pub fn get_first_node_complete_star(quotient: FFPolynomial, dim: usize) -> HGraph<NodeData, ()> {
    let mut bfs_manager = BFSState::new(quotient.clone(), dim, vec![usize::MAX]);
    let lookup = Arc::new(RwLock::new(GaloisField::new(quotient.clone())));
    let mut hg = HGraph::new();
    while bfs_manager.current_bfs_distance < 3 {
        bfs_manager.step(&KTypeSubgroup::new(dim, lookup.clone()), &mut hg);
    }
    hg.star([0])
}

pub fn compute_classical_parity_check(quotient: FFPolynomial, dim: usize) {
    let star = get_first_node_complete_star(quotient.clone(), dim);
    let checks = star.edges_of_size(dim - 1);
    let message_ids = star.edges_of_size(dim);
    let local_code = ReedSolomon::new(quotient.field_mod, (quotient.field_mod - 1) as usize);
    let mut mat = build_matrices_from_checks(checks, message_ids, &star, &local_code);
    println!("Computed hgraph, now reducing matrix.");
    let (rank, mut grassmannian) = mat.rref();

    grassmannian.swap_layout();
    let mut min_nnz = usize::MAX;
    println!(
        "grassmannian num cols: {:}",
        grassmannian.ix_to_section.len()
    );
    if grassmannian.memory_layout == MemoryLayout::ColMajor {
        for (col_ix, col) in grassmannian.ix_to_section.iter() {
            min_nnz = min_nnz.min(col.nnz());
        }
    }
    println!("rank: {rank}");
    println!("Distance: {:}", min_nnz + 1);
}

/// Constructs the interior and border check matrices from the given parity checks.
/// The interior matrix will be in row-echelon form and the border matrix will be as
/// reduced as possible from the interior check pivots.
///
/// Assumes provided inputs satisfy the following:
/// - each parity check, whether interior or border, has a complete check view for the provided
///     local code.
/// - Each message_id has at least one parity check that sees the given message symbol.
fn build_matrices_from_checks(
    checks: Vec<u64>,
    message_ids: Vec<u64>,
    hgraph: &HGraph<u16, ()>,
    local_code: &ReedSolomon,
) -> SparseFFMatrix {
    let num_messages = message_ids.len();
    let message_id_to_ix: HashMap<u64, usize> =
        message_ids.into_iter().zip(0..num_messages).collect();
    let mut message_id_to_checks: HashMap<u64, HashSet<u64>> =
        HashMap::with_capacity(message_id_to_ix.len());
    let check_to_message_ixs: HashMap<_, _> = checks
        .iter()
        .map(|check| {
            let max_edges = hgraph.maximal_edges(&check);
            let mut max_ixs = Vec::new();
            for max_edge in max_edges {
                let max_edge_to_checks = message_id_to_checks.entry(max_edge).or_default();
                max_edge_to_checks.insert(*check);
                max_ixs.push(*message_id_to_ix.get(&max_edge).expect("Interior check was collected that did not have a corresponding message_id accounted for."));
            }
            max_ixs.sort();
            max_ixs.dedup();
            (*check, max_ixs)
        })
        .collect();

    let mut interior_entries: Vec<(usize, usize, FFRep)> = Vec::new();
    let num_rows_per_check = local_code.parity_check_len();
    let mut counter = 0;
    let interior_check_to_row_ix_offset: HashMap<u64, usize> = check_to_message_ixs
        .keys()
        .map(|check| {
            let ret = counter * num_rows_per_check;
            counter += 1;
            (*check, ret)
        })
        .collect();
    for (message_id, message_ix) in message_id_to_ix {
        // These are the checks that can see the current message_id

        let viewing_checks = message_id_to_checks.get(&message_id).unwrap();
        for viewing_check in viewing_checks.iter() {
            let local_view = check_to_message_ixs.get(viewing_check).unwrap();
            let parity_check_input = local_view
                .into_iter()
                .map(|ix| {
                    FiniteField::new(
                        if *ix == message_ix { 1 } else { 0 },
                        local_code.field_mod(),
                    )
                })
                .collect();
            let parity_check = local_code.parity_check(&parity_check_input);
            for ix in 0..parity_check.len() {
                let row_ix = ix + interior_check_to_row_ix_offset.get(viewing_check).unwrap();
                let col_ix = message_ix;
                let entry = parity_check[ix].0;
                interior_entries.push((row_ix, col_ix, entry));
            }
        }
    }
    SparseFFMatrix::new_with_entries(
        0,
        0,
        local_code.field_mod(),
        MemoryLayout::RowMajor,
        interior_entries,
    )
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use crate::{first_node::compute_classical_parity_check, math::polynomial::FFPolynomial};

    #[test]
    fn small_example_first_node() {
        let q = FFPolynomial::from_str("1x^2 + 4x + 2 % 5").unwrap();
        let dim = 3;
        println!("q: {:}", q);
        if q.is_irreducible() == false {
            panic!()
        }
        compute_classical_parity_check(q, dim);
    }
}
