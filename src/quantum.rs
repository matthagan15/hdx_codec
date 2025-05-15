use std::collections::HashMap;

use mhgl::{HGraph, HyperGraph};

use crate::matrices::sparse_ffmatrix::{MemoryLayout, SparseFFMatrix};

/// Computes the linear map from `input_dim` -> `input_dim + 1`.
pub fn boundary_up<N, E>(
    hgraph: &HGraph<N, E>,
    mut input_edges: Vec<u64>,
    input_dim: usize,
) -> SparseFFMatrix {
    //    let mut input_edges = hgraph.edges_of_size(input_dim);
    input_edges.sort();
    let num_input_edges = input_edges.len();

    let input_id_to_ix: HashMap<u64, usize> =
        input_edges.into_iter().zip(0..num_input_edges).collect();
    let mut output_edges = hgraph.edges_of_size(input_dim + 1);
    output_edges.sort();
    let num_output_edges = output_edges.len();
    let output_edge_id_to_ix: HashMap<u64, usize> =
        output_edges.into_iter().zip(0..num_output_edges).collect();
    let mut entries = Vec::new();
    let mut count = 0;
    let mut n_rows = 0;
    for input_edge in input_id_to_ix.keys() {
        let local_boundary_up = hgraph.boundary_up(input_edge);
        if local_boundary_up.len() == 0 {
            println!("no boundary up on this edge: {input_edge}");
            println!(
                "containing_edges: {:?}",
                hgraph.containing_edges(input_edge)
            );
            println!("maximal edges: {:?}", hgraph.maximal_edges(input_edge));
            println!("nodes of edge: {:?}", hgraph.query_edge(input_edge));
        }
        let mut new_entries = local_boundary_up
            .into_iter()
            .map(|output_id| {
                let row_ix = *input_id_to_ix.get(input_edge).unwrap();
                n_rows = n_rows.max(row_ix);
                (
                    row_ix,
                    *output_edge_id_to_ix.get(&output_id).unwrap(),
                    1_u32,
                )
            })
            .collect();
        entries.append(&mut new_entries);
        count += 1;
    }
    let ret = SparseFFMatrix::new_with_entries(0, 0, 2, MemoryLayout::RowMajor, entries);
    ret
}

/// Computes the linear map from `input_dim` -> `input_dim - 1`.
pub fn boundary_down<N, E>(
    hgraph: &HGraph<N, E>,
    mut input_edges: Vec<u64>,
    input_dim: usize,
) -> SparseFFMatrix {
    if input_dim == 0 || input_dim != 2 {
        panic!("Cannot compute boundary down operator of just nodes.")
    }
    // let mut input_edges = hgraph.edges_of_size(input_dim);
    input_edges.sort();
    let num_input_edges = input_edges.len();
    let input_id_to_ix: HashMap<u64, usize> =
        input_edges.into_iter().zip(0..num_input_edges).collect();

    let mut output_edges = hgraph.nodes();
    output_edges.sort();
    let num_output_edges = output_edges.len();
    let output_edge_id_to_ix: HashMap<u32, usize> =
        output_edges.into_iter().zip(0..num_output_edges).collect();
    let entries = Vec::new();
    for input_edge in input_id_to_ix.keys() {
        panic!("there is a bug in the below code.");
        let local_boundary_down = hgraph.boundary_down(input_edge);
        let mut new_entries = local_boundary_down
            .into_iter()
            .map(|output_id| {
                (
                    *input_id_to_ix.get(input_edge).unwrap(),
                    *output_edge_id_to_ix.get(&(output_id as u32)).unwrap(),
                    1_u32,
                )
            })
            .collect();
        entries.append(&mut new_entries);
    }
    SparseFFMatrix::new_with_entries(0, 0, 2, MemoryLayout::RowMajor, entries)
}

#[cfg(test)]
mod test {
    use mhgl::HGraph;

    use super::{boundary_down, boundary_up};

    #[test]
    fn simple_example() {
        let mut hg = HGraph::<(), ()>::new();
        hg.add_nodes(5);
        let e1 = hg.add_edge([0, 1], ());
        let e2 = hg.add_edge([0, 1, 2], ());
        let e3 = hg.add_edge([0, 2, 3, 4], ());
        let e4 = hg.add_edge([0, 1, 2, 3], ());
        let bd = boundary_down(&hg, vec![e2], 3);
        let bu = boundary_up(&hg, vec![e2], 3);
        dbg!(&bd);
        dbg!(&bu);
        println!("boundary up. Dims: {:} x {:}", bu.n_rows, bu.n_cols);
        bu.dense_print();
        println!("boundary down. Dims: {:} x {:}", bd.n_rows, bd.n_cols);
        bd.dense_print();
    }
}
