use std::collections::HashMap;

use uuid::Uuid;

use crate::math::coset_complex::CosetComplex;
use crate::math::finite_field::FiniteField;
use crate::math::polynomial::FiniteFieldPolynomial;
use crate::reed_solomon::ReedSolomon;

/// Data necessary to encode and decode the classical codes from [New Codes on High Dimensional Expanders](https://arxiv.org/abs/2308.15563)
struct HDXCode {
    coset_complex: CosetComplex,
    /// Currently assume the same local code for each edge. In the future could look more like `edge_id_to_code: HashMap<Uuid, ReedSolomon>`
    local_code: ReedSolomon,
    /// Stores the encoded message on the triangles of the HDX. So we need a map from Uuid of the triangle to the message symbol.
    encoded_message: HashMap<Uuid, FiniteField>,
}

impl HDXCode {
    // pub fn new(coset_dir: String, dim: usize, quotient: FiniteFieldPolynomial) -> Self {
    //     let cc = CosetComplex::new(coset_dir, dim, quotient);
    //     let rs = ReedSolomon::new();
    //     HDXCode {
    //         coset_complex: cc,
    //         local_code: rs,
    //     }
    // }

    /// Returns the edges local view of the encoded message.
    fn get_sorted_local_view(&self, edge: (u32, u32)) -> Vec<FiniteField> {
        let mut triangles = self.coset_complex.get_triangles_containing_edge(edge);
        triangles.sort();
        triangles
            .into_iter()
            .map(|id| {
                self.encoded_message
                    .get(&id)
                    .expect("Given triangle has not been assigned a message symbol.")
            })
            .cloned()
            .collect()
    }
    fn compute_parity_check_matrix(&self) {
        // Going to assume a triangle complex for now.
        // The parity checks occur at each edge.
        // The bits are defined along the triangles.
        // I'm going to assume that the parity check acts
        // on column vectors. Therefore the matrix takes
        // in a vector defined over the triangles
        // and returns a vector defined over the
        // edges. So the matrix should be
        // # edges x # triangles.
        // As the parity check is linear
        // Can have entries in F_p or F_p[x] % q(x) ?
    }
}
