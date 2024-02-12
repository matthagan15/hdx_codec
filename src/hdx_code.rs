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

    pub fn is_message_in_code(&self, message: &HashMap<Uuid, FiniteField>) -> bool {
        let edge_ids = self.coset_complex.hgraph_ref().edges_of_size(2);
        let mut pass_all_checks = true;
        for edge_id in edge_ids {
            let edge_nodes = self.coset_complex.hgraph_ref().query_edge_id(&edge_id).expect("That edge better be in there!");
            if edge_nodes.len() != 2 {
                panic!("edge is not of size 2");
            }
            let mut triangles = self.coset_complex.get_triangles_containing_edge((edge_nodes[0], edge_nodes[1]));
            triangles.sort();
            let local_view: Vec<FiniteField> = triangles
                .into_iter()
                .map(|id| {
                    message.get(&id).expect("Triangle is not included in message.")
                })
                .cloned()
                .collect();
            if self.local_code.is_message_in_code(&local_view) == false {
                pass_all_checks = false;
                break;
            }
        }
        pass_all_checks
    }

}
