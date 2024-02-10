use mhgl::HGraph;
use uuid::Uuid;

use crate::math::polynomial::FiniteFieldPolynomial;

fn build_parity_check_matrix(
    hgraph: &HGraph,
    parity_checks: Vec<Uuid>,
    message_spots: Vec<Uuid>,
    field_mod: u32,
    quotient: FiniteFieldPolynomial,
) {
    // Iterate over local codewords input space
    for ix in 0..message_spots.len() {}
    for c in 1..field_mod {}
}
