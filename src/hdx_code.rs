use crate::math::coset_complex::CosetComplex;
use crate::reed_solomon::ReedSolomon;

struct HdxCode {
    coset_complex: CosetComplex,
    local_code: ReedSolomon,
}

impl HdxCode {
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
