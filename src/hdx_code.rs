use crate::math::coset_complex::CosetComplex;
use crate::reed_solomon::ReedSolomon;

struct HdxCode {
    coset_complex: CosetComplex,
    local_code: ReedSolomon
}