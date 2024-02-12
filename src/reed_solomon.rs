use serde::{Deserialize, Serialize};

use crate::math::{
    ffmatrix::{vandermonde, FFMatrix}, finite_field::FiniteField, polynomial::FiniteFieldPolynomial
};

use crate::tanner_code::*;

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct ReedSolomon {
    message_len: usize,
    encoded_len: usize,
    parity_checker: FFMatrix,
    generator_matrix: FFMatrix,
    field_mod: u32,
}

impl ReedSolomon {
    pub fn field_mod(&self) -> u32 { self.field_mod }
    pub fn message_len(&self) -> usize { self.message_len }
    pub fn encoded_len(&self) -> usize { self.encoded_len }


    /// `max_degree` is not inclusive!
    pub fn new(field_mod: u32, max_degree: usize) -> Self {
        let eval_points: Vec<FiniteField> = (0..field_mod).into_iter().map(|c| FiniteField::new(c, field_mod)).collect();
        let mut generator_matrix = vandermonde(&eval_points, max_degree);
        generator_matrix.transpose();
        Self {
            message_len: 1,
            encoded_len: 1,
            parity_checker: get_parity_check_matrix_from(&generator_matrix),
            generator_matrix,
            field_mod,
        }
    }

    /// Encode using matrix-vector multiplication with the generator matrix $G$ associated with the code acting on the message input vector $m$. As the columns of $G$ are guaranteed to be in the kernel of the parity check matrix $H$, we get that any message here has a unique output in the encoded space. Note if your message is not the right size this function will `panic!`.
    pub fn encode(&self, message: Vec<FiniteField>) -> Vec<FiniteField> {
        let parsed_message = if message.len() > self.message_len {
            Vec::from(&message[0..self.message_len])
        } else if message.len() == self.message_len {
            message
        } else {
            let mut v = message.clone();
            for _ in 0..(self.message_len - message.len()) {
                v.push((0, self.field_mod).into());
            }
            v
        };
        let mut deg_coeff_buf = Vec::new();
        for d in 0..parsed_message.len() {
            if parsed_message[d].0 != 0 {
                deg_coeff_buf.push((d, parsed_message[d]));
            }
        }
        if deg_coeff_buf.len() == 0 {
            deg_coeff_buf.push((0, FiniteField(0, self.field_mod)));
        }
        let poly = FiniteFieldPolynomial::from(&deg_coeff_buf[..]);
        (0..self.field_mod)
            .into_iter()
            .map(|alpha| {
                poly.evaluate(&(alpha, self.field_mod).into())
            })
            .collect()
    }

    /// Returns true if the given `message` is in the codespace and false 
    /// otherwise. Works via matrix - vector multiplication, so takes time O()
    pub fn is_message_in_code(&self, message: &Vec<FiniteField>) -> bool {
        let parity_checks = &self.parity_checker * message;
        let mut is_zero = true;
        for pc in parity_checks {
            if pc.0 != 0 {
                is_zero = false;
                break;
            }
        }
        is_zero
    }

    pub fn parity_check(&self, encoded_message: &Vec<FiniteField>) -> Vec<FiniteField> {
        &self.parity_checker * encoded_message
    }

    /// Decodes the given message using the algorithm given by Shuhong Gao in this [paper](http://www.math.clemson.edu/~sgao/papers/RS.pdf). 
    pub fn decode(&self, encoded_message: Vec<FiniteField>) -> Option<Vec<FiniteField>> {
        let mut all_zero = true;
        for e_i in encoded_message.iter() {
            if e_i.0 != 0 {
                all_zero = false;
            }
        }
        if all_zero {
            let mut ret = Vec::new();
            for _ in 0..self.message_len {
                ret.push(FiniteField::from((0, self.field_mod)));
            }
            return Some(ret);
        }

        let mut g0 = FiniteFieldPolynomial::constant(1, self.field_mod);
        for alpha in (0..self.field_mod).into_iter() {
            let tmp_coeffs = vec![(0, FiniteField::from((alpha as i32 * -1_i32, self.field_mod))), (1, FiniteField::new(1, self.field_mod))];
            let tmp_factor = FiniteFieldPolynomial::from(&tmp_coeffs[..]);
            g0 *= tmp_factor;
        }

        let interpolation_points = (0..self.field_mod)
            .map(|k| (FiniteField::new(k, self.field_mod), encoded_message[k as usize]))
            .collect();
        let interpolated_poly = FiniteFieldPolynomial::interpolation(interpolation_points);
        let deg_cutoff = (self.field_mod as usize + self.message_len) / 2;
        let (u, v, g) = g0.partial_gcd(&interpolated_poly, deg_cutoff);
        let (f, r) = g / v;
        if r.is_zero() && f.degree() < self.message_len {
            let mut ret = Vec::new();
            for ix in 0..self.message_len {
                if f.coeffs.contains_key(&ix) {
                    ret.push(f.coeffs.get(&ix).unwrap().clone());
                } else {
                    ret.push((0, self.field_mod).into())
                }
            }
            Some(ret)
        } else {
            None
        }
    }
}

impl crate::code::Code for ReedSolomon {
    fn encode(message: Vec<FiniteField>) -> Vec<FiniteField> {
        todo!()
    }

    fn decode(encrypted: Vec<FiniteField>) -> Vec<FiniteField> {
        todo!()
    }

    fn code_check(word: &Vec<FiniteField>) -> bool {
        todo!()
    }

    fn parity_check(&self, message: &Vec<FiniteField>) -> Vec<FiniteField> {
        todo!()
    }
}

mod tests {
    use crate::math::finite_field::FiniteField;

    use super::ReedSolomon;

    // fn smallest_rs() -> ReedSolomon {
    //     let p = 3u32;
    //     let eval_points = vec![FiniteField(0, p), FiniteField(1, p), FiniteField(2, p)];
    //     ReedSolomon {
    //         evaluation_points: eval_points,
    //         message_len: 2,
    //         field_mod: p,
    //     }
    // }

    #[test]
    fn test_new_new() {
        let rs = ReedSolomon::new(11, 5);
        println!("generator: {:}", rs.generator_matrix);
        println!("parity check matrix: {:}", rs.parity_checker);
        let g = rs.generator_matrix.clone();
        let h = rs.parity_checker.clone();
        let multiplied = &h * &g;
        println!("parity_check * generator: {:}", multiplied);
    }
    #[test]
    fn test_small_small_example() {
        let p = 3u32;
        // let rs = smallest_rs();
        // for a in 0..p {
        //     for b in 0..p {
        //         let input = vec![(a, p).into(), (b, p).into()];
        //         let out = rs.encode(input.clone());
        //         println!("message: ({:}, {:})", a, b);
        //         println!("encoding: ({:}, {:}, {:})", out[0].0, out[1].0, out[2].0);
        //         let decoded = rs.decode(out);
        //         assert_eq!(input, decoded.unwrap());
        //     }
        // }
    }
}
