use serde::{Deserialize, Serialize};

use crate::math::{
    ffmatrix::{vandermonde, FFMatrix},
    finite_field::{FFRep, FiniteField},
    group_ring_field::Ring,
    polynomial::{FiniteFieldPolynomial, PolyDegree},
};

use crate::code::*;
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
    pub fn field_mod(&self) -> u32 {
        self.field_mod
    }
    pub fn encoded_len(&self) -> usize {
        self.encoded_len
    }

    pub fn parity_check_matrix(&self) -> FFMatrix {
        self.parity_checker.clone()
    }

    pub fn print(&self) {
        println!("{:}", "#".repeat(50));
        println!(
            "# {:} Reed-Solomon({:}, {:}) {:} #",
            " ".repeat(13),
            self.encoded_len,
            self.message_len,
            " ".repeat(13)
        );
        println!("{:}", "#".repeat(50));
        println!("PARITY CHECK MATRIX");
        println!("{:}", self.parity_checker);
        println!(
            "n_rows, ncols = {:}, {:}",
            self.parity_checker.n_rows, self.parity_checker.n_cols
        );
        println!("\nGENERATOR MATRIX\n{:}", self.generator_matrix);
    }

    /// `max_degree` is not inclusive!
    pub fn new(field_mod: u32, max_degree: usize) -> Self {
        let eval_points: Vec<FiniteField> = (0..field_mod)
            .into_iter()
            .map(|c| FiniteField::new(c, field_mod))
            .collect();
        let mut generator_matrix = vandermonde(&eval_points, max_degree);
        generator_matrix.transpose();

        let pc = if generator_matrix.is_square() {
            FFMatrix::zero(1, 1, field_mod)
        } else {
            get_parity_check_matrix_from(&generator_matrix)
        };
        Self {
            message_len: generator_matrix.n_cols,
            encoded_len: generator_matrix.n_rows,
            parity_checker: pc,
            generator_matrix,
            field_mod,
        }
    }

    /// ReedSolomon(n, k) where `num_eval_points` is the number of logical symbols, which is determined by the number of points a message polynomial is evaluated on to produce an encoded output. and n is the number of evaluation points.
    pub fn new_with_parity_check_input(
        num_eval_points: usize,
        max_degree: usize,
        field_mod: u32,
    ) -> Self {
        // TODO: random or deterministic?
        // TODO: This generator matrix business right now is pretty sketch
        // when the num_eval_points and degree get low. idk what to do.
        let eval_points: Vec<FiniteField> = (0..num_eval_points)
            .into_iter()
            .map(|c| FiniteField::new(c as FFRep, field_mod))
            .collect();
        let mut generator_matrix = vandermonde(&eval_points, max_degree);
        if generator_matrix.n_cols > generator_matrix.n_rows {
            generator_matrix.transpose();
        }
        // println!("getting parity check from generator matrix\n{:}", generator_matrix);
        let pc = get_parity_check_matrix_from(&generator_matrix);
        // println!("reed solomon constructor.");
        // println!("pc nrows: {:}", pc.n_rows);
        Self {
            message_len: generator_matrix.n_cols,
            encoded_len: generator_matrix.n_rows,
            parity_checker: pc,
            generator_matrix,
            field_mod,
        }
    }

    /// Encode using matrix-vector multiplication with the generator matrix $G$ associated with the code acting on the message input vector $m$. As the columns of $G$ are guaranteed to be in the kernel of the parity check matrix $H$, we get that any message here has a unique output in the encoded space. Note if your message is not the right size this function will `panic!`.
    pub fn manual_encode(&self, message: Vec<FiniteField>) -> Vec<FiniteField> {
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
                deg_coeff_buf.push((d as PolyDegree, parsed_message[d]));
            }
        }
        if deg_coeff_buf.len() == 0 {
            deg_coeff_buf.push((0, FiniteField(0, self.field_mod)));
        }
        let poly = FiniteFieldPolynomial::from(&deg_coeff_buf[..]);
        (0..self.field_mod)
            .into_iter()
            .map(|alpha| poly.evaluate(&(alpha, self.field_mod).into()))
            .collect()
    }

    /// Returns true if the given `message` is in the codespace and false
    /// otherwise. Works via matrix - vector multiplication, so takes time O()
    pub fn is_message_in_code(&self, message: &Vec<FiniteField>) -> bool {
        if (self.parity_checker.n_rows, self.parity_checker.n_cols) == (1, 1) {
            return true;
        }
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
        if (self.parity_checker.n_rows, self.parity_checker.n_cols) == (1, 1) {
            vec![FiniteField::new(0, self.field_mod)]
        } else {
            let ret = &self.parity_checker * encoded_message;
            println!(
                "[ReedSolomon::parity_check] message: {:?}, parity check: {:?}",
                encoded_message, ret
            );
            println!("parity check matrix: {:}", self.parity_checker);
            ret
        }
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
            let tmp_coeffs = vec![
                (
                    0,
                    FiniteField::from((alpha as i32 * -1_i32, self.field_mod)),
                ),
                (1, FiniteField::new(1, self.field_mod)),
            ];
            let tmp_factor = FiniteFieldPolynomial::from(&tmp_coeffs[..]);
            g0 *= tmp_factor;
        }

        let interpolation_points = (0..self.field_mod)
            .map(|k| {
                (
                    FiniteField::new(k, self.field_mod),
                    encoded_message[k as usize],
                )
            })
            .collect();
        let interpolated_poly = FiniteFieldPolynomial::interpolation(interpolation_points);
        let deg_cutoff = (self.field_mod as usize + self.message_len) / 2;
        let (u, v, g) = g0.partial_gcd(&interpolated_poly, deg_cutoff as PolyDegree);
        let (f, r) = g / v;
        if r.is_zero() && (f.degree() as usize) < self.message_len {
            let mut ret = Vec::new();
            for ix in 0..self.message_len {
                let ix = ix as PolyDegree;
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

impl Code for ReedSolomon {
    fn encode(&self, message: &Vec<FiniteField>) -> Vec<FiniteField> {
        if message.len() != self.generator_matrix.n_cols {
            println!(
                "Message to encode must be of length: {:}",
                self.generator_matrix.n_cols
            );
            panic!()
        }
        &self.generator_matrix * message
    }

    fn decode(&self, encrypted: &Vec<FiniteField>) -> Vec<FiniteField> {
        todo!()
    }

    fn is_codeword(&self, word: &Vec<FiniteField>) -> bool {
        let parity_check = self.parity_check(word);
        let mut entries_all_zero = true;
        for e in parity_check {
            if e.0 != 0 {
                entries_all_zero = false;
            }
        }
        entries_all_zero
    }

    fn parity_check(&self, message: &Vec<FiniteField>) -> Vec<FiniteField> {
        &self.parity_checker * message
    }

    fn parity_check_matrix(&self) -> FFMatrix {
        self.parity_checker.clone()
    }

    /// length of the output parity check
    fn parity_check_len(&self) -> usize {
        self.parity_checker.n_rows
    }

    fn message_len(&self) -> usize {
        self.message_len()
    }
}

mod tests {
    use crate::math::finite_field::FiniteField;

    use super::ReedSolomon;

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
        let p = 199;
        let rs = ReedSolomon::new_with_parity_check_input(p as usize, 50, p);
        rs.print();
        let pc = rs.parity_check_matrix();
        println!("rank: {:}", pc.rank());
    }
}
