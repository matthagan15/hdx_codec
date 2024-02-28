use std::collections::HashSet;

use crate::math::{
    polymatrix::{dim_three_det, PolyMatrix},
    polynomial::FiniteFieldPolynomial,
};

pub struct MatrixEnumerator {
    pub indices: Vec<usize>,
    pub dimension: usize,
    pub polys: Vec<FiniteFieldPolynomial>,
    pub quotient: FiniteFieldPolynomial,
}

impl MatrixEnumerator {
    /// returns true if indices can be incremented again, false if it cannot be incremented.
    pub fn increment_indices(&mut self) -> bool {
        let mut carry_increment = true;
        let mut ix = 0;
        while carry_increment {
            let mut tmp = self.indices[ix];
            tmp += 1;
            tmp %= self.polys.len();
            self.indices[ix] = tmp;
            if tmp != 0 {
                carry_increment = false;
            } else {
                if ix == self.indices.len() - 1 {
                    return false;
                }
                carry_increment = true;
                ix += 1;
            }
        }
        true
    }

    pub fn construct_matrix(&self) -> PolyMatrix {
        let entries: Vec<FiniteFieldPolynomial> = self
            .indices
            .clone()
            .into_iter()
            .map(|ix| self.polys[ix].clone())
            .collect();
        let p = entries[0].field_mod;
        PolyMatrix {
            entries,
            n_rows: self.dimension,
            n_cols: self.dimension,
            field_mod: p,
            quotient: self.quotient.clone(),
        }
    }

    pub fn generate_sl3(&mut self) -> HashSet<PolyMatrix> {
        let upper_bound = self.polys.len().pow(self.indices.len() as u32);
        let mut counter = 0;
        let mut can_increment = true;
        let mut sl3 = HashSet::new();
        while can_increment {
            if counter % 10000 == 0 {
                let percent_done = counter as f64 / upper_bound as f64;
                println!("{:.4} % done", percent_done);
            }
            counter += 1;
            let m = self.construct_matrix();
            if dim_three_det(&m).is_one() {
                sl3.insert(m);
            }
            can_increment = self.increment_indices();
        }
        sl3
    }
}
