use crate::math::{polynomial::FiniteFieldPolynomial, finite_field::FiniteField};

struct ReedSolomon {
    evaluation_points: Vec<FiniteField>,
    input_length: usize,
    field_mod: u32,
}

impl ReedSolomon {
    fn encode(&self, message: Vec<FiniteField>) -> Vec<FiniteField> {
        let parsed_message = if message.len() > self.input_length {
            Vec::from(&message[0..self.input_length])
        } else if message.len() == self.input_length {
            message.clone()
        } else {
            let mut v = message.clone();
            for _ in 0..(self.input_length - message.len()) {
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
            deg_coeff_buf.push((0, FiniteField (0, self.field_mod)));
        }
        let poly = FiniteFieldPolynomial::from(&deg_coeff_buf[..]);
        self.evaluation_points.iter().map(|alpha| {poly.evaluate(alpha)}).collect()
    }

    fn decode(&self, encoded_message: Vec<FiniteField>) -> Option<Vec<FiniteField>> {
        let mut all_zero = true;
        for e_i in encoded_message.iter() {
            if e_i.0 != 0 {
                all_zero = false;
            }
        }
        if all_zero {
            let mut ret = Vec::new();
            for _ in 0..self.input_length {
                ret.push(FiniteField::from((0, self.field_mod)));
            }
            return Some(ret);
        }

        let mut g0 = FiniteFieldPolynomial::constant(1, self.field_mod);
        for alpha in self.evaluation_points.iter() {
            let tmp_coeffs = vec![(0, *alpha * -1_i32), (1, (1, self.field_mod).into())];
            let tmp_factor = FiniteFieldPolynomial::from(&tmp_coeffs[..]);
            g0 *= tmp_factor;
        }
        let interpolation_points = (0..self.evaluation_points.len()).map(|k| {
            (self.evaluation_points[k], encoded_message[k])
        }).collect();
        let interpolated_poly = FiniteFieldPolynomial::interpolation(interpolation_points);
        let deg_cutoff = (self.evaluation_points.len() + self.input_length) / 2;
        let (u, v, g) = g0.partial_gcd(&interpolated_poly, deg_cutoff);
        let (f, r) = g / v;
        if r.is_zero() && f.degree() < self.input_length {
            let mut ret = Vec::new();
            for ix in 0..self.input_length {
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

mod tests {
    use crate::math::finite_field::FiniteField;

    use super::ReedSolomon;

    fn smallest_rs() -> ReedSolomon {
        let p = 3u32;
        let eval_points = vec![
            FiniteField(0, p), FiniteField(1, p), FiniteField(2, p)
        ];
        ReedSolomon {
            evaluation_points: eval_points,
            input_length: 2,
            field_mod: p,
        }
    }

    #[test]
    fn test_small_small_example() {
        let p = 3u32;
        let rs = smallest_rs();
        for a in 0..p {
            for b in 0..p {
                let input = vec![(a, p).into(), (b, p).into()];
                let out = rs.encode(input.clone());
                println!("message: ({:}, {:})", a, b);
                println!("encoding: ({:}, {:}, {:})", out[0].0, out[1].0, out[2].0);
                let decoded = rs.decode(out);
                assert_eq!(input, decoded.unwrap());
            }
        }
    }

}
