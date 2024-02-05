use std::{collections::{HashMap, HashSet}, ops::{Add, Mul}};

use super::finite_field::FiniteField;

struct Multinomial {
    degrees_to_coeff: HashMap<Vec<usize>, u32>,
    num_indeterminates: usize,
    max_degree: usize,
    field_mod: u32,
}

impl<'a> Add<&'a Multinomial> for &'a Multinomial {
    type Output = Multinomial;

    fn add(self, rhs: Self) -> Self::Output {
        // currently unchecked.
        let mut new_deg_to_c = HashMap::new();
        let mut max_degree = 0;
        for (k, v) in self.degrees_to_coeff.iter() {
            let e = if rhs.degrees_to_coeff.contains_key(k) {
                *rhs.degrees_to_coeff.get(k).unwrap()
            } else {
                0_u32
            };
            new_deg_to_c.insert(k.clone(), (v + e) % self.field_mod);
        }
        for (k, v) in rhs.degrees_to_coeff.iter() {
            if self.degrees_to_coeff.contains_key(k) {
                continue;
            }
            let e = new_deg_to_c.entry(k.clone()).or_default();
            *e += *v;
        }
        let mut keys_for_deletion = HashSet::new();
        for (k, v) in new_deg_to_c.iter() {
            if *v == 0 {
                keys_for_deletion.insert(k.clone());
                continue;
            }
            let deg = k.iter().fold(0 as usize, |a, b| a + *b);
            max_degree = max_degree.max(deg);
        }
        for k in keys_for_deletion.into_iter() {
            new_deg_to_c.remove(&k);
        }
        Multinomial {
            degrees_to_coeff: new_deg_to_c,
            num_indeterminates: self.num_indeterminates,
            max_degree,
            field_mod: self.field_mod,
        }
    }
}