use core::panic;
use std::{
    collections::HashMap,
    fs::File,
    io::{Read, Write},
    path::PathBuf,
    sync::{Arc, RwLock},
};

use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use serde::{Deserialize, Serialize};

use super::{
    finite_field::FFRep,
    galois_field::GaloisField,
    polynomial::{FFPolynomial, PolyDegree},
};
use crate::matrices::galois_matrix::GaloisMatrix;

/// Helper function for creating subgroups for BFS generation.
pub fn compute_deg(dim: usize, type_ix: usize, row_ix: usize, col_ix: usize) -> PolyDegree {
    if type_ix >= dim {
        panic!("Incorrect type_ix provided for degree computation.")
    }
    let mut how_many_over: HashMap<usize, usize> = HashMap::new();
    let mut row = type_ix;
    let mut hoppable_cols_left = (dim - 1) as i32;
    for _ in 0..dim {
        let e = how_many_over.entry(row).or_default();
        *e = hoppable_cols_left as usize;
        hoppable_cols_left -= 1;
        row += 1;
        row %= dim;
    }
    let hoppable_cols = how_many_over[&row_ix];
    let mut dist_from_diag = col_ix as i32 - row_ix as i32;
    while dist_from_diag < 0 {
        dist_from_diag += dim as i32;
    }
    dist_from_diag %= dim as i32;
    if dist_from_diag > (hoppable_cols as i32) {
        0
    } else {
        dist_from_diag.try_into().unwrap()
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct CosetRep {
    pub rep: GaloisMatrix,
    pub type_ix: u16,
}

impl Serialize for CosetRep {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        let mut s = String::new();
        let rep_serde = serde_json::to_string(&self.rep);
        if let Ok(rep_s) = rep_serde {
            s.push_str(&rep_s);
        } else {
            panic!("cannot serialize matrix rep")
        }
        s.push('@');
        s.push_str(&self.type_ix.to_string());
        serializer.serialize_str(&s)
    }
}

impl<'de> Deserialize<'de> for CosetRep {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let data = <String>::deserialize(deserializer)?;
        let splitted: Vec<&str> = data.split('@').collect();
        if splitted.len() != 2 {
            panic!("improper CosetRep encountered.")
        }
        let matrix = splitted[0];
        let type_ix = splitted[1].parse::<u16>();
        let mat_out = serde_json::from_str::<GaloisMatrix>(matrix);
        if mat_out.is_err() {
            panic!("Could not deserialize matrix");
        }
        let mat = mat_out.unwrap();
        Ok(CosetRep {
            rep: mat,
            type_ix: type_ix.expect("Could not parse CosetRep type"),
        })
    }
}
#[derive(Debug, Clone)]
pub struct KTypeSubgroup {
    generators: Vec<GaloisMatrix>,
    dim: usize,
    lookup: Arc<RwLock<GaloisField>>,
}
impl KTypeSubgroup {
    // TODO: Make this take ownership of an Arc! No need to take a reference to a reference.
    pub fn new(dim: usize, lookup: Arc<RwLock<GaloisField>>) -> Self {
        let mut type_zero = k_type_subgroup(dim, 0, lookup.clone());
        let mut type_one = k_type_subgroup(dim, 1, lookup.clone());
        let mut type_two = k_type_subgroup(dim, 2, lookup.clone());
        let mut gens = Vec::with_capacity(type_zero.len() + type_one.len() + type_two.len());
        gens.append(&mut type_zero);
        gens.append(&mut type_one);
        gens.append(&mut type_two);
        Self {
            generators: gens,
            dim,
            lookup: lookup.clone(),
        }
    }

    pub fn field_mod(&self) -> FFRep {
        self.lookup.read().unwrap().field_mod
    }

    pub fn generate_left_mul(&self, mat: &GaloisMatrix) -> Vec<GaloisMatrix> {
        self.generators
            .par_iter()
            .map(|h| h.mul(mat, self.lookup.clone()))
            .collect()
    }

    pub fn generate_right_mul(&self, mat: &GaloisMatrix) -> Vec<GaloisMatrix> {
        self.generators
            .par_iter()
            .map(|h| mat.mul(h, self.lookup.clone()))
            .collect()
    }

    /// Computes the coset reps given the output of a `generate_right_mul`
    pub fn coset_reps(&self, mats: &[GaloisMatrix]) -> Vec<CosetRep> {
        if mats.len() % self.dim != 0 {
            panic!("Improper number of matrices.")
        }
        let coset_len = mats.len() / self.dim;
        let mut start_ix = 0;
        let mut end_ix = coset_len;
        (0..self.dim)
            .map(|type_ix| {
                let rep = mats[start_ix..end_ix].iter().min().unwrap();
                start_ix += coset_len;
                end_ix += coset_len;
                CosetRep {
                    rep: rep.clone(),
                    type_ix: type_ix as u16,
                }
            })
            .collect()
    }

    pub fn get_coset_reps(&self, mat: &GaloisMatrix) -> Vec<CosetRep> {
        let total_cosets = self.generate_right_mul(mat);
        self.coset_reps(&total_cosets[..])
    }

    pub fn to_disk(&self, fiilename: PathBuf) {
        let mut buf = String::new();
        for mat in self.generators.iter() {
            if let Ok(s) = serde_json::to_string(mat) {
                buf.push_str(&s);
                buf.push(',');
            }
        }
        buf.pop();
        let mut file = File::create(fiilename).expect("Could not make file.");
        file.write_all(buf.as_bytes()).expect("Could not write.");
    }

    pub fn from_file(filename: PathBuf, lookup: &Arc<RwLock<GaloisField>>) -> Self {
        let mut buf = String::new();
        let mut file = File::open(&filename).expect("Cannot open for read.");
        file.read_to_string(&mut buf).expect("Cannot read.");
        let mut generators = Vec::new();
        for serialized_generator in buf.split(',') {
            let generator = serde_json::from_str::<GaloisMatrix>(serialized_generator)
                .expect("Cannot read matrix?");
            generators.push(generator);
        }
        if generators.len() == 0 {
            panic!("Did not find any generator matrices.")
        }
        let (n_rows, n_cols) = (generators[0].n_rows, generators[0].n_cols);
        if n_rows != n_cols {
            panic!("Only works for square generators")
        }
        Self {
            generators,
            dim: n_rows,
            lookup: lookup.clone(),
        }
    }
}

fn k_type_subgroup(
    dim: usize,
    type_ix: usize,
    lookup: Arc<RwLock<GaloisField>>,
) -> Vec<GaloisMatrix> {
    // let dim = 3;
    let p = lookup.read().unwrap().field_mod;
    let id = GaloisMatrix::id(dim);
    let mut ret = vec![id];
    // TODO: this is where the business logic goes
    for row_ix in 0..dim {
        for col_ix in 0..dim {
            let deg = compute_deg(dim, type_ix, row_ix, col_ix);
            if deg == 0 {
                continue;
            }
            let mut new_ret = Vec::new();
            for mut mat in ret.drain(..) {
                for a in 0..p {
                    let new_entry = FFPolynomial::monomial((a, p).into(), deg);
                    mat.set_entry(row_ix, col_ix, new_entry, lookup.clone());
                    new_ret.push(mat.clone());
                }
            }
            ret = new_ret;
        }
    }

    ret
}

#[cfg(test)]
mod tests {
    use std::{
        str::FromStr,
        sync::{Arc, RwLock},
    };

    use crate::{
        math::{
            coset_complex_subgroups::{compute_deg, k_type_subgroup},
            galois_field::GaloisField,
            polynomial::FFPolynomial,
        },
        matrices::ffmatrix::FFMatrix,
    };

    #[test]
    fn compute_degree() {
        // test for dim = 3 for now.
        let dim = 3;
        for type_ix in 0..dim {
            let mut mat = Vec::new();
            for row_ix in 0..dim {
                for col_ix in 0..dim {
                    mat.push((compute_deg(dim, type_ix, row_ix, col_ix), u32::MAX).into());
                }
            }
            let pretty_mat = FFMatrix::new(mat, dim, dim);
            println!("{:}", pretty_mat);
        }
    }

    #[test]
    fn k_type_subgroups() {
        let poly = FFPolynomial::from_str("1*x^3 + 2 * x^1 + 1 * x^0 % 3").unwrap();
        let lookup = Arc::new(RwLock::new(GaloisField::new(poly)));
        let subs = k_type_subgroup(3, 0, lookup.clone());
        println!("subs len: {:}", subs.len());
        for sub in subs {
            println!("{:}", sub.pretty_print(lookup.clone()));
        }
    }
}
