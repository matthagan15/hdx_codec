use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::{Read, Write},
    path::PathBuf,
    sync::Arc,
};

use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use serde::{Deserialize, Serialize};

use super::{
    finite_field::FiniteField,
    galois_field::GaloisField,
    polynomial::{FFPolynomial, PolyDegree},
};
use crate::matrices::galois_matrix::GaloisMatrix;

/// Helper function for creating subgroups for BFS generation.
pub fn compute_deg(dim: usize, type_ix: usize, row_ix: usize, col_ix: usize) -> PolyDegree {
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
    type_zero_end: usize,
    type_one_end: usize,
    lookup: Arc<GaloisField>,
}
impl KTypeSubgroup {
    pub fn new(lookup: &Arc<GaloisField>) -> Self {
        let mut type_zero = k_type_subgroup(0, lookup.clone());
        let mut type_one = k_type_subgroup(1, lookup.clone());
        let mut type_two = k_type_subgroup(2, lookup.clone());
        let mut gens = Vec::with_capacity(type_zero.len() + type_one.len() + type_two.len());
        gens.append(&mut type_zero);
        let type_zero_end = gens.len() - 1;
        gens.append(&mut type_one);
        let type_one_end = gens.len() - 1;
        gens.append(&mut type_two);
        Self {
            generators: gens,
            type_zero_end,
            type_one_end,
            lookup: lookup.clone(),
        }
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
    pub fn coset_reps(&self, mats: &[GaloisMatrix]) -> (CosetRep, CosetRep, CosetRep) {
        let rep0 = mats[..=self.type_zero_end].iter().min().unwrap();
        let rep1 = mats[self.type_zero_end + 1..=self.type_one_end]
            .iter()
            .min()
            .unwrap();
        let rep2 = mats[self.type_one_end + 1..].iter().min().unwrap();
        let c0 = CosetRep {
            rep: rep0.clone(),
            type_ix: 0,
        };
        let c1 = CosetRep {
            rep: rep1.clone(),
            type_ix: 1,
        };
        let c2 = CosetRep {
            rep: rep2.clone(),
            type_ix: 2,
        };
        (c0, c1, c2)
    }

    pub fn get_coset_reps(&self, mat: &GaloisMatrix) -> (CosetRep, CosetRep, CosetRep) {
        let total_cosets = self.generate_right_mul(mat);
        let (coset0, coset1, coset2) = (
            total_cosets[..=self.type_zero_end].to_vec(),
            total_cosets[self.type_zero_end + 1..=self.type_one_end].to_vec(),
            total_cosets[self.type_one_end + 1..].to_vec(),
        );
        let rep0 = coset0.iter().min().unwrap();
        let rep1 = coset1.iter().min().unwrap();
        let rep2 = coset2.iter().min().unwrap();
        let c0 = CosetRep {
            rep: rep0.clone(),
            type_ix: 0,
        };
        let c1 = CosetRep {
            rep: rep1.clone(),
            type_ix: 1,
        };
        let c2 = CosetRep {
            rep: rep2.clone(),
            type_ix: 2,
        };
        (c0, c1, c2)
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

    pub fn from_file(filename: PathBuf, lookup: &Arc<GaloisField>) -> Self {
        let mut buf = String::new();
        let mut file = File::open(&filename).expect("Cannot open for read.");
        file.read_to_string(&mut buf).expect("Cannot read.");
        let mut generators = Vec::new();
        for serialized_generator in buf.split(',') {
            let generator = serde_json::from_str::<GaloisMatrix>(serialized_generator)
                .expect("Cannot read matrix?");
            generators.push(generator);
        }
        if generators.len() % 3 != 0 {
            panic!("I only work for dim 3 matrices.")
        }
        let type_zero_end = generators.len() / 3;
        let type_one_end = 2 * type_zero_end;
        Self {
            generators,
            type_zero_end,
            type_one_end,
            lookup: lookup.clone(),
        }
    }
}

fn k_type_subgroup(type_ix: usize, lookup: Arc<GaloisField>) -> Vec<GaloisMatrix> {
    let dim = 3;
    let p = lookup.field_mod;
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

#[derive(Debug, Clone)]
struct HTypeSubgroup {
    generators: Vec<GaloisMatrix>,
    type_zero_end: usize,
    type_one_end: usize,
    lookup: Arc<GaloisField>,
}

impl HTypeSubgroup {
    pub fn new(lookup: &Arc<GaloisField>) -> Self {
        let mut type_zero = h_type_subgroup(0, lookup.clone());
        let mut type_one = h_type_subgroup(1, lookup.clone());
        let mut type_two = h_type_subgroup(2, lookup.clone());
        let mut gens = Vec::with_capacity(type_zero.len() + type_one.len() + type_two.len());
        gens.append(&mut type_zero);
        let type_zero_end = gens.len() - 1;
        gens.append(&mut type_one);
        let type_one_end = gens.len() - 1;
        gens.append(&mut type_two);
        Self {
            generators: gens,
            type_zero_end,
            type_one_end,
            lookup: lookup.clone(),
        }
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

    pub fn from_file(filename: PathBuf, lookup: &Arc<GaloisField>) -> Self {
        let mut buf = String::new();
        let mut file = File::open(&filename).expect("Cannot open for read.");
        file.read_to_string(&mut buf).expect("Cannot read.");
        let mut generators = Vec::new();
        for serialized_generator in buf.split(',') {
            let generator = serde_json::from_str::<GaloisMatrix>(serialized_generator)
                .expect("Cannot read matrix?");
            generators.push(generator);
        }
        if generators.len() % 3 != 0 {
            panic!("I only work for dim 3 matrices.")
        }
        let type_zero_end = generators.len() / 3;
        let type_one_end = 2 * type_zero_end;
        Self {
            generators,
            type_zero_end,
            type_one_end,
            lookup: lookup.clone(),
        }
    }
}

fn h_type_subgroup(type_ix: usize, lookup: Arc<GaloisField>) -> Vec<GaloisMatrix> {
    let mut ret = Vec::new();
    let dim = 3;
    let p = lookup.field_mod;
    let id = GaloisMatrix::id(dim);
    let mut row_ix = type_ix as i32 - 1;
    while row_ix <= 0 {
        row_ix += dim as i32;
    }
    row_ix %= dim as i32;
    for a in 0..p {
        let mut tmp = id.clone();
        let new_entry = FFPolynomial::monomial(FiniteField::new(a, p), 1);
        tmp.set_entry(row_ix as usize, type_ix, new_entry, lookup.clone());
        ret.push(tmp);
    }
    ret
}

#[derive(Debug, Clone)]
struct ParallelSubgroups {
    generators: Vec<GaloisMatrix>,
    type_zero_end: usize,
    type_one_end: usize,
    lookup: Arc<GaloisField>,
}

impl ParallelSubgroups {
    pub fn print_gens(&self) {
        println!("{:}", "#".repeat(50));
        println!(
            "Have {:} matrices across all subgroups.",
            self.generators.len()
        );
        for g in self.generators.iter() {
            println!("{:}", g.pretty_print(self.lookup.clone()));
        }
        println!("{:}", "#".repeat(50));
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

    pub fn from_file(filename: PathBuf, lookup: &Arc<GaloisField>) -> Self {
        let mut buf = String::new();
        let mut file = File::open(&filename).expect("Cannot open for read.");
        file.read_to_string(&mut buf).expect("Cannot read.");
        let mut generators = Vec::new();
        for serialized_generator in buf.split(',') {
            let generator = serde_json::from_str::<GaloisMatrix>(serialized_generator)
                .expect("Cannot read matrix?");
            generators.push(generator);
        }
        if generators.len() % 3 != 0 {
            panic!("I only work for dim 3 matrices.")
        }
        let type_zero_end = generators.len() / 3;
        let type_one_end = 2 * type_zero_end;
        Self {
            generators,
            type_zero_end,
            type_one_end,
            lookup: lookup.clone(),
        }
    }

    pub fn new(dim: usize, lookup: &Arc<GaloisField>) -> Self {
        let p = lookup.field_mod;
        let id = GaloisMatrix::id(dim);
        let mut type_to_gens = HashMap::new();

        for type_ix in 0..dim {
            let ret_type_i: &mut HashSet<GaloisMatrix> = type_to_gens.entry(type_ix).or_default();
            ret_type_i.insert(id.clone());
            for row_ix in 0..dim {
                for col_ix in 0..dim {
                    if row_ix == col_ix {
                        continue;
                    }
                    let deg = compute_deg(dim, type_ix, row_ix, col_ix);
                    if deg > 0 {
                        let matrices = type_to_gens
                            .get_mut(&type_ix)
                            .expect("couldn't get matrices");
                        let mut new_matrices = HashSet::new();
                        for matrix in matrices.iter() {
                            for c in 1..p {
                                let mut new_matrix = matrix.clone();
                                let monomial = FFPolynomial::monomial((c, p).into(), deg);
                                new_matrix.set_entry(row_ix, col_ix, monomial, lookup.clone());
                                new_matrices.insert(new_matrix);
                            }
                        }
                        for matrix in new_matrices {
                            matrices.insert(matrix);
                        }
                    }
                }
            }
        }
        let mut gens = Vec::new();
        let zero_gens = type_to_gens.remove(&0).expect("No type zero gens");
        let one_gens = type_to_gens.remove(&1).expect("No type zero gens");
        let two_gens = type_to_gens.remove(&2).expect("No type zero gens");
        for m in zero_gens.into_iter() {
            gens.push(m);
        }
        let type_zero_end = gens.len() - 1;
        for m in one_gens.into_iter() {
            gens.push(m);
        }
        let type_one_end = gens.len() - 1;
        for m in two_gens.into_iter() {
            gens.push(m);
        }
        Self {
            generators: gens,
            type_zero_end,
            type_one_end,
            lookup: lookup.clone(),
        }
    }

    /// Returns the sorted coset of the given representative (member of the coset) and the type of the coset.
    pub fn get_coset(&self, coset_rep: &GaloisMatrix, type_ix: usize) -> Vec<GaloisMatrix> {
        if type_ix == 0 {
            let subs = &self.generators[0..=self.type_zero_end];
            let mut coset: Vec<GaloisMatrix> = subs
                .par_iter()
                .map(|k| coset_rep.mul(k, self.lookup.clone()))
                .collect();
            coset.sort();
            coset
        } else if type_ix == 1 {
            let subs = &self.generators[self.type_zero_end + 1..=self.type_one_end];
            let mut coset: Vec<GaloisMatrix> = subs
                .par_iter()
                .map(|k| coset_rep.mul(k, self.lookup.clone()))
                .collect();
            coset.sort();
            coset
        } else if type_ix == 2 {
            let subs = &self.generators[self.type_one_end + 1..self.generators.len()];
            let mut coset: Vec<GaloisMatrix> = subs
                .par_iter()
                .map(|k| coset_rep.mul(k, self.lookup.clone()))
                .collect();
            coset.sort();
            coset
        } else {
            panic!("This type is not implemented yet.")
        }
    }

    /// Returns the sorted cosets of a given matrix of all 3 types.
    /// Returns all 3 at once as this can be more easily parallelized than
    /// getting each coset one at a time.
    pub fn get_all_cosets(
        &self,
        rep: &GaloisMatrix,
    ) -> (Vec<GaloisMatrix>, Vec<GaloisMatrix>, Vec<GaloisMatrix>) {
        let v: Vec<GaloisMatrix> = self
            .generators
            .par_iter()
            .map(|k| rep.mul(k, self.lookup.clone()))
            .collect();
        let (mut v0, mut v1, mut v2) = (
            v[..=self.type_zero_end].to_vec(),
            v[self.type_zero_end + 1..=self.type_one_end].to_vec(),
            v[self.type_one_end + 1..].to_vec(),
        );
        v0.sort();
        v1.sort();
        v2.sort();
        (v0, v1, v2)
    }

    pub fn get_canonical_rep(&self, rep: &GaloisMatrix, type_ix: usize) -> GaloisMatrix {
        let coset = self.get_coset(rep, type_ix);
        coset[0].clone()
    }

    pub fn get_coset_reps(&self, mat: &GaloisMatrix) -> (CosetRep, CosetRep, CosetRep) {
        let (coset0, coset1, coset2) = self.get_all_cosets(mat);
        let c0 = CosetRep {
            rep: coset0[0].clone(),
            type_ix: 0,
        };
        let c1 = CosetRep {
            rep: coset1[0].clone(),
            type_ix: 1,
        };
        let c2 = CosetRep {
            rep: coset2[0].clone(),
            type_ix: 2,
        };
        (c0, c1, c2)
    }
}

mod tests {
    use std::{str::FromStr, sync::Arc};

    use crate::math::{
        coset_complex_subgroups::{compute_deg, k_type_subgroup},
        galois_field::GaloisField,
        polynomial::FFPolynomial,
    };
    use crate::matrices::ffmatrix::FFMatrix;

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
        let poly = FFPolynomial::from_str("1*x^2 + 2 * x^1 + 2 * x^0 % 3").unwrap();
        let lookup = Arc::new(GaloisField::new(poly));
        let subs = k_type_subgroup(0, lookup.clone());
        println!("subs len: {:}", subs.len());
        for sub in subs {
            println!("{:}", sub.pretty_print(lookup.clone()));
        }
    }
}
