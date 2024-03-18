use core::{num, panic};
use std::collections::{HashMap, HashSet};
use std::env;
use std::ffi::OsStr;
use std::fs::File;
use std::io::{BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::rc::Rc;
use std::str::FromStr;

use bitvec::view;
use mhgl::HGraph;
use serde::{Deserialize, Serialize};
use uuid::Uuid;

use crate::code::Code;
use crate::math::ffmatrix::FFMatrix;
use crate::math::finite_field::{self, FFRep, FiniteField as FF, FiniteFieldExt};
use crate::math::polynomial::FiniteFieldPolynomial;
use crate::reed_solomon::ReedSolomon;
use crate::tanner_code::get_generator_from_parity_check;

pub const HDX_CONFIG_FILENAME: &str = "hdx_codec_config.json";

/// Everything needed to compute the `HGraph` required by the HDXCode
/// struct
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HDXCodeConfig {
    pub field_mod: u32,
    pub quotient_poly: FiniteFieldPolynomial,
    pub dim: usize,
    pub reed_solomon_degree: usize,
    pub base_dir: String,
}

impl HDXCodeConfig {
    /// Attempts to read the current working directory for a file "hdx_codec_config.json".
    pub fn new(cur_dir: PathBuf) -> Option<Self> {
        let mut config_path = cur_dir.canonicalize().expect("Cannot find canonical path.");
        let canonical_dir = config_path.clone();
        config_path.push(HDX_CONFIG_FILENAME);
        println!("Checking {:} for a config file.", config_path.display());
        let mut s = String::new();
        if let Ok(mut file) = File::open(config_path) {
            let res = file.read_to_string(&mut s);
            if res.is_err() {
                return None;
            }
            let deserialized = serde_json::from_str::<HDXCodeConfig>(&s);
            if deserialized.is_err() {
                println!(
                    "Received this error while deserializing a HDXCodeConfig: {:?}",
                    deserialized
                );
                None
            } else {
                let mut conf = deserialized.unwrap();
                conf.base_dir = canonical_dir.to_string_lossy().to_string();
                Some(conf)
            }
        } else {
            None
        }
    }

    pub fn save_to_disk(&self) {
        let s = serde_json::to_string(&self).expect("could not serialize config.");
        let mut filepath =
            PathBuf::from_str(&self.base_dir).expect("Cannot make pathbuf to save file.");
        filepath.push(HDX_CONFIG_FILENAME);
        let mut file = File::create(filepath).expect("Cannot create file for writing.");
        file.write_all(s.as_bytes()).expect("Could not write");
    }

    /// Will generate the required HGraph from scratch if worst
    /// comes to worst.
    pub fn get_hgraph_at_all_costs(&self) -> HGraph {
        todo!()
    }
}

pub struct NewHDXCode {
    hg: HGraph,
    view_size_to_code: HashMap<usize, ReedSolomon>,
    lines: Vec<Uuid>,
    /// Maps a triangle Uuid to it's index in the message
    triangles: HashMap<Uuid, usize>,
    field_mod: FFRep,
    quotient: FiniteFieldPolynomial,
}

impl NewHDXCode {
    pub fn new(hg_file: &Path, field_mod: FFRep, quotient: &FiniteFieldPolynomial) -> Self {
        if let Some(hg) = HGraph::from_file(hg_file) {
            // now need to construct local codes for each edge. These are isomorphic
            // to a ReedSolomon over the prime base and the number of triangles that
            // each line can see.
            let lines = hg.edges_of_size(2);
            let triangle_vec = hg.edges_of_size(3);
            let mut counter = 0;
            let triangles = triangle_vec
                .into_iter()
                .map(|id| {
                    let out = (id, counter);
                    counter += 1;
                    out
                })
                .collect();
            let mut view_size_to_code = HashMap::new();
            for line in lines.iter() {
                let star = hg.star_id(line);
                if view_size_to_code.contains_key(&star.len()) == false && star.len() == field_mod as usize {
                    // todo: I don't know if the quotient degree is the right degree
                    // to be using here.
                    let rs = ReedSolomon::new_with_parity_check_input(star.len(), quotient.degree() as usize, field_mod);
                    rs.print();
                    println!(
                        "Created ReedSolomon::new({:}, {:}) yields input message length: {:}",
                        star.len(),
                        quotient.degree(), 
                        rs.encoded_len()
                    );
                    let pcm = rs.parity_check_matrix();
                    println!("Parity check matrix:\n{:}", pcm);
                    println!("nrows, ncols = {:}, {:}", pcm.n_rows, pcm.n_cols);
                    view_size_to_code.insert(star.len(), rs);
                }
            }
            return Self {
                hg,
                view_size_to_code,
                lines,
                triangles,
                field_mod,
                quotient: quotient.clone(),
            };
        } else {
            panic!("No hypergraph file found, idk what to do.")
        }
    }

    // TODO: return a slice that lives as long as the provided message instead
    // cloning the values.
    pub fn get_local_view(&self, local_check: &Uuid, message: &Vec<FF>) -> Vec<FF> {
        let star = self.hg.star_id(local_check);
        star.into_iter()
            .map(|id| {
                let ix = self.triangles.get(&id).expect("Could not find triangle.");
                message.get(*ix).unwrap().clone()
            })
            .collect()
    }

    /// Returns the Uuids of the local codes that report an error, regardless of
    /// the distance from a local codeword.
    pub fn get_failing_checks(&self, message: &Vec<FF>) -> HashSet<Uuid> {
        self.lines
            .iter()
            .filter(|line| {
                let local_view = self.get_local_view(line, message);
                let local_view_formatted: Vec<FF> = local_view
                    .into_iter()
                    .map(|ff| FF::new(ff.0, self.field_mod))
                    .collect();
                let code = self
                    .view_size_to_code
                    .get(&local_view_formatted.len())
                    .expect("Could not find local code.");
                code.code_check(&local_view_formatted)
            })
            .cloned()
            .collect()
    }
    fn get_lines_of_triangle(&self, triangle: &Uuid) -> (Uuid, Uuid, Uuid) {
        let t_nodes = self.hg.query_edge_id(triangle).expect("No edge?");
        let l1 = vec![t_nodes[0], t_nodes[1]];
        let l2 = vec![t_nodes[0], t_nodes[2]];
        let l3 = vec![t_nodes[1], t_nodes[2]];
        let l1_id = self.hg.get_edge_id(&l1).unwrap();
        let l2_id = self.hg.get_edge_id(&l2).unwrap();
        let l3_id = self.hg.get_edge_id(&l3).unwrap();
        (l1_id, l2_id, l3_id)
    }

    fn line_parity_check(&self, line: &Uuid, message: &Vec<FF>) -> Vec<FF> {
        let local_view = self.get_local_view(line, message);
        if let Some(code) = self.view_size_to_code.get(&local_view.len()) {
            println!("local_message: {:?}", local_view);
            let ret = code.parity_check(&local_view);
            println!("local check: {:?}", ret);
            ret
        } else {
            panic!("Could not find the code for the provided line.")
        }
    }

    fn check_line(&self, line: &Uuid, message: &Vec<FF>) -> bool {
        let local_parity = self.line_parity_check(line, message);
        let mut is_zero = true;
        for symbol in local_parity.iter() {
            if symbol.0 != 0 {
                is_zero = false;
                break;
            }
        }
        is_zero
    }

    pub fn check_lines(&self, message: &Vec<FF>) -> HashMap<Uuid, bool> {
        let parity_checks = self.get_line_parity_checks(message);
        parity_checks
            .into_iter()
            .map(|(id, pc)| {
                let mut is_zero = true;
                for symbol in pc.iter() {
                    if symbol.0 != 0 {
                        is_zero = false;
                        break;
                    }
                }
                (id, is_zero)
            })
            .collect()
    }

    fn number_of_failing_checks(&self, triangle: &Uuid, message: &Vec<FF>) -> usize {
        let (l1, l2, l3) = self.get_lines_of_triangle(triangle);
        let mut tot = 0;
        if self.check_line(&l1, message) {
            tot += 1;
        }
        if self.check_line(&l2, message) {
            tot += 1;
        }
        if self.check_line(&l3, message) {
            tot += 1;
        }
        tot
    }

    fn get_line_parity_checks(&self, message: &Vec<FF>) -> HashMap<Uuid, Vec<FF>> {
        self.lines
            .iter()
            .filter(|line_id| {
                let local_view = self.get_local_view(*line_id, message);
                local_view.len() == self.field_mod as usize
            })
            .map(|id| {
                let pc = self.line_parity_check(id, message);
                (id.clone(), pc)
            })
            .collect()
    }

    pub fn parity_check(&self, message: &Vec<FF>) -> Vec<FF> {
        let mut s = String::new();
        for m in message.iter() {
            s.push_str(&m.0.to_string())
        }
        println!("[parity_check] input: {:}", s);

        let line_to_parity_check = self.get_line_parity_checks(message);
        let mut parity_out: Vec<(Uuid, Vec<FF>)> = line_to_parity_check.into_iter().collect();
        parity_out.sort_by(|(id_1, _), (id_2, _)| id_1.cmp(id_2));
        let mut ret = Vec::with_capacity(parity_out.len());
        for (_, mut pc) in parity_out.drain(..) {
            ret.append(&mut pc);
        }
        ret
    }

    pub fn parity_check_matrix(&self) -> FFMatrix {
        let mut zero_vec = finite_field::zero_vec(self.triangles.len(), self.field_mod);
        let mut parity_checks: Vec<Vec<FF>> = Vec::new();
        for ix in 0..zero_vec.len() {
            let e = zero_vec.get_mut(ix).unwrap();
            e.0 = 1;
            let parity_check = self.parity_check(&zero_vec);
            let e = zero_vec.get_mut(ix).unwrap();
            e.0 = 0;
            parity_checks.push(parity_check);
        }
        let output_dim = parity_checks[0].len();
        let mut row_major_entries = Vec::with_capacity(output_dim * parity_checks.len());
        for row_ix in 0..output_dim {
            for col_ix in 0..parity_checks.len() {
                row_major_entries.push(*parity_checks[col_ix].get(row_ix).unwrap())
            }
        }
        FFMatrix::new(row_major_entries, output_dim, parity_checks.len())
    }

    /// Returns the triangles that have edges that report a failure. Triangles are sorted in a decreasing
    /// order by number of edges that report failure.
    pub fn get_failing_triangles(&self, message: Vec<FF>) {
        if message.len() != self.triangles.len() {
            println!("Number of message symbols must match number of triangles.");
            return;
        }

        // Map each message symbol to a triangle.
    }

    pub fn get_encoder_from_parity_check(&self) -> FFMatrix {
        let pc = self.parity_check_matrix();
        println!("Computed parity check matrix:\n{:}", pc);
        get_generator_from_parity_check(&pc)
    }
}

mod tests {
    use std::path::{Path, PathBuf};

    use crate::math::{iterative_bfs::GroupBFS, polynomial::FiniteFieldPolynomial};

    use super::NewHDXCode;

    #[test]
    fn test_whole_shebang() {
        let dir = PathBuf::from("/Users/matt/repos/qec/tmp");
        let p = 3_u32;
        let primitive_coeffs = [(2, (1, p).into()), (1, (2, p).into()), (0, (2, p).into())];
        let q = FiniteFieldPolynomial::from(&primitive_coeffs[..]);
        let mut bfs = GroupBFS::new(&dir, &q);
        bfs.bfs((2 as usize).pow(6));
        let hdx_code = NewHDXCode::new(&bfs.get_hgraph_file_path(), p, &q);
        let gen_mat = hdx_code.get_encoder_from_parity_check();
        println!("Generator:\n{:}", gen_mat);
    }
}
