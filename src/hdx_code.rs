use core::num;
use std::collections::HashMap;
use std::env;
use std::ffi::OsStr;
use std::fs::File;
use std::io::{BufReader, Read, Write};
use std::path::PathBuf;
use std::str::FromStr;

use mhgl::HGraph;
use serde::{Deserialize, Serialize};
use uuid::Uuid;

use crate::math::coset_complex::{self, CosetComplex};
use crate::math::ffmatrix::FFMatrix;
use crate::math::finite_field::{self, FiniteField, FiniteFieldRep};
use crate::math::polynomial::FiniteFieldPolynomial;
use crate::reed_solomon::ReedSolomon;

pub const HDX_CONFIG_FILENAME: &str = "hdx_codec_config.json";

/// Everything needed to compute the `HGraph` required by the HDXCode
/// struct
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HDXCodeConfig {
    pub field_mod: FiniteFieldRep,
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

    /// Checks to see if an hgraph was already computed or not
    /// in the `base_dir`, this checks for intermediate artifacts
    /// that should also be present in order for the hgraph to be
    /// properly computed.
    // pub fn get_hgraph_from_disk(&self) -> Option<HGraph> {
    //     let mut hgraph_path = self.base_dir.clone();
    //     hgraph_path.push(coset_complex::HGRAPH_FILENAME);
    //     let file = File::open(hgraph_path).expect("Cannot open hgraph file.");
    //     let reader = BufReader::new(file);
    //     if let Ok(hg) = serde_json::from_reader(reader) {
    //         Some(hg)
    //     } else {
    //         println!("Could not deserialize hgraph.");
    //         None
    //     }
    // }

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
        let mut cc = CosetComplex::new(self.base_dir.clone(), self.dim, &self.quotient_poly);
        cc.compute_hgraph_at_all_costs()
    }
}

/// Data necessary to encode and decode the classical codes from [New Codes on High Dimensional Expanders](https://arxiv.org/abs/2308.15563)
pub struct HDXCode {
    /// Currently assume the same local code for each edge. In the future could look more like `edge_id_to_code: HashMap<Uuid, ReedSolomon>`
    local_code: ReedSolomon,

    // Unclear if we need the whole CosetComplex structure, but we might be able
    // to get by without it.
    hgraph: HGraph,

    /// Stores the encoded message on the triangles of the HDX. So we need a map from Uuid of the triangle to the message symbol.
    // encoded_message: HashMap<Uuid, FiniteField>,
    field_mod: FiniteFieldRep,
}

impl HDXCode {
    /// Will attempt to read in the necessary data from the
    pub fn new(config: HDXCodeConfig) -> Self {
        let hg = config.get_hgraph_at_all_costs();
        let rs = ReedSolomon::new(config.field_mod, config.reed_solomon_degree);
        Self {
            local_code: rs,
            hgraph: hg,
            field_mod: config.field_mod,
        }
    }

    pub fn generate_random_encoder(&self, num_codewords: usize) -> FFMatrix {
        let mut generator_transpose: Vec<FiniteField> = Vec::new();
        for k in 0..num_codewords {
            // let random_message = finite_field::random_message(self.coset_complex.num_triangles(), self.field_mod);
        }
        FFMatrix::id(1, self.field_mod)
    }

    /// Returns the edges local view of the encoded message.
    fn get_sorted_local_view(
        &self,
        edge: (u32, u32),
        message: &HashMap<Uuid, FiniteField>,
    ) -> Vec<FiniteField> {
        let mut triangles = self.hgraph.get_containing_edges(&[edge.0, edge.1]);
        triangles.sort();
        triangles
            .into_iter()
            .map(|id| {
                message
                    .get(&id)
                    .expect("Given triangle has not been assigned a message symbol.")
            })
            .cloned()
            .collect()
    }

    fn get_failing_checks(&self, message: &HashMap<Uuid, FiniteField>) {
        let edge_ids = self.hgraph.edges_of_size(2);
        
    }

    pub fn is_message_in_code(&self, message: &HashMap<Uuid, FiniteField>) -> bool {
        let edge_ids = self.hgraph.edges_of_size(2);
        let mut pass_all_checks = true;
        for edge_id in edge_ids {
            let edge_nodes = self
                .hgraph
                .query_edge_id(&edge_id)
                .expect("That edge better be in there!");
            if edge_nodes.len() != 2 {
                panic!("edge is not of size 2");
            }
            let mut triangles = self
                .hgraph
                .get_containing_edges(&[edge_nodes[0], edge_nodes[1]]);
            triangles.sort();
            let local_view: Vec<FiniteField> = triangles
                .into_iter()
                .map(|id| {
                    message
                        .get(&id)
                        .expect("Triangle is not included in message.")
                })
                .cloned()
                .collect();
            if self.local_code.is_message_in_code(&local_view) == false {
                pass_all_checks = false;
                break;
            }
        }
        pass_all_checks
    }
}
