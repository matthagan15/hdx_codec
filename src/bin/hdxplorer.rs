use core::panic;
use std::{
    collections::HashMap,
    env,
    fs::File,
    io::{Read, Write},
    path::{Path, PathBuf},
    str::FromStr,
    time::Instant,
};

use clap::*;
use hdx_codec::{
    code::Code,
    hdx_code::HDXCodeConfig,
    math::{
        finite_field::FFRep,
        iterative_bfs_new::GroupBFS,
        lps::{self, compute_lps_graph},
        polynomial::FiniteFieldPolynomial,
    },
    reed_solomon::ReedSolomon,
    tanner_code::{ParityCode, TannerCode},
};
use mhgl::{ConGraph, HGraph, HyperGraph};

use log::{debug, error, info, trace, warn};
use simple_logger::SimpleLogger;

fn get_config_from_current_working_dir() -> Option<HDXCodeConfig> {
    // First check if the current directory contains a config.
    let cur_dir = env::current_dir().expect("Cannot get current working directory.");
    let mut cur_dir_config_path = cur_dir.clone();
    cur_dir_config_path.push(hdx_codec::hdx_code::HDX_CONFIG_FILENAME);
    println!(
        "Checking for a config here: {:}",
        cur_dir_config_path.display()
    );
    if let Ok(b) = cur_dir_config_path.try_exists() {
        if b {
            println!("New?");
            // try to read it
            let ret = HDXCodeConfig::new(cur_dir.clone());
            if ret.is_none() {
                println!("Did not find one.");
            } else {
                println!("Found a config!");
            }
            ret
        } else {
            None
        }
    } else {
        None
    }
}

fn get_hdx_config_from_user() -> HDXCodeConfig {
    println!("No config found, create one now.");
    let mut user_input = String::new();
    println!("Enter prime for base alphabet: ");
    std::io::stdin()
        .read_line(&mut user_input)
        .expect("Could not read user input.");
    let p = user_input.trim().parse::<u32>().expect("Could not parse.");
    println!("Now enter a polynomial to quotient matrix entries by: ");
    user_input.clear();
    std::io::stdin()
        .read_line(&mut user_input)
        .expect("Could not read user input.");
    let q_string = user_input.trim().to_string();
    let q = FiniteFieldPolynomial::from_str(&q_string).expect("Could not parse polynomial.");
    if q.field_mod != p {
        println!("Improper field_mod entered.");
        panic!("Idk what to do.")
    }
    user_input.clear();
    println!("Enter a dimension to use for the coset complex (enter 3): ");
    std::io::stdin()
        .read_line(&mut user_input)
        .expect("Could not read user input.");
    let dim = user_input
        .trim()
        .parse::<usize>()
        .expect("Could not parse.");

    user_input.clear();
    println!("Enter a max degree (non-inclusive) to use for the local Reed-Solomon code: ");
    std::io::stdin()
        .read_line(&mut user_input)
        .expect("Could not read user input.");
    let rs_degree = user_input
        .trim()
        .parse::<usize>()
        .expect("Could not parse.");
    let cur_dir = env::current_dir().expect("Cannot get current working directory.");
    HDXCodeConfig {
        field_mod: p,
        quotient_poly: q,
        dim,
        reed_solomon_degree: rs_degree,
        base_dir: cur_dir.to_string_lossy().to_string(),
    }
}

fn degree_stats(hg: &ConGraph) {
    println!("Checking degrees.");
    let nodes = hg.nodes();
    let mut node_to_node_stats = HashMap::new();
    let mut node_to_edges_stats = HashMap::new();
    let edges = hg.edges_of_size(2);
    let mut edge_stats = HashMap::new();
    println!("nodes.");
    for node in nodes {
        let link = hg.link_of_nodes(&[node]);
        let mut num_nodes = 0;
        let mut num_edges = 0;
        for (id, set) in link.into_iter() {
            if set.len() == 1 {
                num_nodes += 1;
                continue;
            }
            if set.len() == 2 {
                num_edges += 1;
            }
        }
        let e: &mut usize = node_to_node_stats.entry(num_nodes).or_default();
        *e += 1;
        let ee: &mut usize = node_to_edges_stats.entry(num_edges).or_default();
        *ee += 1;
    }
    println!("Edges.");
    for edge in edges {
        if let Some(edge_set) = hg.query_edge(&edge) {
            let mut num_nodes = 0;
            for (id, set) in hg.link_of_nodes(&edge_set[..]) {
                if set.len() == 1 {
                    num_nodes += 1;
                }
            }
            let e: &mut usize = edge_stats.entry(num_nodes).or_default();
            *e += 1;
        }
    }
    println!("Node statistics. The node -> node degree stats are:");
    println!("{:?}", node_to_node_stats);
    println!("node -> edge degree stats are:");
    println!("{:?}", node_to_edges_stats);
    println!("And the edge -> node degree stats are:");
    println!("{:?}", edge_stats);
}

#[derive(Subcommand)]
pub enum HdxCommands {}

#[derive(Parser)]
#[command(version, about, long_about = None)]
enum Cli {
    /// Builds the Coset Complex given by Dinur, Liu, and Zhang, which is a
    /// variant of those given by Kaufman and Oppenheim.
    BuildCosetComplex {
        #[arg(short, long, value_name = "QUOTIENT")]
        quotient: String,
        #[arg(long, value_name = "DIM")]
        dim: usize,
        /// Where to save the file. If you do not specify
        /// then it will be saved to the current working directory
        /// with the  default filename
        /// `p_<QUOTIENT.field_mod>_<QUOTIENT.degree>_<DIM>.hg' and
        /// then whatever extensions we use `.hg`, `.dot`, `.cache`, `.conf`
        #[arg(short, long, value_name = "FILENAME")]
        filename: Option<String>,

        /// Upper limit how many BFS steps the walker will traverse of the graph. If this is `None` then `usize::MAX` will be used.
        #[arg(short, long, value_name = "MAX_BFS_STEPS")]
        max_bfs_steps: Option<usize>,
    },
    /// Builds LPS(p, q) graph.
    BuildLPS {
        p: FFRep,
        q: FFRep,
        filename: Option<String>,
    },

    View {
        /// Currently you write out the entire damn thing.
        #[arg(short, long, value_name = "FILENAME")]
        filename: String,
    },
    CosetCodeRank {
        #[arg(short, long, value_name = "QUOTIENT")]
        quotient: String,
        #[arg(long, value_name = "DIM")]
        dim: usize,
        /// Where to save the file. If you do not specify
        /// then it will be saved to the current working directory
        /// with the  default filename
        /// `p_<QUOTIENT.field_mod>_<QUOTIENT.degree>_<DIM>.hg' and
        /// then whatever extensions we use `.hg`, `.dot`, `.cache`, `.conf`
        #[arg(short, long, value_name = "FILENAME")]
        filename: Option<String>,

        /// Upper limit how many BFS steps the walker will traverse of the graph. If this is `None` then `usize::MAX` will be used.
        #[arg(short, long, value_name = "MAX_BFS_STEPS")]
        max_bfs_steps: Option<usize>,

        #[arg(long, value_name = "RS_FIELD_MOD")]
        rs_field_mod: FFRep,

        #[arg(long, value_name = "RS_DEGREE")]
        rs_degree: usize,
    },
}

fn main() {
    let logger = SimpleLogger::new().init().unwrap();
    log::debug!("Testing logging.");
    let cli = Cli::parse();
    match cli {
        Cli::BuildCosetComplex {
            quotient,
            dim,
            filename,
            max_bfs_steps,
        } => {
            let q = FiniteFieldPolynomial::from_str(&quotient)
                .expect("Could not parse quotient argument.");
            let path_buf = match filename {
                Some(path_string) => PathBuf::from(path_string),
                None => {
                    let mut dir = env::current_dir().expect(
                        "Expected user to be operating in a shell with a valid current directory.",
                    );
                    dir.push(format!(
                        "hdx_coset_dim_{:}_deg_{:}_mod_{:}.hg",
                        dim,
                        q.degree(),
                        q.field_mod
                    ));
                    dir
                }
            };
            let directory = path_buf.parent().expect("No directory stem.");
            let filename = path_buf
                .file_stem()
                .expect("Just made sure filename existed")
                .to_str()
                .unwrap();
            let mut bfs = GroupBFS::new(directory, String::from(filename), &q);
            bfs.bfs(max_bfs_steps.unwrap_or(usize::MAX));
        }
        Cli::View { filename } => {
            let mut pathbuf = PathBuf::from(&filename);
            let hg = ConGraph::from_file(&pathbuf).expect("Could not find hgraph.");
            println!("hg: {:}", hg);
            let good_lines: Vec<_> = hg
                .edges_of_size(2)
                .into_iter()
                .filter(|line| hg.maximal_edges(line).len() == 3)
                .collect();

            degree_stats(&hg);
        }
        Cli::CosetCodeRank {
            quotient,
            dim,
            filename,
            max_bfs_steps,
            rs_field_mod,
            rs_degree,
        } => {
            let q = FiniteFieldPolynomial::from_str(&quotient)
                .expect("Could not parse quotient argument.");
            let path_buf = match filename {
                Some(path_string) => PathBuf::from(path_string),
                None => {
                    let mut dir = env::current_dir().expect(
                        "Expected user to be operating in a shell with a valid current directory.",
                    );
                    dir.push(format!(
                        "hdx_coset_dim_{:}_deg_{:}_mod_{:}.hg",
                        dim,
                        q.degree(),
                        q.field_mod
                    ));
                    dir
                }
            };
            let directory = path_buf.parent().expect("No directory stem.");
            let filename = path_buf
                .file_stem()
                .expect("Just made sure filename existed")
                .to_str()
                .unwrap();
            let mut bfs = GroupBFS::new(directory, String::from(filename), &q);
            bfs.bfs(max_bfs_steps.unwrap_or(usize::MAX));
            let hg = bfs.hgraph();
        }
        Cli::BuildLPS { p, q, filename } => todo!(),
    }
    // let q =
    //     FiniteFieldPolynomial::from_str(&hdx_builder.quotient).expect("Could not parse quotient");
    // let mut hdx_bfs = GroupBFS::new(&hdx_builder.directory, &q);
    // hdx_bfs.bfs((2 as usize).pow(1));
    // let hg_path = hdx_bfs.get_hgraph_file_path();
    // println!("Graph is here: {:?}", hg_path);
    // let start = Instant::now();
    // let hg = HGraph::from_file(hg_path.as_path()).expect("Could not open file");
    // println!(
    //     "Took this long to parse: {:}",
    //     start.elapsed().as_secs_f64()
    // );
    // degree_stats(&hg);
    // println!("{:}", hg);
}
