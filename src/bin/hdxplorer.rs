use core::panic;
use std::{collections::HashMap, env, path::PathBuf, str::FromStr, time::Instant};

use clap::*;
use hdx_codec::{
    hdx_code::HDXCodeConfig,
    iterative_rank_estimator_old::{IterativeRankEstimator, RankEstimatorConfig},
    math::{finite_field::FFRep, iterative_bfs_new::GroupBFS, polynomial::FFPolynomial},
};
use mhgl::{ConGraph, HyperGraph};

use simple_logger::SimpleLogger;

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

        /// Whether to enable caching or not.
        #[arg(short, long, value_name = "CACHE")]
        cache: bool,

        /// Upper limit how many BFS steps the walker will traverse of the graph. If this is `None` then `usize::MAX` will be used.
        #[arg(short, long, value_name = "MAX_BFS_STEPS")]
        max_bfs_steps: Option<usize>,
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
    IterativeRankUpperBound {
        #[arg(short, long, value_name = "QUOTIENT")]
        quotient: String,

        #[arg(long, value_name = "DIM")]
        dim: usize,

        #[arg(long, value_name = "RS_DEG")]
        reed_solomon_degree: usize,

        #[arg(short, long, value_name = "FILENAME")]
        filename: Option<String>,
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
            cache,
            max_bfs_steps,
        } => {
            let q = FFPolynomial::from_str(&quotient).expect("Could not parse quotient argument.");
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
            let mut bfs = GroupBFS::new(directory, String::from(filename), &q, cache);
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
            let q = FFPolynomial::from_str(&quotient).expect("Could not parse quotient argument.");
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
            let mut bfs = GroupBFS::new(directory, String::from(filename), &q, false);
            bfs.bfs(max_bfs_steps.unwrap_or(usize::MAX));
            let hg = bfs.hgraph();
        }
        Cli::IterativeRankUpperBound {
            quotient,
            dim,
            reed_solomon_degree,
            filename,
        } => {
            let start = Instant::now();
            let conf =
                RankEstimatorConfig::new(quotient, dim, reed_solomon_degree, filename.unwrap());
            let mut iterator = IterativeRankEstimator::new(conf);
            let rel_rate_upper_bound = iterator.compute_rate();
            let elapsed = start.elapsed().as_secs_f64();
            // log::trace!("Estimated rate upper bound: {:}", rel_rate_upper_bound);
            log::trace!("Took {:} seconds", elapsed);
        }
    }
}
