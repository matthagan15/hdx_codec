use std::{collections::HashMap, env, io::Write, path::PathBuf, str::FromStr, time::Instant};

use clap::*;
use hdx_codec::{
    math::{coset_complex_bfs::GroupBFS, finite_field::FFRep, polynomial::FFPolynomial},
    rank_estimator_sparse::{IterativeRankEstimator, RankEstimatorConfig},
};
use mhgl::{EdgeSet, HGraph, HyperGraph};

use simple_logger::SimpleLogger;

pub enum HgClientCommand {
    Link,
    ContainingEdges,
    Maximal,
    Quit,
}

impl FromStr for HgClientCommand {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let trimmed = s.trim().to_ascii_lowercase();
        match &trimmed[..] {
            "l" | "link" => Ok(HgClientCommand::Link),
            "contain" => Ok(HgClientCommand::ContainingEdges),
            "max" | "maximal" => Ok(HgClientCommand::Maximal),
            "q" | "quit" => Ok(HgClientCommand::Quit),
            _ => {
                println!("Error parsing input, here it is post trim: {:?}", trimmed);
                println!("What you want to do.");
                Err(())
            }
        }
    }
}

pub fn hgraph_client_loop(hg: HGraph<u16, ()>) {
    let mut input_buf = String::new();
    let nodes = hg.nodes();
    let mut state = nodes[0];
    println!("HGraph Navigator:");
    println!("[q | quit] Quit this sub menu.");
    println!("[max | maximal] Maximal Edges");
    println!("[contain] Containing edges");
    println!("[l | link] Link of current state");
    loop {
        print!("hg {:} > ", state);
        std::io::stdout().flush().unwrap();
        input_buf.clear();
        std::io::stdin()
            .read_line(&mut input_buf)
            .expect("Could not read input.");
        let command = HgClientCommand::from_str(&input_buf[..]);
        if let Ok(c) = command {
            match c {
                HgClientCommand::Link => {
                    let link = hg.link_of_nodes([state]);
                    let mut s = String::new();
                    let mut link_size_to_edge_set: HashMap<usize, Vec<EdgeSet<u32>>> =
                        HashMap::new();
                    for (_, nodes) in link {
                        let e = EdgeSet::from(nodes);
                        link_size_to_edge_set.entry(e.len()).or_default();
                        let e_string = e.to_string();
                        s.push_str(&e_string[..]);
                        s.push_str(", ");
                    }
                    println!("{s}");
                }
                HgClientCommand::ContainingEdges => {
                    let containers = hg.containing_edges_of_nodes([state]);
                    let mut s = String::new();
                    for id in containers {
                        let nodes = hg.query_edge(&id).unwrap();
                        let e = EdgeSet::from(nodes);
                        s.push_str(&e.to_string()[..]);
                        s.push_str(", ");
                    }
                    println!("{s}");
                }
                HgClientCommand::Maximal => {
                    let maximal = hg.maximal_edges_of_nodes([state]);
                    let mut s = String::new();
                    for id in maximal {
                        let nodes = hg.query_edge(&id).unwrap();
                        let e = EdgeSet::from(nodes);
                        s.push_str(&e.to_string()[..]);
                        s.push_str(", ");
                    }
                    println!("{s}");
                }
                HgClientCommand::Quit => {
                    println!("Done.");
                    return;
                }
            }
        }
    }
}

fn degree_stats<N, E>(hg: &HGraph<N, E>) {
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
        let mut num_nodes = 0;
        for (id, set) in hg.link(&edge) {
            if set.len() == 1 {
                num_nodes += 1;
            }
        }
        let e: &mut usize = edge_stats.entry(num_nodes).or_default();
        *e += 1;
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
    IterativeRank {
        #[arg(short, long, value_name = "CONFIG")]
        config: PathBuf,
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
            let hg = HGraph::<u16, ()>::from_file(&pathbuf).expect("Could not find hgraph.");
            // let good_lines: Vec<_> = hg
            //     .edges_of_size(2)
            //     .into_iter()
            //     .filter(|line| hg.maximal_edges(line).len() == 3)
            //     .collect();

            degree_stats(&hg);
            hgraph_client_loop(hg);
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
        Cli::IterativeRank { config } => {
            let start = Instant::now();
            // let conf = RankEstimatorConfig::new(
            //     quotient,
            //     dim,
            //     reed_solomon_degree,
            //     cache_file.clone(),
            //     hgraph_filename,
            //     if num_threads.is_none() {
            //         std::thread::available_parallelism().unwrap().into()
            //     } else {
            //         num_threads.unwrap()
            //     },
            // );
            let conf = RankEstimatorConfig::from_disk(&config);
            let mut iterator = if let Some(s) = conf.cache_file.clone() {
                let p = PathBuf::from_str(&s[..]).unwrap();
                if let Some(i) = IterativeRankEstimator::load_from_cache(p) {
                    i
                } else {
                    log::trace!("Did not find usable cache file, starting from scratch.");
                    IterativeRankEstimator::new(conf)
                }
            } else {
                log::trace!("No cache file provided, starting from scratch.");
                IterativeRankEstimator::new(conf)
            };
            iterator.compute_rate();
            let elapsed = start.elapsed().as_secs_f64();
            // log::trace!("Estimated rate upper bound: {:}", rel_rate_upper_bound);
            log::trace!("Took {:} seconds", elapsed);
        }
    }
}
