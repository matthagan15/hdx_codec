use std::{
    collections::{HashMap, HashSet},
    env,
    io::Write,
    path::PathBuf,
    str::FromStr,
    time::Instant,
};

use clap::*;
use hdx_codec::{
    math::{coset_complex_bfs::bfs, finite_field::FFRep, polynomial::FFPolynomial},
    matrices::sparse_ffmatrix::benchmark_rate,
    quantum::{boundary_down, boundary_up},
    rank_estimator_sparse::{split_hg_into_checks, RankConfig, RankEstimatorConfig},
    reed_solomon::ReedSolomon,
};
use mhgl::{EdgeSet, HGraph, HyperGraph};

use serde_json::json;
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
    Bench {
        #[arg(short)]
        dim: usize,

        #[arg(short)]
        samples: usize,
    },
    Quantum {
        #[arg(short, long, value_name = "QUOTIENT")]
        quotient_poly: String,

        #[arg(short, long, value_name = "DIM")]
        dim: usize,

        #[arg(short, long, value_name = "TRUNCATION")]
        truncation: usize,

        #[arg(short, long, value_name = "OUTPUT")]
        output: Option<PathBuf>,
    },
    /// Builds the Coset Complex given by Dinur, Liu, and Zhang, which is a
    /// variant of those given by Kaufman and Oppenheim.
    Build {
        #[arg(short, long, value_name = "QUOTIENT")]
        quotient: String,
        #[arg(long, value_name = "DIM")]
        dim: usize,
        /// Where to save the hgraph file. If you do not specify a name
        /// then it will be saved to the current working directory
        /// with the  default filename
        /// `coset_complex_<QUOTIENT.field_mod>_<QUOTIENT.degree>_<DIM>.hg'.
        /// Ex: For quotient = "1*x^2 + 2*x^1 + 2*x^0 % 3" with dim = 3, the filename
        /// would be `coset_complex_3_2_3.hg`.
        /// If a caching directory is passed in then the hgraph file will be saved in
        /// that corresponding directory.
        #[arg(short, long, value_name = "FILENAME")]
        filename: Option<String>,

        /// An optional caching directory. If a caching directory is specified
        /// the coset_complex output will be saved in this directory.
        #[arg(short, long, value_name = "CACHE")]
        cache: Option<PathBuf>,

        /// Upper limit how many BFS steps the walker will traverse of the graph. If this is `None` then `usize::MAX` will be used.
        #[arg(short = 't', long = "trunc", value_name = "TRUNCATION")]
        truncation: Option<usize>,
    },

    View {
        /// Currently you write out the entire damn thing.
        #[arg(short, long, value_name = "FILENAME")]
        filename: String,
    },

    Rank {
        #[arg(short, long, value_name = "QUOTIENT")]
        quotient_poly: String,
        #[arg(short, long, value_name = "DIM")]
        dim: usize,
        #[arg(short, long, value_name = "REED_SOLOMON_DEGREE")]
        rs_degree: usize,
        #[arg(short, long, value_name = "CACHE_DIR")]
        cache: PathBuf,
        #[arg(short, long, value_name = "TRUNCATION")]
        truncation: Option<usize>,

        /// (Optional) The number of steps between each cache.
        #[arg(long, value_name = "CACHE_RATE")]
        cache_rate: Option<usize>,

        /// (Optional) The number of steps between logging data about the rank computation.
        #[arg(long, value_name = "LOG_RATE")]
        log_rate: Option<usize>,

        /// (Optional) The number of threads to use for the matrix reduction. Defaults to the
        /// number of available cores.
        #[arg(short = 'j', long, value_name = "NUM_THREADS")]
        num_threads: Option<usize>,
    },
}

fn main() {
    let logger = SimpleLogger::new().init().unwrap();
    let cli = Cli::parse();
    match cli {
        Cli::Build {
            quotient,
            dim,
            filename,
            cache,
            truncation,
        } => {
            let q = FFPolynomial::from_str(&quotient).expect("Could not parse quotient argument.");
            let dir = match cache {
                Some(cache_dir) => {
                    if cache_dir.is_dir() {
                        cache_dir
                    } else {
                        println!("Cache parameter received, but the path is not a directory.");
                        panic!()
                    }
                }
                None => std::env::current_dir().expect(
                    "Do not have access to current working directory and no cache directory found.",
                ),
            };
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
            bfs(q, dim, truncation, Some(dir), Some(5));
            // let mut bfs = GroupBFS::new(directory, String::from(filename), &q, cache);
            // bfs.bfs(max_bfs_steps.unwrap_or(usize::MAX));
        }
        Cli::View { filename } => {
            let mut pathbuf = PathBuf::from(&filename);
            let hg = HGraph::<u16, ()>::from_file(&pathbuf).expect("Could not find hgraph.");
            degree_stats(&hg);
            hgraph_client_loop(hg);
        }
        Cli::Bench { dim, samples } => {
            benchmark_rate(dim, samples);
        }
        Cli::Rank {
            quotient_poly,
            dim,
            rs_degree,
            cache,
            truncation,
            cache_rate,
            log_rate,
            num_threads,
        } => {
            if let Some(mut rank_config) = RankConfig::from_disk(cache.as_path()) {
                rank_config.run(
                    cache.as_path(),
                    std::thread::available_parallelism().unwrap().into(),
                );
            } else {
                let q = FFPolynomial::from_str(&quotient_poly)
                    .expect("Could not parse quotient polynomial.");
                let mut rank_config = RankConfig::new(
                    q,
                    dim,
                    rs_degree,
                    truncation.unwrap(),
                    truncation.unwrap() / 50,
                );
                rank_config.save_to_disk(cache.as_path());
                rank_config.run(
                    cache.as_path(),
                    num_threads.unwrap_or(std::thread::available_parallelism().unwrap().into()),
                );
            }
            // if cache.is_dir() == false {
            //     log::error!("Cache needs to be a directory!");
            //     panic!()
            // }
            // hdx_codec::rank_estimator_sparse::compute_rank_bounds(
            //     q,
            //     dim,
            //     rs_degree,
            //     cache,
            //     truncation,
            //     cache_rate,
            //     log_rate,
            //     num_threads,
            // );
        }
        Cli::Quantum {
            quotient_poly,
            dim,
            truncation,
            output,
        } => {
            let q = FFPolynomial::from_str(&quotient_poly[..]).unwrap();
            println!(
                "In Quantum.\nQuotient: {:?}, dim: {:}, truncation: {:}, output: {:?}",
                q, dim, truncation, output
            );
            let (mut hgraph, _) = bfs(q.clone(), 3, Some(truncation), None, None);
            // println!("hg:\n{:?}", hgraph);
            let local_rs = ReedSolomon::new(q.field_mod, 2);
            let qubits: HashSet<u64> = hgraph
                .edges_of_size(dim)
                .into_iter()
                .filter(|edge_id| hgraph.maximal_edges(edge_id).len() == q.field_mod as usize)
                .collect();
            let mut z_checks = HashSet::new();
            let mut x_checks = HashSet::new();
            for qubit in qubits.iter() {
                let maximal_edges = hgraph.maximal_edges(qubit);
                for z_check in maximal_edges {
                    z_checks.insert(z_check);
                }
                for node in hgraph.query_edge(qubit).unwrap() {
                    x_checks.insert(node);
                }
            }
            let total_num_nodes = hgraph.nodes().len();
            let total_num_qubits = hgraph.edges_of_size(2).len();
            let total_num_triangles = hgraph.edges_of_size(3).len();
            println!(
                "Number of nodes:     {:} / {:}",
                x_checks.len(),
                total_num_nodes
            );
            println!(
                "Number of edges:     {:} / {:}",
                qubits.len(),
                total_num_qubits
            );
            println!(
                "Number of triangles: {:} / {:}",
                z_checks.len(),
                total_num_triangles
            );
            let bd = boundary_down(&hgraph, qubits.iter().cloned().collect(), 2);
            let bu = boundary_up(&hgraph, qubits.iter().cloned().collect(), 2);
            println!("boundary up. Dims: {:} x {:}", bu.n_rows, bu.n_cols);
            // bu.dense_print();
            println!("boundary down. Dims: {:} x {:}", bd.n_rows, bd.n_cols);
            // bd.dense_print();

            if let Some(filename) = output {
                let boundary_up = bu.to_dense();
                let boundary_up_entries: Vec<u32> =
                    boundary_up.entries.iter().map(|ff| ff.0).collect();
                let boundary_down = bd.to_dense();
                let boundary_down_entries: Vec<u32> =
                    boundary_down.entries.iter().map(|ff| ff.0).collect();
                let data = json!({
                    "memory_layout": "row_major",
                    "boundary_up_n_rows": boundary_up.n_rows,
                    "boundary_up_n_cols": boundary_up.n_cols,
                    "boundary_down_n_rows": boundary_down.n_rows,
                    "boundary_down_n_cols": boundary_down.n_cols,
                    "boundary_up": boundary_up_entries,
                    "boundary_down": boundary_down_entries,
                });
                if let Ok(serialized_data) = serde_json::to_string(&data) {
                    std::fs::write(filename, serialized_data)
                        .expect("Could not write data to disk");
                }
            }
        }
    }
}
