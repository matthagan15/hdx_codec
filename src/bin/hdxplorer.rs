use std::{
    collections::{HashMap, HashSet},
    path::PathBuf,
    str::FromStr,
    time::Instant,
};

use clap::*;
use hdx_codec::{
    first_node::get_first_node_complete_star,
    math::{coset_complex_bfs::bfs, polynomial::FFPolynomial},
    matrices::sparse_ffmatrix::{MemoryLayout, SparseFFMatrix},
    quantum::{boundary_down, boundary_up},
    rank_estimator_sparse::RankConfig,
    rate_and_dist_estimator::RateAndDistConfig,
    reed_solomon::ReedSolomon,
};
use mhgl::{HGraph, HyperGraph};

use serde_json::json;
use simple_logger::SimpleLogger;

const NUMBER_OF_CHECKPOINTS: usize = 10;

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
        for (_id, set) in link.into_iter() {
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
        for (_id, set) in hg.link(&edge) {
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

    RateAndDistance {
        #[arg(short, long, value_name = "QUOTIENT")]
        quotient_poly: String,

        #[arg(short, long, value_name = "DIM")]
        dim: usize,

        #[arg(long, value_name = "DIRECTORY")]
        directory: PathBuf,

        #[arg(short, long, value_name = "TRUNCATION")]
        /// The number of steps to perform the BFS search for during
        /// the initial hypergraph discovery phase.
        truncation: Option<usize>,

        #[arg(short = 'l', long = "log_rate", value_name = "HYPERGRAPH_LOG_RATE")]
        /// The number of steps to report statistics during the BFS
        hgraph_log_rate: Option<usize>,

        #[arg(long, value_name = "NUM_CHECKPOINTS")]
        /// The number of checkpoints to report the error correcting code parameters
        /// n, k, and d.
        num_checkpoints: usize,
    },
    /// Builds the Coset Complex given by Dinur, Liu, and Zhang, which is a
    /// variant of those given by Kaufman and Oppenheim.
    Build {
        #[arg(short, long, value_name = "QUOTIENT")]
        quotient: String,
        #[arg(short, long, value_name = "DIM")]
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

        /// The number of steps between logging stats of the coset complex
        /// build process.
        #[arg(short, long = "log", value_name = "LOG_RATE")]
        log_rate: Option<usize>,

        /// Upper limit how many BFS steps the walker will traverse of the graph. If this is `None` then `usize::MAX` will be used.
        #[arg(short = 't', long = "trunc", value_name = "TRUNCATION")]
        truncation: Option<usize>,
    },

    View {
        /// Currently you write out the entire damn thing.
        #[arg(short, long, value_name = "FILENAME")]
        filename: String,
    },

    /// Compute the rank of a coset complex code. The size of the complex can be
    /// limited by setting `TRUNCATION` and the number of intermediate steps to report
    /// can be set with `CACHE_RATE`. Uses multiple threads to compute the rank of the resulting
    /// parity check matrices, the default is the number of cores on the device but it can be set
    /// with `NUM_THREADS`.
    Rank {
        /// The polynomial that is used to compute the base field for the matrix group
        /// that is used to compute the coset complex.
        #[arg(short, long, value_name = "QUOTIENT")]
        quotient_poly: String,

        /// Dimension of the matrix group $GL_(dim) (F_q)$
        #[arg(short, long, value_name = "DIM")]
        dim: usize,

        /// Where to store caches, configs, and other outputs.
        #[arg(short = 'o', long = "dir", value_name = "DIRECTORY")]
        output_directory: PathBuf,

        /// The amount of steps to perform on the BFS. Directly corresponds
        /// to the number of triangles in the computed simplicial complex. Defaults
        /// to the entire complex.
        #[arg(short, long, value_name = "TRUNCATION")]
        truncation: Option<usize>,

        /// (Optional) The number of threads to use for the matrix reduction. Defaults to the
        /// number of available cores.
        #[arg(short = 'j', long, value_name = "NUM_THREADS")]
        num_threads: Option<usize>,
    },
}

fn main() {
    let cli = Cli::parse();
    match cli {
        Cli::Build {
            quotient,
            dim,
            truncation,
            log_rate,
            filename,
        } => {
            if log_rate.is_some() {
                let _logger = SimpleLogger::new().init().unwrap();
            }
            let q = FFPolynomial::from_str(&quotient).expect("Could not parse quotient argument.");
            let hg = bfs(
                q,
                dim,
                truncation,
                filename.map(|s| PathBuf::from_str(&s[..]).ok()).flatten(),
                log_rate,
            );
        }
        Cli::View { filename } => {
            let pathbuf = PathBuf::from(&filename);
            let hg = HGraph::<u16, ()>::from_file(&pathbuf).expect("Could not find hgraph.");
            degree_stats(&hg);
        }
        Cli::Rank {
            quotient_poly,
            dim,
            output_directory,
            truncation,
            num_threads,
        } => {
            let _logger = SimpleLogger::new().init().unwrap();
            if let Some(mut rank_config) = RankConfig::from_disk(output_directory.as_path()) {
                rank_config.run(
                    output_directory.as_path(),
                    std::thread::available_parallelism().unwrap().into(),
                );
            } else {
                let q = FFPolynomial::from_str(&quotient_poly)
                    .expect("Could not parse quotient polynomial.");
                let mut rank_config = RankConfig::new(
                    q.clone(),
                    dim,
                    (q.field_mod - 1) as usize,
                    truncation.unwrap(),
                    truncation.unwrap() / NUMBER_OF_CHECKPOINTS,
                );
                rank_config.save_to_disk(output_directory.as_path());
                rank_config.run(
                    output_directory.as_path(),
                    num_threads.unwrap_or(std::thread::available_parallelism().unwrap().into()),
                );
            }
        }
        Cli::Quantum {
            quotient_poly,
            dim,
            truncation,
            output,
        } => {
            let start = Instant::now();
            let q = FFPolynomial::from_str(&quotient_poly[..]).unwrap();
            println!(
                "In Quantum.\nQuotient: {:?}, dim: {:}, truncation: {:}, output: {:?}",
                q, dim, truncation, output
            );
            let hgraph = bfs(q.clone(), 3, Some(truncation), None, None);
            let _local_rs: ReedSolomon = ReedSolomon::new(q.field_mod, 2);
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
            let mut z_parity_check = bu.clone();
            z_parity_check.transpose();

            let mut boundary_up_tester = bu.clone();
            boundary_up_tester.transpose();
            boundary_up_tester.swap_layout();

            let mut clean_boundary_down = bd.clone();
            clean_boundary_down.transpose();
            clean_boundary_down.swap_layout();

            let mut good_rows = Vec::new();
            for row_ix in 0..clean_boundary_down.n_rows {
                let row = clean_boundary_down.get_row(row_ix);
                if (&boundary_up_tester * &row).is_zero() {
                    good_rows.push(row_ix);
                }
            }
            let mut new_entries = Vec::new();
            for good_row in good_rows {
                let row = clean_boundary_down.get_row(good_row);
                let mut entries = row
                    .0
                    .into_iter()
                    .map(|(col_ix, entry)| (good_row, col_ix, entry))
                    .collect();
                new_entries.append(&mut entries);
            }
            let x_parity_check = SparseFFMatrix::new_with_entries(
                0,
                0,
                2,
                hdx_codec::matrices::sparse_ffmatrix::MemoryLayout::RowMajor,
                new_entries,
            );

            let mut x_rank_mat = x_parity_check.clone();
            let mut z_rank_mat = z_parity_check.clone();
            if x_rank_mat.memory_layout == MemoryLayout::ColMajor {
                x_rank_mat.swap_layout();
            }
            if z_rank_mat.memory_layout == MemoryLayout::ColMajor {
                z_rank_mat.swap_layout();
            }
            let x_rank = x_rank_mat.rank();
            let z_rank = z_rank_mat.rank();

            let mut column_checker = x_parity_check.clone();
            column_checker.transpose();
            column_checker.swap_layout();
            let mut num_zero_cols = 0;
            for row_ix in 0..column_checker.n_rows {
                let row = column_checker.get_row(row_ix);
                if row.is_zero() {
                    num_zero_cols += 1;
                }
            }
            println!(
                "x_parity_check shape: {:} x {:}, rank: {:}, num_zero_cols: {:}",
                x_parity_check.n_rows, x_parity_check.n_cols, x_rank, num_zero_cols
            );
            println!(
                "z_parity_check shape: {:} x {:}, rank: {:}, memory layout: {:?}",
                z_parity_check.n_rows, z_parity_check.n_cols, z_rank, z_parity_check.memory_layout
            );
            let k_x = x_parity_check.n_cols - x_rank;
            let k_z = z_parity_check.n_cols - z_rank;
            let k = k_x as i64 + k_z as i64 - qubits.len() as i64;
            println!("number logical qubits? = {:}", k);

            if let Some(filename) = output {
                let boundary_up = bu.to_dense();
                let boundary_down = bd.to_dense();
                let data = json!({
                    "memory_layout": "row_major",
                    "boundary_up_n_rows": boundary_up.n_rows,
                    "boundary_up_n_cols": boundary_up.n_cols,
                    "boundary_down_n_rows": boundary_down.n_rows,
                    "boundary_down_n_cols": boundary_down.n_cols,
                    "boundary_up": boundary_up.entries,
                    "boundary_down": boundary_down.entries,
                });
                if let Ok(serialized_data) = serde_json::to_string(&data) {
                    std::fs::write(filename, serialized_data)
                        .expect("Could not write data to disk");
                }
            }
            let tot_time = start.elapsed().as_secs_f64();
            println!("took this many seconds: {:}", tot_time);
        }
        Cli::RateAndDistance {
            quotient_poly,
            dim,
            directory,
            truncation,
            hgraph_log_rate,
            num_checkpoints,
        } => {
            let _logger = SimpleLogger::new().init().unwrap();
            let q = FFPolynomial::from_str(&quotient_poly[..]).unwrap();
            let mut rate_estimator = RateAndDistConfig::new(
                q,
                dim,
                truncation,
                hgraph_log_rate,
                Some(1),
                directory,
                num_checkpoints,
            );
            rate_estimator.run();
            let (i, b) = rate_estimator.quit();
        }
    }
}
