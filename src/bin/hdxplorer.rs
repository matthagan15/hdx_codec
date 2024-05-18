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
    hdx_code::HDXCodeConfig,
    math::{
        finite_field::FFRep,
        iterative_bfs_new::GroupBFS,
        lps::{self, compute_lps_graph},
        polynomial::FiniteFieldPolynomial,
    },
    reed_solomon::ReedSolomon,
    tanner_code::TannerCode,
};
use mhgl::HGraph;

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

fn degree_stats(hg: &HGraph) {
    println!("Checking degrees.");
    let nodes = hg.nodes();
    let mut node_to_node_stats = HashMap::new();
    let mut node_to_edges_stats = HashMap::new();
    let edges = hg.edges_of_size(2);
    let mut edge_stats = HashMap::new();
    println!("nodes.");
    for node in nodes {
        let link = hg.link_as_vec(&[node]);
        let mut num_nodes = 0;
        let mut num_edges = 0;
        for (set, weight) in link.into_iter() {
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
        if let Some(edge_set) = hg.query_edge_id(&edge) {
            let mut num_nodes = 0;
            for (set, weight) in hg.link_as_vec(&edge_set[..]) {
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
pub enum HdxCommands {
    Build {
        #[arg(short, long, value_name = "DIRECTORY")]
        directory: PathBuf,
        #[arg(short, long, value_name = "QUOTIENT")]
        quotient: String,
        #[arg(long, value_name = "DIM")]
        dim: usize,
        /// If you want to label the output files.
        /// The default filename has the format `p_<QUOTIENT.field_mod>_<QUOTIENT.degree>_<DIM>' and
        /// then whatever extensions we use `.hg`, `.dot`, `.cache`, `.conf`
        #[arg(short, long, value_name = "FILENAME")]
        filename: Option<String>,

        /// Upper limit how many BFS steps the walker will traverse of the graph. Recorded in the conf file.
        #[arg(short, long, value_name = "MAX_BFS_STEPS")]
        max_bfs_steps: Option<usize>,
    },
    View {
        /// Currently you write out the entire damn thing.
        #[arg(short, long, value_name = "FILENAME")]
        filename: String,
    },
    Compute {
        #[arg(short, long, value_name = "DIRECTORY")]
        directory: PathBuf,
    },

    Ranks,
}

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    /// The config for the specified complex to use. If none is provided we ask the user via a wizard to figure out a config.
    #[arg(short, long, value_name = "FILE")]
    config: Option<PathBuf>,

    #[command(subcommand)]
    command: HdxCommands,
}

fn main() {
    println!("testing, 1,2,3.");
    // let conf_from_cur_dir = get_config_from_current_working_dir();
    // let hdx_conf = if conf_from_cur_dir.is_some() {
    //     conf_from_cur_dir.unwrap()
    // } else {
    //     get_hdx_config_from_user()
    // };
    // hdx_conf.save_to_disk();
    // let hdx_code = HDXCode::new(hdx_conf);
    // let cli = Cli::parse();

    let cli = Cli::parse();
    match cli.command {
        HdxCommands::Build {
            directory,
            quotient,
            dim,
            filename,
            max_bfs_steps,
        } => {
            let q = FiniteFieldPolynomial::from_str(&quotient)
                .expect("Could not parse quotient argument.");
            let mut bfs = GroupBFS::new(&directory, &q);
            bfs.bfs(max_bfs_steps.unwrap_or(usize::MAX));
        }
        HdxCommands::View { filename } => {
            let mut pathbuf = PathBuf::from(&filename);
            let hg = HGraph::from_file(&pathbuf).expect("Could not find hgraph.");
            println!("hg: {:}", hg);
            degree_stats(&hg);
            dbg!(hg);
        }
        HdxCommands::Compute { directory } => {
            let p = 3_u32;
            let primitive_coeffs = [(2, (1, p).into()), (1, (2, p).into()), (0, (2, p).into())];
            let q = FiniteFieldPolynomial::from(&primitive_coeffs[..]);
            let mut bfs_manager = hdx_codec::math::iterative_bfs_new::GroupBFS::new(&directory, &q);
            bfs_manager.bfs(usize::MAX);
        }
        HdxCommands::Ranks => {
            let directory = "/u/hagan/hdx_outputs/big_one/";
            let names = vec![
                "bfs_100000.hg",
                "bfs_1000000.hg",
                "bfs_1100000.hg",
                "bfs_1200000.hg",
                "bfs_1300000.hg",
                "bfs_1400000.hg",
                "bfs_1500000.hg",
                "bfs_1600000.hg",
                "bfs_1700000.hg",
                "bfs_1800000.hg",
                "bfs_1900000.hg",
                "bfs_200000.hg",
                "bfs_2000000.hg",
                "bfs_2100000.hg",
                "bfs_2200000.hg",
                "bfs_2300000.hg",
                "bfs_2400000.hg",
                "bfs_2500000.hg",
                "bfs_2600000.hg",
                "bfs_2700000.hg",
                "bfs_2800000.hg",
                "bfs_2900000.hg",
                "bfs_300000.hg",
                "bfs_3000000.hg",
                "bfs_3100000.hg",
                "bfs_3200000.hg",
                "bfs_3300000.hg",
                "bfs_3400000.hg",
                "bfs_3500000.hg",
                "bfs_3600000.hg",
                "bfs_3700000.hg",
                "bfs_3800000.hg",
                "bfs_3900000.hg",
                "bfs_400000.hg",
                "bfs_4000000.hg",
                "bfs_4100000.hg",
                "bfs_4200000.hg",
                "bfs_4300000.hg",
                "bfs_4400000.hg",
                "bfs_4500000.hg",
                "bfs_4600000.hg",
                "bfs_4700000.hg",
                "bfs_4800000.hg",
                "bfs_4900000.hg",
                "bfs_500000.hg",
                "bfs_5000000.hg",
                "bfs_5100000.hg",
                "bfs_5200000.hg",
                "bfs_5300000.hg",
                "bfs_5400000.hg",
                "bfs_5500000.hg",
                "bfs_5600000.hg",
                "bfs_5700000.hg",
                "bfs_5800000.hg",
                "bfs_5900000.hg",
                "bfs_600000.hg",
                "bfs_6000000.hg",
                "bfs_6100000.hg",
                "bfs_6200000.hg",
                "bfs_6300000.hg",
                "bfs_6400000.hg",
                "bfs_6500000.hg",
                "bfs_6600000.hg",
                "bfs_6700000.hg",
                "bfs_6800000.hg",
                "bfs_6900000.hg",
                "bfs_700000.hg",
                "bfs_7000000.hg",
                "bfs_7100000.hg",
                "bfs_7200000.hg",
                "bfs_7300000.hg",
                "bfs_7400000.hg",
                "bfs_7500000.hg",
                "bfs_7600000.hg",
                "bfs_7700000.hg",
                "bfs_7800000.hg",
                "bfs_7900000.hg",
                "bfs_800000.hg",
                "bfs_8000000.hg",
                "bfs_8100000.hg",
                "bfs_8200000.hg",
                "bfs_8300000.hg",
                "bfs_8400000.hg",
                "bfs_8500000.hg",
                "bfs_8600000.hg",
                "bfs_8700000.hg",
                "bfs_8800000.hg",
                "bfs_8900000.hg",
                "bfs_900000.hg",
            ];
            for name in names {
                let mut hg_path = PathBuf::from(directory);
                hg_path.push(name);
                let hg = HGraph::from_file(&hg_path).expect("Could not read hgraph.");
                let tc = TannerCode::<ReedSolomon>::new(hg, 2, 3, 2);
                println!("Code created, computing matrix.");
                let mut mat = tc.sparse_parity_check_matrix();
                println!("Computed matrix, size: {:} x {:}", mat.n_rows, mat.n_cols);
                println!("Changing memory layout to row layout for future use.");
                mat.to_row_layout();
                println!("computing rank.");
                println!("rank per dim: {:}", mat.rank() as f64 / mat.n_cols as f64);
            }
        }
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
