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
        iterative_bfs::GroupBFS,
        lps::{self, compute_lps_graph},
        polynomial::FiniteFieldPolynomial,
    },
};
use mhgl::{HGraph};

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

#[derive(Debug)]
enum LpsCommand {
    ComputeGraph,
    SaveGraph,
    NavigateGraph,
    PrintGraph,
    Quit,
}

struct LpsError {}
impl FromStr for LpsCommand {
    type Err = LpsError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let trimmed = s.trim().to_ascii_lowercase();
        match &trimmed[..] {
            "compute" | "c" => Ok(LpsCommand::ComputeGraph),
            "save" | "s" => Ok(LpsCommand::SaveGraph),
            "nav" | "navigate" => Ok(LpsCommand::NavigateGraph),
            "print" | "p" => Ok(LpsCommand::PrintGraph),
            "q" | "quit" => Ok(LpsCommand::Quit),
            _ => {
                println!("Error parsing input, here it is post trim: {:?}", trimmed);
                println!("What you want to do.");
                Err(LpsError {})
            }
        }
    }
}

fn lps_menu() {
    let mut input_buf = String::new();
    let mut lps_graph: Option<HGraph> = None;
    loop {
        println!("LPS specific commands:");
        println!("[compute | c] Compute an LPS graph.");
        println!("[save | s] Save the graph to disk.");
        println!("[nav | navigate] Do a manual walk on the graph");
        println!("[p | print] Print the graph");
        println!("[q | quit] Quit this sub menu.");
        print!("lps> ");
        std::io::stdout().flush().unwrap();
        input_buf.clear();
        std::io::stdin()
            .read_line(&mut input_buf)
            .expect("Could not read input.");
        let command = LpsCommand::from_str(&input_buf[..]);
        if let Ok(cmd) = &command {
            match cmd {
                LpsCommand::ComputeGraph => {
                    println!("Enter a prime p: ");
                    let mut prime_buf = String::new();
                    std::io::stdin()
                        .read_line(&mut prime_buf)
                        .expect("Could not read input.");
                    if let Ok(p) = prime_buf.trim_end().parse::<u32>() {
                        let mut q_buf = String::new();
                        println!("Enter second number:");
                        std::io::stdin()
                            .read_line(&mut q_buf)
                            .expect("Enter second number.");
                        if let Ok(q) = q_buf.trim_end().parse::<u32>() {
                            if let Some(g) = compute_lps_graph(p, q) {
                                lps_graph = Some(g);
                                println!("computed graph successfully.");
                            } else {
                                println!("didn't work, try again.");
                            }
                        }
                    } else {
                        println!("Could not parse input prime, try again.");
                    }
                }
                LpsCommand::SaveGraph => {
                    println!("Enter filename to save graph in json.");
                    let mut filename = String::new();
                    std::io::stdin()
                        .read_to_string(&mut filename)
                        .expect("could not read filename input.");
                    if let Some(hg) = &lps_graph {
                        let hg_string =
                            serde_json::to_string(hg).expect("could not serialize graph.");
                        let mut lps_file = std::fs::File::create(filename)
                            .expect("could not open file for writing.");
                        lps_file
                            .write_all(hg_string.as_bytes())
                            .expect("Could not write file to disk.");
                        println!("Done writing graph to disk.");
                    } else {
                        println!("No graph exists to save. Try again.")
                    }
                }
                LpsCommand::NavigateGraph => {
                    if let Some(hg) = &lps_graph {
                        println!("Entering walk loop.");
                        println!("Not Implemented yet.");
                    } else {
                        println!("No graph to walk on. Try computing one first.");
                    }
                }
                LpsCommand::PrintGraph => {
                    if let Some(hg) = &lps_graph {
                        println!("{:}", hg);
                    }
                }
                LpsCommand::Quit => {
                    println!("Returning to main menu.");
                    return;
                }
            }
        }
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
pub enum Commands {
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
}

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    /// The config for the specified complex to use. If none is provided we ask the user via a wizard to figure out a config.
    #[arg(short, long, value_name = "FILE")]
    config: Option<PathBuf>,

    #[command(subcommand)]
    command: Commands,
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
        Commands::Build {
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
        Commands::View { filename } => 
        {
            let mut pathbuf = PathBuf::from(&filename);
            let hg = HGraph::from_file(&pathbuf).expect("Could not find hgraph.");
            println!("hg: {:}", hg);
            degree_stats(&hg);
            dbg!(hg);
        },
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
