use core::panic;
use std::{
    env,
    fs::File,
    io::{Read, Write},
    path::{Path, PathBuf},
    str::FromStr,
};

use clap::Parser;
use hdx_codec::{
    hdx_code::{HDXCode, HDXCodeConfig},
    math::{
        coset_complex::CosetComplex,
        lps::{self, compute_graph},
        polynomial::FiniteFieldPolynomial,
    },
};
use mhgl::{HGraph, SparseBasis};

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

#[derive(Debug, Clone)]
enum CosetComplexCommands {
    NewComplex,
    SaveComplex,
    LoadComplex,
    ComputeGroup,
    ComputeVertices,
    ComputeEdges,
    ComputeTriangles,
    DegreeStatistics,
    Navigate,
    Lps,
    Print,
    PrintSubgroups,
    Quit,
}

#[derive(Debug)]
struct UserCommandParseError {}

impl FromStr for CosetComplexCommands {
    type Err = UserCommandParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match &s.trim().to_ascii_lowercase()[..] {
            "create" | "c" => Ok(CosetComplexCommands::NewComplex),
            "save" | "s" => Ok(CosetComplexCommands::SaveComplex),
            "load" | "l" => Ok(CosetComplexCommands::LoadComplex),
            "group" => Ok(CosetComplexCommands::ComputeGroup),
            "vertices" => Ok(CosetComplexCommands::ComputeVertices),
            "edges" => Ok(CosetComplexCommands::ComputeEdges),
            "triangles" => Ok(CosetComplexCommands::ComputeTriangles),
            "degrees" | "deg" => Ok(CosetComplexCommands::DegreeStatistics),
            "navigate" | "nav" => Ok(CosetComplexCommands::Navigate),
            "print" => Ok(CosetComplexCommands::Print),
            "print_subs" => Ok(CosetComplexCommands::PrintSubgroups),
            "quit" | "q" => Ok(CosetComplexCommands::Quit),
            "lps" => Ok(CosetComplexCommands::Lps),
            _ => Err(UserCommandParseError {}),
        }
    }
}

fn main_loop_help() {
    println!("Enter command from the following list:");
    println!("[create / c] to create a complex");
    println!("[save / s] to save a complex");
    println!("[load / l] to load a complex");
    println!("[group] to compute group");
    println!("[vertices] to compute vertices");
    println!("[edges] to compute edges");
    println!("[triangles] to compute triangles.");
    println!("[degrees | deg] to compute degree statistics");
    println!("[navigation | nav] to navigate the coset complex.");
    println!("[lps] enter the lps menu.");
    println!("[print] to print the hypergraph.");
    println!("[print_subs] to print subgroups.");
    println!("[quit / q] to quit.");
}

fn input_loop() {
    let mut input_buf = String::new();
    let mut coset_complex: Option<CosetComplex> = None;
    loop {
        main_loop_help();
        print!("> ");
        std::io::stdout().flush().unwrap();
        std::io::stdin()
            .read_line(&mut input_buf)
            .expect("Could not read input.");
        let command = CosetComplexCommands::from_str(&input_buf[..]);
        println!("entered command {:?}", command);
        if let Ok(cmd) = command {
            match cmd {
                CosetComplexCommands::NewComplex => {
                    if coset_complex.is_none() {
                        let mut file_base = String::new();
                        let mut polynomial_string = String::new();
                        let mut dimension_string = String::new();

                        println!("Enter a base name for data files:");
                        std::io::stdin()
                            .read_line(&mut file_base)
                            .expect("could not read file base input");

                        println!("Enter a polynomial to quotient by. Format is 'c_d * x^d + ... + c_0 * x ^ 0 % p', where the degree of each term must be specified, terms must be separated by a '+', and there must be a % p at the end where p is a prime that the coefficients are modulo'd by.");
                        std::io::stdin().read_line(&mut polynomial_string);

                        println!("Enter a dimension of matrices to use.");
                        std::io::stdin().read_line(&mut dimension_string);
                        let quotient = FiniteFieldPolynomial::from_str(&polynomial_string[..])
                            .expect("could not parse quotient polynomial.");
                        let dim: usize = dimension_string
                            .trim()
                            .parse()
                            .expect("could not parse dimension.");
                        let dm = CosetComplex::new(file_base, dim, &quotient);
                        coset_complex = Some(dm);
                        dbg!(&coset_complex);
                    }
                }
                CosetComplexCommands::SaveComplex => {
                    if let Some(dm) = &coset_complex {
                        dm.save_to_disk();
                    }
                }
                CosetComplexCommands::LoadComplex => {
                    if let Some(dm) = &mut coset_complex {
                        dm.load_from_disk();
                    }
                }
                CosetComplexCommands::ComputeGroup => {
                    if let Some(dm) = &mut coset_complex {
                        dm.generate_group();
                    }
                }
                CosetComplexCommands::ComputeVertices => {
                    if let Some(dm) = &mut coset_complex {
                        dm.compute_vertices()
                    }
                }
                CosetComplexCommands::ComputeEdges => {
                    if let Some(dm) = &mut coset_complex {
                        dm.compute_edges();
                    }
                }
                CosetComplexCommands::ComputeTriangles => {
                    if let Some(dm) = &mut coset_complex {
                        dm.compute_triangles();
                    }
                }
                CosetComplexCommands::DegreeStatistics => {
                    if let Some(dm) = &mut coset_complex {
                        dm.check_degrees();
                    }
                }
                CosetComplexCommands::Navigate => {
                    if let Some(cc) = &coset_complex {
                        graph_walk_loop(cc.hgraph_ref());
                    }
                }
                CosetComplexCommands::Lps => {
                    lps_menu();
                }
                CosetComplexCommands::Print => {
                    if let Some(dm) = &coset_complex {
                        dm.print_hgraph();
                    }
                }
                CosetComplexCommands::PrintSubgroups => {
                    if let Some(dm) = &coset_complex {
                        dm.print_subgroups();
                    } else {
                        println!("No CosetComplex, create one first.");
                    }
                }
                CosetComplexCommands::Quit => {
                    return;
                }
            }
        }
        input_buf.clear();
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
                            if let Some(g) = compute_graph(p, q) {
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
                        graph_walk_loop(hg);
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

fn graph_walk_loop(hgraph: &HGraph) {
    println!("Nodes in the graph:");
    println!("{:?}", hgraph.nodes());
    println!("enter a subset to start, format should be like [1, 2, 3]:");
    let mut start_set_buf = String::from("\"");
    std::io::stdin()
        .read_line(&mut start_set_buf)
        .expect("Could not read input.");
    let mut trimmed = start_set_buf.trim().to_string();
    trimmed.push('\"');

    let mut walker_location: SparseBasis<u32> =
        serde_json::from_str(&trimmed[..]).expect("Could not parse input.");
    println!("This is what was entered: {:}", walker_location);
    println!("What kind of walk to perform: link, star, up-down, or down-up?");
    let mut walk_buf = String::new();
    std::io::stdin()
        .read_line(&mut walk_buf)
        .expect("Could not read input.");
    match &walk_buf.trim().to_ascii_lowercase()[..] {
        "link" | "l" => loop {
            println!("Currently here: {:}", walker_location);
            let link = hgraph.link_as_vec(&walker_location.node_vec()[..]);
            println!("Here is where you can go:");
            for ix in 0..link.len() {
                println!("({:}) - {:?}", ix, link[ix].0);
            }
            walk_buf.clear();
            println!("Which one?");
            print!("lps> ");
            std::io::stdout().flush().expect("could not flush");
            std::io::stdin()
                .read_line(&mut walk_buf)
                .expect("Could not read input");
            if let Ok(choice) = walk_buf.trim().parse::<usize>() {
                if let Some(new_loc) = link.get(choice) {
                    let new_loc_vec: Vec<u32> = new_loc.0.clone().into_iter().collect();
                    walker_location = SparseBasis::from_slice(&new_loc_vec[..]);
                }
            } else {
                if walk_buf.starts_with("q") {
                    println!("Leaving navigation.");
                    break;
                }
            }
        },
        "star" | "s" => {
            println!("Not implemented.");
        }
        "up-down" | "ud" => {
            println!("Not implemented.");
        }
        "down-up" | "du" => {
            println!("Not implemented.");
        }

        _ => {
            println!("idk what you said, exiting.");
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

#[derive(Parser)]
struct HDXBuilder {
    directory: PathBuf,
    quotient: String,
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
    input_loop();
}
