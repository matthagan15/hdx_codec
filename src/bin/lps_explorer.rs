use std::{
    io::{Read, Write},
    path::PathBuf,
    str::FromStr,
};

use clap::{Parser, Subcommand};
use hdx_codec::{lps::compute_lps_graph, math::finite_field::FFRep};
use mhgl::HGraph;
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

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    /// The config for the specified complex to use. If none is provided we ask the user via a wizard to figure out a config.
    #[arg(short, long, value_name = "FILE")]
    config: Option<PathBuf>,

    #[command(subcommand)]
    lps_commands: LpsCommands,
}

#[derive(Subcommand)]
pub enum LpsCommands {
    /// Computes and saves the graph LPS(p, q) to
    /// disk with the filename lps_p_q.hg
    Build {
        #[arg(short, long, value_name = "DIRECTORY")]
        directory: PathBuf,

        #[arg(short, value_name = "PRIME_1")]
        /// The "p" prime used in LPS(p, q)
        p: FFRep,

        #[arg(short, value_name = "PRIME_2")]
        /// The "q" prime used in LPS(p, q)
        q: FFRep,
    },
    View {
        #[arg(short, long, value_name = "FILENAME")]
        filename: PathBuf,
    },
    Navigate {},
    /// Compute the rank and distance of the
    ParityCode {
        #[arg(short, long, value_name = "FILENAME")]
        filename: PathBuf,
    },
    ReedSolomonCode {
        #[arg(short, long, value_name = "FILENAME")]
        filename: PathBuf,

        #[arg(short, long, value_name = "FIELD_MOD")]
        field_mod: FFRep,

        #[arg(short, long, value_name = "MAX_DEGREE")]
        /// The max degree (non-inclusive) of the polynomial in the ReedSolomon Code
        degree: usize,
    },
}

fn main() {
    let cli = Cli::parse();
    match cli.lps_commands {
        LpsCommands::Build { directory, p, q } => {
            if let Some(graph) = compute_lps_graph(p, q) {
                let mut directory = directory.clone();
                directory.push("lps_");
                directory.push(p.to_string());
                directory.push("_");
                directory.push(q.to_string());
                directory.push(".hg");
                graph.to_disk(&directory);
            } else {
                println!("Could not compute graph.");
            }
        }
        LpsCommands::View { filename } => todo!(),
        LpsCommands::Navigate {} => todo!(),
        LpsCommands::ParityCode { filename } => todo!(),
        LpsCommands::ReedSolomonCode {
            filename,
            field_mod,
            degree,
        } => todo!(),
    }
}
