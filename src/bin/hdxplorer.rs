use std::str::FromStr;

use qec::math::{coset_complex::DiskManager, polynomial::FiniteFieldPolynomial};

#[derive(Debug, Clone)]
enum UserCommand {
    NewComplex,
    SaveComplex,
    LoadComplex,
    ComputeGroup,
    ComputeVertices,
    ComputeEdges,
    ComputeTriangles,
    Print,
    Quit,
}

#[derive(Debug)]
struct UserCommandParseError {}

impl FromStr for UserCommand {
    type Err = UserCommandParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match &s.trim().to_ascii_lowercase()[..] {
            "create" | "c" => Ok(UserCommand::NewComplex),
            "save" | "s" => Ok(UserCommand::SaveComplex),
            "load" | "l" => Ok(UserCommand::LoadComplex),
            "group" => Ok(UserCommand::ComputeGroup),
            "vertices" => Ok(UserCommand::ComputeVertices),
            "edges" => Ok(UserCommand::ComputeEdges),
            "triangles" => Ok(UserCommand::ComputeTriangles),
            "print" => Ok(UserCommand::Print),
            "quit" | "q" => Ok(UserCommand::Quit),
            _ => Err(UserCommandParseError {}),
        }
    }
}

fn input_loop() {
    let mut input_buf = String::new();
    let mut disk_manager: Option<DiskManager> = None;
    loop {
        println!("Enter command from the following list:");
        println!("[create / c] to create a complex");
        println!("[save / s] to save a complex");
        println!("[load / l] to load a complex");
        println!("[group] to compute group");
        println!("[vertices] to compute vertices");
        println!("[edges] to compute edges");
        println!("[triangles] to compute triangles.");
        println!("[print] to print the hypergraph.");
        println!("[quit / q] to quit.");
        std::io::stdin()
            .read_line(&mut input_buf)
            .expect("Could not read input.");
        let command = UserCommand::from_str(&input_buf[..]);
        println!("entered command {:?}", command);
        if let Ok(cmd) = command {
            match cmd {
                UserCommand::NewComplex => {
                    if disk_manager.is_none() {
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
                        let dm = DiskManager::new(file_base, dim, &quotient);
                        disk_manager = Some(dm);
                        dbg!(&disk_manager);
                    }
                }
                UserCommand::SaveComplex => {
                    if let Some(dm) = &disk_manager {
                        dm.save_to_disk();
                    }
                }
                UserCommand::LoadComplex => {
                    if let Some(dm) = &mut disk_manager {
                        dm.load_from_disk();
                    }
                }
                UserCommand::ComputeGroup => {
                    if let Some(dm) = &mut disk_manager {
                        dm.generate_group();
                    }
                }
                UserCommand::ComputeVertices => {
                    if let Some(dm) = &mut disk_manager {
                        dm.compute_vertices()
                    }
                }
                UserCommand::ComputeEdges => {
                    if let Some(dm) = &mut disk_manager {
                        dm.compute_edges();
                    }
                }
                UserCommand::ComputeTriangles => {
                    if let Some(dm) = &mut disk_manager {
                        dm.compute_triangles();
                    }
                }
                UserCommand::Print => {
                    if let Some(dm) = &disk_manager {
                        dm.print_hgraph();
                    }
                }
                UserCommand::Quit => {
                    return;
                }
            }
        }
        input_buf.clear();
    }
}

fn main() {
    println!("testing, 1,2,3.");
    input_loop();
}
