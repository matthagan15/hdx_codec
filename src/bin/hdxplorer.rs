use std::str::FromStr;

use qec::math::coset_complex::DiskManager;

#[derive(Debug, Clone)]
enum UserCommand {
    CreateComplex,
    SaveComplex,
    LoadComplex,
}

#[derive(Debug)]
struct UserCommandParseError {}

impl FromStr for UserCommand {
    type Err = UserCommandParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match &s.trim().to_ascii_lowercase()[..] {
            "create" | "c" => {Ok(UserCommand::CreateComplex)},
            "save" | "s" => {Ok(UserCommand::SaveComplex)},
            "load" | "l" => {Ok(UserCommand::LoadComplex)},
            _ => {Err(UserCommandParseError {  })}
        }
    }
}

fn input_loop() {
    let mut input_buf = String::new();
    loop {
        println!("Enter command from the following list:");
        println!("[create / c] to create a complex");
        println!("[save / s] to save a complex");
        println!("[load / l] to load a complex");
        std::io::stdin().read_line(&mut input_buf).expect("Could not read input.");
        let command = UserCommand::from_str(&input_buf[..]);
        println!("entered command {:?}", command);
        input_buf.clear();
    }
}

fn main() {
    println!("testing, 1,2,3.");
    input_loop();
}