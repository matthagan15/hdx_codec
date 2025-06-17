use clap::Parser;

#[derive(Parser)]
#[command(version, about, long_about = None)]
enum Cli {
    RowReduction {},
    BFS,
}

fn main() {
    use hdx_codec::row_reduction_benchmark;
    let cli = Cli::parse();
    match cli {
        Cli::RowReduction {} => {
            row_reduction_benchmark("/Users/matt/repos/qec/bench/galois_bench_data.json".into());
        }
        Cli::BFS => {
            hdx_codec::math::coset_complex_bfs::bfs_benchmark();
        }
    }
}
