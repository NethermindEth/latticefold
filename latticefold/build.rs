use serde::Deserialize;
use std::fmt::format;
use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::{env, fs};
#[derive(Debug, Deserialize)]
pub struct BenchmarkRecord {
    pub id: u32,
    pub X_LEN: usize,   // Matches "X_LEN" in TOML
    pub C: usize,       // Matches "C" in TOML
    pub W: usize,       // Matches "W" in TOML
    pub B: u64,         // Matches "B" in TOML
    pub L: usize,       // Matches "L" in TOML
    pub B_small: usize, // Matches "B_small" in TOML
    pub K: usize,       // Matches "K" in TOML
}

#[derive(Debug, Deserialize)]
pub struct Benchmarks {
    pub single_non_scalar_goldilocks: Vec<BenchmarkRecord>, // Matches "single_non_scalar_goldilocks" in TOML
}

#[derive(Debug, Deserialize)]
pub struct BenchmarkConfig {
    pub benchmarks: Benchmarks, // Matches the "[benchmarks]" table in TOML
}

fn single_benchmark(name: &str, cs: &str, ring: &str, br: BenchmarkRecord) -> String {
    format!(
        r#"
    {{
        const X_LEN: usize = {};
        const C: usize = {};
        const W: usize = {};
        const WIT_LEN: usize = {};
        #[derive(Clone)]
        struct DP {{}}
        impl DecompositionParams for DP {{
            const B: u128 = {};
            const L: usize = {};
            const B_SMALL: usize = {};
            const K: usize = {};
        }}

        {name}::<X_LEN, C, W, WIT_LEN, {cs}, {ring}, DP>(group);
    }}
    "#,
        br.X_LEN,
        br.C,
        br.W,
        br.W * br.L,
        br.B,
        br.L,
        br.B_small,
        br.K
    )
}

fn fn_call_as_str(name: &str, br: BenchmarkRecord) -> String {
    format!("fn {name}::<>")
}

fn main() -> Result<(), String> {
    // Read the TOML configuration file
    let toml_content = fs::read_to_string("benches/config.toml").expect("Failed to read TOML file");
    let config: BenchmarkConfig =
        toml::from_str(&toml_content).expect("Failed to deserialize benches/config.toml");

    // Get the output directory
    let out_dir = env::var("OUT_DIR").expect("OUT_DIR is not set");
    let dest_path = Path::new(&out_dir).join("generated_benchmarks.rs");

    // Open the file for writing
    let mut file = File::create(&dest_path).expect("Failed to create generated file");

    // Write the header for the function
    writeln!(file, "fn folding_single_goldilocks_non_scalar_benchmark(group: &mut BenchmarkGroup<WallTime>) {{")
        .expect("Failed to write to file");

    // Write benchmarks iteratively
    for benchmark in config.benchmarks.single_non_scalar_goldilocks {
        writeln!(
            file,
            "{}",
            single_benchmark(
                "folding_benchmarks_non_scalar",
                "GoldilocksChallengeSet",
                "GoldilocksRingNTT",
                benchmark
            )
        )
        .expect("Failed to write to file");
    }

    // Write the closing brace for the function
    writeln!(file, "}}").expect("Failed to write to file");

    Ok(())
}
