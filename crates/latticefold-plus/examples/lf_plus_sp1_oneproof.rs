//! LF+ one-proof harness for SP1 shrink verifier R1LF (research).
//!
//! This is the intended entrypoint for provers: consume SP1’s `.r1lf` + `{path}.chunks` cache
//! (fast random-access, padded dimensions) and run LF+ experiments on top of that statement.
//!
//! Usage:
//!   SP1_R1LF=/path/to/shrink_verifier.r1lf \
//!     cargo run -p latticefold-plus --example lf_plus_sp1_oneproof --features we_gate --release

use stark_rings::cyclotomic_ring::models::frog_ring::RqPoly as R;
use std::time::Instant;

fn main() {
    let r1lf_path = std::env::var("SP1_R1LF").expect("Set SP1_R1LF=/path/to/shrink.r1lf");
    let chunk_size: usize = std::env::var("CHUNK_SIZE")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(1 << 20);
    let pad_cols_to_multiple_of: usize = std::env::var("PAD_COLS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(256);

    println!("=========================================================");
    println!("LF+ SP1 One-Proof (R1LF loader + cache sanity)");
    println!("=========================================================");
    println!("  CHUNK_SIZE={chunk_size} PAD_COLS={pad_cols_to_multiple_of}");

    let t0 = Instant::now();
    let cache =
        latticefold_plus::sp1_r1lf::open_sp1_r1lf_chunk_cache::<R>(&r1lf_path, chunk_size, pad_cols_to_multiple_of)
            .expect("open_sp1_r1lf_chunk_cache");
    println!("  cache open: {:?}", t0.elapsed());
    println!("  chunks={} ncols={}", cache.num_chunks, cache.ncols);
    println!(
        "  stats: num_vars={} num_constraints={} num_public={} p_bb={} total_nonzeros={}",
        cache.stats.num_vars,
        cache.stats.num_constraints,
        cache.stats.num_public,
        cache.stats.p_bb,
        cache.stats.total_nonzeros
    );

    // Read chunk 0 as a basic sanity/perf check.
    let t1 = Instant::now();
    let [a, b, c] = cache.read_chunk(0).expect("read_chunk(0)");
    println!(
        "  read chunk0: {:?} (nrows={}, ncols={})",
        t1.elapsed(),
        a.nrows,
        a.ncols
    );
    let _ = (b, c);
}

