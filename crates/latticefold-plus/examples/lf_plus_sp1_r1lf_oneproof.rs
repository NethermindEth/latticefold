//! LF+ harness for reading SP1 "R1LF" lifted R1CS files (research).
//!
//! Usage:
//!   SP1_R1LF=/path/to/shrink_verifier.r1lf \
//!     cargo run -p latticefold-plus --example lf_plus_sp1_r1lf_oneproof --features we_gate --release

#![cfg(feature = "we_gate")]

use stark_rings::cyclotomic_ring::models::frog_ring::RqPoly as R;
use stark_rings::PolyRing;
use ark_ff::Field;

fn is_const_coeff_ring(x: &R) -> bool {
    x.coeffs()
        .iter()
        .skip(1)
        .all(|c| *c == <R as PolyRing>::BaseRing::ZERO)
}

fn is_const_coeff_sparse_matrix(m: &stark_rings_linalg::SparseMatrix<R>) -> bool {
    for row in &m.coeffs {
        for (c, _j) in row {
            if !is_const_coeff_ring(c) {
                return false;
            }
        }
    }
    true
}

fn main() {
    let path = std::env::var("SP1_R1LF").expect("Set SP1_R1LF=/path/to/file.r1lf");
    let chunk_size: usize = std::env::var("CHUNK_SIZE")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(1 << 20);
    let pad_cols_to_multiple_of: usize = std::env::var("PAD_COLS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(256);
    let scan_aux: bool = std::env::var("SCAN_AUX")
        .ok()
        .as_deref()
        .map(|s| s == "1")
        .unwrap_or(false);

    println!("=========================================================");
    println!("LF+ SP1 R1LF One-Proof (loader + const-coeff inspection)");
    println!("=========================================================");
    println!("  CHUNK_SIZE={chunk_size} PAD_COLS={pad_cols_to_multiple_of} SCAN_AUX={scan_aux}");

    let t0 = std::time::Instant::now();
    let cache =
        latticefold_plus::sp1_r1lf::open_sp1_r1lf_chunk_cache::<R>(&path, chunk_size, pad_cols_to_multiple_of)
            .expect("open R1LF chunk cache");
    println!("  open/build cache: {:?}", t0.elapsed());
    println!(
        "  header: num_vars={} num_constraints={} num_public={} p_bb={} chunks={} ncols={}",
        cache.stats.num_vars,
        cache.stats.num_constraints,
        cache.stats.num_public,
        cache.stats.p_bb,
        cache.num_chunks,
        cache.ncols
    );
    let nvars = latticefold_plus::sp1_r1lf::nvars_from_ncols_pow2(cache.ncols)
        .expect("cache.ncols should be pow2");
    println!("  derived: nvars_cm=nvars_setchk={} (from ncols=2^{})", nvars, nvars);
    println!("  digest={:02x?}...", &cache.stats.digest[..8]);

    if scan_aux {
        let t = std::time::Instant::now();
        let (_kinds, summary) =
            latticefold_plus::sp1_lift_witness::scan_aux_vars_from_r1lf_chunks(&cache)
                .expect("scan_aux_vars_from_r1lf_chunks");
        println!(
            "  aux vars: total={} carry={} quotient={}  (scan {:?})",
            summary.num_aux_vars,
            summary.num_carry,
            summary.num_quotient,
            t.elapsed()
        );
    }

    let t1 = std::time::Instant::now();
    let [a, b, c] = cache.read_chunk(0).expect("read_chunk(0)");
    println!("  read chunk0: {:?}", t1.elapsed());
    println!("  chunk0 dims: nrows={} ncols={}", a.nrows, a.ncols);
    println!(
        "  const-coeff chunk0: A={} B={} C={}",
        is_const_coeff_sparse_matrix(&a),
        is_const_coeff_sparse_matrix(&b),
        is_const_coeff_sparse_matrix(&c)
    );
}

