//! LF+ one-proof harness for SP1 shrink verifier R1CS (research).
//!
//! Goal: apples-to-apples with Symphony’s `symphony_sp1_oneproof`:
//! - load the same SP1 R1CS chunk cache
//! - inspect whether the **const-coeff** fast paths are applicable
//! - (optionally) run a tiny streaming sumcheck round over a `SparseMatVecConstCoeff` MLE
//!
//! Usage:
//!   SP1_R1CS=/path/to/shrink_verifier.r1cs \
//!     cargo run -p latticefold-plus --example lf_plus_sp1_oneproof --features sp1_import --release

use std::sync::Arc;
use std::time::Instant;

use stark_rings::cyclotomic_ring::models::frog_ring::RqPoly as R;
use stark_rings::{CoeffRing, PolyRing, Ring};

use latticefold_plus::streaming_sumcheck::StreamingMleEnum;

// Reuse Symphony’s on-disk loader format.
use symphony::sp1_r1cs_loader::FieldFromU64;

/// BabyBear field element for loading R1CS.
#[derive(Debug, Clone, Copy, Default)]
struct BabyBear(u64);

const BABYBEAR_P: u64 = 0x78000001; // 2013265921

impl FieldFromU64 for BabyBear {
    fn from_canonical_u64(val: u64) -> Self {
        BabyBear(val % BABYBEAR_P)
    }
    fn as_canonical_u64(&self) -> u64 {
        self.0
    }
}

#[inline]
fn is_const_coeff_ring(x: &R) -> bool {
    x.coeffs().iter().skip(1).all(|c| *c == <R as PolyRing>::BaseRing::ZERO)
}

#[inline]
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
    let r1cs_path = std::env::var("SP1_R1CS").expect("Set SP1_R1CS=/path/to/shrink.r1cs");
    let chunk_size: usize = std::env::var("CHUNK_SIZE")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(1 << 20);
    let pad_cols_to_multiple_of: usize = std::env::var("PAD_COLS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(256);
    let scan_all: bool = std::env::var("SCAN_ALL")
        .ok()
        .as_deref()
        .map(|s| s == "1")
        .unwrap_or(false);
    let alloc_witness0: bool = std::env::var("ALLOC_WITNESS0")
        .ok()
        .as_deref()
        .map(|s| s == "1")
        .unwrap_or(false);

    println!("=========================================================");
    println!("LF+ SP1 One-Proof (loader + const-coeff inspection)");
    println!("=========================================================");
    println!(
        "  CHUNK_SIZE={chunk_size} PAD_COLS={pad_cols_to_multiple_of} SCAN_ALL={} ALLOC_WITNESS0={}",
        scan_all, alloc_witness0
    );

    let t0 = Instant::now();
    let cache = latticefold_plus::sp1_r1cs::open_sp1_cache::<R, BabyBear>(
        &r1cs_path,
        chunk_size,
        pad_cols_to_multiple_of,
    )
    .expect("open_sp1_cache");
    println!("  cache open: {:?}", t0.elapsed());
    println!("  chunks={} ncols={}", cache.num_chunks, cache.ncols);
    println!(
        "  stats: num_vars={} num_constraints={} num_public={} total_nonzeros={}",
        cache.stats.num_vars, cache.stats.num_constraints, cache.stats.num_public, cache.stats.total_nonzeros
    );

    // Read chunk 0 and inspect const-coeff property.
    let t1 = Instant::now();
    let [a, b, c] = cache.read_chunk(0).expect("read_chunk(0)");
    println!(
        "  read chunk0: {:?} (nrows={}, ncols={})",
        t1.elapsed(),
        a.nrows,
        a.ncols
    );

    let a_const = is_const_coeff_sparse_matrix(&a);
    let b_const = is_const_coeff_sparse_matrix(&b);
    let c_const = is_const_coeff_sparse_matrix(&c);
    println!("  const-coeff matrices: A={} B={} C={}", a_const, b_const, c_const);

    if scan_all {
        let t_scan = Instant::now();
        let mut ok_all = true;
        for i in 0..cache.num_chunks {
            let [aa, bb, cc] = cache.read_chunk(i).expect("read_chunk(i)");
            let ok = is_const_coeff_sparse_matrix(&aa)
                && is_const_coeff_sparse_matrix(&bb)
                && is_const_coeff_sparse_matrix(&cc);
            if !ok {
                println!("  scan_all: chunk {i} is NOT const-coeff");
                ok_all = false;
                break;
            }
        }
        println!("  scan_all: all_chunks_const_coeff={} ({:?})", ok_all, t_scan.elapsed());
    }

    // Optional demo: if A is const-coeff and you allow allocating an O(ncols) witness0, evaluate a few rows
    // using `SparseMatVecConstCoeff` exactly like Symphony’s fast path.
    //
    // NOTE: for SP1 shrink, ncols can be huge (e.g. 67M). This allocation can be hundreds of MB.
    if a_const {
        if alloc_witness0 {
            let a_arc = Arc::new(a);
            let ncols = a_arc.ncols;
            let t_w = Instant::now();
            let mut witness0 = vec![<R as PolyRing>::BaseRing::ZERO; ncols];
            witness0[0] = <R as PolyRing>::BaseRing::ONE;
            let witness0 = Arc::new(witness0);
            println!("  alloc witness0: {:?} (len={})", t_w.elapsed(), ncols);

            let mle = StreamingMleEnum::<R>::SparseMatVecConstCoeff {
                matrix: a_arc.clone(),
                witness0,
                num_vars: (chunk_size.trailing_zeros() as usize), // chunk nrows is a power-of-two in this setup
            };

            let t_eval = Instant::now();
            let v0 = mle.eval_at_index(0);
            let v1 = mle.eval_at_index(1);
            println!(
                "  demo: SparseMatVecConstCoeff evals: v[0].ct={:?} v[1].ct={:?} (eval_time={:?})",
                v0.ct(),
                v1.ct(),
                t_eval.elapsed()
            );
        } else {
            println!("  demo: skipped witness0 allocation (set ALLOC_WITNESS0=1 to run)");
        }
    } else {
        println!("  demo: skipped (A not const-coeff in this representation)");
    }
}

