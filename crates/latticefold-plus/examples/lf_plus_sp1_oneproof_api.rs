//! LF+ one-proof harness for SP1 shrink verifier R1LF (API-driven).
//!
//! This example is intentionally a **thin wrapper** over:
//! `latticefold_plus::sp1_oneproof_api::run_sp1_oneproof_we_gate_from_files`.
//!
//! Usage:
//!   SP1_R1LF=/path/to/shrink_verifier.r1lf \
//!   SP1_WITNESS=/path/to/shrink_verifier.witness.bundle \
//!     cargo run -p latticefold-plus --example lf_plus_sp1_oneproof --features we_gate --release

#![cfg(feature = "we_gate")]

fn hex32(bytes: [u8; 32]) -> String {
    const HEX: &[u8; 16] = b"0123456789abcdef";
    let mut out = String::with_capacity(64);
    for &b in &bytes {
        out.push(HEX[(b >> 4) as usize] as char);
        out.push(HEX[(b & 0x0f) as usize] as char);
    }
    out
}

fn main() {
    let r1lf_path = std::env::var("SP1_R1LF").expect("Set SP1_R1LF=/path/to/shrink.r1lf");
    let witness_path =
        std::env::var("SP1_WITNESS").expect("Set SP1_WITNESS=/path/to/shrink_verifier.witness.bundle");

    let out = latticefold_plus::sp1_oneproof_api::run_sp1_oneproof_we_gate_from_files(
        &r1lf_path,
        &witness_path,
    )
    .expect("run_sp1_oneproof_we_gate_from_files");

    println!("stmt_digest=0x{}", hex32(out.stmt_digest));
    println!("lock_coin_seed=0x{}", hex32(out.lock_coin_seed));
}

