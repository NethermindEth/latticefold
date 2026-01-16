//! Cheat / soundness regression tests for the SP1 BabyBear→Frog lift.
//!
//! Goal: exercise the exact failure mode:
//! - If lift vars (`t_i` = carry/quotient) are *not* bounded, the lifted constraint is vacuous:
//!     t := (A*B - C) * p_bb^{-1}  (mod q_frog)
//!   always satisfies in the host field.
//! - With LF+ boundedness enforced via **balanced base-16 digits with k=8** (SP1/Frog64 setting),
//!   that cheating `t`
//!   is (overwhelmingly) out of range and cannot be decomposed -> rejected.

#![cfg(all(test, feature = "we_gate"))]

use ark_ff::Field;
use stark_rings::{
    balanced_decomposition::Decompose,
    cyclotomic_ring::models::frog_ring::Fq,
    Zq,
};

const P_BB: u64 = 2013265921;
const BOUND_BASE: u128 = 16;
const BOUND_K: usize = 8;

fn centered_i64(x: Fq) -> i64 {
    let mag = x
        .center()
        .to_u64()
        .expect("centered magnitude should fit in u64") as i64;
    let neg = x.sign() == -Fq::ONE;
    if neg { -mag } else { mag }
}

fn decompose_fits_sp1_frog64_bound(x: Fq) -> bool {
    // If x needs more than k digits in balanced base B, Decompose will panic (slice too short).
    // We treat that as "does not fit bound", which is exactly what the boundedness layer should enforce.
    let mut out = [Fq::ZERO; BOUND_K];
    std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        x.decompose_to(BOUND_BASE, &mut out);
    }))
    .is_ok()
}

#[test]
fn test_lift_vacuity_exists_but_boundedness_rejects_random() {
    // Show: for random a,b,c in [0,p_bb), we can always solve t in Frog field so that
    //   a*b = c + p_bb*t (mod q_frog)
    // but that t is (almost always) not representable under b=2^16,k=2 boundedness.
    use rand::Rng;

    let mut rng = ark_std::test_rng();
    let inv_p = Fq::from(P_BB).inverse().unwrap();

    // Find a sample where the computed t does NOT fit in the SP1/Frog64 boundedness regime.
    for _ in 0..10_000 {
        let a_u: u64 = rng.gen_range(0..P_BB);
        let b_u: u64 = rng.gen_range(0..P_BB);
        let c_u: u64 = rng.gen_range(0..P_BB);
        let a = Fq::from(a_u);
        let b = Fq::from(b_u);
        let c = Fq::from(c_u);

        // Cheating choice of t in the host field:
        //   t := (a*b - c) / p_bb  (mod q_frog)
        let t = (a * b - c) * inv_p;

        // This ALWAYS satisfies in the host field.
        let lhs = a * b;
        let rhs = c + Fq::from(P_BB) * t;
        assert_eq!(lhs, rhs);

        // Under boundedness (base=16,k=8), this t should almost never fit.
        if !decompose_fits_sp1_frog64_bound(t) {
            // Nice to sanity-print magnitude if needed while debugging.
            let _mag = centered_i64(t);
            return;
        }
    }

    panic!("could not find an out-of-range cheating t in 10k trials; unexpected");
}

#[test]
fn test_lift_nonvacuous_accepts_small_valid_case() {
    // Construct a small valid BabyBear-style mul where quotient is 0 (so boundedness trivially holds).
    // Choose a,b < sqrt(p_bb) so a*b < p_bb.
    let a_u: u64 = 12345;
    let b_u: u64 = 23456;
    let prod = (a_u as u128) * (b_u as u128);
    assert!(prod < (P_BB as u128));

    let a = Fq::from(a_u);
    let b = Fq::from(b_u);
    let c = Fq::from(prod as u64);
    let t = Fq::ZERO; // quotient/carry is 0

    assert_eq!(a * b, c + Fq::from(P_BB) * t);
    assert!(decompose_fits_sp1_frog64_bound(t));
}

#[test]
fn test_lift_add_carry_example_is_bounded() {
    // BabyBear add overflow example:
    //   (p-1) + 5 = 4 + p*1
    // carry=1 must be small/bounded.
    let a_u: u64 = P_BB - 1;
    let b_u: u64 = 5;
    let c_u: u64 = 4;
    let carry_u: u64 = 1;

    let a = Fq::from(a_u);
    let b = Fq::from(b_u);
    let c = Fq::from(c_u);
    let carry = Fq::from(carry_u);

    // Model linear constraint as (A=1, B = a+b, C = c + p*carry) or equivalently just check integer equality:
    assert_eq!(a + b, c + Fq::from(P_BB) * carry);
    assert!(decompose_fits_sp1_frog64_bound(carry));
}

#[test]
fn test_lift_mul_cheat_t_is_rejected_by_boundedness() {
    // Construct a single lifted mul constraint:
    //   a*b = c + p*t
    // where a,b,c are in [0,p_bb). Choose random values and set
    //   t := (a*b - c)/p  (mod q_frog)
    // which always satisfies in Frog, but almost never fits b=2^16,k=2 boundedness.
    use rand::Rng;

    let mut rng = ark_std::test_rng();
    let inv_p = Fq::from(P_BB).inverse().unwrap();
    for _ in 0..50_000 {
        let a_u: u64 = rng.gen_range(0..P_BB);
        let b_u: u64 = rng.gen_range(0..P_BB);
        let c_u: u64 = rng.gen_range(0..P_BB);
        let a = Fq::from(a_u);
        let b = Fq::from(b_u);
        let c = Fq::from(c_u);
        let t = (a * b - c) * inv_p;
        assert_eq!(a * b, c + Fq::from(P_BB) * t);
        if !decompose_fits_sp1_frog64_bound(t) {
            return;
        }
    }
    panic!("could not find a cheating t rejected by boundedness in 50k trials; unexpected");
}

