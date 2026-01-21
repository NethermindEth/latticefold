use ark_ff::{BigInteger, Field, Fp64, MontBackend, MontConfig, PrimeField, Zero};
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;

use dpp::accepting_set::accepting_set_for_packed_query_sparse;
use dpp::boolean_proof::BooleanProofFlpcpSparse;
use dpp::dr1cs_flpcp::{Dr1csInstanceSparse, RsDr1csNpFlpcpSparse};
use dpp::pipeline::build_rev2_dpp_sparse_boolean_auto;
use dpp::{EmbeddingParams, SparseVec};

#[derive(MontConfig)]
#[modulus = "5"]
#[generator = "2"]
pub struct F5Config;
type F5 = Fp64<MontBackend<F5Config, 1>>;

#[derive(MontConfig)]
// Goldilocks prime (2^64 - 2^32 + 1)
#[modulus = "18446744069414584321"]
#[generator = "7"]
pub struct GoldilocksConfig;
type Goldilocks = Fp64<MontBackend<GoldilocksConfig, 1>>;

fn lift_small_to_gold(x: F5) -> Goldilocks {
    Goldilocks::from_le_bytes_mod_order(&x.into_bigint().to_bytes_le())
}

#[test]
fn test_tiny_field_accepting_set_membership_roundtrip() {
    // Tiny dR1CS over F5 with one constraint: z0 * z1 = z2.
    // Public: z0. Witness: (z1, z2).
    let n_total = 3usize;
    let a_row = SparseVec::new(vec![(F5::ONE, 0)]);
    let b_row = SparseVec::new(vec![(F5::ONE, 1)]);
    let c_row = SparseVec::new(vec![(F5::ONE, 2)]);
    let inst = Dr1csInstanceSparse::<F5> {
        n: n_total,
        a: vec![a_row],
        b: vec![b_row],
        c: vec![c_row],
    };

    let l_public = 1usize;
    let k_rows = inst.k();
    let ell = 2 * k_rows;
    assert!(ell <= 5, "F5 requires ell <= |F|; keep this test tiny");

    let flpcp = RsDr1csNpFlpcpSparse::<F5>::new(inst, l_public, ell);

    // Satisfying assignment in F5.
    let z0 = F5::from(2u64);
    let z1 = F5::from(3u64);
    let z2 = z0 * z1; // = 1 in F5
    let x_small = vec![z0];
    let z_w_small = vec![z1, z2];

    // Underlying NP proof π_field = (z_w || w).
    let pi_field_small = flpcp.prove(&x_small, &z_w_small);

    // Booleanize π_field into π_bits (Rev2 uses this proof space).
    let boolized = BooleanProofFlpcpSparse::<F5, _>::new(flpcp.clone());
    let pi_bits_small = boolized.encode_proof_bits(&pi_field_small);

    // Embed instance/proof into Goldilocks.
    let x_gold = x_small.iter().copied().map(lift_small_to_gold).collect::<Vec<_>>();
    let pi_bits_gold = pi_bits_small
        .iter()
        .map(|b| if b.is_zero() { Goldilocks::ZERO } else { Goldilocks::ONE })
        .collect::<Vec<_>>();

    let dppv = build_rev2_dpp_sparse_boolean_auto::<F5, Goldilocks, _>(
        flpcp.clone(),
        EmbeddingParams {
            gamma: 2,
            assume_boolean_proof: true,
            k_prime: 0,
        },
    )
    .expect("build_rev2_dpp_sparse_boolean_auto");

    // Sample a packed query and enumerate the induced accepting set A.
    let mut rng = ChaCha20Rng::seed_from_u64(12345);
    let query = dppv.sample_query(&mut rng, &x_gold).expect("sample_query");
    let a_set = accepting_set_for_packed_query_sparse::<Goldilocks>(&query, 1_000_000)
        .expect("accepting_set_for_packed_query_sparse");
    assert!(!a_set.is_empty());

    // Compute packed answer a = <q, (x||π)> and check membership in A.
    let a = query.q.dot_two_slices(&x_gold, &pi_bits_gold);
    assert!(
        a_set.binary_search(&a).is_ok(),
        "packed answer must lie in explicit accepting set"
    );
    let ok = dppv
        .verify_packed_answer(&a, &query)
        .expect("verify_packed_answer");
    assert!(ok);

    // Negative: flip a proof bit; with overwhelming probability it should leave A.
    let mut bad_pi = pi_bits_gold.clone();
    bad_pi[0] = if bad_pi[0].is_zero() { Goldilocks::ONE } else { Goldilocks::ZERO };
    let a_bad = query.q.dot_two_slices(&x_gold, &bad_pi);
    let in_a = a_set.binary_search(&a_bad).is_ok();
    let ok_bad = dppv
        .verify_packed_answer(&a_bad, &query)
        .expect("verify_packed_answer");
    assert!(!ok_bad);
    // Membership failure is the stronger check; allow rare collision but still require verifier reject.
    if in_a {
        // If this happens, it means the accepting set is too coarse for this tiny instance;
        // keep the soundness check above as the authoritative condition.
        eprintln!("warning: flipped proof still landed in A (rare for tiny instances)");
    }
}

