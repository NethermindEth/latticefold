use ark_ff::{Field, Fp64, MontBackend, MontConfig};
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;

use dpp::accepting_set::{accepting_set_for_packed_query_sparse, AcceptingSetError};
use dpp::dr1cs_flpcp::{Dr1csInstanceSparse, RsDr1csNpFlpcpSparse};
use dpp::pipeline::build_rev2_dpp_sparse_boolean_auto;
use dpp::{EmbeddingParams, SparseVec};

#[derive(MontConfig)]
// Goldilocks prime (2^64 - 2^32 + 1)
#[modulus = "18446744069414584321"]
#[generator = "7"]
pub struct GoldilocksConfig;
type F = Fp64<MontBackend<GoldilocksConfig, 1>>;

#[derive(MontConfig)]
// NIST P-384 prime (same as `we_dpp` bench uses as FBig)
#[modulus = "39402006196394479212279040100143613805079739270465446667948293404245721771496870329047266088258938001861606973112319"]
#[generator = "2"]
pub struct Secp384r1Config;
type FBig = ark_ff::Fp384<MontBackend<Secp384r1Config, 6>>;

/// This is the “realistic next step” test: show that for non-tiny instances, the induced
/// accepting set of the packed-predicate DPP is astronomically large, so explicit enumeration
/// (Case A) is impossible.
///
/// We keep the instance size moderate so this test runs fast: it should fail early with
/// `DomainTooLarge` before doing any heavy work.
#[test]
fn test_packed_predicate_accepting_set_is_too_large_for_realistic_sizes() {
    // Choose k_rows large enough that the bound b explodes, but still small enough that
    // constructing the instance is cheap.
    let k_rows = 64usize;
    let n_total = 256usize;
    let l_public = 8usize;
    assert!(l_public < n_total);

    // Use a sparse dr1cs instance with trivial structure:
    // (z_i) * (z_{i+1}) = z_{i+2} on each row (wrapping indices).
    let mut a = Vec::with_capacity(k_rows);
    let mut b = Vec::with_capacity(k_rows);
    let mut c = Vec::with_capacity(k_rows);
    for i in 0..k_rows {
        let i0 = i % n_total;
        let i1 = (i + 1) % n_total;
        let i2 = (i + 2) % n_total;
        a.push(SparseVec::new(vec![(F::ONE, i0)]));
        b.push(SparseVec::new(vec![(F::ONE, i1)]));
        c.push(SparseVec::new(vec![(F::ONE, i2)]));
    }
    let inst = Dr1csInstanceSparse::<F> { n: n_total, a, b, c };
    let ell = 2 * k_rows;

    // Dummy statement x (public prefix of z). It does not need to be satisfying for this test.
    let x = vec![FBig::ZERO; l_public];

    let flpcp = RsDr1csNpFlpcpSparse::<F>::new(inst, l_public, ell);
    let dppv = build_rev2_dpp_sparse_boolean_auto::<F, FBig, _>(
        flpcp,
        EmbeddingParams {
            gamma: 2,
            assume_boolean_proof: true,
            k_prime: 0,
        },
    )
    .expect("build_rev2_dpp_sparse_boolean_auto");

    // Sample a packed query (cheap relative to any enumeration).
    let mut rng = ChaCha20Rng::seed_from_u64(1);
    let query = dppv.sample_query(&mut rng, &x).expect("sample_query");

    // Attempt to enumerate with a very small limit: this must fail with DomainTooLarge.
    let res = accepting_set_for_packed_query_sparse::<FBig>(&query, 1_000_000);
    assert!(matches!(res, Err(AcceptingSetError::DomainTooLarge { .. })));
}

