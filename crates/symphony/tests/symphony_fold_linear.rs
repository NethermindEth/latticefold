

use symphony::symphony_fold::{fold_batchlin, fold_instances, SymphonyBatchLin, SymphonyInstance};
use stark_rings::{cyclotomic_ring::models::frog_ring::RqPoly as R, PolyRing, Ring};

#[test]
fn test_linear_folding_matches_manual_combination() {
    // This is a unit test for the *linear folding operators* (Figure 4, Step 5),
    // independent of Π_rg internals. We just need two consistent instances with shared randomness.
    let r = vec![
        <R as PolyRing>::BaseRing::from(11u128),
        <R as PolyRing>::BaseRing::from(22u128),
    ];
    let r_prime = vec![
        <R as PolyRing>::BaseRing::from(11u128),
        <R as PolyRing>::BaseRing::from(22u128),
        <R as PolyRing>::BaseRing::from(33u128),
    ];
    let kappa = 3;
    let k_g = 4;

    let inst0 = SymphonyInstance { c: vec![R::from(1u128), R::from(2u128)], r: r.clone(), v: R::from(7u128) };
    let inst1 = SymphonyInstance { c: vec![R::from(3u128), R::from(5u128)], r: r.clone(), v: R::from(9u128) };

    let bat0 = SymphonyBatchLin {
        r_prime: r_prime.clone(),
        c_g: (0..k_g).map(|i| vec![R::from((10 + i) as u128); kappa]).collect(),
        u: (0..k_g).map(|i| R::from((100 + i) as u128)).collect(),
    };
    let bat1 = SymphonyBatchLin {
        r_prime: r_prime.clone(),
        c_g: (0..k_g).map(|i| vec![R::from((20 + i) as u128); kappa]).collect(),
        u: (0..k_g).map(|i| R::from((200 + i) as u128)).collect(),
    };

    // In Symphony, β is sampled from a low-norm subset S ⊆ R_q.
    // For this unit test we just use small *constant* ring elements.
    let beta = vec![R::from(3u128), R::from(5u128)];

    let folded_inst = fold_instances(&beta, &[inst0.clone(), inst1.clone()]);
    let folded_bat = fold_batchlin(&beta, &[bat0.clone(), bat1.clone()]);

    // Manual: (c,r,v) fold.
    let b0_r = beta[0];
    let b1_r = beta[1];

    let mut c = vec![R::ZERO; inst0.c.len()];
    for i in 0..c.len() {
        c[i] = b0_r * inst0.c[i] + b1_r * inst1.c[i];
    }
    let r = inst0.r.clone(); // shared randomness (must match)
    let v = b0_r * inst0.v + b1_r * inst1.v;

    assert_eq!(folded_inst, SymphonyInstance { c, r, v });

    // Manual: batchlin fold (both c_g and u).
    assert_eq!(folded_bat.r_prime, bat0.r_prime);
    let mut u = vec![R::ZERO; bat0.u.len()];
    for i in 0..u.len() {
        u[i] = b0_r * bat0.u[i] + b1_r * bat1.u[i];
    }
    let mut c_g = vec![vec![R::ZERO; kappa]; k_g];
    for dig in 0..k_g {
        for j in 0..kappa {
            c_g[dig][j] = b0_r * bat0.c_g[dig][j] + b1_r * bat1.c_g[dig][j];
        }
    }
    assert_eq!(folded_bat, SymphonyBatchLin { r_prime: bat0.r_prime, c_g, u });
}

