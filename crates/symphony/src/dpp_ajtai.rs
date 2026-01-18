//! Ajtai opening-check arithmetization (prime-field sparse dR1CS).
//!
//! This module provides a tiny building block needed for the WE/DPP gate:
//! encode `AjtaiOpen(commitment, message)` as *linear* constraints over the Poseidon
//! sponge's base prime field.
//!
//! Note: we intentionally target **small** committed messages (e.g. `cfs_*` transcript
//! message commitments). Arithmetizing `cm_f` directly would require including the full
//! witness `f` in the DPP witness, which is not the production Architecture‑T shape.

use ark_ff::Field;
use latticefold::commitment::AjtaiCommitmentScheme;
use rayon::prelude::*;
use stark_rings::{OverField, PolyRing, Ring};

use crate::dpp_poseidon::{Constraint, SparseDr1csInstance};

/// Build a sparse dR1CS instance that enforces `AjtaiCommitmentScheme::commit(message) == commitment`.
///
/// - Variables are **just** the scalar coefficients `message[j].coeffs()[0]` (plus the constant 1 slot).
/// - Constraints are linear: for each commitment row `i` and ring coefficient lane `ℓ`,
///   \sum_j A[i,j,ℓ] * msg[j] = commitment[i].coeffs()[ℓ]
///
/// Returns `(dr1cs, assignment)` where `assignment[0]=1` and `assignment[1+j]=msg[j]`.
pub fn ajtai_open_dr1cs_from_scheme<R>(
    scheme: &AjtaiCommitmentScheme<R>,
    message: &[R],
    commitment: &[R],
) -> Result<
    (
        SparseDr1csInstance<<<R as PolyRing>::BaseRing as Field>::BasePrimeField>,
        Vec<<<R as PolyRing>::BaseRing as Field>::BasePrimeField>,
    ),
    String,
>
where
    R: OverField + Ring,
    R::BaseRing: Field,
{
    let n = scheme.width();
    let kappa = scheme.kappa();
    if message.len() != n {
        return Err(format!("Ajtai dR1CS: message length mismatch (got {}, expected {n})", message.len()));
    }
    if commitment.len() != kappa {
        return Err(format!(
            "Ajtai dR1CS: commitment length mismatch (got {}, expected {kappa})",
            commitment.len()
        ));
    }

    // We only support BaseRing extension degree 1 for now (true for Frog base field regime).
    if R::BaseRing::extension_degree() != 1 {
        return Err("Ajtai dR1CS: BaseRing extension_degree != 1 not supported".to_string());
    }

    let d = R::dimension();

    // Assignment variables: [1, msg_0, ..., msg_{n-1}] in BF.
    let mut assignment: Vec<<<R as PolyRing>::BaseRing as Field>::BasePrimeField> =
        Vec::with_capacity(1 + n);
    assignment.push(<<<R as PolyRing>::BaseRing as Field>::BasePrimeField as Field>::ONE); // const slot
    for x in message {
        let s = x.coeffs()[0];
        let fp = s
            .to_base_prime_field_elements()
            .into_iter()
            .next()
            .ok_or_else(|| "Ajtai dR1CS: BaseRing to_base_prime_field_elements empty".to_string())?;
        assignment.push(fp);
    }

    // Precompute A columns by committing to basis vectors.
    // Parallelized: each column is independent.
    let cols: Vec<Vec<R>> = (0..n)
        .into_par_iter()
        .map(|j| {
        let mut basis = vec![R::ZERO; n];
        basis[j] = R::ONE;
        let col = scheme
            .commit(&basis)
                .expect("Ajtai dR1CS: commit(basis) failed")
            .as_ref()
            .to_vec();
        debug_assert_eq!(col.len(), kappa);
            col
        })
        .collect();

    // Build linear constraints.
    // Variable indices: 0 is 1, and 1+j is msg[j].
    let mut constraints: Vec<Constraint<<<R as PolyRing>::BaseRing as Field>::BasePrimeField>> =
        Vec::with_capacity(kappa * d);
    for i in 0..kappa {
        for lane in 0..d {
            let mut a_lc: Vec<(<<R as PolyRing>::BaseRing as Field>::BasePrimeField, usize)> =
                Vec::with_capacity(n + 1);
            for j in 0..n {
                let aij_lane = cols[j][i].coeffs()[lane];
                let fp = aij_lane
                    .to_base_prime_field_elements()
                    .into_iter()
                    .next()
                    .ok_or_else(|| "Ajtai dR1CS: aij to_base_prime_field_elements empty".to_string())?;
                let coeff = fp;
                if coeff != <<R as PolyRing>::BaseRing as Field>::BasePrimeField::ZERO {
                    a_lc.push((coeff, 1 + j));
                }
            }

            let rhs = commitment[i].coeffs()[lane];
            let fp = rhs
                .to_base_prime_field_elements()
                .into_iter()
                .next()
                .ok_or_else(|| "Ajtai dR1CS: rhs to_base_prime_field_elements empty".to_string())?;
            let rhs_bf = fp;
            // Move RHS to LHS: ... - rhs = 0
            if rhs_bf != <<R as PolyRing>::BaseRing as Field>::BasePrimeField::ZERO {
                a_lc.push((-rhs_bf, 0));
            }

            constraints.push(Constraint {
                a: a_lc,
                b: vec![(<<R as PolyRing>::BaseRing as Field>::BasePrimeField::ONE, 0)], // * 1
                c: vec![(<<R as PolyRing>::BaseRing as Field>::BasePrimeField::ZERO, 0)], // = 0
            });
        }
    }

    Ok((SparseDr1csInstance { nvars: assignment.len(), constraints }, assignment))
}

/// Same as `ajtai_open_dr1cs_from_scheme`, but supports **general ring elements** in `message`.
///
/// This treats every coefficient of every message ring element as a separate BF variable.
/// The Ajtai commitment is linear over the base field, because multiplication by a fixed
/// ring element is a linear map on coefficients.
///
/// This is still intended only for **small** messages (e.g. `cfs_mon_b`), since it requires
/// `n * d` basis commitments.
pub fn ajtai_open_dr1cs_from_scheme_full<R>(
    scheme: &AjtaiCommitmentScheme<R>,
    message: &[R],
    commitment: &[R],
) -> Result<
    (
        SparseDr1csInstance<<<R as PolyRing>::BaseRing as Field>::BasePrimeField>,
        Vec<<<R as PolyRing>::BaseRing as Field>::BasePrimeField>,
    ),
    String,
>
where
    R: OverField + Ring,
    R::BaseRing: Field,
{
    let n = scheme.width();
    let kappa = scheme.kappa();
    if message.len() != n {
        return Err(format!("Ajtai dR1CS(full): message length mismatch (got {}, expected {n})", message.len()));
    }
    if commitment.len() != kappa {
        return Err(format!(
            "Ajtai dR1CS(full): commitment length mismatch (got {}, expected {kappa})",
            commitment.len()
        ));
    }
    if R::BaseRing::extension_degree() != 1 {
        return Err("Ajtai dR1CS(full): BaseRing extension_degree != 1 not supported".to_string());
    }

    let d = R::dimension();

    // Assignment: [1] || flattened message coeffs (n * d).
    let mut assignment: Vec<<<R as PolyRing>::BaseRing as Field>::BasePrimeField> =
        Vec::with_capacity(1 + n * d);
    assignment.push(<<R as PolyRing>::BaseRing as Field>::BasePrimeField::ONE);
    for x in message {
        for lane in 0..d {
            let s = x.coeffs()[lane];
            let fp = s
                .to_base_prime_field_elements()
                .into_iter()
                .next()
                .ok_or_else(|| "Ajtai dR1CS(full): BaseRing to_base_prime_field_elements empty".to_string())?;
            assignment.push(fp);
        }
    }

    // Precompute columns for each (j,lane) basis coefficient.
    // Parallelized: each (j,lane) combination is independent.
    let cols: Vec<Vec<R>> = (0..(n * d))
        .into_par_iter()
        .map(|idx| {
            let j = idx / d;
            let lane = idx % d;
            let mut coeffs = vec![<<R as PolyRing>::BaseRing as Field>::ZERO; d];
            coeffs[lane] = <<R as PolyRing>::BaseRing as Field>::ONE;
            let mono_vec = R::promote_from_coeffs(coeffs)
                .expect("Ajtai dR1CS(full): promote_from_coeffs failed");
            assert_eq!(mono_vec.len(), 1, "Ajtai dR1CS(full): expected single ring element");
            let mono = mono_vec[0];
            let mut basis = vec![R::ZERO; n];
            basis[j] = mono;
            let col = scheme
                .commit(&basis)
                .expect("Ajtai dR1CS(full): commit(basis) failed")
                .as_ref()
                .to_vec();
            debug_assert_eq!(col.len(), kappa);
            col
        })
        .collect();

    // Constraints: for each commitment row i and output coeff lane_out:
    // Σ_{j,lane_in} coeff(j,lane_in -> lane_out) * msg[j][lane_in] = commitment[i][lane_out]
    let mut constraints: Vec<Constraint<<<R as PolyRing>::BaseRing as Field>::BasePrimeField>> =
        Vec::with_capacity(kappa * d);
    for i in 0..kappa {
        for lane_out in 0..d {
            let mut a_lc: Vec<(<<R as PolyRing>::BaseRing as Field>::BasePrimeField, usize)> = Vec::new();
            for j in 0..n {
                for lane_in in 0..d {
                    let col_idx = j * d + lane_in;
                    let aij_lane = cols[col_idx][i].coeffs()[lane_out];
                    let fp = aij_lane
                        .to_base_prime_field_elements()
                        .into_iter()
                        .next()
                        .ok_or_else(|| "Ajtai dR1CS(full): aij to_base_prime_field_elements empty".to_string())?;
                    if fp != <<R as PolyRing>::BaseRing as Field>::BasePrimeField::ZERO {
                        let var_idx = 1 + j * d + lane_in;
                        a_lc.push((fp, var_idx));
                    }
                }
            }
            let rhs = commitment[i].coeffs()[lane_out];
            let fp = rhs
                .to_base_prime_field_elements()
                .into_iter()
                .next()
                .ok_or_else(|| "Ajtai dR1CS(full): rhs to_base_prime_field_elements empty".to_string())?;
            if fp != <<R as PolyRing>::BaseRing as Field>::BasePrimeField::ZERO {
                a_lc.push((-fp, 0));
            }

            constraints.push(Constraint {
                a: a_lc,
                b: vec![(<<R as PolyRing>::BaseRing as Field>::BasePrimeField::ONE, 0)],
                c: vec![(<<R as PolyRing>::BaseRing as Field>::BasePrimeField::ZERO, 0)],
            });
        }
    }

    Ok((SparseDr1csInstance { nvars: assignment.len(), constraints }, assignment))
}

#[cfg(test)]
mod tests {
    use super::*;
    use latticefold::commitment::AjtaiCommitmentScheme;
    use stark_rings::cyclotomic_ring::models::frog_ring::RqPoly as R;

    #[test]
    fn test_ajtai_open_dr1cs_satisfiable_small() {
        // Tiny scheme: kappa=4, n=6.
        const MASTER_SEED: [u8; 32] = *b"SYMPHONY_AJTAI_SEED_V1_000000000";
        let scheme = AjtaiCommitmentScheme::<R>::seeded(b"test", MASTER_SEED, 4, 6);

        let msg = (0..6)
            .map(|i| R::from(<R as PolyRing>::BaseRing::from((i + 1) as u128)))
            .collect::<Vec<_>>();
        let cm = scheme.commit(&msg).unwrap().as_ref().to_vec();

        let (inst, asg) = ajtai_open_dr1cs_from_scheme::<R>(&scheme, &msg, &cm).unwrap();
        inst.check(&asg).unwrap();

        // Sanity: wrong commitment should fail at least one constraint.
        let mut bad = cm.clone();
        bad[0] += R::ONE;
        let (inst2, asg2) = ajtai_open_dr1cs_from_scheme::<R>(&scheme, &msg, &bad).unwrap();
        assert!(inst2.check(&asg2).is_err());
    }
}

