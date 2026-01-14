use cyclotomic_rings::rings::SuitableRing;
use stark_rings::{
    balanced_decomposition::DecomposeToVec,
    cyclotomic_ring::{CRT, ICRT},
    PolyRing,
    Ring,
};
use stark_rings_linalg::Matrix;
use rand_chacha::ChaCha20Rng;
use rand::SeedableRng;
use sha2::{Digest, Sha256};
#[cfg(feature = "parallel")]
use rayon::prelude::*;
#[cfg(feature = "parallel")]
use ark_std::cfg_into_iter;

use super::homomorphic_commitment::Commitment;
use crate::{
    ark_base::*, commitment::CommitmentError, decomposition_parameters::DecompositionParams,
};

#[derive(Clone, Debug)]
enum AjtaiMatrix<R> {
    Explicit(Matrix<R>),
    /// A pseudorandom Ajtai matrix defined implicitly by a seed.
    ///
    /// This avoids materializing a `kappa x n` dense matrix in memory.
    Seeded {
        seed: [u8; 32],
        domain: Vec<u8>,
        kappa: usize,
        n: usize,
    },
}

/// A concrete instantiation of the Ajtai commitment scheme.
/// Contains either an explicit Ajtai matrix (kappa x n) or a seeded implicit matrix.
#[derive(Clone, Debug)]
pub struct AjtaiCommitmentScheme<R> {
    matrix: AjtaiMatrix<R>,
}

impl<R> AjtaiCommitmentScheme<R> {
    /// Create a new scheme using the provided Ajtai matrix
    pub fn new(matrix: Matrix<R>) -> Self {
        Self { matrix: AjtaiMatrix::Explicit(matrix) }
    }
}

impl<R: Ring> AjtaiCommitmentScheme<R> {
    /// Returns a random Ajtai commitment matrix
    pub fn rand<Rng: rand::Rng + ?Sized>(kappa: usize, n: usize, rng: &mut Rng) -> Self {
        Self::new(vec![vec![R::rand(rng); n]; kappa].into())
    }

    /// Create a scheme with an implicitly-defined pseudorandom Ajtai matrix.
    ///
    /// The matrix is derived deterministically from `(domain, seed)`; different domains should use
    /// different labels and/or different seeds.
    pub fn seeded(domain: impl AsRef<[u8]>, seed: [u8; 32], kappa: usize, n: usize) -> Self {
        Self {
            matrix: AjtaiMatrix::Seeded {
                seed,
                domain: domain.as_ref().to_vec(),
                kappa,
                n,
            },
        }
    }

    /// Derive the seed for generating a specific Ajtai matrix column.
    ///
    /// This is a deterministic function of `(domain, seed, col)` used to generate
    /// column `col` of the Ajtai matrix on-the-fly.
    pub fn derive_col_seed(domain: &[u8], seed: &[u8; 32], col: u64) -> [u8; 32] {
        let mut h = Sha256::new();
        h.update(b"AJTAI_COL_V1");
        h.update((domain.len() as u64).to_le_bytes());
        h.update(domain);
        h.update(seed);
        h.update(col.to_le_bytes());
        let out = h.finalize();
        let mut s = [0u8; 32];
        s.copy_from_slice(&out);
        s
    }
}

impl<R: Ring> AjtaiCommitmentScheme<R> {
    /// Commit to a witness
    pub fn commit(&self, f: &[R]) -> Result<Commitment<R>, CommitmentError> {
        match &self.matrix {
            AjtaiMatrix::Explicit(matrix) => {
                if f.len() != matrix.ncols {
                    return Err(CommitmentError::WrongWitnessLength(f.len(), matrix.ncols));
                }
                let commitment = matrix
                    .checked_mul_vec(f)
                    .ok_or(CommitmentError::WrongWitnessLength(f.len(), matrix.ncols))?;
                Ok(Commitment::from_vec_raw(commitment))
            }
            AjtaiMatrix::Seeded { seed, domain, kappa, n } => {
                if f.len() != *n {
                    return Err(CommitmentError::WrongWitnessLength(f.len(), *n));
                }
                // Speed: pre-hash the (domain, seed) prefix once, then only hash `col` per column.
                // This preserves *exactly* the same `derive_col_seed` outputs.
                let prefix_hasher = {
                    let mut h = Sha256::new();
                    h.update(b"AJTAI_COL_V1");
                    h.update((domain.len() as u64).to_le_bytes());
                    h.update(domain);
                    h.update(seed);
                    h
                };
                #[cfg(feature = "parallel")]
                {
                    let kappa = *kappa;
                    let domain = domain.as_slice();
                    let seed = *seed;
                    let acc = cfg_into_iter!(0..f.len())
                        .fold(
                            || (vec![R::ZERO; kappa], prefix_hasher.clone()),
                            |(mut local, prefix), j| {
                                let fj = f[j];
                                if fj == R::ZERO {
                                    return (local, prefix);
                                }
                                // Equivalent to `derive_col_seed(domain, &seed, col=j)`, but avoids
                                // re-hashing the prefix for every column.
                                let col_seed = {
                                    let mut h = prefix.clone();
                                    h.update((j as u64).to_le_bytes());
                                    let out = h.finalize();
                                    let mut s = [0u8; 32];
                                    s.copy_from_slice(&out);
                                    s
                                };
                                let mut rng = ChaCha20Rng::from_seed(col_seed);
                                for i in 0..kappa {
                                    let aij = R::rand(&mut rng);
                                    local[i] += aij * fj;
                                }
                                (local, prefix)
                            },
                        )
                        .reduce(
                            || (vec![R::ZERO; kappa], prefix_hasher.clone()),
                            |(mut a, prefix_a), (b, _prefix_b)| {
                                for i in 0..kappa {
                                    a[i] += b[i];
                                }
                                (a, prefix_a)
                            },
                        );
                    Ok(Commitment::from_vec_raw(acc.0))
                }
                #[cfg(not(feature = "parallel"))]
                {
                    let mut acc = vec![R::ZERO; *kappa];
                    let prefix = prefix_hasher;
                    for (j, fj) in f.iter().enumerate() {
                        if *fj == R::ZERO {
                            continue;
                        }
                        let col_seed = {
                            let mut h = prefix.clone();
                            h.update((j as u64).to_le_bytes());
                            let out = h.finalize();
                            let mut s = [0u8; 32];
                            s.copy_from_slice(&out);
                            s
                        };
                        let mut rng = ChaCha20Rng::from_seed(col_seed);
                        for i in 0..*kappa {
                            let aij = R::rand(&mut rng);
                            acc[i] += aij * *fj;
                        }
                    }
                    Ok(Commitment::from_vec_raw(acc))
                }
            }
        }
    }

    
    /// Commit to `t` witnesses of length `n` provided implicitly.
    ///
    /// This computes `t` commitments under the same Ajtai matrix, reusing the per-column randomness
    /// generation for all `t` vectors.
    pub fn commit_many_with<F>(
        &self,
        n: usize,
        t: usize,
        fill_values_at: F,
    ) -> Result<Vec<Commitment<R>>, CommitmentError>
    where
        F: Fn(usize, &mut [R]) + Sync,
        R: Send + Sync + Copy,
    {
        if t == 0 {
            return Ok(vec![]);
        }
        match &self.matrix {
            AjtaiMatrix::Explicit(matrix) => {
                if n != matrix.ncols {
                    return Err(CommitmentError::WrongWitnessLength(n, matrix.ncols));
                }
                // Fallback: materialize each vector and call `commit()`.
                let mut outs = Vec::with_capacity(t);
                let mut tmp = vec![R::ZERO; n];
                for which in 0..t {
                    for j in 0..n {
                        // Fill a length-`t` scratch and take the `which` slot.
                        let mut scratch = vec![R::ZERO; t];
                        fill_values_at(j, &mut scratch);
                        tmp[j] = scratch[which];
                    }
                    outs.push(self.commit(&tmp)?);
                }
                Ok(outs)
            }
            AjtaiMatrix::Seeded { seed, domain, kappa, n: n_expected } => {
                if n != *n_expected {
                    return Err(CommitmentError::WrongWitnessLength(n, *n_expected));
                }
                let prefix_hasher = {
                    let mut h = Sha256::new();
                    h.update(b"AJTAI_COL_V1");
                    h.update((domain.len() as u64).to_le_bytes());
                    h.update(domain);
                    h.update(seed);
                    h
                };
                #[cfg(feature = "parallel")]
                {
                    let kappa = *kappa;
                    let domain = domain.as_slice();
                    let seed = *seed;
                    let acc = cfg_into_iter!(0..n)
                        .fold(
                            || (vec![vec![R::ZERO; kappa]; t], vec![R::ZERO; t], prefix_hasher.clone()),
                            |(mut local, mut scratch, prefix), j| {
                                scratch.fill(R::ZERO);
                                fill_values_at(j, &mut scratch);
                                if scratch.iter().all(|x| *x == R::ZERO) {
                                    return (local, scratch, prefix);
                                }
                                let col_seed = {
                                    let mut h = prefix.clone();
                                    h.update((j as u64).to_le_bytes());
                                    let out = h.finalize();
                                    let mut s = [0u8; 32];
                                    s.copy_from_slice(&out);
                                    s
                                };
                                let mut rng = ChaCha20Rng::from_seed(col_seed);
                                for i in 0..kappa {
                                    let aij = R::rand(&mut rng);
                                    for which in 0..t {
                                        let fj = scratch[which];
                                        if fj != R::ZERO {
                                            local[which][i] += aij * fj;
                                        }
                                    }
                                }
                                (local, scratch, prefix)
                            },
                        )
                        .reduce(
                            || (vec![vec![R::ZERO; kappa]; t], vec![R::ZERO; t], prefix_hasher.clone()),
                            |(mut a, scratch_a, prefix_a), (b, _scratch_b, _prefix_b)| {
                                for which in 0..t {
                                    for i in 0..kappa {
                                        a[which][i] += b[which][i];
                                    }
                                }
                                (a, scratch_a, prefix_a)
                            },
                        )
                        .0;
                    Ok(acc.into_iter().map(Commitment::from_vec_raw).collect())
                }
                #[cfg(not(feature = "parallel"))]
                {
                    let mut acc = vec![vec![R::ZERO; *kappa]; t];
                    let mut scratch = vec![R::ZERO; t];
                    let prefix = prefix_hasher;
                    for j in 0..n {
                        scratch.fill(R::ZERO);
                        fill_values_at(j, &mut scratch);
                        if scratch.iter().all(|x| *x == R::ZERO) {
                            continue;
                        }
                        let col_seed = {
                            let mut h = prefix.clone();
                            h.update((j as u64).to_le_bytes());
                            let out = h.finalize();
                            let mut s = [0u8; 32];
                            s.copy_from_slice(&out);
                            s
                        };
                        let mut rng = ChaCha20Rng::from_seed(col_seed);
                        for i in 0..*kappa {
                            let aij = R::rand(&mut rng);
                            for which in 0..t {
                                let fj = scratch[which];
                                if fj != R::ZERO {
                                    acc[which][i] += aij * fj;
                                }
                            }
                        }
                    }
                    Ok(acc.into_iter().map(Commitment::from_vec_raw).collect())
                }
            }
        }
    }

    /// Like `commit_many_with`, but only processes indices `j ∈ [start, end)`.
    ///
    /// This is useful when the witness is mostly constant/zero and the caller can represent the
    /// nontrivial part as a contiguous range.
    pub fn commit_many_with_range<F>(
        &self,
        n: usize,
        t: usize,
        start: usize,
        end: usize,
        fill_values_at: F,
    ) -> Result<Vec<Commitment<R>>, CommitmentError>
    where
        F: Fn(usize, &mut [R]) + Sync,
        R: Send + Sync + Copy,
    {
        if t == 0 {
            return Ok(vec![]);
        }
        if start > end || end > n {
            return Err(CommitmentError::WrongWitnessLength(end, n));
        }
        match &self.matrix {
            AjtaiMatrix::Explicit(_) => self.commit_many_with(n, t, fill_values_at),
            AjtaiMatrix::Seeded { seed, domain, kappa, n: n_expected } => {
                if n != *n_expected {
                    return Err(CommitmentError::WrongWitnessLength(n, *n_expected));
                }
                #[cfg(feature = "parallel")]
                {
                    let kappa = *kappa;
                    let domain = domain.as_slice();
                    let seed = *seed;
                    let acc = cfg_into_iter!(start..end)
                        .fold(
                            || (vec![vec![R::ZERO; kappa]; t], vec![R::ZERO; t]),
                            |(mut local, mut scratch), j| {
                                scratch.fill(R::ZERO);
                                fill_values_at(j, &mut scratch);
                                if scratch.iter().all(|x| *x == R::ZERO) {
                                    return (local, scratch);
                                }
                                let col_seed = Self::derive_col_seed(domain, &seed, j as u64);
                                let mut rng = ChaCha20Rng::from_seed(col_seed);
                                for i in 0..kappa {
                                    let aij = R::rand(&mut rng);
                                    for which in 0..t {
                                        let fj = scratch[which];
                                        if fj != R::ZERO {
                                            local[which][i] += aij * fj;
                                        }
                                    }
                                }
                                (local, scratch)
                            },
                        )
                        .reduce(
                            || (vec![vec![R::ZERO; kappa]; t], vec![R::ZERO; t]),
                            |(mut a, scratch_a), (b, _scratch_b)| {
                                for which in 0..t {
                                    for i in 0..kappa {
                                        a[which][i] += b[which][i];
                                    }
                                }
                                (a, scratch_a)
                            },
                        )
                        .0;
                    Ok(acc.into_iter().map(Commitment::from_vec_raw).collect())
                }
                #[cfg(not(feature = "parallel"))]
                {
                    let mut acc = vec![vec![R::ZERO; *kappa]; t];
                    let mut scratch = vec![R::ZERO; t];
                    for j in start..end {
                        scratch.fill(R::ZERO);
                        fill_values_at(j, &mut scratch);
                        if scratch.iter().all(|x| *x == R::ZERO) {
                            continue;
                        }
                        let col_seed = Self::derive_col_seed(domain, &seed, j as u64);
                        let mut rng = ChaCha20Rng::from_seed(col_seed);
                        for i in 0..*kappa {
                            let aij = R::rand(&mut rng);
                            for which in 0..t {
                                let fj = scratch[which];
                                if fj != R::ZERO {
                                    acc[which][i] += aij * fj;
                                }
                            }
                        }
                    }
                    Ok(acc.into_iter().map(Commitment::from_vec_raw).collect())
                }
            }
        }
    }

    /// Commit to a witness, using a fast-path when all witness entries are *constant-coefficient*
    /// ring elements (i.e., elements embedded from the base ring).
    ///
    /// This is a pure optimization: if any entry is not constant-coefficient, we fall back to
    /// `commit()` with identical outputs.
    pub fn commit_const_coeff_fast(&self, f: &[R]) -> Result<Commitment<R>, CommitmentError>
    where
        R: PolyRing,
        R::BaseRing: Ring,
        R: core::ops::Mul<R::BaseRing, Output = R>,
    {
        // Quick check: if any nonzero entry has a nonzero non-constant coefficient, fall back.
        for fj in f {
            if *fj == R::ZERO {
                continue;
            }
            let coeffs = fj.coeffs();
            if coeffs.iter().skip(1).any(|c| *c != R::BaseRing::ZERO) {
                return self.commit(f);
            }
        }

        match &self.matrix {
            AjtaiMatrix::Explicit(_) => self.commit(f),
            AjtaiMatrix::Seeded { seed, domain, kappa, n } => {
                if f.len() != *n {
                    return Err(CommitmentError::WrongWitnessLength(f.len(), *n));
                }
                #[cfg(feature = "parallel")]
                {
                    let kappa = *kappa;
                    let domain = domain.as_slice();
                    let seed = *seed;
                    let acc = cfg_into_iter!(0..f.len())
                        .fold(
                            || vec![R::ZERO; kappa],
                            |mut local, j| {
                                let fj0 = f[j].coeffs()[0];
                                if fj0 == R::BaseRing::ZERO {
                                    return local;
                                }
                                let col_seed = Self::derive_col_seed(domain, &seed, j as u64);
                                let mut rng = ChaCha20Rng::from_seed(col_seed);
                                for i in 0..kappa {
                                    let aij = R::rand(&mut rng);
                                    local[i] += aij * fj0;
                                }
                                local
                            },
                        )
                        .reduce(
                            || vec![R::ZERO; kappa],
                            |mut a, b| {
                                for i in 0..kappa {
                                    a[i] += b[i];
                                }
                                a
                            },
                        );
                    Ok(Commitment::from_vec_raw(acc))
                }
                #[cfg(not(feature = "parallel"))]
                {
                    let mut acc = vec![R::ZERO; *kappa];
                    for (j, fj) in f.iter().enumerate() {
                        let fj0 = fj.coeffs()[0];
                        if fj0 == R::BaseRing::ZERO {
                            continue;
                        }
                        let col_seed = Self::derive_col_seed(domain, &seed, j as u64);
                        let mut rng = ChaCha20Rng::from_seed(col_seed);
                        for i in 0..*kappa {
                            let aij = R::rand(&mut rng);
                            acc[i] += aij * fj0;
                        }
                    }
                    Ok(Commitment::from_vec_raw(acc))
                }
            }
        }
    }

    /// Ajtai matrix number of rows
    ///
    /// This value affects the security of the scheme.
    pub fn kappa(&self) -> usize {
        match &self.matrix {
            AjtaiMatrix::Explicit(m) => m.nrows,
            AjtaiMatrix::Seeded { kappa, .. } => *kappa,
        }
    }

    /// Ajtai matrix number of columns
    ///
    /// The size of the witness must be equal to this value.
    pub fn width(&self) -> usize {
        match &self.matrix {
            AjtaiMatrix::Explicit(m) => m.ncols,
            AjtaiMatrix::Seeded { n, .. } => *n,
        }
    }

    /// Commit to a **constant-coefficient** witness given directly as base scalars.
    ///
    /// This is equivalent to calling `commit_const_coeff_fast` on `f_ring[j] = R::from(f0[j])`,
    /// but avoids allocating the `Vec<R>` witness (critical for SP1-scale n=2^27).
    ///
    /// Semantics:
    /// - Treats the witness as ring elements embedded from `R::BaseRing` (i.e. only the constant term).
    /// - For explicit Ajtai matrices, falls back to `commit()` via embedding (since the explicit path is
    ///   not the performance bottleneck in our large-n seeded setting).
    pub fn commit_const_coeff_base_fast(
        &self,
        f0: &[<R as PolyRing>::BaseRing],
    ) -> Result<Commitment<R>, CommitmentError>
    where
        R: PolyRing,
        R::BaseRing: Ring,
        R: core::ops::Mul<R::BaseRing, Output = R>,
    {
        match &self.matrix {
            AjtaiMatrix::Explicit(m) => {
                // Implement the explicit matrix case without requiring `R: From<R::BaseRing>`:
                // treat each base scalar as a constant-coeff ring element and multiply via `aij * fj0`.
                if f0.len() != m.ncols {
                    return Err(CommitmentError::WrongWitnessLength(f0.len(), m.ncols));
                }
                let kappa = m.nrows;
                let mut acc = vec![R::ZERO; kappa];
                for j in 0..m.ncols {
                    let fj0 = f0[j];
                    if fj0 == R::BaseRing::ZERO {
                        continue;
                    }
                    for i in 0..kappa {
                        acc[i] += m.vals[i][j] * fj0;
                    }
                }
                Ok(Commitment::from_vec_raw(acc))
            }
            AjtaiMatrix::Seeded { seed, domain, kappa, n } => {
                if f0.len() != *n {
                    return Err(CommitmentError::WrongWitnessLength(f0.len(), *n));
                }
                #[cfg(feature = "parallel")]
                {
                    let kappa = *kappa;
                    let domain = domain.as_slice();
                    let seed = *seed;
                    let acc = cfg_into_iter!(0..f0.len())
                        .fold(
                            || vec![R::ZERO; kappa],
                            |mut local, j| {
                                let fj0 = f0[j];
                                if fj0 == R::BaseRing::ZERO {
                                    return local;
                                }
                                let col_seed = Self::derive_col_seed(domain, &seed, j as u64);
                                let mut rng = ChaCha20Rng::from_seed(col_seed);
                                for i in 0..kappa {
                                    let aij = R::rand(&mut rng);
                                    local[i] += aij * fj0;
                                }
                                local
                            },
                        )
                        .reduce(
                            || vec![R::ZERO; kappa],
                            |mut a, b| {
                                for i in 0..kappa {
                                    a[i] += b[i];
                                }
                                a
                            },
                        );
                    Ok(Commitment::from_vec_raw(acc))
                }
                #[cfg(not(feature = "parallel"))]
                {
                    let mut acc = vec![R::ZERO; *kappa];
                    for (j, &fj0) in f0.iter().enumerate() {
                        if fj0 == R::BaseRing::ZERO {
                            continue;
                        }
                        let col_seed = Self::derive_col_seed(domain, &seed, j as u64);
                        let mut rng = ChaCha20Rng::from_seed(col_seed);
                        for i in 0..*kappa {
                            let aij = R::rand(&mut rng);
                            acc[i] += aij * fj0;
                        }
                    }
                    Ok(Commitment::from_vec_raw(acc))
                }
            }
        }
    }

    /// Batched variant of `commit_const_coeff_base_fast`.
    ///
    /// Computes `t` Ajtai commitments to *constant-coefficient* witnesses of length `n`,
    /// provided on-the-fly as base scalars by `fill_values_at(j, out)`, where `out.len() == t`
    /// and `out[w]` is the base scalar at column `j` for witness `w`.
    ///
    /// This matches Symphony’s streaming pattern: one pass over columns, reuse the same Ajtai column
    /// RNG sample `a_ij` across many outputs. It avoids allocating `t` separate `Vec<R>` witnesses.
    pub fn commit_many_const_coeff_base_fast<F>(
        &self,
        n: usize,
        t: usize,
        fill_values_at: F,
    ) -> Result<Vec<Commitment<R>>, CommitmentError>
    where
        R: PolyRing,
        R::BaseRing: Ring,
        R: core::ops::Mul<R::BaseRing, Output = R>,
        F: Fn(usize, &mut [R::BaseRing]) + Send + Sync,
    {
        if t == 0 {
            return Ok(Vec::new());
        }
        match &self.matrix {
            AjtaiMatrix::Explicit(m) => {
                if n != m.ncols {
                    return Err(CommitmentError::WrongWitnessLength(n, m.ncols));
                }
                let kappa = m.nrows;
                let mut acc = vec![vec![R::ZERO; kappa]; t];
                let mut scratch = vec![R::BaseRing::ZERO; t];
                for j in 0..n {
                    scratch.fill(R::BaseRing::ZERO);
                    fill_values_at(j, &mut scratch);
                    if scratch.iter().all(|x| *x == R::BaseRing::ZERO) {
                        continue;
                    }
                    for i in 0..kappa {
                        let aij = m.vals[i][j];
                        for which in 0..t {
                            let fj0 = scratch[which];
                            if fj0 != R::BaseRing::ZERO {
                                acc[which][i] += aij * fj0;
                            }
                        }
                    }
                }
                Ok(acc.into_iter().map(Commitment::from_vec_raw).collect())
            }
            AjtaiMatrix::Seeded { seed, domain, kappa, n: nn } => {
                if n != *nn {
                    return Err(CommitmentError::WrongWitnessLength(n, *nn));
                }
                let kappa = *kappa;
                let domain = domain.as_slice();
                let seed = *seed;
                #[cfg(feature = "parallel")]
                {
                    let acc = cfg_into_iter!(0..n)
                        .fold(
                            || (vec![vec![R::ZERO; kappa]; t], vec![R::BaseRing::ZERO; t]),
                            |(mut local, mut scratch), j| {
                                scratch.fill(R::BaseRing::ZERO);
                                fill_values_at(j, &mut scratch);
                                if scratch.iter().all(|x| *x == R::BaseRing::ZERO) {
                                    return (local, scratch);
                                }
                                let col_seed = Self::derive_col_seed(domain, &seed, j as u64);
                                let mut rng = ChaCha20Rng::from_seed(col_seed);
                                for i in 0..kappa {
                                    let aij = R::rand(&mut rng);
                                    for which in 0..t {
                                        let fj0 = scratch[which];
                                        if fj0 != R::BaseRing::ZERO {
                                            local[which][i] += aij * fj0;
                                        }
                                    }
                                }
                                (local, scratch)
                            },
                        )
                        .reduce(
                            || (vec![vec![R::ZERO; kappa]; t], vec![R::BaseRing::ZERO; t]),
                            |(mut a, scratch_a), (b, _scratch_b)| {
                                for which in 0..t {
                                    for i in 0..kappa {
                                        a[which][i] += b[which][i];
                                    }
                                }
                                (a, scratch_a)
                            },
                        )
                        .0;
                    Ok(acc.into_iter().map(Commitment::from_vec_raw).collect())
                }
                #[cfg(not(feature = "parallel"))]
                {
                    let mut acc = vec![vec![R::ZERO; kappa]; t];
                    let mut scratch = vec![R::BaseRing::ZERO; t];
                    for j in 0..n {
                        scratch.fill(R::BaseRing::ZERO);
                        fill_values_at(j, &mut scratch);
                        if scratch.iter().all(|x| *x == R::BaseRing::ZERO) {
                            continue;
                        }
                        let col_seed = Self::derive_col_seed(domain, &seed, j as u64);
                        let mut rng = ChaCha20Rng::from_seed(col_seed);
                        for i in 0..kappa {
                            let aij = R::rand(&mut rng);
                            for which in 0..t {
                                let fj0 = scratch[which];
                                if fj0 != R::BaseRing::ZERO {
                                    acc[which][i] += aij * fj0;
                                }
                            }
                        }
                    }
                    Ok(acc.into_iter().map(Commitment::from_vec_raw).collect())
                }
            }
        }
    }
}

// SuitableRing helpers
impl<NTT: SuitableRing> AjtaiCommitmentScheme<NTT> {
    /// Commit to a witness in the NTT form.
    /// The most basic one just multiplies by the matrix.
    pub fn commit_ntt(&self, f: &[NTT]) -> Result<Commitment<NTT>, CommitmentError> {
        self.commit(f)
    }

    /// Commit to a witness in the coefficient form.
    /// Performs NTT on each component of the witness and then does Ajtai commitment.
    pub fn commit_coeff<P: DecompositionParams>(
        &self,
        f: Vec<NTT::CoefficientRepresentation>,
    ) -> Result<Commitment<NTT>, CommitmentError> {
        self.commit_ntt(&CRT::elementwise_crt(f))
    }

    /// Takes a coefficient form witness, decomposes it vertically in radix-B,
    /// i.e. computes a preimage G_B^{-1}(w), and Ajtai commits to the result.
    pub fn decompose_and_commit_coeff<P: DecompositionParams>(
        &self,
        f: &[NTT::CoefficientRepresentation],
    ) -> Result<Commitment<NTT>, CommitmentError> {
        let f = f
            .decompose_to_vec(P::B, P::L)
            .into_iter()
            .flatten()
            .collect::<Vec<_>>();

        self.commit_coeff::<P>(f)
    }

    /// Takes an NTT form witness, transforms it into the coefficient form,
    /// decomposes it vertically in radix-B, i.e.
    /// computes a preimage G_B^{-1}(w), and Ajtai commits to the result.
    pub fn decompose_and_commit_ntt<P: DecompositionParams>(
        &self,
        w: Vec<NTT>,
    ) -> Result<Commitment<NTT>, CommitmentError> {
        let coeff: Vec<NTT::CoefficientRepresentation> = ICRT::elementwise_icrt(w);

        self.decompose_and_commit_coeff::<P>(&coeff)
    }
}

#[cfg(test)]
mod tests {
    use cyclotomic_rings::rings::GoldilocksRingNTT;
    use stark_rings::OverField;

    use super::{AjtaiCommitmentScheme, CommitmentError};
    use crate::ark_base::*;

    pub(crate) fn generate_ajtai<NTT: OverField>(
        kappa: usize,
        n: usize,
    ) -> AjtaiCommitmentScheme<NTT> {
        let mut matrix = Vec::<Vec<NTT>>::new();

        for i in 0..kappa {
            let mut row = Vec::<NTT>::new();
            for j in 0..n {
                row.push(NTT::from((i * n + j) as u128));
            }
            matrix.push(row)
        }

        AjtaiCommitmentScheme::new(matrix.into())
    }

    #[test]
    fn test_commit_ntt() -> Result<(), CommitmentError> {
        const WITNESS_SIZE: usize = 1 << 15;
        const OUTPUT_SIZE: usize = 9;

        let ajtai_data: AjtaiCommitmentScheme<GoldilocksRingNTT> =
            generate_ajtai(OUTPUT_SIZE, WITNESS_SIZE);
        let witness: Vec<_> = (0..(1 << 15)).map(|_| 2_u128.into()).collect();

        let committed = ajtai_data.commit_ntt(&witness)?;

        for (i, &x) in committed.as_ref().iter().enumerate() {
            let expected: u128 =
                ((WITNESS_SIZE) * (2 * i * WITNESS_SIZE + (WITNESS_SIZE - 1))) as u128;
            assert_eq!(x, expected.into());
        }

        Ok(())
    }
}
