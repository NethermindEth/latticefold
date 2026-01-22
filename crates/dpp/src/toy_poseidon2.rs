//! Toy Poseidon2-like permutation/sponge for *cost de-risking* in tiny-field prototypes.
//!
//! IMPORTANT:
//! - This is **not** a standard, audited Poseidon2 instantiation.
//! - The only goal is to have a realistic-ish algebraic permutation (S-box + linear layer)
//!   whose *constraint cost* is comparable in shape to Poseidon2.
//!
//! We keep it deterministic, parameterized, and cheap to embed in toy gates.

use ark_ff::PrimeField;

/// Parameters for the toy permutation.
#[derive(Clone, Debug)]
pub struct ToyPoseidon2Params {
    /// State width `t`.
    pub width: usize,
    /// Number of full rounds.
    pub full_rounds: usize,
    /// Number of partial rounds.
    pub partial_rounds: usize,
}

/// Simple op-count estimate (field ops) for the permutation.
#[derive(Clone, Debug, Default)]
pub struct ToyPoseidon2Cost {
    pub muls: u64,
    pub adds: u64,
}

/// A very small sponge built from the toy permutation.
#[derive(Clone, Debug)]
pub struct ToyPoseidon2Sponge<F: PrimeField> {
    params: ToyPoseidon2Params,
    state: Vec<F>,
    pos: usize,
    permutes: u64,
}

impl<F: PrimeField> ToyPoseidon2Sponge<F> {
    pub fn new(params: ToyPoseidon2Params) -> Self {
        assert!(params.width >= 2);
        Self {
            state: vec![F::ZERO; params.width],
            params,
            pos: 0,
            permutes: 0,
        }
    }

    /// Absorb a sequence of field elements (rate = width-1).
    pub fn absorb(&mut self, inputs: &[F]) {
        let rate = self.params.width - 1;
        for x in inputs {
            if self.pos == rate {
                permute_in_place(&self.params, &mut self.state);
                self.permutes += 1;
                self.pos = 0;
            }
            self.state[self.pos] += *x;
            self.pos += 1;
        }
    }

    /// Squeeze `n` field elements (rate = width-1).
    pub fn squeeze(&mut self, n: usize) -> Vec<F> {
        let rate = self.params.width - 1;
        let mut out = Vec::with_capacity(n);
        for _ in 0..n {
            if self.pos == rate {
                permute_in_place(&self.params, &mut self.state);
                self.permutes += 1;
                self.pos = 0;
            }
            out.push(self.state[self.pos]);
            self.pos += 1;
        }
        out
    }

    /// Number of permutation calls used so far.
    pub fn permute_count(&self) -> u64 {
        self.permutes
    }
}

/// Apply the toy permutation in place.
pub fn permute_in_place<F: PrimeField>(params: &ToyPoseidon2Params, state: &mut [F]) {
    assert_eq!(state.len(), params.width);
    let t = params.width;

    // Deterministic “round constants”: (round_index + lane_index + 1).
    let rc = |round: usize, lane: usize| -> F { F::from((round + lane + 1) as u64) };

    let rf = params.full_rounds;
    let rp = params.partial_rounds;
    let total = rf + rp;

    for r in 0..total {
        // Add round constants.
        for i in 0..t {
            state[i] += rc(r, i);
        }

        // S-box layer: x -> x^5.
        if r < rf / 2 || r >= rf / 2 + rp {
            // Full round: all lanes.
            for i in 0..t {
                state[i] = pow5(state[i]);
            }
        } else {
            // Partial round: first lane only (Poseidon2-like).
            state[0] = pow5(state[0]);
        }

        // Linear layer (toy MDS): y_i = Σ_j M_{i,j} * x_j with M_{i,j} = (i+1) + 2*(j+1).
        let old = state.to_vec();
        for i in 0..t {
            let mut acc = F::ZERO;
            for j in 0..t {
                let mij = F::from(((i + 1) + 2 * (j + 1)) as u64);
                acc += mij * old[j];
            }
            state[i] = acc;
        }
    }
}

/// Return a rough cost estimate (field adds/muls) per permutation call.
pub fn estimate_cost(params: &ToyPoseidon2Params) -> ToyPoseidon2Cost {
    let t = params.width as u64;
    let rf = params.full_rounds as u64;
    let rp = params.partial_rounds as u64;
    let total = rf + rp;

    // Add round constants: t adds per round.
    let adds_rc = total * t;

    // S-box x^5: 3 muls (x2=x*x, x4=x2*x2, x5=x4*x).
    // Full rounds: t S-boxes; partial: 1 S-box.
    let full_sboxes = rf * t;
    let partial_sboxes = rp * 1;
    let muls_sbox = 3 * (full_sboxes + partial_sboxes);

    // Linear layer: t^2 muls and t*(t-1) adds per round (naive matmul).
    let muls_lin = total * t * t;
    let adds_lin = total * t * (t - 1);

    ToyPoseidon2Cost {
        muls: muls_sbox + muls_lin,
        adds: adds_rc + adds_lin,
    }
}

#[inline]
fn pow5<F: PrimeField>(x: F) -> F {
    let x2 = x * x;
    let x4 = x2 * x2;
    x4 * x
}

