use crate::{
    poly_utils::{ MultiPoly, UnivPoly },
    transcript::{ poseidon::PoseidonTranscript, Transcript },
};
use lattirust_arithmetic::polynomials::VirtualPolynomial;
use ark_crypto_primitives::sponge::{ poseidon::PoseidonConfig, Absorb };
use ark_ff::PrimeField;
use lattirust_arithmetic::challenge_set::latticefold_challenge_set::{
    LatticefoldChallengeSet,
    OverField,
};

pub struct SumCheckIP<F: PrimeField, R: OverField<F>, CS: LatticefoldChallengeSet<F, R>>
    where F: Absorb {
    pub claimed_sum: R,
    pub polynomial: VirtualPolynomial<R>,
    pub rounds: Vec<SumCheckRound<F, R>>,
    pub hasher: PoseidonTranscript<F, R, CS>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SumCheckRound<F: PrimeField, R: OverField<F>> {
    _marker: std::marker::PhantomData<F>,
    pub challenge: R,
    var_index: usize,
    pub unipoly: VirtualPolynomial<R>,
}

impl<F: PrimeField, R: OverField<F>, CS: LatticefoldChallengeSet<F, R>> SumCheckIP<F, R, CS>
    where F: Absorb
{
    pub fn new(
        claimed_sum: R,
        polynomial: VirtualPolynomial<R>,
        num_rounds: &usize
    ) -> SumCheckIP<F, R, CS> {
        let config = PoseidonConfig {
            full_rounds: 8, // Example values, adjust according to your needs
            partial_rounds: 57,
            alpha: 5,
            ark: vec![vec![F::zero(); 3]; 8 + 57], // Adjust to actual ark parameters
            mds: vec![vec![F::zero(); 3]; 3], // Adjust to actual MDS matrix parameters
            rate: 2,
            capacity: 1,
        };

        SumCheckIP {
            claimed_sum,
            polynomial,
            rounds: Vec::with_capacity(*num_rounds),
            hasher: PoseidonTranscript::new(&config),
        }
    }

    pub fn add_round(&mut self, challenge: R, var_index: usize, unipoly: VirtualPolynomial<R>) {
        let absorbtion_vec: Vec<R> = unipoly.flattened_ml_extensions
            .iter()
            .map(|mle| mle.evaluations[0])
            .collect();
        self.hasher.absorb_ring_vec(&absorbtion_vec);
        let round = SumCheckRound {
            challenge,
            var_index,
            unipoly,
            _marker: std::marker::PhantomData,
        };

        self.rounds.push(round);
    }
}

impl<F: PrimeField, R: OverField<F>> SumCheckRound<F, R> {
    pub fn new(
        unipoly: VirtualPolynomial<R>,
        var_index: usize,
        challenge: R
    ) -> SumCheckRound<F, R> {
        SumCheckRound {
            challenge,
            var_index,
            unipoly,
            _marker: std::marker::PhantomData,
        }
    }
}
