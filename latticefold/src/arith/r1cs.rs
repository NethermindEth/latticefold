//! Defines behaviour of R1CS, a degree two constraint system

use ark_std::collections::{hash_map, HashMap};
use cyclotomic_rings::rings::SuitableRing;
use stark_rings::Ring;
use stark_rings_linalg::{sparse_matrix::dense_matrix_u64_to_sparse, SparseMatrix};

use super::{
    error::CSError as Error,
    utils::{mat_vec_mul, vec_add, vec_scalar_mul},
};
use crate::{arith::hadamard, ark_base::*};

/// a representation of a R1CS instance
#[derive(Debug, Clone, PartialEq)]
pub struct R1CS<R: Ring> {
    /// Length of public input
    pub l: usize,
    /// First constraint matrix
    pub A: SparseMatrix<R>,
    /// Second constraint matrix
    pub B: SparseMatrix<R>,
    /// Third constraint matrix
    pub C: SparseMatrix<R>,
}

impl<R: Ring> R1CS<R> {
    /// check that a R1CS structure is satisfied by a z vector. Only for testing.
    pub fn check_relation(&self, z: &[R]) -> Result<(), Error> {
        let Az = mat_vec_mul(&self.A, z)?;
        let Bz = mat_vec_mul(&self.B, z)?;

        let Cz = mat_vec_mul(&self.C, z)?;
        let AzBz = hadamard(&Az, &Bz)?;

        if AzBz != Cz {
            Err(Error::NotSatisfied)
        } else {
            Ok(())
        }
    }

    /// Converts the R1CS instance into a RelaxedR1CS as described in
    /// [Nova](https://eprint.iacr.org/2021/370.pdf#page=14)
    pub fn relax(self) -> RelaxedR1CS<R> {
        RelaxedR1CS::<R> {
            l: self.l,
            E: vec![R::zero(); self.A.nrows()],
            A: self.A,
            B: self.B,
            C: self.C,
            u: R::one(),
        }
    }

    /// Create an R1CS from a constraint system
    pub fn from_constraint_system(cs: ConstraintSystem<R>) -> Self {
        cs.to_r1cs()
    }
}

/// A RelaxedR1CS instance as described in
/// [Nova](https://eprint.iacr.org/2021/370.pdf#page=14).
///
/// A witness $z$ is satisfying if $(A \cdot z) \circ (B \cdot z) = u \cdot (C \cdot z) + E$.
#[derive(Debug, Clone, PartialEq)]
pub struct RelaxedR1CS<R: Ring> {
    /// Public input length
    pub l: usize,
    /// First constraint matrix
    pub A: SparseMatrix<R>,
    /// Second constraint matrix
    pub B: SparseMatrix<R>,
    /// Third constraint matrix
    pub C: SparseMatrix<R>,
    /// Scalar coefficient of $(C \cdot z)$
    pub u: R,
    /// The error matrix
    pub E: Vec<R>,
}

impl<R: Ring> RelaxedR1CS<R> {
    /// check that a RelaxedR1CS structure is satisfied by a z vector.
    pub fn check_relation(&self, z: &[R]) -> Result<(), Error> {
        let Az = mat_vec_mul(&self.A, z)?;
        let Bz = mat_vec_mul(&self.B, z)?;
        let Cz = mat_vec_mul(&self.C, z)?;

        let uCz = vec_scalar_mul(&Cz, &self.u);
        let uCzE = vec_add(&uCz, &self.E)?;
        let AzBz = hadamard(&Az, &Bz)?;
        if AzBz != uCzE {
            Err(Error::NotSatisfied)
        } else {
            Ok(())
        }
    }
}

/// Returns a matrix of ring elements given a matrix of unsigned ints
pub fn to_F_matrix<R: Ring>(M: Vec<Vec<usize>>) -> SparseMatrix<R> {
    // dense_matrix_to_sparse(to_F_dense_matrix::<R>(M))
    let M_u64: Vec<Vec<u64>> = M
        .iter()
        .map(|m| m.iter().map(|r| *r as u64).collect())
        .collect();
    dense_matrix_u64_to_sparse(M_u64)
}

/// Returns a dense matrix of ring elements given a matrix of unsigned ints
pub fn to_F_dense_matrix<R: Ring>(M: Vec<Vec<usize>>) -> Vec<Vec<R>> {
    M.iter()
        .map(|m| m.iter().map(|r| R::from(*r as u64)).collect())
        .collect()
}

/// Returns a vector of ring elements given a vector of unsigned ints
pub fn to_F_vec<R: Ring>(z: Vec<usize>) -> Vec<R> {
    z.iter().map(|c| R::from(*c as u64)).collect()
}

#[cfg(test)]
pub(crate) fn get_test_r1cs<R: Ring>() -> R1CS<R> {
    // R1CS for: x^3 + x + 5 = y (example from article
    // https://www.vitalik.ca/general/2016/12/10/qap.html )
    let A = to_F_matrix::<R>(vec![
        vec![1, 0, 0, 0, 0, 0],
        vec![0, 0, 0, 1, 0, 0],
        vec![1, 0, 0, 0, 1, 0],
        vec![0, 5, 0, 0, 0, 1],
    ]);
    let B = to_F_matrix::<R>(vec![
        vec![1, 0, 0, 0, 0, 0],
        vec![1, 0, 0, 0, 0, 0],
        vec![0, 1, 0, 0, 0, 0],
        vec![0, 1, 0, 0, 0, 0],
    ]);
    let C = to_F_matrix::<R>(vec![
        vec![0, 0, 0, 1, 0, 0],
        vec![0, 0, 0, 0, 1, 0],
        vec![0, 0, 0, 0, 0, 1],
        vec![0, 0, 1, 0, 0, 0],
    ]);

    R1CS::<R> { l: 1, A, B, C }
}

/// Return a R1CS instance of arbitrary size, useful for benching.
/// Only works when z vector consists of multiplicative identities.
pub fn get_test_dummy_r1cs<R: Ring, const X_LEN: usize, const WIT_LEN: usize>(
    rows: usize,
) -> R1CS<R> {
    let R1CS_A = create_dummy_identity_sparse_matrix(rows, X_LEN + WIT_LEN + 1);
    let R1CS_B = R1CS_A.clone();
    let R1CS_C = R1CS_A.clone();

    R1CS::<R> {
        l: 1,
        A: R1CS_A,
        B: R1CS_B,
        C: R1CS_C,
    }
}

/// Return a R1CS instance of arbitrary size, useful for benching.
/// Works for arbitrary z vector.
pub fn get_test_dummy_r1cs_non_scalar<R: Ring, const X_LEN: usize, const WIT_LEN: usize>(
    rows: usize,
    witness: &[R],
) -> R1CS<R> {
    let R1CS_A = create_dummy_identity_sparse_matrix(rows, X_LEN + WIT_LEN + 1);
    let R1CS_B = R1CS_A.clone();
    let R1CS_C = create_dummy_squaring_sparse_matrix(rows, X_LEN + WIT_LEN + 1, witness);

    R1CS::<R> {
        l: 1,
        A: R1CS_A,
        B: R1CS_B,
        C: R1CS_C,
    }
}

pub(crate) fn create_dummy_identity_sparse_matrix<R: Ring>(
    rows: usize,
    columns: usize,
) -> SparseMatrix<R> {
    let mut matrix = SparseMatrix {
        n_rows: rows,
        n_cols: columns,
        coeffs: vec![vec![]; rows],
    };
    for (i, row) in matrix.coeffs.iter_mut().enumerate() {
        row.push((R::one(), i));
    }
    matrix
}

// Takes a vector and returns a matrix that will square the vector
pub(crate) fn create_dummy_squaring_sparse_matrix<R: Ring>(
    rows: usize,
    columns: usize,
    witness: &[R],
) -> SparseMatrix<R> {
    assert_eq!(
        rows,
        witness.len(),
        "Length of witness vector must be equal to ccs width"
    );
    let mut matrix = SparseMatrix {
        n_rows: rows,
        n_cols: columns,
        coeffs: vec![vec![]; rows],
    };
    for (i, row) in matrix.coeffs.iter_mut().enumerate() {
        row.push((witness[i], i));
    }
    matrix
}

pub(crate) fn get_test_z<R: Ring>(input: usize) -> Vec<R> {
    // z = (io, 1, w)
    to_F_vec(vec![
        input, // io
        1,
        input * input * input + input + 5, // x^3 + x + 5
        input * input,                     // x^2
        input * input * input,             // x^2 * x
        input * input * input + input,     // x^3 + x
    ])
}

pub(crate) fn get_test_z_ntt<R: SuitableRing>() -> Vec<R> {
    let mut res = Vec::new();
    for input in 0..R::dimension() {
        // z = (io, 1, w)
        res.push(to_F_vec::<R::BaseRing>(vec![
            input, // io
            1,
            input * input * input + input + 5, // x^3 + x + 5
            input * input,                     // x^2
            input * input * input,             // x^2 * x
            input * input * input + input,     // x^3 + x
        ]))
    }

    let mut ret: Vec<R> = Vec::new();
    for j in 0..res[0].len() {
        let mut vec = Vec::new();
        for witness in &res {
            vec.push(witness[j]);
        }
        ret.push(R::from(vec));
    }

    ret
}

/// Return scalar z vector for Vitalik's [R1CS example](https://medium.com/@VitalikButerin/quadratic-arithmetic-programs-from-zero-to-hero-f6d558cea649#81e4),
/// split into statement, constant, and witness.
pub fn get_test_z_split<R: Ring>(input: usize) -> (R, Vec<R>, Vec<R>) {
    let z = get_test_z(input);
    (z[1], vec![z[0]], z[2..].to_vec())
}

/// Return non-scalar z vector for Vitalik's [R1CS example](https://medium.com/@VitalikButerin/quadratic-arithmetic-programs-from-zero-to-hero-f6d558cea649#81e4),
/// split into statement, constant, and witness.
pub fn get_test_z_ntt_split<R: SuitableRing>() -> (R, Vec<R>, Vec<R>) {
    let z = get_test_z_ntt();
    (z[1], vec![z[0]], z[2..].to_vec())
}

/// Return z vector consisting only of multiplicative identities,
/// split into statement, constant, and witness.
pub fn get_test_dummy_z_split<R: Ring, const X_LEN: usize, const WIT_LEN: usize>(
) -> (R, Vec<R>, Vec<R>) {
    (
        R::one(),
        to_F_vec(vec![1; X_LEN]),
        to_F_vec(vec![1; WIT_LEN]),
    )
}

/// Return z vector consisting of non scalar ring elements,
/// split into statement, constant, and witness.
pub fn get_test_dummy_z_split_ntt<R: SuitableRing, const X_LEN: usize, const WIT_LEN: usize>(
) -> (R, Vec<R>, Vec<R>) {
    let statement_vec = (0..X_LEN).map(|_| R::one()).collect();

    let witness_vec = (0..WIT_LEN)
        .map(|_| {
            R::from(
                (0..R::dimension())
                    .map(|i| R::BaseRing::from(i as u128))
                    .collect::<Vec<R::BaseRing>>(),
            )
        })
        .collect();

    (R::one(), statement_vec, witness_vec)
}

/// A linear combination of variables
#[derive(Debug, Clone, PartialEq)]
pub struct LinearCombination<R: Ring> {
    /// The terms in the linear combination
    /// For each element, `.0` is the coefficient, `.1` is the variable index
    pub terms: Vec<(R, usize)>,
}

impl<R: Ring> LinearCombination<R> {
    /// Create a new empty linear combination
    pub fn new() -> Self {
        Self { terms: Vec::new() }
    }

    /// Initializes a new linear combination with a single term
    pub fn single_term(coeff: impl Into<R>, index: usize) -> Self {
        Self {
            terms: vec![(coeff.into(), index)],
        }
    }

    /// Add a term to the linear combination
    pub fn add_term(mut self, coeff: impl Into<R>, index: usize) -> Self {
        self.terms.push((coeff.into(), index));
        self
    }

    /// Add a term to the linear combination
    pub fn add_terms(mut self, terms: &[(R, usize)]) -> Self {
        terms.iter().for_each(|&term| self.terms.push(term));
        self
    }

    /// Evaluate the linear combination given a variable assignment
    pub fn evaluate(&self, assignment: &[R]) -> R {
        self.terms
            .iter()
            .fold(R::zero(), |acc, term| acc + term.0 * assignment[term.1])
    }

    /// Check if the linear combination is valid given the number of variables
    pub fn is_valid(&self, nvars: usize) -> bool {
        self.terms.iter().all(|term| term.1 < nvars)
    }
}

/// A single R1CS constraint of the form A * B = C
#[derive(Debug, Clone, PartialEq)]
pub struct Constraint<R: Ring> {
    /// The A linear combination
    pub a: LinearCombination<R>,
    /// The B linear combination
    pub b: LinearCombination<R>,
    /// The C linear combination
    pub c: LinearCombination<R>,
}

impl<R: Ring> Constraint<R> {
    /// Create a new constraint from three linear combinations
    pub fn new(a: LinearCombination<R>, b: LinearCombination<R>, c: LinearCombination<R>) -> Self {
        Self { a, b, c }
    }
}

/// A system of R1CS constraints
#[derive(Debug, Clone, PartialEq)]
pub struct ConstraintSystem<R: Ring> {
    /// The number of public inputs
    pub ninputs: usize,
    /// The number of private inputs (auxiliary variables)
    pub nauxs: usize,
    /// The constraints in the system
    pub constraints: Vec<Constraint<R>>,
    /// The variable map
    pub vars: VariableMap,
}

impl<R: Ring> ConstraintSystem<R> {
    /// Create a new empty constraint system
    pub fn new() -> Self {
        Self {
            ninputs: 0,
            nauxs: 0,
            constraints: Vec::new(),
            vars: VariableMap::new(),
        }
    }

    /// Get the total number of variables
    pub fn nvars(&self) -> usize {
        self.ninputs + self.nauxs
    }

    /// Get the number of constraints
    pub fn nconstraints(&self) -> usize {
        self.constraints.len()
    }

    /// Add a constraint to the system
    pub fn add_constraint(&mut self, constraint: Constraint<R>) -> &mut Self {
        self.constraints.push(constraint);
        self
    }

    /// Check if the constraint system is valid
    pub fn is_valid(&self) -> bool {
        if self.ninputs > self.nvars() {
            return false;
        }

        for constraint in &self.constraints {
            if !(constraint.a.is_valid(self.nvars())
                && constraint.b.is_valid(self.nvars())
                && constraint.c.is_valid(self.nvars()))
            {
                return false;
            }
        }

        true
    }

    /// Check if the constraint system is satisfied by the given inputs
    pub fn is_satisfied(&self, primary_input: &[R], auxiliary_input: &[R]) -> Result<(), Error> {
        if primary_input.len() != self.ninputs {
            return Err(Error::LengthsNotEqual(
                "primary_input".to_string(),
                "num_inputs".to_string(),
                primary_input.len(),
                self.ninputs,
            ));
        }

        if primary_input.len() + auxiliary_input.len() != self.nvars() {
            return Err(Error::LengthsNotEqual(
                "primary_input + auxiliary_input".to_string(),
                "num_variables".to_string(),
                primary_input.len() + auxiliary_input.len(),
                self.nvars(),
            ));
        }

        // Combine primary and auxiliary inputs into a full assignment
        let mut full_assignment = primary_input.to_vec();
        full_assignment.extend_from_slice(auxiliary_input);

        // Check each constraint
        self.constraints.iter().try_for_each(|constraint| {
            let a_res = constraint.a.evaluate(&full_assignment);
            let b_res = constraint.b.evaluate(&full_assignment);
            let c_res = constraint.c.evaluate(&full_assignment);

            if a_res * b_res != c_res {
                return Err(Error::NotSatisfied);
            }

            Ok(())
        })
    }

    /// Swap A and B matrices if it would be beneficial for performance
    pub fn swap_AB_if_beneficial(&mut self) {
        // Count non-zero entries in A and B
        let mut touched_by_A = vec![false; self.nvars() + 1];
        let mut touched_by_B = vec![false; self.nvars() + 1];

        for constraint in &self.constraints {
            for term in &constraint.a.terms {
                touched_by_A[term.1] = true;
            }
            for term in &constraint.b.terms {
                touched_by_B[term.1] = true;
            }
        }

        let non_zero_A_count = touched_by_A.iter().filter(|&&x| x).count();
        let non_zero_B_count = touched_by_B.iter().filter(|&&x| x).count();

        // If B has more non-zero entries than A, swap them
        if non_zero_B_count > non_zero_A_count {
            for constraint in &mut self.constraints {
                ark_std::mem::swap(&mut constraint.a, &mut constraint.b);
            }
        }
    }

    /// Convert to a sparse matrix representation (R1CS)
    pub fn to_r1cs(&self) -> R1CS<R> {
        let nconstraints = self.nconstraints();
        let nvars = self.nvars();

        // Create empty sparse matrices
        let mut A = SparseMatrix {
            n_rows: nconstraints,
            n_cols: nvars,
            coeffs: vec![vec![]; nconstraints],
        };

        let mut B = SparseMatrix {
            n_rows: nconstraints,
            n_cols: nvars,
            coeffs: vec![vec![]; nconstraints],
        };

        let mut C = SparseMatrix {
            n_rows: nconstraints,
            n_cols: nvars,
            coeffs: vec![vec![]; nconstraints],
        };

        // Fill the matrices
        for (i, constraint) in self.constraints.iter().enumerate() {
            for &term in &constraint.a.terms {
                A.coeffs[i].push(term);
            }
            for &term in &constraint.b.terms {
                B.coeffs[i].push(term);
            }
            for &term in &constraint.c.terms {
                C.coeffs[i].push(term);
            }
        }

        R1CS {
            l: self.ninputs,
            A,
            B,
            C,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct VariableMap {
    map: HashMap<String, (usize, usize)>,
    one: usize,
}

impl VariableMap {
    pub fn new() -> Self {
        Self {
            map: HashMap::new(),
            one: 0,
        }
    }

    pub fn add(&mut self, name: impl Into<String>, index: usize, len: usize) {
        self.map.insert(name.into(), (index, len));
    }

    pub fn get(&self, name: &str) -> Option<(usize, usize)> {
        self.map.get(name).copied()
    }

    pub fn add_one(&mut self, index: usize) -> usize {
        self.one = index;
        index
    }

    pub fn get_one(&self) -> usize {
        self.one
    }

    pub fn vars(&self) -> hash_map::Iter<String, (usize, usize)> {
        self.map.iter()
    }
}

#[cfg(test)]
mod tests {
    use cyclotomic_rings::rings::GoldilocksRingNTT;

    use super::*;

    type RqNTT = GoldilocksRingNTT;

    #[test]
    fn test_r1cs_check_relation() {
        let r1cs = get_test_r1cs::<RqNTT>();
        let z = get_test_z(5);

        r1cs.check_relation(&z).unwrap();
        r1cs.relax().check_relation(&z).unwrap();
    }

    #[test]
    fn test_r1cs_linear_combination() {
        let lc = LinearCombination::<RqNTT>::new()
            .add_term(2u64, 0)
            .add_term(3u64, 1);
        let assignment = vec![RqNTT::from(5u64), RqNTT::from(7u64)];
        let result = lc.evaluate(&assignment);
        assert_eq!(result, RqNTT::from(31u64)); // 2*5 + 3*7 = 31
    }

    #[test]
    fn test_r1cs_constraint_system() {
        // Create a constraint system for x * y = z
        let mut cs = ConstraintSystem::<RqNTT>::new();
        cs.ninputs = 2; // x and y are public inputs
        cs.nauxs = 1; // z is an auxiliary variable

        // Create the constraint x * y = z
        let a = LinearCombination::new().add_term(1u64, 0); // x
        let b = LinearCombination::new().add_term(1u64, 1); // y
        let c = LinearCombination::new().add_term(1u64, 2); // z

        let constraint = Constraint::new(a, b, c);
        cs.add_constraint(constraint);

        let primary_input = vec![RqNTT::from(3u64), RqNTT::from(4u64)]; // x=3, y=4
        let aux_input = vec![RqNTT::from(12u64)]; // z=12
        cs.is_satisfied(&primary_input, &aux_input).unwrap();
    }

    #[test]
    fn test_r1cs_from_constraint_system() {
        // Create a constraint system for x * y = z
        let mut cs = ConstraintSystem::<RqNTT>::new();
        cs.ninputs = 2; // x and y are public inputs
        cs.nauxs = 1; // z is an auxiliary variable

        // Create the constraint x * y = z
        let a = LinearCombination::new().add_term(1u64, 0); // x
        let b = LinearCombination::new().add_term(1u64, 1); // y
        let c = LinearCombination::new().add_term(1u64, 2); // z

        let constraint = Constraint::new(a, b, c);
        cs.add_constraint(constraint);

        // Convert to R1CS
        let r1cs = R1CS::from_constraint_system(cs);

        // Check that the R1CS is correct
        assert_eq!(r1cs.l, 2);
        assert_eq!(r1cs.A.n_rows, 1);
        assert_eq!(r1cs.A.n_cols, 3);
        assert_eq!(r1cs.B.n_rows, 1);
        assert_eq!(r1cs.B.n_cols, 3);
        assert_eq!(r1cs.C.n_rows, 1);
        assert_eq!(r1cs.C.n_cols, 3);

        // Test with a valid assignment
        let z = vec![RqNTT::from(3u64), RqNTT::from(4u64), RqNTT::from(12u64)]; // x=3, y=4, z=12
        r1cs.check_relation(&z).unwrap();
    }

    #[test]
    fn test_r1cs_example_from_constraint_system() {
        // x^3 + x + 5 = y
        let mut cs = ConstraintSystem::<RqNTT>::new();

        // 1 public input (x) and 5 auxiliary variables
        cs.ninputs = 1;
        cs.nauxs = 5;

        // Variables:
        // 0: x (public input)
        // 1: 1 (constant)
        // 2: y (output)
        // 3: x^2
        // 4: x^3
        // 5: x^3 + x

        // Constraint 1: x * x = x^2
        let a1 = LinearCombination::new().add_term(1u64, 0); // x
        let b1 = LinearCombination::new().add_term(1u64, 0); // x
        let c1 = LinearCombination::new().add_term(1u64, 3); // x^2
        cs.add_constraint(Constraint::new(a1, b1, c1));

        // Constraint 2: x^2 * x = x^3
        let a2 = LinearCombination::new().add_term(1u64, 3); // x^2
        let b2 = LinearCombination::new().add_term(1u64, 0); // x
        let c2 = LinearCombination::new().add_term(1u64, 4); // x^3
        cs.add_constraint(Constraint::new(a2, b2, c2));

        // Constraint 3: (x + x^3) * 1 = x^3 + x
        let a3 = LinearCombination::new().add_term(1u64, 0).add_term(1u64, 4); // x + x^3
        let b3 = LinearCombination::new().add_term(1u64, 1); // 1
        let c3 = LinearCombination::new().add_term(1u64, 5); // x^3 + x
        cs.add_constraint(Constraint::new(a3, b3, c3));

        // Constraint 4: (5*1 + x^3 + x) * 1 = y
        let a4 = LinearCombination::new().add_term(5u64, 1).add_term(1u64, 5);
        let b4 = LinearCombination::new().add_term(1u64, 1); // 1
        let c4 = LinearCombination::new().add_term(1u64, 2); // y
        cs.add_constraint(Constraint::new(a4, b4, c4));

        let r1cs_from_cs = R1CS::from_constraint_system(cs);
        let r1cs_original = get_test_r1cs::<RqNTT>();

        // Test with x = 5
        let x = 5;
        let z = get_test_z(x);

        r1cs_original.check_relation(&z).unwrap();
        r1cs_from_cs.check_relation(&z).unwrap();

        let y = x * x * x + x + 5;
        assert_eq!(z[2], RqNTT::from(y as u64));
    }
}
