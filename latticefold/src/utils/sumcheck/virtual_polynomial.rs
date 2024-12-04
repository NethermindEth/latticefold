// Copyright (c) 2023 Espresso Systems (espressosys.com)
// This file is part of the HyperPlonk library.

// Adapted for rings by Nethermind

//! This module defines our main mathematical object `VirtualPolynomial`; and
//! various functions associated with it.

use ark_serialize::CanonicalSerialize;
use ark_std::{
    cfg_iter_mut, end_timer,
    rand::{Rng, RngCore},
    start_timer,
    string::ToString,
    vec::*,
};
use ark_std::{cmp::max, marker::PhantomData, ops::Add};
use lattirust_poly::mle::DenseMultilinearExtension;
use lattirust_poly::polynomials::{random_mle_list, random_zero_mle_list, ArithErrors, RefCounter};

#[cfg(feature = "std")]
use ark_std::collections::HashSet;
#[cfg(not(feature = "std"))]
use hashbrown::HashSet;

use lattirust_ring::Ring;
#[cfg(feature = "parallel")]
use rayon::prelude::*;

#[rustfmt::skip]
/// A virtual polynomial is a sum of products of multilinear polynomials;
/// where the multilinear polynomials are stored via their multilinear
/// extensions:  `(coefficient, DenseMultilinearExtension)`
///
/// * Number of products n = `polynomial.products.len()`,
/// * Number of multiplicands of ith product m_i =
///   `polynomial.products[i].1.len()`,
/// * Coefficient of ith product c_i = `polynomial.products[i].0`
///
/// The resulting polynomial is
///
/// $$ \sum_{i=0}^{n} c_i \cdot \prod_{j=0}^{m_i} P_{ij} $$
///
/// Example:
///  f = c0 * f0 * f1 * f2 + c1 * f3 * f4
/// where f0 ... f4 are multilinear polynomials
///
/// - flattened_ml_extensions stores the multilinear extension representation of
///   f0, f1, f2, f3 and f4
/// - products is
///     \[
///         (c0, \[0, 1, 2\]),
///         (c1, \[3, 4\])
///     \]
/// - raw_pointers_lookup_table maps fi to i
///
#[derive(Clone, Debug, Default, PartialEq)]
pub struct VirtualPolynomial<R: Ring> {
    /// Aux information about the multilinear polynomial
    pub aux_info: VPAuxInfo<R>,
    /// Stores multilinear extensions in which product multiplicand can refer
    /// to.
    pub flattened_ml_extensions: Vec<RefCounter<DenseMultilinearExtension<R>>>,
    /// Pointers to the above poly extensions
    raw_pointers_lookup_table: HashSet<*const DenseMultilinearExtension<R>>,
}

#[derive(Clone, Debug, Default, PartialEq, Eq, CanonicalSerialize)]
/// Auxiliary information about the multilinear polynomial
pub struct VPAuxInfo<R: Ring> {
    /// max number of multiplicands in each product
    pub max_degree: usize,
    /// number of variables of the polynomial
    pub num_variables: usize,
    /// Associated field
    #[doc(hidden)]
    pub phantom: PhantomData<R>,
}

// TODO: convert this into a trait
impl<R: Ring> VirtualPolynomial<R> {
    /// Creates an empty virtual polynomial with `num_variables`.
    pub fn new(num_variables: usize) -> Self {
        VirtualPolynomial {
            aux_info: VPAuxInfo {
                max_degree: 0,
                num_variables,
                phantom: PhantomData,
            },
            flattened_ml_extensions: Vec::new(),
            raw_pointers_lookup_table: HashSet::new(),
        }
    }

    /// Creates an new virtual polynomial from a MLE and its coefficient.
    pub fn new_from_mle(mle: &RefCounter<DenseMultilinearExtension<R>>, coefficient: R) -> Self {
        let mle_ptr: *const DenseMultilinearExtension<R> = RefCounter::as_ptr(mle);
        let mut hm = HashSet::new();
        hm.insert(mle_ptr);

        VirtualPolynomial {
            aux_info: VPAuxInfo {
                // The max degree is the max degree of any individual variable
                max_degree: 1,
                num_variables: mle.num_vars,
                phantom: PhantomData,
            },
            // here `0` points to the first polynomial of `flattened_ml_extensions`
            flattened_ml_extensions: vec![mle.clone()],
            raw_pointers_lookup_table: hm,
        }
    }

    /// Add a product of list of multilinear extensions to self
    /// Returns an error if the list is empty, or the MLE has a different
    /// `num_vars` from self.
    ///
    /// The MLEs will be multiplied together, and then multiplied by the scalar
    /// `coefficient`.
    pub fn add_mle_list(
        &mut self,
        mle_list: impl IntoIterator<Item = RefCounter<DenseMultilinearExtension<R>>>,
    ) -> Result<(), ArithErrors> {
        let mle_list: Vec<RefCounter<DenseMultilinearExtension<R>>> =
            mle_list.into_iter().collect();

        if mle_list.is_empty() {
            return Err(ArithErrors::InvalidParameters(
                "input mle_list is empty".to_string(),
            ));
        }

        self.aux_info.max_degree = max(self.aux_info.max_degree, mle_list.len());

        for mle in mle_list {
            if mle.num_vars != self.aux_info.num_variables {
                return Err(ArithErrors::InvalidParameters(format!(
                    "product has a multiplicand with wrong number of variables {} vs {}",
                    mle.num_vars, self.aux_info.num_variables
                )));
            }

            let mle_ptr: *const DenseMultilinearExtension<R> = RefCounter::as_ptr(&mle);
            if self.raw_pointers_lookup_table.get(&mle_ptr).is_none() {
                let curr_index = self.flattened_ml_extensions.len();
                self.flattened_ml_extensions.push(mle.clone());
                self.raw_pointers_lookup_table.insert(mle_ptr);
            }
        }
        Ok(())
    }

    /// Multiple the current VirtualPolynomial by an MLE:
    /// - add the MLE to the MLE list;
    /// - multiple each product by MLE and its coefficient.
    ///
    /// Returns an error if the MLE has a different `num_vars` from self.
    pub fn mul_by_mle(
        &mut self,
        mle: RefCounter<DenseMultilinearExtension<R>>,
    ) -> Result<(), ArithErrors> {
        let start = start_timer!(|| "mul by mle");

        if mle.num_vars != self.aux_info.num_variables {
            return Err(ArithErrors::InvalidParameters(format!(
                "product has a multiplicand with wrong number of variables {} vs {}",
                mle.num_vars, self.aux_info.num_variables
            )));
        }

        let mle_ptr: *const DenseMultilinearExtension<R> = RefCounter::as_ptr(&mle);

        // check if this mle already exists in the virtual polynomial
        if self.raw_pointers_lookup_table.get(&mle_ptr).is_none() {
            self.raw_pointers_lookup_table.insert(mle_ptr);
            self.flattened_ml_extensions.push(mle);
        };

        // increase the max degree by one as the MLE has degree 1.
        self.aux_info.max_degree += 1;
        end_timer!(start);
        Ok(())
    }

    /// Evaluate the virtual polynomial at point `point`.
    /// Returns an error is point.len() does not match `num_variables`.
    pub fn evaluate(&self, point: &[R]) -> Result<R, ArithErrors> {
        let start = start_timer!(|| "evaluation");

        if self.aux_info.num_variables != point.len() {
            return Err(ArithErrors::InvalidParameters(format!(
                "wrong number of variables {} vs {}",
                self.aux_info.num_variables,
                point.len()
            )));
        }

        let evals: Vec<R> = self
            .flattened_ml_extensions
            .iter()
            .map(|x| {
                x.evaluate(point).unwrap() // safe unwrap here since we have
                                           // already checked that num_var
                                           // matches
            })
            .collect();
        let res = R::zero();
        //let res = self
        //    .products
        //    .iter()
        //    .map(|(c, p)| *c * p.iter().map(|&i| evals[i]).product::<R>())
        //    .sum();

        end_timer!(start);
        Ok(res)
    }

    /// Sample a random virtual polynomial, return the polynomial and its sum.
    pub fn rand<Rn: RngCore>(
        nv: usize,
        num_multiplicands_range: (usize, usize),
        num_products: usize,
        rng: &mut Rn,
    ) -> Result<(Self, R), ArithErrors> {
        let start = start_timer!(|| "sample random virtual polynomial");

        let mut sum = R::zero();
        let mut poly = VirtualPolynomial::new(nv);
        for _ in 0..num_products {
            let num_multiplicands =
                rng.gen_range(num_multiplicands_range.0..num_multiplicands_range.1);
            let (product, product_sum) = random_mle_list(nv, num_multiplicands, rng);
            //let coefficient = R::rand(rng);
            //poly.add_mle_list(product.into_iter(), coefficient)?;
            poly.add_mle_list(product.into_iter())?;
            //sum += product_sum * coefficient;
        }

        end_timer!(start);
        Ok((poly, sum))
    }

    /// Sample a random virtual polynomial that evaluates to zero everywhere
    /// over the boolean hypercube.
    pub fn rand_zero<Rn: RngCore>(
        nv: usize,
        num_multiplicands_range: (usize, usize),
        num_products: usize,
        rng: &mut Rn,
    ) -> Result<Self, ArithErrors> {
        let mut poly = VirtualPolynomial::new(nv);
        for _ in 0..num_products {
            let num_multiplicands =
                rng.gen_range(num_multiplicands_range.0..num_multiplicands_range.1);
            let product = random_zero_mle_list(nv, num_multiplicands, rng);
            //let coefficient = R::rand(rng);
            //poly.add_mle_list(product.into_iter(), coefficient)?;
            poly.add_mle_list(product.into_iter())?;
        }

        Ok(poly)
    }

    // Input poly f(x) and a random vector r, output
    //      \hat f(x) = \sum_{x_i \in eval_x} f(x_i) eq(x, r)
    // where
    //      eq(x,y) = \prod_i=1^num_var (x_i * y_i + (1-x_i)*(1-y_i))
    //
    // This function is used in ZeroCheck.
    pub fn build_f_hat(&self, r: &[R]) -> Result<Self, ArithErrors> {
        let start = start_timer!(|| "zero check build hat f");

        if self.aux_info.num_variables != r.len() {
            return Err(ArithErrors::InvalidParameters(format!(
                "r.len() is different from number of variables: {} vs {}",
                r.len(),
                self.aux_info.num_variables
            )));
        }

        let eq_x_r = build_eq_x_r(r)?;
        let mut res = self.clone();
        res.mul_by_mle(eq_x_r)?;

        end_timer!(start);
        Ok(res)
    }

    /// Print out the evaluation map for testing. Panic if the num_vars > 5.
    #[cfg(feature = "std")]
    pub fn print_evals(&self) {
        if self.aux_info.num_variables > 5 {
            panic!("this function is used for testing only. cannot print more than 5 num_vars");
        }
        for i in 0..1 << self.aux_info.num_variables {
            let point = bit_decompose(i, self.aux_info.num_variables);
            let point_fr: Vec<R> = point.iter().map(|&x| R::from(x)).collect();
            println!("{} {}", i, self.evaluate(point_fr.as_ref()).unwrap());
        }
        println!()
    }
}

impl<R: Ring> VPAuxInfo<R> {
    pub fn new(num_variables: usize, max_degree: usize) -> Self {
        VPAuxInfo {
            max_degree,
            num_variables,
            phantom: PhantomData,
        }
    }
}

/// Evaluate eq polynomial.
pub fn eq_eval<R: Ring>(x: &[R], y: &[R]) -> Result<R, ArithErrors> {
    if x.len() != y.len() {
        return Err(ArithErrors::InvalidParameters(
            "x and y have different length".to_string(),
        ));
    }
    let start = start_timer!(|| "eq_eval");
    let mut res = R::one();
    for (&xi, &yi) in x.iter().zip(y.iter()) {
        let xi_yi = xi * yi;
        res *= xi_yi + xi_yi - xi - yi + R::one();
    }
    end_timer!(start);
    Ok(res)
}

/// This function build the eq(x, r) polynomial for any given r.
///
/// Evaluate
///      eq(x,y) = \prod_i=1^num_var (x_i * y_i + (1-x_i)*(1-y_i))
/// over r, which is
///      eq(x,y) = \prod_i=1^num_var (x_i * r_i + (1-x_i)*(1-r_i))
pub fn build_eq_x_r<R: Ring>(
    r: &[R],
) -> Result<RefCounter<DenseMultilinearExtension<R>>, ArithErrors> {
    let evals = build_eq_x_r_vec(r)?;
    let mle = DenseMultilinearExtension::from_evaluations_vec(r.len(), evals);

    Ok(RefCounter::new(mle))
}
/// This function build the eq(x, r) polynomial for any given r, and output the
/// evaluation of eq(x, r) in its vector form.
///
/// Evaluate
///      eq(x,y) = \prod_i=1^num_var (x_i * y_i + (1-x_i)*(1-y_i))
/// over r, which is
///      eq(x,y) = \prod_i=1^num_var (x_i * r_i + (1-x_i)*(1-r_i))
pub fn build_eq_x_r_vec<R: Ring>(r: &[R]) -> Result<Vec<R>, ArithErrors> {
    // we build eq(x,r) from its evaluations
    // we want to evaluate eq(x,r) over x \in {0, 1}^num_vars
    // for example, with num_vars = 4, x is a binary vector of 4, then
    //  0 0 0 0 -> (1-r0)   * (1-r1)    * (1-r2)    * (1-r3)
    //  1 0 0 0 -> r0       * (1-r1)    * (1-r2)    * (1-r3)
    //  0 1 0 0 -> (1-r0)   * r1        * (1-r2)    * (1-r3)
    //  1 1 0 0 -> r0       * r1        * (1-r2)    * (1-r3)
    //  ....
    //  1 1 1 1 -> r0       * r1        * r2        * r3
    // we will need 2^num_var evaluations

    let mut eval = Vec::new();
    build_eq_x_r_helper(r, &mut eval)?;

    Ok(eval)
}

/// A helper function to build eq(x, r) recursively.
/// This function takes `r.len()` steps, and for each step it requires a maximum
/// `r.len()-1` multiplications.
fn build_eq_x_r_helper<R: Ring>(r: &[R], buf: &mut Vec<R>) -> Result<(), ArithErrors> {
    if r.is_empty() {
        return Err(ArithErrors::InvalidParameters("r length is 0".to_string()));
    } else if r.len() == 1 {
        // initializing the buffer with [1-r_0, r_0]
        buf.push(R::one() - r[0]);
        buf.push(r[0]);
    } else {
        build_eq_x_r_helper(&r[1..], buf)?;

        // suppose at the previous step we received [b_1, ..., b_k]
        // for the current step we will need
        // if x_0 = 0:   (1-r0) * [b_1, ..., b_k]
        // if x_0 = 1:   r0 * [b_1, ..., b_k]
        // let mut res = vec![];
        // for &b_i in buf.iter() {
        //     let tmp = r[0] * b_i;
        //     res.push(b_i - tmp);
        //     res.push(tmp);
        // }
        // *buf = res;

        let mut res = vec![R::zero(); buf.len() << 1];
        cfg_iter_mut!(res).enumerate().for_each(|(i, val)| {
            let bi = buf[i >> 1];
            let tmp = r[0] * bi;
            if (i & 1) == 0 {
                *val = bi - tmp;
            } else {
                *val = tmp;
            }
        });
        *buf = res;
    }

    Ok(())
}

/// Decompose an integer into a binary vector in little endian.
#[cfg(feature = "std")]
pub fn bit_decompose(input: u64, num_var: usize) -> Vec<bool> {
    let mut res = Vec::with_capacity(num_var);
    let mut i = input;
    for _ in 0..num_var {
        res.push((i & 1) == 1);
        i >>= 1;
    }
    res
}
