use std::{
    fmt::{Debug, Display},
    iter::{Product, Sum},
    ops::{Mul, MulAssign, Neg},
};

use num_traits::{NumAssign, One};

mod mersenne31;

// pub trait FieldOps<F: Field>: ColumnOps<F> {
//     // TODO(Ohad): change to use a mutable slice.
//     fn batch_inverse(column: &Self::Column, dst: &mut Self::Column);
// }

pub trait FieldExpOps: Mul<Output = Self> + MulAssign + Sized + One + Copy {
    fn square(&self) -> Self {
        (*self) * (*self)
    }

    fn pow(&self, exp: u128) -> Self {
        let mut res = Self::one();
        let mut base = *self;
        let mut exp = exp;
        while exp > 0 {
            if exp & 1 == 1 {
                res *= base;
            }
            base = base.square();
            exp >>= 1;
        }
        res
    }

    fn inverse(&self) -> Self;

    // /// Inverts a batch of elements using Montgomery's trick.
    // fn batch_inverse(column: &[Self], dst: &mut [Self]) {
    //     const WIDTH: usize = 4;
    //     let n = column.len();
    //     debug_assert!(dst.len() >= n);
    //
    //     if n <= WIDTH || n % WIDTH != 0 {
    //         batch_inverse_classic(column, dst);
    //         return;
    //     }
    //
    //     // First pass. Compute 'WIDTH' cumulative products in an interleaving fashion, reducing
    //     // instruction dependency and allowing better pipelining.
    //     let mut cum_prod: [Self; WIDTH] = [Self::one(); WIDTH];
    //     dst[..WIDTH].copy_from_slice(&cum_prod);
    //     for i in 0..n {
    //         cum_prod[i % WIDTH] *= column[i];
    //         dst[i] = cum_prod[i % WIDTH];
    //     }
    //
    //     // Inverse cumulative products.
    //     // Use classic batch inversion.
    //     let mut tail_inverses = [Self::one(); WIDTH];
    //     batch_inverse_classic(&dst[n - WIDTH..], &mut tail_inverses);
    //
    //     // Second pass.
    //     for i in (WIDTH..n).rev() {
    //         dst[i] = dst[i - WIDTH] * tail_inverses[i % WIDTH];
    //         tail_inverses[i % WIDTH] *= column[i];
    //     }
    //     dst[0..WIDTH].copy_from_slice(&tail_inverses);
    // }
}

pub trait Field:
    NumAssign
    + Neg<Output = Self>
    + ComplexConjugate
    + Copy
    + Default
    + Debug
    + Display
    + PartialOrd
    + Ord
    + Send
    + Sync
    + Sized
    + FieldExpOps
    + Product
    + for<'a> Product<&'a Self>
    + Sum
    + for<'a> Sum<&'a Self>
{
    fn double(&self) -> Self {
        (*self) + (*self)
    }
}

pub trait ComplexConjugate {}

#[macro_export]
macro_rules! impl_field {
    ($field_name: ty, $field_size: ident) => {
        use std::iter::{Product, Sum};

        use num_traits::{Num, One, Zero};
        use $crate::fields::Field;

        impl Num for $field_name {
            type FromStrRadixErr = Box<dyn std::error::Error>;

            fn from_str_radix(_str: &str, _radix: u32) -> Result<Self, Self::FromStrRadixErr> {
                unimplemented!(
                    "Num::from_str_radix is not implemented for {}",
                    stringify!($field_name)
                );
            }
        }

        impl Field for $field_name {}

        impl AddAssign for $field_name {
            fn add_assign(&mut self, rhs: Self) {
                *self = *self + rhs;
            }
        }

        impl SubAssign for $field_name {
            fn sub_assign(&mut self, rhs: Self) {
                *self = *self - rhs;
            }
        }

        impl MulAssign for $field_name {
            fn mul_assign(&mut self, rhs: Self) {
                *self = *self * rhs;
            }
        }

        impl Div for $field_name {
            type Output = Self;

            #[allow(clippy::suspicious_arithmetic_impl)]
            fn div(self, rhs: Self) -> Self::Output {
                self * rhs.inverse()
            }
        }

        impl DivAssign for $field_name {
            fn div_assign(&mut self, rhs: Self) {
                *self = *self / rhs;
            }
        }

        impl Rem for $field_name {
            type Output = Self;

            fn rem(self, _rhs: Self) -> Self::Output {
                unimplemented!("Rem is not implemented for {}", stringify!($field_name));
            }
        }

        impl RemAssign for $field_name {
            fn rem_assign(&mut self, _rhs: Self) {
                unimplemented!(
                    "RemAssign is not implemented for {}",
                    stringify!($field_name)
                );
            }
        }

        impl FieldExpOps for $field_name {
            fn inverse(&self) -> Self {
                assert!(!self.is_zero(), "0 has no inverse");
                self.pow(($field_size - 2) as u128)
            }
        }

        impl Product for $field_name {
            fn product<I>(mut iter: I) -> Self
            where
                I: Iterator<Item = Self>,
            {
                let first = iter.next().unwrap_or_else(Self::one);
                iter.fold(first, |a, b| a * b)
            }
        }

        impl<'a> Product<&'a Self> for $field_name {
            fn product<I>(iter: I) -> Self
            where
                I: Iterator<Item = &'a Self>,
            {
                iter.map(|&v| v).product()
            }
        }

        impl Sum for $field_name {
            fn sum<I>(mut iter: I) -> Self
            where
                I: Iterator<Item = Self>,
            {
                let first = iter.next().unwrap_or_else(Self::zero);
                iter.fold(first, |a, b| a + b)
            }
        }

        impl<'a> Sum<&'a Self> for $field_name {
            fn sum<I>(iter: I) -> Self
            where
                I: Iterator<Item = &'a Self>,
            {
                iter.map(|&v| v).sum()
            }
        }
    };
}
