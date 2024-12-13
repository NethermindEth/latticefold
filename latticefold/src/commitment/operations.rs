//! Provide macros to expand the implementation of commitment operations

/// Given the additive operation for two references of a type,
/// implement the additive operations for non-references.
#[macro_export]
macro_rules! impl_additive_ops_from_ref {
    ($type: ident, $params: ident, $constant: ident) => {
        #[allow(unused_qualifications)]
        impl<'a, const C: $constant, R: $params> Add<$type<C, R>> for &'a $type<C, R> {
            type Output = $type<C, R>;

            fn add(self, rhs: $type<C, R>) -> Self::Output {
                self + &rhs
            }
        }
        #[allow(unused_qualifications)]
        impl<'a, const C: usize, R: Ring> Add<&'a $type<C, R>> for $type<C, R> {
            type Output = $type<C, R>;

            fn add(self, rhs: &'a $type<C, R>) -> Self::Output {
                &self + rhs
            }
        }
        #[allow(unused_qualifications)]
        impl<const C: usize, R: Ring> Add<$type<C, R>> for $type<C, R> {
            type Output = $type<C, R>;

            fn add(self, rhs: $type<C, R>) -> Self::Output {
                &self + &rhs
            }
        }

        #[allow(unused_qualifications)]
        impl<const C: $constant, R: $params> core::ops::AddAssign<Self> for $type<C, R> {
            fn add_assign(&mut self, other: Self) {
                *self = &*self + &other;
            }
        }

        #[allow(unused_qualifications)]
        impl<'a, const C: $constant, R: $params> core::ops::AddAssign<&'a mut Self>
            for $type<C, R>
        {
            fn add_assign(&mut self, other: &'a mut Self) {
                *self = &*self + &*other;
            }
        }

        #[allow(unused_qualifications)]
        impl<const C: $constant, R: $params> core::iter::Sum<Self> for $type<C, R> {
            fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
                iter.fold(Self::zero(), core::ops::Add::add)
            }
        }

        #[allow(unused_qualifications)]
        impl<'a, const C: $constant, R: $params> core::iter::Sum<&'a Self> for $type<C, R> {
            fn sum<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
                iter.fold(Self::zero(), core::ops::Add::add)
            }
        }
    };
}

/// Given the subtractive operation for two references of a type,
/// implement the additive operations for non-references.
#[macro_export]
macro_rules! impl_subtractive_ops_from_ref {
    ($type: ident, $params: ident, $constant: ident) => {
        #[allow(unused_qualifications)]
        impl<'a, const C: usize, R: Ring> Sub<$type<C, R>> for &'a $type<C, R> {
            type Output = $type<C, R>;

            fn sub(self, rhs: $type<C, R>) -> Self::Output {
                self - &rhs
            }
        }
        #[allow(unused_qualifications)]
        impl<'a, const C: usize, R: Ring> Sub<&'a $type<C, R>> for $type<C, R> {
            type Output = $type<C, R>;

            fn sub(self, rhs: &'a $type<C, R>) -> Self::Output {
                &self - rhs
            }
        }
        #[allow(unused_qualifications)]
        impl<const C: usize, R: Ring> Sub<$type<C, R>> for $type<C, R> {
            type Output = $type<C, R>;

            fn sub(self, rhs: $type<C, R>) -> Self::Output {
                &self - &rhs
            }
        }
        #[allow(unused_qualifications)]
        impl<const C: $constant, R: $params> core::ops::SubAssign<Self> for $type<C, R> {
            fn sub_assign(&mut self, other: Self) {
                *self = &*self - &other;
            }
        }

        #[allow(unused_qualifications)]
        impl<'a, const C: $constant, R: $params> core::ops::SubAssign<&'a mut Self>
            for $type<C, R>
        {
            fn sub_assign(&mut self, other: &'a mut Self) {
                *self = &*self - &*other;
            }
        }
    };
}

/// Given the multiplicative operation for two references of a type,
/// implement the additive operations for non-references.
#[macro_export]
macro_rules! impl_multiplicative_ops_from_ref {
    ($type: ident, $params: ident, $constant: ident) => {
        #[allow(unused_qualifications)]
        impl<'a, const C: $constant, R: $params> Mul<R> for &'a $type<C, R> {
            type Output = $type<C, R>;

            fn mul(self, rhs: R) -> Self::Output {
                self * &rhs
            }
        }

        #[allow(unused_qualifications)]
        impl<'a, const C: usize, R: Ring> Mul<&'a R> for $type<C, R> {
            type Output = $type<C, R>;

            fn mul(self, rhs: &'a R) -> Self::Output {
                &self * rhs
            }
        }

        #[allow(unused_qualifications)]
        impl<const C: usize, R: Ring> Mul<R> for $type<C, R> {
            type Output = $type<C, R>;

            fn mul(self, rhs: R) -> Self::Output {
                &self * &rhs
            }
        }

        #[allow(unused_qualifications)]
        impl<const C: $constant, R: $params> core::ops::MulAssign<R> for $type<C, R> {
            fn mul_assign(&mut self, other: R) {
                *self = &*self * &other;
            }
        }

        #[allow(unused_qualifications)]
        impl<'a, const C: $constant, R: $params> core::ops::MulAssign<&'a mut R> for $type<C, R> {
            fn mul_assign(&mut self, other: &'a mut R) {
                *self = &*self * &*other;
            }
        }
    };
}
