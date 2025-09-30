//! Provide macros to expand the implementation of commitment operations

/// Given the additive operation for two references of a type,
/// implement the additive operations for non-references.
#[macro_export]
macro_rules! impl_additive_ops {
    ($type: ident, $params: ident, $constant: ident) => {
        #[allow(unused_qualifications)]
        impl<'a, R: $params> core::ops::Add<&'a Self> for $type<R> {
            type Output = $type<R>;
            fn add(self, other: &'a Self) -> Self::Output {
                let mut res = self;
                res += other;
                res
            }
        }

        #[allow(unused_qualifications)]
        impl<'a, R: $params> core::ops::Add<$type<R>> for &'a $type<R> {
            type Output = $type<R>;
            fn add(self, rhs: $type<R>) -> Self::Output {
                self.clone() + rhs
            }
        }

        #[allow(unused_qualifications)]
        impl<R: $params> core::ops::Add<$type<R>> for $type<R> {
            type Output = $type<R>;
            fn add(self, rhs: $type<R>) -> Self::Output {
                self + &rhs
            }
        }

        #[allow(unused_qualifications)]
        impl<R: $params> core::ops::AddAssign<Self> for $type<R> {
            fn add_assign(&mut self, other: Self) {
                *self += &other;
            }
        }
    };
}

/// Given the subtractive operation for two references of a type,
/// implement the subtractive operations for non-references.
#[macro_export]
macro_rules! impl_subtractive_ops {
    ($type: ident, $params: ident, $constant: ident) => {
        #[allow(unused_qualifications)]
        impl<'a, R: $params> core::ops::Sub<&'a Self> for $type<R> {
            type Output = $type<R>;
            fn sub(self, other: &'a Self) -> Self::Output {
                let mut res = self;
                res -= other;
                res
            }
        }

        #[allow(unused_qualifications)]
        impl<'a, R: $params> core::ops::Sub<$type<R>> for &'a $type<R> {
            type Output = $type<R>;
            fn sub(self, rhs: $type<R>) -> Self::Output {
                self.clone() - &rhs
            }
        }

        #[allow(unused_qualifications)]
        impl<R: $params> core::ops::Sub<$type<R>> for $type<R> {
            type Output = $type<R>;
            fn sub(self, rhs: $type<R>) -> Self::Output {
                self - &rhs
            }
        }

        #[allow(unused_qualifications)]
        impl<R: $params> core::ops::SubAssign<Self> for $type<R> {
            fn sub_assign(&mut self, other: Self) {
                *self -= &other;
            }
        }
    };
}

/// Given the multiplicative operation for two references of a type,
/// implement the multiplicative operations for non-references.
#[macro_export]
macro_rules! impl_multiplicative_ops {
    ($type: ident, $params: ident, $constant: ident) => {
        #[allow(unused_qualifications)]
        impl<'a, R: $params> core::ops::Mul<&'a R> for $type<R> {
            type Output = $type<R>;
            fn mul(self, other: &'a R) -> Self::Output {
                let mut res = self;
                res *= other;
                res
            }
        }

        #[allow(unused_qualifications)]
        impl<'a, R: $params> core::ops::Mul<R> for &'a $type<R> {
            type Output = $type<R>;
            fn mul(self, rhs: R) -> Self::Output {
                self.clone() * &rhs
            }
        }

        #[allow(unused_qualifications)]
        impl<R: $params> core::ops::Mul<R> for $type<R> {
            type Output = $type<R>;
            fn mul(self, rhs: R) -> Self::Output {
                self * &rhs
            }
        }

        #[allow(unused_qualifications)]
        impl<R: $params> core::ops::MulAssign<R> for $type<R> {
            fn mul_assign(&mut self, other: R) {
                *self *= &other;
            }
        }
    };
}
