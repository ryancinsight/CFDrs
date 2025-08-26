//! Basic quadrature rules for numerical integration

use crate::integration::traits::Quadrature;
use cfd_core::error::{Error, NumericalErrorKind, Result};
use nalgebra::RealField;
use num_traits::cast::FromPrimitive;

/// Trapezoidal rule for numerical integration
pub struct TrapezoidalRule;

impl<T: RealField + From<f64> + FromPrimitive + Copy> Quadrature<T> for TrapezoidalRule {
    fn integrate<F>(&self, f: F, a: T, b: T) -> T
    where
        F: Fn(T) -> T,
    {
        let two = T::from_f64(2.0).unwrap_or_else(|| T::zero());
        (b - a) * (f(a) + f(b)) / two
    }

    fn order(&self) -> usize {
        2
    }

    fn num_points(&self) -> usize {
        2
    }
}

/// Simpson's rule for numerical integration
pub struct SimpsonsRule;

impl<T: RealField + From<f64> + FromPrimitive + Copy> Quadrature<T> for SimpsonsRule {
    fn integrate<F>(&self, f: F, a: T, b: T) -> T
    where
        F: Fn(T) -> T,
    {
        let two = T::from_f64(2.0).unwrap_or_else(|| T::zero());
        let four = T::from_f64(4.0).unwrap_or_else(|| T::zero());
        let six = T::from_f64(6.0).unwrap_or_else(|| T::zero());

        let mid = (a + b) / two;
        (b - a) * (f(a) + four * f(mid) + f(b)) / six
    }

    fn order(&self) -> usize {
        4
    }

    fn num_points(&self) -> usize {
        3
    }
}

/// Gauss-Legendre quadrature
pub struct GaussQuadrature<T: RealField + Copy> {
    points: Vec<T>,
    weights: Vec<T>,
    order: usize,
}

impl<T: RealField + Copy + FromPrimitive> GaussQuadrature<T> {
    /// Create a new Gauss quadrature rule with specified order
    pub fn new(order: usize) -> Result<Self> {
        if order == 0 || order > 5 {
            return Err(Error::InvalidInput(format!(
                "Gauss quadrature order must be between 1 and 5, got {}",
                order
            )));
        }

        let (points, weights) = Self::get_points_and_weights(order)?;

        Ok(Self {
            points,
            weights,
            order,
        })
    }

    /// Get Gauss-Legendre points and weights for given order
    fn get_points_and_weights(order: usize) -> Result<(Vec<T>, Vec<T>)> {
        match order {
            1 => Self::gauss_1_point(),
            2 => Self::gauss_2_point(),
            3 => Self::gauss_3_point(),
            4 => Self::gauss_4_point(),
            5 => Self::gauss_5_point(),
            _ => Err(Error::InvalidInput(format!(
                "Unsupported Gauss quadrature order: {}",
                order
            ))),
        }
    }

    fn gauss_1_point() -> Result<(Vec<T>, Vec<T>)> {
        Ok((
            vec![T::zero()],
            vec![T::from_f64(2.0).ok_or_else(|| {
                Error::Numerical(NumericalErrorKind::InvalidValue {
                    value: "Cannot convert weight 2.0".to_string(),
                })
            })?],
        ))
    }

    fn gauss_2_point() -> Result<(Vec<T>, Vec<T>)> {
        let sqrt3_inv = T::from_f64(1.0 / 3.0_f64.sqrt()).ok_or_else(|| {
            Error::Numerical(NumericalErrorKind::InvalidValue {
                value: "Cannot convert sqrt(1/3)".to_string(),
            })
        })?;
        Ok((vec![-sqrt3_inv, sqrt3_inv], vec![T::one(), T::one()]))
    }

    fn gauss_3_point() -> Result<(Vec<T>, Vec<T>)> {
        let sqrt15 = T::from_f64(15.0_f64.sqrt()).ok_or_else(|| {
            Error::Numerical(NumericalErrorKind::InvalidValue {
                value: "Cannot convert sqrt(15)".to_string(),
            })
        })?;
        let five = T::from_f64(5.0).ok_or_else(|| {
            Error::Numerical(NumericalErrorKind::InvalidValue {
                value: "Cannot convert 5".to_string(),
            })
        })?;
        let nine = T::from_f64(9.0).ok_or_else(|| {
            Error::Numerical(NumericalErrorKind::InvalidValue {
                value: "Cannot convert 9".to_string(),
            })
        })?;
        let eight = T::from_f64(8.0).ok_or_else(|| {
            Error::Numerical(NumericalErrorKind::InvalidValue {
                value: "Cannot convert 8".to_string(),
            })
        })?;

        Ok((
            vec![-sqrt15 / five, T::zero(), sqrt15 / five],
            vec![five / nine, eight / nine, five / nine],
        ))
    }

    fn gauss_4_point() -> Result<(Vec<T>, Vec<T>)> {
        let sqrt6_5 = (6.0 / 5.0_f64).sqrt();
        let term1 = T::from_f64((3.0 - 2.0 * sqrt6_5) / 7.0)
            .ok_or_else(|| {
                Error::Numerical(NumericalErrorKind::InvalidValue {
                    value: "Cannot convert term1".to_string(),
                })
            })?
            .sqrt();
        let term2 = T::from_f64((3.0 + 2.0 * sqrt6_5) / 7.0)
            .ok_or_else(|| {
                Error::Numerical(NumericalErrorKind::InvalidValue {
                    value: "Cannot convert term2".to_string(),
                })
            })?
            .sqrt();

        let sqrt30 = 30.0_f64.sqrt();
        let w1 = T::from_f64((18.0 + sqrt30) / 36.0).ok_or_else(|| {
            Error::Numerical(NumericalErrorKind::InvalidValue {
                value: "Cannot convert weight1".to_string(),
            })
        })?;
        let w2 = T::from_f64((18.0 - sqrt30) / 36.0).ok_or_else(|| {
            Error::Numerical(NumericalErrorKind::InvalidValue {
                value: "Cannot convert weight2".to_string(),
            })
        })?;

        Ok((vec![-term2, -term1, term1, term2], vec![w2, w1, w1, w2]))
    }

    fn gauss_5_point() -> Result<(Vec<T>, Vec<T>)> {
        // 5-point Gauss-Legendre quadrature
        let sqrt10_7 = T::from_f64((10.0 / 7.0_f64).sqrt()).ok_or_else(|| {
            Error::Numerical(NumericalErrorKind::InvalidValue {
                value: "Cannot convert sqrt(10/7)".to_string(),
            })
        })?;
        let seven = T::from_f64(7.0).ok_or_else(|| {
            Error::Numerical(NumericalErrorKind::InvalidValue {
                value: "Cannot convert 7".to_string(),
            })
        })?;
        let five = T::from_f64(5.0).ok_or_else(|| {
            Error::Numerical(NumericalErrorKind::InvalidValue {
                value: "Cannot convert 5".to_string(),
            })
        })?;
        let two = T::from_f64(2.0).ok_or_else(|| {
            Error::Numerical(NumericalErrorKind::InvalidValue {
                value: "Cannot convert 2".to_string(),
            })
        })?;
        let one_third = T::from_f64(1.0 / 3.0).ok_or_else(|| {
            Error::Numerical(NumericalErrorKind::InvalidValue {
                value: "Cannot convert 1/3".to_string(),
            })
        })?;

        let term = two * sqrt10_7 / seven;

        let x1 = one_third * (five - term).sqrt();
        let x2 = one_third * (five + term).sqrt();

        let w0 = T::from_f64(128.0 / 225.0).ok_or_else(|| {
            Error::Numerical(NumericalErrorKind::InvalidValue {
                value: "Cannot convert 128/225".to_string(),
            })
        })?;
        let w1 = T::from_f64((322.0 + 13.0 * 70.0_f64.sqrt()) / 900.0).ok_or_else(|| {
            Error::Numerical(NumericalErrorKind::InvalidValue {
                value: "Cannot convert weight".to_string(),
            })
        })?;
        let w2 = T::from_f64((322.0 - 13.0 * 70.0_f64.sqrt()) / 900.0).ok_or_else(|| {
            Error::Numerical(NumericalErrorKind::InvalidValue {
                value: "Cannot convert weight".to_string(),
            })
        })?;

        Ok((vec![-x2, -x1, T::zero(), x1, x2], vec![w2, w1, w0, w1, w2]))
    }

    /// Create default 2-point Gauss quadrature
    pub fn default() -> Result<Self> {
        Self::new(2)
    }
}

impl<T: RealField + From<f64> + FromPrimitive + Copy> Quadrature<T> for GaussQuadrature<T> {
    fn integrate<F>(&self, f: F, a: T, b: T) -> T
    where
        F: Fn(T) -> T,
    {
        let two = T::from_f64(2.0).unwrap_or_else(|| T::zero());
        let half_interval = (b - a) / two;
        let mid_point = (a + b) / two;

        // Use iterator combinators for zero-copy optimization
        let result = self
            .points
            .iter()
            .zip(self.weights.iter())
            .map(|(point, weight)| {
                let x = mid_point + half_interval * *point;
                *weight * f(x)
            })
            .fold(T::zero(), |acc, term| acc + term);

        result * half_interval
    }

    fn order(&self) -> usize {
        2 * self.order
    }

    fn num_points(&self) -> usize {
        self.points.len()
    }
}
