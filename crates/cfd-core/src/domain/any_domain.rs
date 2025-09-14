//! AnyDomain enum and conversions

use nalgebra::{Point1, Point2, Point3, RealField};
use serde::{Deserialize, Serialize};
use super::common::Domain;
use super::domain_1d::Domain1D;
use super::domain_2d::Domain2D;
use super::domain_3d::Domain3D;

/// Generic domain that can be 1D, 2D, or 3D
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum AnyDomain<T: RealField + Copy> {
    /// 1D domain
    D1(Domain1D<T>),
    /// 2D domain
    D2(Domain2D<T>),
    /// 3D domain
    D3(Domain3D<T>),
}

impl<T: RealField + Copy> Domain<T> for AnyDomain<T> {
    fn dimension(&self) -> usize {
        match self {
            Self::D1(d) => d.dimension(),
            Self::D2(d) => d.dimension(),
            Self::D3(d) => d.dimension(),
        }
    }

    fn contains_1d(&self, point: &Point1<T>) -> Option<bool> {
        match self {
            Self::D1(d) => d.contains_1d(point),
            _ => None,
        }
    }

    fn contains_2d(&self, point: &Point2<T>) -> Option<bool> {
        match self {
            Self::D2(d) => d.contains_2d(point),
            _ => None,
        }
    }

    fn contains_3d(&self, point: &Point3<T>) -> Option<bool> {
        match self {
            Self::D3(d) => d.contains_3d(point),
            _ => None,
        }
    }

    fn volume(&self) -> T {
        match self {
            Self::D1(d) => d.volume(),
            Self::D2(d) => d.volume(),
            Self::D3(d) => d.volume(),
        }
    }
}

// Ergonomic From implementations for AnyDomain conversions
impl<T: RealField + Copy> From<Domain1D<T>> for AnyDomain<T> {
    fn from(domain: Domain1D<T>) -> Self {
        AnyDomain::D1(domain)
    }
}

impl<T: RealField + Copy> From<Domain2D<T>> for AnyDomain<T> {
    fn from(domain: Domain2D<T>) -> Self {
        AnyDomain::D2(domain)
    }
}

impl<T: RealField + Copy> From<Domain3D<T>> for AnyDomain<T> {
    fn from(domain: Domain3D<T>) -> Self {
        AnyDomain::D3(domain)
    }
}