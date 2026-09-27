//! Arc channel generation strategy.

mod curvature;
mod figure_eight;
mod geometry;
mod path;
mod strategy;
mod symmetry;

use crate::config::{ArcConfig, ConstantsRegistry, GeometryConfig};
use crate::geometry::{ChannelType, Point2D};

use super::ChannelTypeStrategy;

pub use strategy::ArcChannelStrategy;
