//! Wave phase direction calculation for bilateral mirror symmetry.

use super::super::ChannelGenerationContext;
use super::SerpentineChannelStrategy;
use crate::geometry::Point2D;
use crate::state_management::bilateral_symmetry::{
    BilateralPhaseDirectionCalculator, BilateralSymmetryConfig, SymmetryContext,
};

impl SerpentineChannelStrategy {
    /// Calculate wave phase direction for perfect bilateral mirror symmetry
    pub(super) fn calculate_wave_phase_direction(
        &self,
        p1: Point2D,
        p2: Point2D,
        context: &ChannelGenerationContext,
    ) -> f64 {
        if self.config.wave_phase_direction.abs() > 1e-6 {
            return self.config.wave_phase_direction;
        }

        let symmetry_config = BilateralSymmetryConfig {
            enable_vertical_symmetry: true,
            enable_horizontal_symmetry: true,
            symmetry_tolerance: 1e-6,
            enable_adaptive_symmetry: true,
            enforcement_strength: 1.0,
        };

        let temp_context = crate::state_management::adaptive::ChannelGenerationContext::new(
            *context.geometry_config,
            context.box_dims,
            context.total_branches,
            context.neighbor_info,
        )
        .with_endpoints(p1, p2);

        let symmetry_context = SymmetryContext::new(temp_context, symmetry_config);

        let phase_calculator = BilateralPhaseDirectionCalculator::default();

        match phase_calculator.calculate_phase_direction(&symmetry_context) {
            Ok(phase_direction) => phase_direction,
            Err(_) => self.calculate_wave_phase_direction_legacy(p1, p2, context.box_dims),
        }
    }

    /// Legacy fallback for phase direction that preserves centerline mirror behavior.
    fn calculate_wave_phase_direction_legacy(
        &self,
        p1: Point2D,
        p2: Point2D,
        box_dims: (f64, f64),
    ) -> f64 {
        let center_y = box_dims.1 / 2.0;
        let channel_center_y = f64::midpoint(p1.1, p2.1);
        let tolerance = box_dims.1 * 0.01;

        if (channel_center_y - center_y).abs() <= tolerance {
            0.0
        } else if channel_center_y > center_y {
            1.0
        } else {
            -1.0
        }
    }
}
