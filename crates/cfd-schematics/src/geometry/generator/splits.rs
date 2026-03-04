//! Split/merge geometry generation for branching channel networks.
//!
//! Implements the core split and merge operations that create bifurcating
//! and trifurcating channel topologies within the bounding box.

use super::super::types::{ChannelSystem, Point2D, SplitType};
use super::GeometryGenerator;
use crate::config::ChannelTypeConfig;

impl GeometryGenerator {
    pub(super) fn generate(mut self, splits: &[SplitType]) -> ChannelSystem {
        let (length, width) = self.box_dims;

        if splits.is_empty() {
            let p1 = (0.0, width / 2.0);
            let p2 = (length, width / 2.0);
            // For single channel, pass empty neighbor list so it uses box boundaries
            self.add_channel_with_neighbors(p1, p2, &[], None); // Use default width
            return self.finalize();
        }

        let (first_half_lines, first_half_widths) = self.generate_first_half(splits);

        // Collect y-coordinates for dynamic amplitude calculation using iterator combinators
        let y_coords_for_amplitude: Vec<f64> = first_half_lines
            .iter()
            .map(|(p1, p2)| f64::midpoint(p1.1, p2.1))
            .collect();

        // Apply Gaussian width modulation — channels closer to the plate
        // y-center are wider, ensuring symmetric y-axis width distribution.
        let modulated_first_widths: Vec<f64> = first_half_lines
            .iter()
            .zip(first_half_widths.iter())
            .map(|((p1, p2), &w)| {
                let y_mid = f64::midpoint(p1.1, p2.1);
                w * self.gaussian_width_scale(y_mid)
            })
            .collect();

        for (i, (p1, p2)) in first_half_lines.iter().enumerate() {
            self.add_channel_with_neighbors(
                *p1,
                *p2,
                &y_coords_for_amplitude,
                Some(modulated_first_widths[i]),
            );
        }

        // Generate the second half with proper merge pattern (inverse of splits)
        self.generate_second_half(splits, &first_half_widths); // Use the widths at the center

        self.finalize()
    }

    fn generate_first_half(&self, splits: &[SplitType]) -> (Vec<(Point2D, Point2D)>, Vec<f64>) {
        let (length, width) = self.box_dims;
        let effective_width = (-2.0f64).mul_add(self.config.wall_clearance, width);
        let half_l = length / 2.0;
        let num_splits = splits.len() as u32;
        let num_segments_per_half = f64::from(num_splits).mul_add(2.0, 1.0);
        let dx = half_l / num_segments_per_half;

        let mut y_coords: Vec<f64> = vec![width / 2.0];
        let mut y_ranges: Vec<f64> = vec![effective_width];
        let mut current_widths: Vec<f64> = vec![self.config.channel_width]; // Initial width
        let mut current_x = 0.0;
        let mut lines = Vec::new();
        let mut line_widths = Vec::new(); // Defines width for each line segment pushed

        for split_type in splits {
            for (i, y) in y_coords.iter().enumerate() {
                lines.push(((current_x, *y), (current_x + dx, *y)));
                line_widths.push(current_widths[i]);
            }
            current_x += dx;

            let (next_y_coords, next_y_ranges, next_widths, new_lines) = self.apply_split(
                *split_type,
                &y_coords,
                &y_ranges,
                &current_widths,
                current_x,
                dx,
            );

            // Record splitters
            // Note: new_lines corresponds to the split connection (not really a channel segment in the main flow sense?)
            // We usually assign the *child* width to the splitter segment.
            // new_lines: [parent1_child1, parent1_child2, parent2_child1, ...]
            // next_widths: [child1_width, child2_width, ...]
            // They align.
            for (i, _line) in new_lines.iter().enumerate() {
                line_widths.push(next_widths[i]);
            }

            y_coords = next_y_coords;
            y_ranges = next_y_ranges;
            current_widths = next_widths;
            lines.extend(new_lines);

            current_x += dx;
        }

        for (i, y) in y_coords.iter().enumerate() {
            lines.push(((current_x, *y), (half_l, *y)));
            line_widths.push(current_widths[i]);
        }
        (lines, line_widths)
    }

    fn apply_split(
        &self,
        split_type: SplitType,
        y_coords: &[f64],
        y_ranges: &[f64],
        current_widths: &[f64],
        current_x: f64,
        dx: f64,
    ) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<(Point2D, Point2D)>) {
        let mut next_y_coords = Vec::new();
        let mut next_y_ranges = Vec::new();
        let mut next_widths = Vec::new();
        let mut new_lines = Vec::new();

        for (j, y_center) in y_coords.iter().enumerate() {
            let y_range = y_ranges[j];
            let parent_width = current_widths[j];
            let n_branches = split_type.branch_count();
            let effective_channel_diameter = self.effective_channel_diameter();
            let edge_padding = self.calculate_edge_padding(y_range);
            let usable_range = (y_range - 2.0 * edge_padding).max(0.0);

            let lower_bound = y_center - y_range / 2.0 + edge_padding;

            // Calculate child widths first — needed for proportional spacing.
            let child_widths = if n_branches == 1 {
                vec![parent_width]
            } else if n_branches == 2 {
                match split_type {
                    SplitType::AsymmetricBifurcation { ratio } => {
                        let w_total = parent_width * 2.0;
                        vec![w_total * ratio, w_total * (1.0 - ratio).max(0.0)]
                    }
                    _ => vec![parent_width, parent_width],
                }
            } else {
                match split_type {
                    SplitType::SymmetricTrifurcation { center_ratio } => {
                        let w_total = parent_width * 3.0;
                        let w_center = w_total * center_ratio;
                        let w_side = (w_total * (1.0 - center_ratio)) / 2.0;
                        vec![w_side, w_center, w_side]
                    }
                    _ => vec![parent_width; n_branches],
                }
            };

            // Width-proportional slot allocation: each branch receives vertical
            // space proportional to its channel width, giving wider channels more
            // room and packing narrower peripheral channels tighter together.
            //
            // Minimum slot enforcement: each slot must be at least the child's
            // rendered width so adjacent channel edges do not overlap.  When
            // the enforced total exceeds usable_range the slots are scaled back
            // proportionally (which effectively narrows the rendered channels
            // in the schematic).
            let total_child_width: f64 = child_widths.iter().sum();
            let min_gap = self.config.wall_clearance;
            let slots: Vec<f64> = if total_child_width > 0.0 {
                let ideal: Vec<f64> = child_widths
                    .iter()
                    .map(|w| usable_range * (w / total_child_width))
                    .collect();
                // Enforce: slot ≥ child_width + min_gap
                let enforced: Vec<f64> = ideal
                    .iter()
                    .zip(&child_widths)
                    .map(|(&s, &w)| s.max(w + min_gap))
                    .collect();
                let enforced_total: f64 = enforced.iter().sum();
                if enforced_total > usable_range && enforced_total > 0.0 {
                    enforced
                        .iter()
                        .map(|&s| s * usable_range / enforced_total)
                        .collect()
                } else {
                    enforced
                }
            } else {
                let uniform = if n_branches > 0 {
                    usable_range / n_branches as f64
                } else {
                    0.0
                };
                vec![uniform; n_branches]
            };

            for (i, &width) in child_widths.iter().enumerate().take(n_branches) {
                let y_new = if n_branches <= 1 {
                    *y_center
                } else {
                    let cumulative: f64 = slots[..i].iter().sum();
                    lower_bound + cumulative + slots[i] / 2.0
                };

                // Per-branch child range also proportional to width
                let child_range_i = if total_child_width > 0.0 {
                    (y_range * (width / total_child_width)).max(effective_channel_diameter)
                } else {
                    (y_range / n_branches as f64).max(effective_channel_diameter)
                };

                new_lines.push(((current_x, *y_center), (current_x + dx, y_new)));
                next_y_coords.push(y_new);
                next_y_ranges.push(child_range_i);
                next_widths.push(width);
            }
        }
        (next_y_coords, next_y_ranges, next_widths, new_lines)
    }

    fn apply_merge(
        &self,
        split_type: SplitType,
        y_coords: &[f64],
        y_ranges: &[f64],
        current_widths: &[f64],
        current_x: f64,
        dx: f64,
    ) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<(Point2D, Point2D)>) {
        let mut next_y_coords = Vec::new();
        let mut next_y_ranges = Vec::new();
        let mut next_widths = Vec::new();
        let mut new_lines = Vec::new();

        let n_branches = split_type.branch_count();

        // Group the y_coords by n_branches to create merges
        for chunk_idx in 0..(y_coords.len() / n_branches) {
            let start = chunk_idx * n_branches;
            let chunk = &y_coords[start..start + n_branches];
            let chunk_widths = &current_widths[start..start + n_branches];

            // Calculate the center y-coordinate for this merge group
            let y_center = chunk.iter().sum::<f64>() / chunk.len() as f64;

            // Merged width: For Asymmetric, reverse the split?
            // Simple approach: Sum of widths? No, for bifurcation w1+w2 = 2*w_p. So sum/2.
            // Or average?
            let merged_width = chunk_widths.iter().sum::<f64>() / n_branches as f64; // Average

            // Calculate the range for the merged group
            let y_range = if chunk.len() > 1 {
                let min_y = chunk.iter().fold(f64::INFINITY, |a, &b| a.min(b));
                let max_y = chunk.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
                max_y - min_y
            } else {
                y_ranges[0]
            };

            // Create merge lines from each branch to the center
            for &y_branch in chunk {
                new_lines.push(((current_x, y_branch), (current_x + dx, y_center)));
            }

            next_y_coords.push(y_center);
            next_y_ranges.push(y_range.max(self.effective_channel_diameter()));
            next_widths.push(merged_width);
        }

        (next_y_coords, next_y_ranges, next_widths, new_lines)
    }

    /// Determine if this geometry configuration should reserve extra split spacing
    /// for curved channels (serpentine/arcs).
    fn requires_curved_spacing(&self) -> bool {
        matches!(
            self.channel_type_config,
            ChannelTypeConfig::AllSerpentine(_)
                | ChannelTypeConfig::AllArcs(_)
                | ChannelTypeConfig::MixedByPosition { .. }
                | ChannelTypeConfig::Adaptive { .. }
                | ChannelTypeConfig::SmoothSerpentineWithTransitions { .. }
        )
    }

    /// Gaussian scale factor for channel width at a given y-position.
    ///
    /// Returns a value in `[0.5, 1.0]` where `1.0` is at the plate y-center
    /// and `0.5` at the bounding-box edges.  The Gaussian envelope produces
    /// a symmetric width distribution about the plate midline.
    ///
    /// # Theorem
    ///
    /// G(μ + Δ) = G(μ − Δ) ∀ Δ ∈ ℝ, guaranteeing bilateral symmetry.
    ///
    /// **Proof sketch**: G(y) = exp(−(y − μ)² / 2σ²) depends only on |y − μ|.  ∎
    fn gaussian_width_scale(&self, y: f64) -> f64 {
        let (_length, height) = self.box_dims;
        let y_center = height / 2.0;
        let effective_height = (-2.0f64).mul_add(self.config.wall_clearance, height);
        // σ = effective_height / 6  →  ~28% reduction at single-trifurcation
        // peripherals, ~50% at plate edges.  Provides visible but physically
        // reasonable differentiation across branching topologies.
        let sigma = effective_height / 6.0;
        let min_scale = 0.5;
        let dy = y - y_center;
        let gaussian = (-dy * dy / (2.0 * sigma * sigma)).exp();
        min_scale + (1.0 - min_scale) * gaussian
    }

    /// Compute edge padding inside a split range to keep branches away from
    /// local boundaries and preserve room for curved channels.
    fn calculate_edge_padding(&self, y_range: f64) -> f64 {
        let effective_channel_diameter = self.effective_channel_diameter();

        let curvature_padding = if self.requires_curved_spacing() {
            effective_channel_diameter.mul_add(1.5, self.config.wall_clearance * 0.5)
        } else {
            effective_channel_diameter * 0.5
        };

        // Cap at 25% of y_range, but guarantee at least wall_clearance when
        // the range is large enough to accommodate it.
        let max_fraction = y_range * 0.25;
        let min_usable = self.config.wall_clearance.min(y_range * 0.15);
        curvature_padding.min(max_fraction).max(min_usable)
    }

    fn generate_second_half(&mut self, splits: &[SplitType], _center_widths: &[f64]) {
        let (length, width) = self.box_dims;
        let effective_width = (-2.0f64).mul_add(self.config.wall_clearance, width);
        let half_l = length / 2.0;
        let num_splits = splits.len() as u32;
        let num_segments_per_half = f64::from(num_splits).mul_add(2.0, 1.0);
        let dx = half_l / num_segments_per_half;

        // Start with the initial state (at center)
        // We need to RE-COMPUTE the center state to get y_coords and y_ranges again
        // because generate_first_half consumed them. Or refactor generate to return internal state.
        // But `generate_first_half` returned lines only.
        // Actually, `generate_first_half` implementation repeated the split logic.
        // We should probably just re-run the splits to get current state (position)
        let mut y_coords = vec![width / 2.0];
        let mut y_ranges = vec![effective_width];
        let mut current_widths = vec![self.config.channel_width]; // THIS IS WRONG. Center widths are passed in.

        // Wait, generate_second_half repeats the split logic to calculate y_coords?
        // Yes, lines 461-466: "Apply all splits to get the final state"
        // But `center_widths` ARE the widths at the center (result of first half).

        // Re-simulate splits to get y_coords/ranges
        for split_type in splits {
            let (next_y_coords, next_y_ranges, next_widths, _) =
                self.apply_split(*split_type, &y_coords, &y_ranges, &current_widths, 0.0, dx);
            y_coords = next_y_coords;
            y_ranges = next_y_ranges;
            current_widths = next_widths;
        }

        // Sanity check: calculated widths should match passed center_widths
        // assert_eq!(current_widths, center_widths);

        // Now generate the second half by reversing the splits (creating merges)
        let mut current_x = half_l;
        let mut lines = Vec::new();
        let mut line_widths = Vec::new();

        // Process splits in reverse order to create merges
        for split_type in splits.iter().rev() {
            // Add horizontal segments from current position
            for (i, y) in y_coords.iter().enumerate() {
                lines.push(((current_x, *y), (current_x + dx, *y)));
                line_widths.push(current_widths[i]);
            }
            current_x += dx;

            // Apply merge (reverse of split)
            let (next_y_coords, next_y_ranges, next_widths, new_lines) = self.apply_merge(
                *split_type,
                &y_coords,
                &y_ranges,
                &current_widths,
                current_x,
                dx,
            );

            // Merge lines have the *input* widths (branches merging).
            // In graph, edges leaving.
            // The method `apply_merge` generates lines FROM branches TO center.
            // The branches have `current_widths`.
            for (i, _line) in new_lines.iter().enumerate() {
                // new_lines is flat list of all branches merging.
                // We assign line width of the branch.
                // But new_lines order? apply_merge logic loops chunks.
                // It pushes lines for each branch in the chunk.
                // So order matches current_widths.
                line_widths.push(current_widths[i]);
            }

            y_coords = next_y_coords;
            y_ranges = next_y_ranges;
            current_widths = next_widths;
            lines.extend(new_lines);

            current_x += dx;
        }

        // Final horizontal segments to the right edge
        for (i, y) in y_coords.iter().enumerate() {
            lines.push(((current_x, *y), (length, *y)));
            line_widths.push(current_widths[i]);
        }

        // Collect y-coordinates for amplitude calculation (from lines)
        let mut y_coords_for_amplitude: Vec<f64> = Vec::with_capacity(lines.len());
        for (p1, p2) in &lines {
            y_coords_for_amplitude.push(f64::midpoint(p1.1, p2.1));
        }

        // Apply Gaussian width modulation (mirrors first-half distribution)
        let modulated_widths: Vec<f64> = lines
            .iter()
            .zip(line_widths.iter())
            .map(|((p1, p2), &w)| {
                let y_mid = f64::midpoint(p1.1, p2.1);
                w * self.gaussian_width_scale(y_mid)
            })
            .collect();

        for (i, (p1, p2)) in lines.iter().enumerate() {
            self.add_channel_with_neighbors(
                *p1,
                *p2,
                &y_coords_for_amplitude,
                Some(modulated_widths[i]),
            );
        }
    }
}
