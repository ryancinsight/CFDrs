//! Bifurcation → Serpentine channel swap example.
//!
//! Uses `create_geometry` with three `SplitType::Bifurcation` levels
//! (identical to the double/triple bifurcation in `comprehensive_split_patterns`)
//! to produce a 3-level binary tree with 8 horizontal leaf channels.
//!
//! Produces 9 SVG schematics: one base bifurcation and eight variants where
//! a single leaf channel has been swapped to serpentine.
//!
//! # Run
//! ```bash
//! cargo run -p cfd-schematics --example bifurcation_serpentine_swap --no-default-features
//! ```
//!
//! Outputs are written to `report/serpentine_swap/`.

use cfd_schematics::config::{ChannelTypeConfig, GeometryConfig};
use cfd_schematics::domain::model::{ChannelShape, NetworkBlueprint};
use cfd_schematics::geometry::generator::create_geometry;
use cfd_schematics::geometry::SplitType;
use std::path::Path;

#[path = "../shared/mod.rs"]
mod shared;

const SERP_SEGMENTS: usize = 6;
const SERP_BEND_RADIUS_M: f64 = 1.2e-3;

fn main() {
    let output_dir = Path::new("report/serpentine_swap");
    std::fs::create_dir_all(output_dir).expect("failed to create output directory");

    // Generate a 3-level bifurcation (8 leaf channels) using the same
    // code path as comprehensive_split_patterns.
    let box_dims = (400.0, 200.0);
    let splits = vec![
        SplitType::Bifurcation,
        SplitType::Bifurcation,
        SplitType::Bifurcation,
    ];
    let base = create_geometry(
        box_dims,
        &splits,
        &GeometryConfig::default(),
        &ChannelTypeConfig::default(),
    );

    println!(
        "Base blueprint: {} nodes, {} channels",
        base.nodes.len(),
        base.channels.len(),
    );

    // Identify the 8 horizontal leaf channels at the center of the tree.
    // Each leaf has two segments (left half + right half) sharing the same y.
    let leaf_pairs = find_leaf_channel_pairs(&base);
    println!("Found {} leaf channel pairs to swap", leaf_pairs.len());


    // Render the base (all straight).
    shared::save_example_output_with_name(
        &base,
        "bifurcation_serpentine_swap",
        "00_base_bifurcation",
    );

    // For each leaf pair, swap both segments to serpentine and render.
    for (i, &(left_idx, right_idx)) in leaf_pairs.iter().enumerate() {
        let mut mutated = base.clone();
        let serp = ChannelShape::Serpentine {
            segments: SERP_SEGMENTS,
            bend_radius_m: SERP_BEND_RADIUS_M,
            wave_type: cfd_schematics::SerpentineWaveType::default(),
        };
        mutated.channels[left_idx].channel_shape = serp;
        mutated.channels[right_idx].channel_shape = serp;

        let filename = format!("{:02}_leaf{}_serpentine", i + 1, i + 1);
        shared::save_example_output_with_name(&mutated, "bifurcation_serpentine_swap", &filename);
    }

    println!(
        "\nDone — {} schematics in {}",
        leaf_pairs.len() + 1,
        output_dir.display(),
    );
}

/// Find the horizontal leaf channels (treatment zone) in a generated bifurcation.
///
/// Leaf channels are the innermost horizontal segments that touch the center
/// of the box (x = `box_width` / 2). Each leaf has two segments: one on the
/// left half and one on the right half, sharing the same y-coordinate.
/// Returns pairs of channel indices (left, right) for each leaf, sorted
/// by y-position.
fn find_leaf_channel_pairs(bp: &NetworkBlueprint) -> Vec<(usize, usize)> {
    use std::collections::HashMap;
    let node_pos: HashMap<&str, (f64, f64)> =
        bp.nodes.iter().map(|n| (n.id.as_str(), n.point)).collect();

    let half_x = bp.box_dims.0 * 0.5;
    let tol = 1e-3;

    // Find all horizontal channels that touch the center line.
    // Key by rounded y-coordinate to group left/right pairs.
    let mut by_y: HashMap<i64, Vec<usize>> = HashMap::new();

    for (i, ch) in bp.channels.iter().enumerate() {
        let from = node_pos.get(ch.from.as_str());
        let to = node_pos.get(ch.to.as_str());
        if let (Some(&(fx, fy)), Some(&(tx, ty))) = (from, to) {
            let horizontal = (fy - ty).abs() < tol;
            let touches_center = (fx - half_x).abs() < tol || (tx - half_x).abs() < tol;
            if horizontal && touches_center {
                #[allow(clippy::cast_possible_truncation)]
                let y_key = (fy * 1e6) as i64;
                by_y.entry(y_key).or_default().push(i);
            }
        }
    }

    // Each leaf y-level should have exactly 2 channels (left + right).
    let mut pairs: Vec<(i64, usize, usize)> = by_y
        .into_iter()
        .filter(|(_, indices)| indices.len() == 2)
        .map(|(y_key, indices)| {
            let (a, b) = (indices[0], indices[1]);
            // Left channel has from.x < half_x, right has from.x >= half_x.
            let fa = node_pos.get(bp.channels[a].from.as_str()).unwrap().0;
            if fa < half_x {
                (y_key, a, b)
            } else {
                (y_key, b, a)
            }
        })
        .collect();

    pairs.sort_by_key(|(y_key, _, _)| *y_key);
    pairs.into_iter().map(|(_, l, r)| (l, r)).collect()
}
