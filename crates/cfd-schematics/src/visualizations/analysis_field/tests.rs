//! Value-semantic overlay and color-law integration tests.

use std::{borrow::Cow, collections::HashMap, ptr};

use iris::color::NamedColorMap;

use super::{colorize, AnalysisField, AnalysisOverlay};
use crate::visualizations::traits::Color;

#[test]
fn iris_blue_red_endpoints_cross_the_byte_boundary_exactly() {
    assert_eq!(
        colorize(0.0, NamedColorMap::BlueRed).expect("finite endpoint is valid"),
        Color::rgb(0, 0, 255)
    );
    assert_eq!(
        colorize(1.0, NamedColorMap::BlueRed).expect("finite endpoint is valid"),
        Color::rgb(255, 0, 0)
    );
}

#[test]
fn finite_out_of_range_values_clamp_and_non_finite_values_fail() {
    assert_eq!(
        colorize(-1.0, NamedColorMap::BlueRed).expect("finite values clamp"),
        colorize(0.0, NamedColorMap::BlueRed).expect("endpoint is valid")
    );
    assert_eq!(
        colorize(2.0, NamedColorMap::BlueRed).expect("finite values clamp"),
        colorize(1.0, NamedColorMap::BlueRed).expect("endpoint is valid")
    );
    assert!(colorize(f64::NAN, NamedColorMap::BlueRed).is_err());
    assert!(colorize(f64::INFINITY, NamedColorMap::BlueRed).is_err());
}

#[test]
fn grayscale_midpoint_uses_iris_nearest_byte_quantization() {
    assert_eq!(
        colorize(0.5, NamedColorMap::Grayscale).expect("midpoint is valid"),
        Color::rgb(128, 128, 128)
    );
}

#[test]
fn borrowed_maps_retain_identity_and_precompute_ranges() {
    let edge_data = HashMap::from([(3, -2.0), (5, 6.0)]);
    let overlay = AnalysisOverlay::new(AnalysisField::Pressure, NamedColorMap::BlueRed)
        .with_edge_data(Cow::Borrowed(&edge_data))
        .expect("finite edge data is valid");

    assert!(ptr::eq(
        ptr::from_ref(overlay.edge_data()),
        ptr::from_ref(&edge_data)
    ));
    assert_eq!(overlay.edge_range(), (-2.0, 6.0));
    assert_eq!(overlay.edge_color(3), Some(Color::rgb(0, 0, 255)));
    assert_eq!(overlay.edge_color(5), Some(Color::rgb(255, 0, 0)));
}

#[test]
fn constant_fields_map_to_the_color_law_midpoint() {
    let edge_data = HashMap::from([(0, 100.0), (1, 100.0)]);
    let overlay = AnalysisOverlay::new(AnalysisField::Pressure, NamedColorMap::BlueRed)
        .with_edge_data(Cow::Owned(edge_data))
        .expect("finite edge data is valid");

    assert_eq!(overlay.edge_range(), (100.0, 100.0));
    assert_eq!(overlay.edge_color(0), Some(Color::rgb(128, 0, 128)));
}

#[test]
fn non_finite_field_data_is_rejected_at_the_ownership_boundary() {
    let edge_data = HashMap::from([(0, f64::NAN)]);
    let result = AnalysisOverlay::new(AnalysisField::FlowRate, NamedColorMap::Viridis)
        .with_edge_data(Cow::Owned(edge_data));
    assert!(result.is_err());
}

#[test]
fn missing_identifiers_have_no_color() {
    let overlay = AnalysisOverlay::new(AnalysisField::FlowRate, NamedColorMap::Viridis);
    assert_eq!(overlay.edge_color(999), None);
}
