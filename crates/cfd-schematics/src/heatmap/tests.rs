use super::palette::cancer_cav_color;
use super::svg::build_svg;

use super::*;

fn sample_candidates() -> Vec<CandidateZoneData> {
    vec![
        CandidateZoneData {
            label: "CCT-lv3".into(),
            cancer_cav: 0.78,
            lysis_risk: 0.001,
            therapy_frac: 0.88,
        },
        CandidateZoneData {
            label: "TBV-cf55".into(),
            cancer_cav: 0.65,
            lysis_risk: 0.002,
            therapy_frac: 0.72,
        },
    ]
}

#[test]
fn svg_starts_with_xml_declaration() {
    let svg = build_svg(&sample_candidates());
    assert!(
        svg.starts_with("<?xml"),
        "expected XML declaration at start"
    );
}

#[test]
fn svg_contains_svg_element() {
    let svg = build_svg(&sample_candidates());
    assert!(svg.contains("<svg "), "expected <svg element");
    assert!(svg.contains("</svg>"), "expected closing </svg>");
    assert!(svg.contains("preserveAspectRatio=\"xMidYMin meet\""));
    assert!(svg.contains("width:min(100%, 100vw, calc(100vh * "));
    assert!(svg.contains("height:auto;display:block;margin:0 auto;"));
}

#[test]
fn svg_uses_split_title_for_plate_heading() {
    let svg = build_svg(&sample_candidates());
    assert!(svg.contains("SDT Millifluidic Device"));
    assert!(svg.contains("96-Well Plate Treatment Zone"));
}

#[test]
fn svg_contains_treatment_zone_rect() {
    let svg = build_svg(&sample_candidates());
    // The dashed treatment zone rect must be present
    assert!(
        svg.contains("stroke-dasharray"),
        "expected dashed border for treatment zone"
    );
}

#[test]
fn svg_treatment_zone_rect_covers_full_six_by_six_envelope() {
    let svg = build_svg(&sample_candidates());
    assert!(
        svg.contains(r#"<rect x="283.3" y="164.4" width="324.0" height="324.0""#),
        "expected treatment zone rect to span the full 6x6 well envelope"
    );
}

#[test]
fn svg_contains_candidate_labels() {
    let svg = build_svg(&sample_candidates());
    assert!(
        svg.contains("CCT-lv3"),
        "expected first candidate label in SVG"
    );
    assert!(
        svg.contains("TBV-cf55"),
        "expected second candidate label in SVG"
    );
}

#[test]
fn svg_empty_candidates_still_valid() {
    let svg = build_svg(&[]);
    assert!(svg.contains("<svg "));
    assert!(svg.contains("</svg>"));
    assert!(!svg.contains("CCT")); // no candidate content
}

#[test]
fn cancer_cav_color_extremes() {
    let low = cancer_cav_color(0.0);
    let high = cancer_cav_color(1.0);
    // Low should be yellowish (#FFD700 area), high should be reddish (#CC0000)
    assert!(low.starts_with('#'), "expected hex color");
    assert!(high.starts_with('#'), "expected hex color");
    assert_ne!(low, high, "colors at 0 and 1 should differ");
}

#[test]
fn writer_serializes_the_candidate_diagram() {
    // Unique per process: nextest runs tests in parallel processes, and a
    // fixed name under the shared temp directory is shared mutable state
    // between them -- one process' cleanup deletes another's output.
    let process_id = std::process::id();
    let output_path = std::env::temp_dir().join(format!("cfd-well-plate-{process_id}-primary.svg"));
    let alternate_output_path =
        std::env::temp_dir().join(format!("cfd-well-plate-{process_id}-alternate.svg"));

    let candidates = sample_candidates();
    let mut alternate_candidates = candidates.clone();
    alternate_candidates
        .first_mut()
        .expect("invariant: sample_candidates returns at least one candidate")
        .label = "ALT-1".into();

    write_well_plate_diagram_svg(&candidates, &output_path)
        .expect("the well-plate diagram must be written");
    write_well_plate_diagram_svg(&alternate_candidates, &alternate_output_path)
        .expect("the alternate well-plate diagram must be written");

    let written_svg =
        std::fs::read_to_string(&output_path).expect("the written diagram must be readable");
    let alternate_written_svg = std::fs::read_to_string(&alternate_output_path)
        .expect("the alternate written diagram must be readable");

    assert_eq!(written_svg, build_svg(&candidates));
    assert_eq!(alternate_written_svg, build_svg(&alternate_candidates));
    assert_ne!(
        written_svg, alternate_written_svg,
        "different candidate inputs must produce different serialized diagrams"
    );

    std::fs::remove_file(&output_path).expect("the temporary diagram must be removable");
    std::fs::remove_file(&alternate_output_path)
        .expect("the alternate temporary diagram must be removable");
}
