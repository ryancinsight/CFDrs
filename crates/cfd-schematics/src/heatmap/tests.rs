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
fn write_to_file_writes_a_closed_svg_document() {
    // Unique per process: nextest runs tests in parallel processes, and a
    // fixed name under the shared temp directory is shared mutable state
    // between them -- one process' cleanup deletes another's output.
    let tmp = std::env::temp_dir().join(format!("cfd-well-plate-{}.svg", std::process::id()));
    write_well_plate_diagram_svg(&sample_candidates(), &tmp)
        .expect("the well-plate diagram must be written");
    let svg = std::fs::read_to_string(&tmp).expect("the written diagram must be readable");
    assert!(
        svg.starts_with("<?xml") || svg.starts_with("<svg"),
        "the diagram must be an SVG document, got: {:?}",
        &svg[..svg.len().min(40)]
    );
    assert!(
        svg.contains("</svg>"),
        "the diagram must be a closed SVG document"
    );
    std::fs::remove_file(&tmp).expect("the temporary diagram must be removable");
}
