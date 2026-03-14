//! Large process-flow SVG figure for the Milestone 12 narrative.

use std::fmt::Write as _;
use std::path::Path;

use super::primitives::{arrow, escape_xml, process_box, svg_end, svg_start, svg_title};

pub(super) fn write_placeholder(
    path: &Path,
    title: &str,
    message: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut svg = String::new();
    svg_start(&mut svg, 1100.0, 420.0);
    svg_title(&mut svg, title);
    let _ = write!(
        svg,
        r##"<rect x="70" y="120" width="960" height="220" fill="#f5f6fa" stroke="#bdc3c7" stroke-width="2"/>"##
    );
    let _ = write!(
        svg,
        r##"<text x="110" y="240" font-size="24" fill="#2c3e50">{}</text>"##,
        escape_xml(message)
    );
    svg_end(&mut svg);
    std::fs::write(path, svg)?;
    Ok(())
}

pub(super) fn write_creation_optimization_process_figure(
    path: &Path,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut svg = String::new();
    let w = 1700.0;
    let h = 780.0;
    svg_start(&mut svg, w, h);
    svg_title(
        &mut svg,
        "Milestone 12 Design Creation & Optimization Process",
    );
    let _ = write!(
        svg,
        r##"<text x="70" y="95" font-size="18" fill="#34495e">Primitive selective-routing pipeline for bifurcation/trifurcation treatment trees.</text>"##
    );

    // Top-row arrows (between boxes 1-5)
    for (x1, y1, x2, y2) in [
        (320.0, 205.0, 370.0, 205.0),
        (650.0, 205.0, 700.0, 205.0),
        (980.0, 205.0, 1030.0, 205.0),
        (1310.0, 205.0, 1360.0, 205.0),
        // Vertical: top row center → bottom row
        (840.0, 290.0, 840.0, 410.0),
        // Bottom row: Option 2 → Option 1
        (710.0, 495.0, 635.0, 495.0),
        // Bottom row: Option 1 → Report
        (375.0, 495.0, 310.0, 495.0),
    ] {
        arrow(&mut svg, x1, y1, x2, y2);
    }

    process_box(
        &mut svg,
        40.0,
        120.0,
        280.0,
        170.0,
        "#eef6ff",
        "#2e86de",
        "1. Primitive candidate generation",
        &[
            "cfd-schematics::create_geometry",
            "Mirrored Bi/Tri split trees",
            "Branch-diameter biasing for",
            "RBC/WBC/CTC routing",
        ],
    );
    process_box(
        &mut svg,
        370.0,
        120.0,
        280.0,
        170.0,
        "#eefcf3",
        "#27ae60",
        "2. ChannelSystem topology SSOT",
        &[
            "Wall ports and true split/merge",
            "Treatment/bypass lanes,",
            "serpentines, venturi throats",
            "No report-local topology recreation",
        ],
    );
    process_box(
        &mut svg,
        700.0,
        120.0,
        280.0,
        170.0,
        "#fff7ea",
        "#d68910",
        "3. NetworkBlueprint projection",
        &[
            "Lossless ChannelSystem ->",
            "NetworkBlueprint",
            "Downstream solver/export only",
            "No independent blueprint authoring",
        ],
    );
    process_box(
        &mut svg,
        1030.0,
        120.0,
        280.0,
        170.0,
        "#fef0f3",
        "#c0392b",
        "4. Cached cfd-1d evaluation",
        &[
            "MetricsCache avoids repeat solves",
            "Pressure, shear, residence, cavitation",
            "Three-pop separation and safety",
        ],
    );
    process_box(
        &mut svg,
        1360.0,
        120.0,
        140.0,
        170.0,
        "#f4f1ff",
        "#7d3c98",
        "5. Ranked outputs",
        &["Option 1", "Option 2", "GA compare"],
    );

    // Bottom-row boxes
    process_box(
        &mut svg,
        710.0,
        410.0,
        260.0,
        170.0,
        "#fff6fb",
        "#b03a78",
        "Option 2 combined ranking",
        &[
            "AsymmetricSplitVenturi",
            "Selective routing + venturi",
            "cavitation + healthy-cell shielding",
        ],
    );
    process_box(
        &mut svg,
        375.0,
        410.0,
        260.0,
        170.0,
        "#eefcf3",
        "#1e8449",
        "Option 1 acoustic ranking",
        &[
            "AsymmetricSplitResidenceSeparation",
            "Center-lane ultrasound/light",
            "RBC-biased peripheral bypass",
        ],
    );
    process_box(
        &mut svg,
        50.0,
        410.0,
        260.0,
        170.0,
        "#f4f6f7",
        "#566573",
        "Report and artifact export",
        &[
            "Schematics from cfd-schematics",
            "Canonical markdown, narrative,",
            "figure manifest and JSON",
        ],
    );

    let _ = write!(
        svg,
        r##"<text x="1520" y="435" font-size="18" fill="#34495e" font-weight="600">Option 1</text>
<text x="1520" y="460" font-size="14" fill="#34495e">Selective acoustic shortlist</text>
<text x="1520" y="495" font-size="18" fill="#34495e" font-weight="600">Option 2</text>
<text x="1520" y="520" font-size="14" fill="#34495e">Combined selective venturi shortlist</text>
<text x="1520" y="555" font-size="18" fill="#34495e" font-weight="600">Figure 6</text>
<text x="1520" y="580" font-size="14" fill="#34495e">HydroSDT GA comparison branch</text>"##
    );
    let _ = write!(
        svg,
        r##"<text x="70" y="660" font-size="14" fill="#566573">Selective routing remains the modeled concept: branch-diameter biasing pushes RBCs peripheral, keeps WBC/CTC enrichment in treatment branches, and changes only the treatment mechanism between Option 1 and Option 2.</text>"##
    );
    svg_end(&mut svg);
    std::fs::write(path, svg)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::write_creation_optimization_process_figure;
    use std::time::{SystemTime, UNIX_EPOCH};

    #[test]
    fn creation_optimization_process_figure_is_generated_not_placeholder() {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("system clock before unix epoch")
            .as_nanos();
        let path = std::env::temp_dir().join(format!("m12_process_{nanos}.svg"));
        write_creation_optimization_process_figure(&path)
            .expect("process figure generation must succeed");
        let svg = std::fs::read_to_string(path).expect("must read generated process svg");
        assert!(svg.contains("Primitive candidate"));
        assert!(svg.contains("Cached cfd-1d"));
        assert!(svg.contains("Option 2 combined ranking"));
        assert!(!svg.contains("placeholder"));
    }
}
