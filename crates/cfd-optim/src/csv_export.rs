//! CSV export for per-channel hemolysis decomposition.
//!
//! [`save_per_channel_csv`] writes a comma-separated file with one row per
//! channel per design, suitable for import into spreadsheets and report
//! appendices.

use std::io::Write;
use std::path::Path;

use crate::optimizer::RankedDesign;

/// Write per-channel hemolysis data for a set of ranked designs to a CSV file.
///
/// Columns:
///
/// | Column            | Description                                      |
/// |-------------------|--------------------------------------------------|
/// | `design_rank`     | 1-based rank within the optimisation batch        |
/// | `design_id`       | Candidate identifier string                      |
/// | `channel_id`      | Channel segment name (e.g. `"center_lv0"`)       |
/// | `is_venturi`      | `true` if this segment is a venturi throat        |
/// | `hi_contribution` | Flow-weighted Giersiepen HI from this channel     |
/// | `wall_shear_pa`   | Wall shear stress \[Pa\]                          |
/// | `transit_time_s`  | Transit time through the segment \[s\]            |
/// | `flow_fraction`   | Fraction of inlet flow carried by this channel    |
///
/// # Errors
/// Returns an error if the file cannot be created or written.
pub fn save_per_channel_csv(
    designs: &[RankedDesign],
    path: &Path,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut file = std::fs::File::create(path)?;
    writeln!(
        file,
        "design_rank,design_id,channel_id,is_venturi,hi_contribution,wall_shear_pa,transit_time_s,flow_fraction"
    )?;

    for d in designs {
        for ch in &d.metrics.per_channel_hemolysis {
            writeln!(
                file,
                "{},{},{},{},{:.6e},{:.4},{:.6e},{:.6}",
                d.rank,
                d.candidate.id,
                ch.channel_id,
                ch.is_venturi_throat,
                ch.hi_contribution,
                ch.wall_shear_pa,
                ch.transit_time_s,
                ch.flow_fraction,
            )?;
        }
    }
    file.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::design::DesignCandidate;
    use crate::design::DesignTopology;
    use crate::metrics::{ChannelHemolysis, SdtMetrics};

    fn test_candidate(id: &str) -> DesignCandidate {
        DesignCandidate {
            id: id.to_owned(),
            topology: DesignTopology::SingleVenturi,
            flow_rate_m3_s: 1e-7,
            inlet_gauge_pa: 5000.0,
            throat_diameter_m: 200e-6,
            inlet_diameter_m: 1e-3,
            throat_length_m: 400e-6,
            channel_width_m: 1e-3,
            channel_height_m: 1e-3,
            serpentine_segments: 6,
            segment_length_m: 10e-3,
            bend_radius_m: 2e-3,
            feed_hematocrit: 0.45,
            trifurcation_center_frac: 1.0 / 3.0,
            cif_pretri_center_frac: 1.0 / 3.0,
            cif_terminal_tri_center_frac: 1.0 / 3.0,
            cif_terminal_bi_treat_frac: 0.68,
            asymmetric_narrow_frac: 0.5,
            trifurcation_left_frac: 1.0 / 3.0,
            cross_section_shape: Default::default(),
        }
    }

    #[test]
    fn csv_round_trip_header_and_rows() {
        let designs = vec![RankedDesign {
            rank: 1,
            candidate: test_candidate("test-001"),
            metrics: SdtMetrics {
                per_channel_hemolysis: vec![
                    ChannelHemolysis {
                        channel_id: "center_lv0".to_owned(),
                        is_venturi_throat: false,
                        hi_contribution: 1.2e-6,
                        wall_shear_pa: 42.5,
                        transit_time_s: 3.1e-3,
                        flow_fraction: 0.33,
                    },
                    ChannelHemolysis {
                        channel_id: "throat".to_owned(),
                        is_venturi_throat: true,
                        hi_contribution: 8.9e-5,
                        wall_shear_pa: 1200.0,
                        transit_time_s: 5.0e-6,
                        flow_fraction: 0.33,
                    },
                ],
                ..Default::default()
            },
            score: 0.85,
        }];

        let dir = std::env::temp_dir().join("cfd_optim_csv_test");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("test_hemolysis.csv");
        save_per_channel_csv(&designs, &path).unwrap();

        let content = std::fs::read_to_string(&path).unwrap();
        let lines: Vec<&str> = content.lines().collect();
        assert_eq!(lines.len(), 3, "header + 2 data rows");
        assert!(lines[0].starts_with("design_rank,"));
        assert!(lines[1].starts_with("1,test-001,center_lv0,false,"));
        assert!(lines[2].starts_with("1,test-001,throat,true,"));

        std::fs::remove_dir_all(&dir).ok();
    }
}
