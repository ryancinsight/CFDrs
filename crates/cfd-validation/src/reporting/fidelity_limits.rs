//! Fidelity-level usage limits and validated accuracy summary.
//!
//! [`FidelityLimit`] captures resolution requirements, validated accuracy,
//! and known limitations for a single solver fidelity level.
//! [`fidelity_limits_table`] returns the complete table suitable for M12
//! report appendices.

use std::fmt;

/// Summary of validated accuracy and usage constraints for one solver fidelity.
#[derive(Debug, Clone)]
pub struct FidelityLimit {
    /// Human-readable fidelity label (e.g. `"1D"`, `"2D"`, `"3D"`).
    pub fidelity: &'static str,
    /// Solver family name (e.g. `"Hagen-Poiseuille network"`).
    pub solver: &'static str,
    /// Analytical or literature benchmarks used for validation.
    pub validated_against: &'static str,
    /// Convergence order on smooth problems.
    pub accuracy_order: &'static str,
    /// Minimum recommended mesh / node count.
    pub min_resolution: &'static str,
    /// Physical phenomena captured at this fidelity.
    pub phenomena: &'static str,
    /// Known limitations or accuracy caveats.
    pub limitations: &'static str,
}

/// Return the complete fidelity-limit table for all solver levels.
///
/// Entries follow the hierarchy: 1D → 2D → 3D.
#[must_use]
pub fn fidelity_limits_table() -> Vec<FidelityLimit> {
    vec![
        FidelityLimit {
            fidelity: "1D",
            solver: "Hagen-Poiseuille lumped network",
            validated_against: "Poiseuille ΔP (exact), Murray's law, Womersley pulsatile",
            accuracy_order: "Exact (analytic)",
            min_resolution: "N/A (lumped)",
            phenomena: "Pressure drop, flow split, Giersiepen hemolysis, Zweifach-Fung routing",
            limitations: "No inertial effects; assumes fully-developed flow; rectangular duct \
                          approximation (wide-duct limit w >> h)",
        },
        FidelityLimit {
            fidelity: "2D",
            solver: "SIMPLE / PISO FVM; LBM BGK D2Q9",
            validated_against: "Ghia cavity (Re 100–1000), Couette, Poiseuille, Taylor-Green",
            accuracy_order: "2nd (FVM); 2nd (LBM, O(Ma²) compressibility error)",
            min_resolution: "40×40 (cavity); 20×N (channel per hydraulic diameter)",
            phenomena: "Viscous flow, pressure-velocity coupling, recirculation, shear fields, \
                        venturi acceleration, Dean vortices (mapped 2D only)",
            limitations: "No 3D secondary flow; assumes 2D plane / depth-averaged; \
                          turbulence models limited to k-ε / k-ω at moderate Re",
        },
        FidelityLimit {
            fidelity: "3D",
            solver: "Stokes FEM (P1-P1 SUPG/PSPG); Picard non-Newtonian iteration",
            validated_against: "Poiseuille ΔP rectangular duct; 1D flow-split conservation",
            accuracy_order: "2nd (linear elements); sub-optimal pressure at coarse resolution",
            min_resolution: "40×8×8 per channel (CIF cascade); 20×5×5 (qualitative only)",
            phenomena: "Full 3D velocity/pressure fields, non-Newtonian (Casson) blood, \
                        wall shear extraction, venturi constriction, hematocrit-viscosity coupling",
            limitations: "Stokes assumption (Re < 1 required); coarse mesh ΔP magnitude \
                          unreliable (factor 10–100× vs analytical); Picard convergence \
                          sensitive to initial viscosity at extreme HCT ratios; \
                          no turbulence, no multiphase",
        },
    ]
}

impl fmt::Display for FidelityLimit {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}: {} | Order: {} | Min res: {} | Validated: {}",
            self.fidelity,
            self.solver,
            self.accuracy_order,
            self.min_resolution,
            self.validated_against,
        )
    }
}

/// Render the fidelity-limit table as a Markdown string.
#[must_use]
pub fn fidelity_limits_markdown() -> String {
    let limits = fidelity_limits_table();
    let mut md = String::with_capacity(2048);
    md.push_str("| Fidelity | Solver | Accuracy | Min Resolution | Validated Against | Phenomena | Limitations |\n");
    md.push_str("|----------|--------|----------|----------------|-------------------|-----------|-------------|\n");
    for l in &limits {
        md.push_str(&format!(
            "| {} | {} | {} | {} | {} | {} | {} |\n",
            l.fidelity,
            l.solver,
            l.accuracy_order,
            l.min_resolution,
            l.validated_against,
            l.phenomena,
            l.limitations,
        ));
    }
    md
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn table_has_three_fidelity_levels() {
        let t = fidelity_limits_table();
        assert_eq!(t.len(), 3);
        assert_eq!(t[0].fidelity, "1D");
        assert_eq!(t[1].fidelity, "2D");
        assert_eq!(t[2].fidelity, "3D");
    }

    #[test]
    fn markdown_contains_header_row() {
        let md = fidelity_limits_markdown();
        assert!(md.starts_with("| Fidelity |"));
        assert!(md.contains("| 1D |"));
        assert!(md.contains("| 2D |"));
        assert!(md.contains("| 3D |"));
    }
}
