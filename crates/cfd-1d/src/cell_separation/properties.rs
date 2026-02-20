//! Cell physical properties for inertial microfluidic separation.
//!
//! # Physical basis
//!
//! Inertial focusing in microchannels depends on the ratio of cell diameter `a`
//! to channel hydraulic diameter `D_h`.  The confinement ratio `κ = a / D_h`
//! governs the strength of inertial lift forces.  Cells with `κ > 0.07` exhibit
//! measurable focusing; cells with `κ > 0.2` focus strongly to equilibrium
//! positions (Di Carlo, 2009, *Lab Chip* 9, 3038–3046).
//!
//! # Published cell data
//!
//! | Cell type | Diameter (µm) | Deformability index | Source |
//! |-----------|--------------|---------------------|--------|
//! | MCF-7 (breast cancer) | 17–20 | 0.15 (stiff) | Suresh (2007) |
//! | MDA-MB-231 (metastatic) | 14–17 | 0.35 (intermediate) | Suresh (2007) |
//! | Red blood cell (RBC) | 6–8 | 0.85 (very deformable) | Hochmuth (2000) |
//! | White blood cell (WBC) | 8–12 | 0.45 (moderate) | Lim (2006) |
//! | Platelet | 2–3 | 0.70 | Hochmuth (2000) |
//!
//! Deformability index (DI) is defined as the ratio of the cell's ability to
//! deform under shear relative to a rigid sphere of the same size.
//! DI = 0 → perfectly rigid; DI = 1 → perfectly deformable membrane.
//!
//! # References
//! - Di Carlo, D. (2009). Inertial microfluidics. *Lab Chip*, 9, 3038–3046.
//! - Suresh, S. (2007). Biomechanics and biophysics of cancer cells. *Acta Biomater.*, 3, 413–438.
//! - Hochmuth, R. M. (2000). Micropipette aspiration of living cells. *J. Biomech.*, 33, 15–22.

use serde::{Deserialize, Serialize};

/// Physical properties of a cell population for inertial focusing analysis.
///
/// # Invariants
/// - `diameter_m > 0`
/// - `deformability_index ∈ [0, 1]`
/// - `density_kg_m3 > 0`
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CellProperties {
    /// Human-readable cell type name (e.g. "MCF-7 breast cancer").
    pub name: &'static str,

    /// Mean cell diameter [m].
    ///
    /// For non-spherical cells (e.g. RBCs), use the effective hydrodynamic
    /// diameter measured in flow (typically the biconcave disc diameter).
    pub diameter_m: f64,

    /// Cell deformability index (DI) ∈ [0, 1].
    ///
    /// DI = 0 → rigid sphere (no membrane deformation under shear).
    /// DI = 1 → perfectly deformable (membrane offers no resistance).
    ///
    /// Rigid cells (low DI) experience stronger inertial lift and focus more
    /// sharply to equilibrium positions.  Deformable cells (high DI) migrate
    /// toward the channel wall due to wall-induced lift reduction.
    pub deformability_index: f64,

    /// Cell density [kg/m³].
    ///
    /// Cancer cells are typically slightly denser than RBCs due to higher
    /// nuclear-to-cytoplasm ratio.
    pub density_kg_m3: f64,
}

impl CellProperties {
    /// MCF-7 breast cancer cell (epithelial, stiff, large).
    ///
    /// Diameter: 17.5 µm (mean of 17–20 µm range, Suresh 2007).
    /// DI: 0.15 (stiff relative to RBCs, Suresh 2007).
    /// Density: 1050 kg/m³ (slightly denser than plasma).
    #[must_use]
    pub fn mcf7_breast_cancer() -> Self {
        Self {
            name: "MCF-7 breast cancer",
            diameter_m: 17.5e-6,
            deformability_index: 0.15,
            density_kg_m3: 1050.0,
        }
    }

    /// MDA-MB-231 metastatic breast cancer cell (mesenchymal, intermediate stiffness).
    ///
    /// Diameter: 15.5 µm (mean of 14–17 µm range, Suresh 2007).
    /// DI: 0.35 (more deformable than MCF-7 due to cytoskeletal remodelling).
    /// Density: 1045 kg/m³.
    #[must_use]
    pub fn mda_mb231_metastatic() -> Self {
        Self {
            name: "MDA-MB-231 metastatic",
            diameter_m: 15.5e-6,
            deformability_index: 0.35,
            density_kg_m3: 1045.0,
        }
    }

    /// Red blood cell (RBC) — highly deformable, small.
    ///
    /// Diameter: 7.0 µm (effective hydrodynamic diameter, Hochmuth 2000).
    /// DI: 0.85 (very deformable biconcave disc).
    /// Density: 1090 kg/m³ (haemoglobin-rich cytoplasm).
    #[must_use]
    pub fn red_blood_cell() -> Self {
        Self {
            name: "Red blood cell (RBC)",
            diameter_m: 7.0e-6,
            deformability_index: 0.85,
            density_kg_m3: 1090.0,
        }
    }

    /// White blood cell (WBC / leukocyte) — moderately deformable, intermediate size.
    ///
    /// Diameter: 10.0 µm (mean of 8–12 µm range, Lim 2006).
    /// DI: 0.45 (moderate deformability, nucleus limits deformation).
    /// Density: 1060 kg/m³.
    #[must_use]
    pub fn white_blood_cell() -> Self {
        Self {
            name: "White blood cell (WBC)",
            diameter_m: 10.0e-6,
            deformability_index: 0.45,
            density_kg_m3: 1060.0,
        }
    }

    /// Platelet (thrombocyte) — small, moderately deformable.
    ///
    /// Diameter: 2.5 µm (mean of 2–3 µm range, Hochmuth 2000).
    /// DI: 0.70 (disc-shaped, deformable under shear).
    /// Density: 1040 kg/m³.
    #[must_use]
    pub fn platelet() -> Self {
        Self {
            name: "Platelet",
            diameter_m: 2.5e-6,
            deformability_index: 0.70,
            density_kg_m3: 1040.0,
        }
    }

    /// Confinement ratio κ = a / D_h.
    ///
    /// Inertial focusing is significant when κ > 0.07 (Di Carlo 2009).
    /// Strong focusing occurs when κ > 0.2.
    #[inline]
    #[must_use]
    pub fn confinement_ratio(&self, hydraulic_diameter_m: f64) -> f64 {
        self.diameter_m / hydraulic_diameter_m
    }

    /// Returns `true` if this cell type will exhibit measurable inertial focusing
    /// in a channel with the given hydraulic diameter.
    ///
    /// Threshold: κ > 0.07 (Di Carlo 2009).
    #[inline]
    #[must_use]
    pub fn will_focus(&self, hydraulic_diameter_m: f64) -> bool {
        self.confinement_ratio(hydraulic_diameter_m) > 0.07
    }

    /// Returns `true` if this cell type will exhibit strong inertial focusing
    /// (κ > 0.2, Di Carlo 2009).
    #[inline]
    #[must_use]
    pub fn will_focus_strongly(&self, hydraulic_diameter_m: f64) -> bool {
        self.confinement_ratio(hydraulic_diameter_m) > 0.2
    }
}
