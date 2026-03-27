//! Pulsatility index computation for Womersley flow.
//!
//! ## Theorem — Gosling-King Pulsatility Index (Gosling & King 1974)
//!
//! The pulsatility index (PI) is a dimensionless measure of flow waveform
//! pulsatility, widely used in Doppler ultrasound assessment of arterial
//! hemodynamics. It quantifies the oscillatory component of flow relative
//! to the mean:
//!
//! ```text
//! PI = (V_peak_systolic - V_end_diastolic) / V_mean
//! ```
//!
//! For Womersley flow decomposed into a steady (Poiseuille) component plus
//! a single oscillatory harmonic:
//!
//! ```text
//! V(t) = V_mean + V_osc · cos(omega·t + phi)
//! ```
//!
//! The peak and trough velocities are:
//!
//! ```text
//! V_peak   = V_mean + |V_osc_max|
//! V_trough = V_mean - |V_osc_max|
//! ```
//!
//! For a single harmonic, `|V_osc_max| = Q_amplitude / A`, giving:
//!
//! ```text
//! PI = 2 · Q_amplitude / Q_mean
//! ```
//!
//! ## Clinical significance
//!
//! | Vessel           | Typical PI  |
//! |------------------|-------------|
//! | Internal carotid | 0.8 - 1.2   |
//! | Middle cerebral  | 0.6 - 1.0   |
//! | Femoral artery   | 4.0 - 12.0  |
//! | Uterine artery   | 0.6 - 1.2   |
//!
//! **Reference**: Gosling, R.G. & King, D.H. (1974). "Arterial Assessment
//! by Doppler-shift Ultrasound", *Proc. R. Soc. Med.* 67:447-449.

/// Compute the Gosling-King pulsatility index from mean and oscillatory flow rates.
///
/// # Arguments
/// * `q_mean` - Mean volumetric flow rate [m³/s] (must be positive)
/// * `q_amplitude` - Oscillatory flow rate amplitude [m³/s] (non-negative)
/// * `cross_section_area` - Cross-sectional area of the vessel [m²] (must be positive)
///
/// # Returns
/// Tuple of `(pulsatility_index, v_peak, v_trough)`:
/// - `pulsatility_index` — PI = (v_peak - v_trough) / v_mean
/// - `v_peak` — peak systolic velocity [m/s]
/// - `v_trough` — end-diastolic (minimum) velocity [m/s]
///
/// Returns `(0.0, v_mean, v_mean)` when `q_mean` is zero or negative.
#[must_use]
pub fn womersley_pulsatility_index(
    q_mean: f64,
    q_amplitude: f64,
    cross_section_area: f64,
) -> (f64, f64, f64) {
    let area = cross_section_area.max(1e-30); // guard against zero area
    let v_mean = q_mean / area;
    let v_osc = q_amplitude.abs() / area;

    let v_peak = v_mean + v_osc;
    let v_trough = v_mean - v_osc;

    let pi = if q_mean.abs() < 1e-30 {
        // Steady flow or zero flow: PI is undefined / zero
        0.0
    } else {
        (v_peak - v_trough) / v_mean.abs()
    };

    (pi, v_peak, v_trough)
}
