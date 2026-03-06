//! Extensible metadata system for channels and nodes
//!
//! This module provides a flexible metadata system that allows for easy addition
//! of new tracking variables without requiring changes to core data structures.
//! It uses trait-based extensibility with type-safe metadata storage.

use std::any::{Any, TypeId};
use std::collections::HashMap;
use serde::{Deserialize, Serialize};
use std::fmt::Debug;

/// Base trait for all metadata types
///
/// This trait provides the foundation for type-safe metadata storage.
/// All metadata types must implement this trait to be stored in the system.
pub trait Metadata: Any + Debug + Send + Sync {
    /// Returns a unique name for this metadata type
    fn metadata_type_name(&self) -> &'static str;

    /// Clone the metadata as a boxed trait object
    fn clone_metadata(&self) -> Box<dyn Metadata>;

    /// Convert to Any for downcasting
    fn as_any(&self) -> &dyn Any;

    /// Convert to mutable Any for downcasting
    fn as_any_mut(&mut self) -> &mut dyn Any;
}

/// Metadata storage container
///
/// This container provides type-safe storage and retrieval of metadata
/// using `TypeId` as keys for efficient lookup.
#[derive(Debug)]
pub struct MetadataContainer {
    data: HashMap<TypeId, Box<dyn Metadata>>,
}

impl MetadataContainer {
    /// Create a new empty metadata container
    #[must_use]
    pub fn new() -> Self {
        Self {
            data: HashMap::new(),
        }
    }

    /// Insert metadata of a specific type
    pub fn insert<T: Metadata + Clone + 'static>(&mut self, metadata: T) {
        self.data.insert(TypeId::of::<T>(), Box::new(metadata));
    }

    /// Get metadata of a specific type
    #[must_use]
    pub fn get<T: Metadata + 'static>(&self) -> Option<&T> {
        self.data
            .get(&TypeId::of::<T>())
            .and_then(|boxed| boxed.as_any().downcast_ref::<T>())
    }

    /// Get mutable metadata of a specific type
    pub fn get_mut<T: Metadata + 'static>(&mut self) -> Option<&mut T> {
        self.data
            .get_mut(&TypeId::of::<T>())
            .and_then(|boxed| boxed.as_any_mut().downcast_mut::<T>())
    }

    /// Remove metadata of a specific type
    pub fn remove<T: Metadata + 'static>(&mut self) -> Option<Box<dyn Metadata>> {
        self.data.remove(&TypeId::of::<T>())
    }

    /// Check if metadata of a specific type exists
    #[must_use]
    pub fn contains<T: Metadata + 'static>(&self) -> bool {
        self.data.contains_key(&TypeId::of::<T>())
    }

    /// Get all metadata type names (for debugging)
    #[must_use]
    pub fn metadata_types(&self) -> Vec<&'static str> {
        self.data
            .values()
            .map(|metadata| metadata.metadata_type_name())
            .collect()
    }

    /// Check if container is empty
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Get number of metadata entries
    #[must_use]
    pub fn len(&self) -> usize {
        self.data.len()
    }
}

impl Clone for MetadataContainer {
    fn clone(&self) -> Self {
        let mut new_container = Self::new();
        for (type_id, metadata) in &self.data {
            new_container
                .data
                .insert(*type_id, metadata.clone_metadata());
        }
        new_container
    }
}

impl Default for MetadataContainer {
    fn default() -> Self {
        Self::new()
    }
}

/// Flow-related metadata for channels
#[derive(Debug, Clone, PartialEq)]
pub struct FlowMetadata {
    /// Flow rate in μL/min
    pub flow_rate: f64,
    /// Pressure drop in Pa
    pub pressure_drop: f64,
    /// Reynolds number
    pub reynolds_number: f64,
    /// Velocity in m/s
    pub velocity: f64,
}

impl Metadata for FlowMetadata {
    fn metadata_type_name(&self) -> &'static str {
        "FlowMetadata"
    }

    fn clone_metadata(&self) -> Box<dyn Metadata> {
        Box::new(self.clone())
    }

    fn as_any(&self) -> &dyn Any {
        self
    }

    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }
}

/// Thermal metadata for channels
#[derive(Debug, Clone, PartialEq)]
pub struct ThermalMetadata {
    /// Temperature in Celsius
    pub temperature: f64,
    /// Heat transfer coefficient in W/(m²·K)
    pub heat_transfer_coefficient: f64,
    /// Thermal conductivity in W/(m·K)
    pub thermal_conductivity: f64,
}

impl Metadata for ThermalMetadata {
    fn metadata_type_name(&self) -> &'static str {
        "ThermalMetadata"
    }

    fn clone_metadata(&self) -> Box<dyn Metadata> {
        Box::new(self.clone())
    }

    fn as_any(&self) -> &dyn Any {
        self
    }

    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }
}

/// Manufacturing tolerance metadata
#[derive(Debug, Clone, PartialEq)]
pub struct ManufacturingMetadata {
    /// Width tolerance in micrometers
    pub width_tolerance: f64,
    /// Height tolerance in micrometers
    pub height_tolerance: f64,
    /// Surface roughness in micrometers
    pub surface_roughness: f64,
    /// Manufacturing method
    pub manufacturing_method: String,
}

impl Metadata for ManufacturingMetadata {
    fn metadata_type_name(&self) -> &'static str {
        "ManufacturingMetadata"
    }

    fn clone_metadata(&self) -> Box<dyn Metadata> {
        Box::new(self.clone())
    }

    fn as_any(&self) -> &dyn Any {
        self
    }

    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }
}

/// Channel geometry metadata for downstream 3D-aware workflows.
#[derive(Debug, Clone, PartialEq)]
pub struct ChannelGeometryMetadata {
    /// Design channel diameter in millimeters used for spacing and reconstruction.
    pub channel_diameter_mm: f64,
}

impl Metadata for ChannelGeometryMetadata {
    fn metadata_type_name(&self) -> &'static str {
        "ChannelGeometryMetadata"
    }

    fn clone_metadata(&self) -> Box<dyn Metadata> {
        Box::new(self.clone())
    }

    fn as_any(&self) -> &dyn Any {
        self
    }

    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }
}

/// Optimization history metadata
#[derive(Debug, Clone, PartialEq)]
pub struct OptimizationMetadata {
    /// Original channel length before optimization
    pub original_length: f64,
    /// Optimized channel length
    pub optimized_length: f64,
    /// Length improvement percentage
    pub improvement_percentage: f64,
    /// Optimization iterations used
    pub iterations: usize,
    /// Optimization time in milliseconds
    pub optimization_time_ms: u64,
    /// Optimization profile used
    pub optimization_profile: String,
}

impl Metadata for OptimizationMetadata {
    fn metadata_type_name(&self) -> &'static str {
        "OptimizationMetadata"
    }

    fn clone_metadata(&self) -> Box<dyn Metadata> {
        Box::new(self.clone())
    }

    fn as_any(&self) -> &dyn Any {
        self
    }

    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }
}

/// Runtime performance metadata
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PerformanceMetadata {
    /// Generation time in microseconds
    pub generation_time_us: u64,
    /// Memory usage in bytes
    pub memory_usage_bytes: usize,
    /// Number of path points generated
    pub path_points_count: usize,
}

impl Metadata for PerformanceMetadata {
    fn metadata_type_name(&self) -> &'static str {
        "PerformanceMetadata"
    }

    fn clone_metadata(&self) -> Box<dyn Metadata> {
        Box::new(self.clone())
    }

    fn as_any(&self) -> &dyn Any {
        self
    }

    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }
}

/// Absolute node position in the schematic layout plane [mm].
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct NodeLayoutMetadata {
    pub x_mm: f64,
    pub y_mm: f64,
}

impl Metadata for NodeLayoutMetadata {
    fn metadata_type_name(&self) -> &'static str {
        "NodeLayoutMetadata"
    }

    fn clone_metadata(&self) -> Box<dyn Metadata> {
        Box::new(*self)
    }

    fn as_any(&self) -> &dyn Any {
        self
    }

    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }
}

/// Rendering role for a channel path in the schematic.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ChannelVisualRole {
    Trunk,
    CenterTreatment,
    PeripheralBypass,
    MergeCollector,
    VenturiThroat,
    Diffuser,
    InternalLink,
}

/// Explicit polyline path for blueprint-native schematic rendering [mm].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ChannelPathMetadata {
    pub polyline_mm: Vec<(f64, f64)>,
    pub visual_role: ChannelVisualRole,
}

impl Metadata for ChannelPathMetadata {
    fn metadata_type_name(&self) -> &'static str {
        "ChannelPathMetadata"
    }

    fn clone_metadata(&self) -> Box<dyn Metadata> {
        Box::new(self.clone())
    }

    fn as_any(&self) -> &dyn Any {
        self
    }

    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }
}

/// Junction family metadata used by rendering and 1D minor-loss models.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum JunctionFamily {
    Bifurcation,
    Trifurcation,
    Tee,
    Cross,
    Merge,
}

/// Geometry metadata for split/merge junctions.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct JunctionGeometryMetadata {
    pub junction_family: JunctionFamily,
    pub branch_angles_deg: Vec<f64>,
    pub merge_angles_deg: Vec<f64>,
}

impl Metadata for JunctionGeometryMetadata {
    fn metadata_type_name(&self) -> &'static str {
        "JunctionGeometryMetadata"
    }

    fn clone_metadata(&self) -> Box<dyn Metadata> {
        Box::new(self.clone())
    }

    fn as_any(&self) -> &dyn Any {
        self
    }

    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }
}

// ── Therapy geometry metadata types ─────────────────────────────────────────

/// Geometry parameters for a venturi constriction channel.
///
/// Attached to the `throat_section` channel in venturi preset factories so
/// that downstream consumers (cfd-optim, cfd-mesh) can query exact throat
/// dimensions without pattern-matching channel IDs.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct VenturiGeometryMetadata {
    /// Throat channel width [m] — the constriction width.
    pub throat_width_m: f64,
    /// Throat channel height [m] — same as inlet height for planar chips.
    pub throat_height_m: f64,
    /// Throat channel length [m].
    pub throat_length_m: f64,
    /// Inlet/outlet channel width [m] upstream and downstream of the throat.
    pub inlet_width_m: f64,
    /// Outlet channel width [m] downstream of the throat.
    pub outlet_width_m: f64,
    /// Convergent half-angle [deg].
    pub convergent_half_angle_deg: f64,
    /// Divergent half-angle [deg].
    pub divergent_half_angle_deg: f64,
}

impl Metadata for VenturiGeometryMetadata {
    fn metadata_type_name(&self) -> &'static str {
        "VenturiGeometryMetadata"
    }
    fn clone_metadata(&self) -> Box<dyn Metadata> {
        Box::new(self.clone())
    }
    fn as_any(&self) -> &dyn Any {
        self
    }
    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }
}

/// Parameters for a cascade-center trifurcation separator.
///
/// Attached to the inlet junction of CCT blueprints so consumers can
/// reconstruct the Zweifach-Fung routing fractions without re-parsing
/// channel names.
#[derive(Debug, Clone, PartialEq)]
pub struct CascadeParams {
    /// Number of trifurcation cascade levels (typically 1–3).
    pub n_levels: u8,
    /// Center-arm width fraction ∈ [0.25, 0.65].
    pub center_frac: f64,
}

impl Metadata for CascadeParams {
    fn metadata_type_name(&self) -> &'static str {
        "CascadeParams"
    }
    fn clone_metadata(&self) -> Box<dyn Metadata> {
        Box::new(self.clone())
    }
    fn as_any(&self) -> &dyn Any {
        self
    }
    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }
}

/// Parameters for an incremental filtration tri-bi separator.
///
/// Attached to the inlet junction of CIF blueprints.
#[derive(Debug, Clone, PartialEq)]
pub struct IncrementalFiltrationParams {
    /// Number of pre-trifurcation stages (typically 1–3).
    pub n_pretri: u8,
    /// Legacy center-arm width fraction ∈ [0.25, 0.65].
    ///
    /// Preserved for backward compatibility with older CIF metadata readers.
    pub center_frac: f64,
    /// Pre-trifurcation center-arm width fraction ∈ [0.25, 0.65].
    pub pretri_center_frac: f64,
    /// Terminal-trifurcation center-arm width fraction ∈ [0.25, 0.65].
    pub terminal_tri_center_frac: f64,
    /// Terminal-bifurcation treatment-arm fraction ∈ [0.50, 0.85].
    pub bi_treat_frac: f64,
    /// Outlet-tail channel length from `outlet_merge` to `outlet` [m].
    ///
    /// Shorter tails represent "remerge near outlet" layouts where treated and
    /// bypass streams converge immediately before exiting the device.
    pub outlet_tail_length_m: f64,
}

impl Metadata for IncrementalFiltrationParams {
    fn metadata_type_name(&self) -> &'static str {
        "IncrementalFiltrationParams"
    }
    fn clone_metadata(&self) -> Box<dyn Metadata> {
        Box::new(self.clone())
    }
    fn as_any(&self) -> &dyn Any {
        self
    }
    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }
}

/// Parameters for an asymmetric 3-stream trifurcation venturi blueprint.
///
/// Attached to the inlet junction of asymmetric trifurcation blueprints.
/// The three arms have independent width fractions; only the center arm
/// receives a venturi throat for selective SDT treatment.
#[derive(Debug, Clone, PartialEq)]
pub struct AsymmetricTrifurcationParams {
    /// Center arm width fraction ∈ [0.20, 0.60] — receives venturi treatment.
    pub center_frac: f64,
    /// Left arm width fraction ∈ [0.15, 0.50] — WBC collection port.
    pub left_frac: f64,
    /// Right arm width fraction = 1 - center_frac - left_frac — RBC waste port.
    pub right_frac: f64,
}

impl Metadata for AsymmetricTrifurcationParams {
    fn metadata_type_name(&self) -> &'static str {
        "AsymmetricTrifurcationParams"
    }
    fn clone_metadata(&self) -> Box<dyn Metadata> {
        Box::new(self.clone())
    }
    fn as_any(&self) -> &dyn Any {
        self
    }
    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }
}

// ── Per-channel differential venturi throat specification ─────────────────────

/// Per-channel venturi throat specification for differential SDT treatment.
///
/// Encodes the number of serial venturi throats placed on a single channel
/// segment and whether the channel is in the CTC-enriched treatment stream.
/// Attached to individual channel segments in selective-separation blueprints
/// so that `cfd-optim` metrics and `cfd-1d` resistance solvers can faithfully
/// model cumulative cavitation exposure without confounding bypass-arm RBC
/// hemolysis.
///
/// # Physical Justification
///
/// The SDT trifurcation figure shows center CTC-enriched channels traversing
/// the sonication zone with multiple venturi constrictions (N_throats ≥ 2),
/// while outer bypass channels carrying the RBC-depleted fraction bypass the
/// zone entirely (N_throats = 0).  This asymmetry is critical for:
///
/// - **Cumulative cavitation dose**: each throat adds one Rayleigh-collapse
///   event per pass: `dose_n = 1 − (1 − cav_potential)^N_throats`.
/// - **Hemolysis partitioning**: bypass channels accumulate wall-shear HI
///   only, zero cavitation HI; center channels accumulate both.
/// - **FDA MI compliance**: the mechanical index at each throat must satisfy
///   `MI = P_neg / √f_us < 1.9` (FDA 510(k) guidance, 2019).
///   For hydrodynamic cavitation the equivalent is `σ < σ_crit` with
///   controlled pressure so `ΔP_throat ≤ P_FDA_max`.
///
/// # Invariants
///
/// - `n_throats ≥ 0`.  `0` means the channel is a bypass (no cavitation).
/// - `n_throats = 1` is a standard single-throat venturi.
/// - `n_throats ≥ 2` applies to the center CTC-enriched channel only.
/// - `is_ctc_stream = true` iff this channel carries the cancer-enriched
///   center fraction from the upstream trifurcation/CIF separator.
#[derive(Debug, Clone, PartialEq)]
pub struct ChannelVenturiSpec {
    /// Number of serial venturi throats on this channel segment.
    ///
    /// - `0` → bypass channel; no cavitation dose, no venturi HI.
    /// - `1` → standard single-throat venturi.
    /// - `≥ 2` → multi-throat CTC treatment channel (center stream only).
    pub n_throats: u8,

    /// True when this channel segment carries the CTC-enriched center stream.
    ///
    /// Used by `compute_metrics` to gate selective-HCT correction and
    /// per-channel bypass hemolysis accounting.
    pub is_ctc_stream: bool,

    /// Throat width [m] used at each serial throat on this channel.
    ///
    /// All serial throats on one channel share the same width to preserve
    /// the downstream pressure budget (Idelchik 1994, §4-9: adding identical
    /// constrictions in series multiplies the pressure drop linearly).
    pub throat_width_m: f64,

    /// Channel height [m] — same as device height for planar PDMS/SU-8 chips.
    pub height_m: f64,

    /// Physical inter-throat spacing [m]: distance between successive
    /// throat exit planes, must satisfy `l_spacing > 10 × D_h` so that
    /// the flow re-develops between throats (Shah & London 1978, §2-3).
    pub inter_throat_spacing_m: f64,
}

impl ChannelVenturiSpec {
    /// Cumulative cavitation dose factor for `N` serial throats.
    ///
    /// # Theorem — Serial Cavitation Accumulation
    ///
    /// Each throat generates an independent Rayleigh-collapse event with
    /// probability `cav_potential`.  Assuming statistically independent
    /// collapse events (valid when `l_spacing > 10 D_h`, so bubbles from
    /// one throat fully collapse before the next):
    ///
    /// ```text
    /// dose_N = 1 − (1 − p)^N
    /// ```
    ///
    /// where `p = cav_potential ∈ [0, 1]` and `N = n_throats`.
    ///
    /// **Proof**: By independence, the probability that a given bubble
    /// *escapes* all N throats without collapsing is `(1 − p)^N`.
    /// Complementing gives the dose fraction.  At `N = 1` this reduces to
    /// the single-throat result `p`.  As `N → ∞` the dose → 1 regardless
    /// of `p > 0`.
    #[must_use]
    pub fn cumulative_cavitation_dose(&self, cav_potential: f64) -> f64 {
        let p = cav_potential.clamp(0.0, 1.0);
        if p <= 0.0 || self.n_throats == 0 {
            return 0.0;
        }
        (1.0 - (1.0 - p).powi(i32::from(self.n_throats))).clamp(0.0, 1.0)
    }

    /// Total throat pressure drop for N serial throats [Pa].
    ///
    /// Each identical throat contributes the same Bernoulli drop
    /// `ΔP_throat = ½ρ(v_throat² − v_inlet²)`.  For N throats in series:
    ///
    /// ```text
    /// ΔP_total = N × ΔP_single × (1 − C_D)
    /// ```
    ///
    /// where `C_D` is the diffuser recovery coefficient (Idelchik 1994).
    #[must_use]
    pub fn total_throat_pressure_drop_pa(
        &self,
        flow_m3_s: f64,
        blood_density_kg_m3: f64,
        diffuser_coeff: f64,
    ) -> f64 {
        if self.n_throats == 0 {
            return 0.0;
        }
        let inlet_area = (self.throat_width_m * 2.0) * self.height_m; // 2× throat for inlet width
        let throat_area = self.throat_width_m * self.height_m;
        let v_in = flow_m3_s / inlet_area.max(1e-18);
        let v_throat = flow_m3_s / throat_area.max(1e-18);
        let single_dp = 0.5
            * blood_density_kg_m3
            * (v_throat * v_throat - v_in * v_in).max(0.0)
            * (1.0 - diffuser_coeff);
        single_dp * f64::from(self.n_throats)
    }
}

impl Metadata for ChannelVenturiSpec {
    fn metadata_type_name(&self) -> &'static str {
        "ChannelVenturiSpec"
    }
    fn clone_metadata(&self) -> Box<dyn Metadata> {
        Box::new(self.clone())
    }
    fn as_any(&self) -> &dyn Any {
        self
    }
    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }
}

/// FDA Mechanical Index compliance record for a venturi throat channel.
///
/// For hydrodynamic cavitation in millifluidic SDT devices the equivalent
/// mechanical index is derived from the throat pressure drop:
///
/// ```text
/// MI_equiv = sqrt(2 × ΔP_throat / ρ_blood) / v_sound_blood
/// ```
///
/// FDA 510(k) guidance (2019) mandates `MI < 1.9` for diagnostic ultrasound;
/// therapeutic SDT typically targets `0.3 < MI_equiv < 1.5` for controlled
/// inertial cavitation without uncontrolled bubble collapse chains.
#[derive(Debug, Clone, PartialEq)]
pub struct FdaCavitationCompliance {
    /// Equivalent mechanical index at this throat (dimensionless).
    pub mi_equiv: f64,

    /// Whether this throat is FDA-compliant: `mi_equiv < FDA_MI_LIMIT`.
    pub fda_mi_compliant: bool,

    /// Target MI for therapeutic efficacy: lower bound for detectable SDT.
    /// Defaults to 0.3 (threshold for inertial cavitation nucleation in blood).
    pub therapeutic_mi_lower: f64,

    /// Whether this throat is in the therapeutic window:
    /// `therapeutic_mi_lower ≤ mi_equiv < FDA_MI_LIMIT`.
    pub in_therapeutic_window: bool,
}

impl FdaCavitationCompliance {
    /// FDA upper MI limit for therapeutic devices (FDA 510(k), 2019).
    pub const FDA_MI_LIMIT: f64 = 1.9;

    /// Lower MI threshold for inertial cavitation nucleation in blood.
    /// Brennen (1995) _Cavitation and Bubble Dynamics_, §1.6.
    pub const INERTIAL_CAV_THRESHOLD_MI: f64 = 0.3;

    /// Speed of sound in blood [m/s] (Shung 2006).
    pub const SOUND_SPEED_BLOOD_M_S: f64 = 1540.0;

    /// Compute compliance from throat velocity and blood properties.
    ///
    /// # Arguments
    ///
    /// * `dp_throat_pa` — pressure drop at the vena contracta [Pa]
    /// * `blood_density` — blood density [kg/m³]
    #[must_use]
    pub fn from_throat_dp(dp_throat_pa: f64, blood_density: f64) -> Self {
        // Equivalent velocity fluctuation: v_rms ≈ sqrt(2 ΔP / ρ)
        let v_rms = (2.0 * dp_throat_pa.max(0.0) / blood_density.max(1.0)).sqrt();
        // Dimensionless MI analogue: v_rms / c_sound
        let mi_equiv = v_rms / Self::SOUND_SPEED_BLOOD_M_S;
        Self {
            mi_equiv,
            fda_mi_compliant: mi_equiv < Self::FDA_MI_LIMIT,
            therapeutic_mi_lower: Self::INERTIAL_CAV_THRESHOLD_MI,
            in_therapeutic_window: mi_equiv >= Self::INERTIAL_CAV_THRESHOLD_MI
                && mi_equiv < Self::FDA_MI_LIMIT,
        }
    }
}

impl Metadata for FdaCavitationCompliance {
    fn metadata_type_name(&self) -> &'static str {
        "FdaCavitationCompliance"
    }
    fn clone_metadata(&self) -> Box<dyn Metadata> {
        Box::new(self.clone())
    }
    fn as_any(&self) -> &dyn Any {
        self
    }
    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }
}

/// Blueprint-level rendering hints written by `cfd-optim` into a
/// [`NetworkBlueprint`] so that a generic renderer can produce fully-annotated
/// schematics without any `cfd-optim`-specific logic in the render path.
///
/// Attach to a blueprint via
/// [`NetworkBlueprint::with_render_hints`][crate::domain::model::NetworkBlueprint::with_render_hints]
/// and read back with
/// [`NetworkBlueprint::render_hints`][crate::domain::model::NetworkBlueprint::render_hints].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BlueprintRenderHints {
    /// Human-readable split-stage sequence, e.g. `"Bi→Tri"` or `"Tri→Bi→Bi"`.
    pub stage_sequence: String,
    /// Number of visible split layers in the topology tree.
    pub split_layers: usize,
    /// Fallback throat count when channel-level metadata yields zero
    /// (e.g. topologies where throats are implicit rather than tagged).
    pub throat_count_hint: usize,
    /// Treatment zone label: `"venturi"` for hydrodynamic SDT, `"ultrasound"`
    /// for acoustic-only designs.
    pub treatment_label: String,
}

impl Metadata for BlueprintRenderHints {
    fn metadata_type_name(&self) -> &'static str {
        "BlueprintRenderHints"
    }

    fn clone_metadata(&self) -> Box<dyn Metadata> {
        Box::new(self.clone())
    }

    fn as_any(&self) -> &dyn Any {
        self
    }

    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }
}

/// Convenience macro for implementing Metadata trait
#[macro_export]
macro_rules! impl_metadata {
    ($type:ty, $name:expr) => {
        impl Metadata for $type {
            fn metadata_type_name(&self) -> &'static str {
                $name
            }

            fn clone_metadata(&self) -> Box<dyn Metadata> {
                Box::new(self.clone())
            }

            fn as_any(&self) -> &dyn Any {
                self
            }

            fn as_any_mut(&mut self) -> &mut dyn Any {
                self
            }
        }
    };
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_metadata_container_basic_operations() {
        let mut container = MetadataContainer::new();

        // Test insertion and retrieval
        let flow_data = FlowMetadata {
            flow_rate: 10.0,
            pressure_drop: 1000.0,
            reynolds_number: 0.1,
            velocity: 0.001,
        };

        container.insert(flow_data.clone());

        let retrieved = container.get::<FlowMetadata>().unwrap();
        assert_eq!(retrieved, &flow_data);

        // Test contains
        assert!(container.contains::<FlowMetadata>());
        assert!(!container.contains::<ThermalMetadata>());

        // Test removal
        let removed = container.remove::<FlowMetadata>();
        assert!(removed.is_some());
        assert!(!container.contains::<FlowMetadata>());
    }

    #[test]
    fn test_multiple_metadata_types() {
        let mut container = MetadataContainer::new();

        let flow_data = FlowMetadata {
            flow_rate: 10.0,
            pressure_drop: 1000.0,
            reynolds_number: 0.1,
            velocity: 0.001,
        };

        let thermal_data = ThermalMetadata {
            temperature: 25.0,
            heat_transfer_coefficient: 100.0,
            thermal_conductivity: 0.6,
        };

        container.insert(flow_data.clone());
        container.insert(thermal_data.clone());

        assert_eq!(container.len(), 2);
        assert!(container.contains::<FlowMetadata>());
        assert!(container.contains::<ThermalMetadata>());

        let retrieved_flow = container.get::<FlowMetadata>().unwrap();
        let retrieved_thermal = container.get::<ThermalMetadata>().unwrap();

        assert_eq!(retrieved_flow, &flow_data);
        assert_eq!(retrieved_thermal, &thermal_data);
    }
}
