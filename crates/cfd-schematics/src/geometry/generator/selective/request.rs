use super::geometry::SelectiveTreeGeometry;
use aequitas::systems::si::quantities::Length;

#[derive(Debug, Clone, Copy, PartialEq)]
/// Describes the center serpentine overlay used by a selective tree.
pub struct CenterSerpentinePathSpec {
    /// Number of segments in the serpentine path.
    pub segments: usize,
    /// Bend radius of the path in metres.
    pub bend_radius_m: Length<f64>,
    /// Waveform type for the serpentine path.  Defaults to Sine when
    /// constructed from legacy topology specs that omit this field.
    pub wave_type: crate::topology::SerpentineWaveType,
}

#[derive(Debug, Clone, PartialEq)]
/// Selects the topology-specific construction parameters for a selective tree.
pub enum SelectiveTreeTopology {
    /// Builds a cascade of center trifurcations.
    CascadeCenterTrifurcation {
        /// Number of trifurcation levels in the cascade.
        n_levels: usize,
        /// Fraction of each split assigned to the center branch.
        center_frac: f64,
        /// Enables Venturi treatment geometry on treatment branches.
        venturi_treatment_enabled: bool,
        /// Optional center serpentine overlay.
        center_serpentine: Option<CenterSerpentinePathSpec>,
    },
    /// Builds incremental filtration stages with terminal tri-bifurcation.
    IncrementalFiltrationTriBi {
        /// Number of preliminary trifurcation stages.
        n_pretri: usize,
        /// Center fraction used by preliminary trifurcations.
        pretri_center_frac: f64,
        /// Center fraction used by the terminal trifurcation.
        terminal_tri_center_frac: f64,
        /// Fraction assigned to the treatment branch at the bifurcation.
        bi_treat_frac: f64,
        /// Enables Venturi treatment geometry on treatment branches.
        venturi_treatment_enabled: bool,
        /// Optional center serpentine overlay.
        center_serpentine: Option<CenterSerpentinePathSpec>,
        /// Length of the outlet tail in metres.
        outlet_tail_length_m: Length<f64>,
    },
    /// Builds a tri-bifurcation-bifurcation-trifurcation selective topology.
    TriBiTriSelective {
        /// Center fraction at the first trifurcation.
        first_center_frac: f64,
        /// Fraction assigned to the treatment branch at the bifurcation.
        bi_treat_frac: f64,
        /// Center fraction at the second trifurcation.
        second_center_frac: f64,
    },
    /// Builds a double-trifurcation CIF topology with a center Venturi train.
    DoubleTrifurcationCif {
        /// Center fraction at the first trifurcation.
        split1_center_frac: f64,
        /// Center fraction at the second trifurcation.
        split2_center_frac: f64,
        /// Number of serial center throats.
        center_throat_count: u8,
        /// Spacing between center throats in metres.
        inter_throat_spacing_m: Length<f64>,
    },
}

#[derive(Debug, Clone, PartialEq)]
/// Defines the physical dimensions and topology of a selective-tree blueprint.
pub struct SelectiveTreeRequest {
    /// Human-readable name assigned to the generated blueprint.
    pub name: String,
    /// Overall rectangular envelope dimensions in metres.
    pub box_dims_m: (Length<f64>, Length<f64>),
    /// Length of the inlet trunk in metres.
    pub trunk_length_m: Length<f64>,
    /// Length of ordinary branches in metres.
    pub branch_length_m: Length<f64>,
    /// Length of hybrid branches in metres.
    pub hybrid_branch_length_m: Length<f64>,
    /// Width of the main channels in metres.
    pub main_width_m: Length<f64>,
    /// Width of Venturi throats in metres.
    pub throat_width_m: Length<f64>,
    /// Length of Venturi throats in metres.
    pub throat_length_m: Length<f64>,
    /// Channel height in metres.
    pub channel_height_m: Length<f64>,
    /// Topology-specific construction parameters.
    pub topology: SelectiveTreeTopology,
}

impl SelectiveTreeRequest {
    fn box_dims_mm(&self) -> (f64, f64) {
        (
            self.box_dims_m.0.into_base() * 1.0e3,
            self.box_dims_m.1.into_base() * 1.0e3,
        )
    }

    pub(super) fn geometry(&self) -> SelectiveTreeGeometry {
        SelectiveTreeGeometry {
            box_dims_mm: self.box_dims_mm(),
            trunk_length_m: self.trunk_length_m.into_base(),
            branch_length_m: self.branch_length_m.into_base(),
            hybrid_branch_length_m: self.hybrid_branch_length_m.into_base(),
            main_width_m: self.main_width_m.into_base(),
            throat_width_m: self.throat_width_m.into_base(),
            throat_length_m: self.throat_length_m.into_base(),
            channel_height_m: self.channel_height_m.into_base(),
        }
    }
}

