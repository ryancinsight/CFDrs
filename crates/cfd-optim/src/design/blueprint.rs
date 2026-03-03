//! Blueprint and channel-system conversion methods for `DesignCandidate`.

use super::candidate::DesignCandidate;
use super::topology::DesignTopology;
use crate::constraints::{PLATE_HEIGHT_MM, PLATE_WIDTH_MM, TREATMENT_HEIGHT_MM};
use cfd_schematics::{
    interface::presets::{
        asymmetric_bifurcation_serpentine_rect, asymmetric_trifurcation_venturi_rect,
        bifurcation_serpentine_rect, bifurcation_trifurcation_venturi_rect,
        bifurcation_venturi_rect, cascade_center_trifurcation_rect,
        cascade_tri_bi_tri_selective_rect, cell_separation_rect, constriction_expansion_array_rect,
        double_bifurcation_venturi_rect, double_trifurcation_venturi_rect,
        incremental_filtration_tri_bi_rect_staged_remerge, parallel_microchannel_array_rect,
        quad_trifurcation_venturi_rect, serial_double_venturi_rect, serpentine_rect,
        spiral_channel_rect, trifurcation_bifurcation_bifurcation_venturi_rect,
        trifurcation_bifurcation_venturi_rect, trifurcation_serpentine_rect,
        trifurcation_venturi_rect, triple_bifurcation_venturi_rect,
        triple_trifurcation_venturi_rect, venturi_rect, venturi_serpentine_rect,
    },
    NetworkBlueprint,
};

impl DesignCandidate {
    /// Build a [`NetworkBlueprint`] for this candidate using `cfd-schematics` presets.
    ///
    /// The blueprint encodes the channel cross-sections in `ChannelSpec.cross_section`
    /// so that downstream consumers (e.g. `compute_metrics`) can read geometry via
    /// `CrossSectionSpec::area()` / `hydraulic_diameter()` rather than recomputing it
    /// from raw parameter fields.
    ///
    /// **Naming convention** (relied on by `compute_metrics`):
    /// - `"throat_section"` — venturi throat (`throat_diameter_m × channel_height_m`)
    /// - `"inlet_section"` — venturi inlet / main channel approach
    /// - `"segment_1"`, `"segment_2"`, … — serpentine straight segments
    /// - `"parent"` — bifurcation / trifurcation trunk
    #[must_use]
    pub fn to_blueprint(&self) -> NetworkBlueprint {
        let w = self.channel_width_m;
        let h = self.channel_height_m;
        let dt = self.throat_diameter_m;
        let tl = self.throat_length_m.max(dt * 2.0);
        let n = self.serpentine_segments;
        let sl = self.segment_length_m;

        match self.topology {
            // ── Single venturi ──
            DesignTopology::SingleVenturi => venturi_rect(&self.id, w, dt, h, tl),

            // ── Pure serpentine grid ──
            DesignTopology::SerpentineGrid => serpentine_rect(&self.id, n, sl, w, h),

            // ── Venturi → serpentine (series, closed-loop) ──
            DesignTopology::VenturiSerpentine => {
                venturi_serpentine_rect(&self.id, w, dt, h, tl, n, sl)
            }

            // ── Bifurcation + venturi in each branch (closed-loop) ──
            DesignTopology::BifurcationVenturi => {
                let trunk_len = TREATMENT_HEIGHT_MM * 0.25e-3;
                bifurcation_venturi_rect(&self.id, trunk_len, w, dt, h, tl)
            }

            // ── Trifurcation + venturi in each branch (closed-loop) ──
            DesignTopology::TrifurcationVenturi => {
                let trunk_len = TREATMENT_HEIGHT_MM * 0.50e-3;
                trifurcation_venturi_rect(&self.id, trunk_len, w, dt, h, tl)
            }

            // ── Cell separation: center venturi + peripheral bypass (closed-loop) ──
            DesignTopology::CellSeparationVenturi | DesignTopology::WbcCancerSeparationVenturi => {
                let approach_len = TREATMENT_HEIGHT_MM * 0.5e-3;
                cell_separation_rect(&self.id, approach_len, w, dt, h, tl, approach_len)
            }

            // ── DoubleBifurcation [Bi,Bi] → 4 parallel venturis (closed-loop) ──
            DesignTopology::DoubleBifurcationVenturi => {
                let trunk_len = TREATMENT_HEIGHT_MM * 0.20e-3;
                let branch1_len = TREATMENT_HEIGHT_MM * 0.15e-3;
                double_bifurcation_venturi_rect(&self.id, trunk_len, branch1_len, w, dt, h, tl)
            }

            // ── TripleBifurcation [Bi,Bi,Bi] → 8 parallel venturis (closed-loop) ──
            DesignTopology::TripleBifurcationVenturi => {
                let trunk_len = TREATMENT_HEIGHT_MM * 0.15e-3;
                let branch1_len = TREATMENT_HEIGHT_MM * 0.12e-3;
                let branch2_len = TREATMENT_HEIGHT_MM * 0.10e-3;
                triple_bifurcation_venturi_rect(
                    &self.id,
                    trunk_len,
                    branch1_len,
                    branch2_len,
                    w,
                    dt,
                    h,
                    tl,
                )
            }

            // ── DoubleTrifurcation [Tri,Tri] → 9 parallel venturis (closed-loop) ──
            DesignTopology::DoubleTrifurcationVenturi => {
                let trunk_len = TREATMENT_HEIGHT_MM * 0.20e-3;
                let branch_len = TREATMENT_HEIGHT_MM * 0.15e-3;
                double_trifurcation_venturi_rect(&self.id, trunk_len, branch_len, w, dt, h, tl)
            }

            // ── BifurcationTrifurcation [Bi,Tri] → 6 parallel venturis (closed-loop) ──
            DesignTopology::BifurcationTrifurcationVenturi => {
                let trunk_len = TREATMENT_HEIGHT_MM * 0.25e-3;
                let branch1_len = TREATMENT_HEIGHT_MM * 0.20e-3;
                bifurcation_trifurcation_venturi_rect(
                    &self.id,
                    trunk_len,
                    branch1_len,
                    w,
                    dt,
                    h,
                    tl,
                )
            }

            // ── TrifurcationBifurcation [Tri,Bi] → 6 parallel venturis (closed-loop) ──
            DesignTopology::TrifurcationBifurcationVenturi => {
                let trunk_len = TREATMENT_HEIGHT_MM * 0.25e-3;
                let branch1_len = TREATMENT_HEIGHT_MM * 0.20e-3;
                trifurcation_bifurcation_venturi_rect(
                    &self.id,
                    trunk_len,
                    branch1_len,
                    w,
                    dt,
                    h,
                    tl,
                )
            }

            // ── TripleTrifurcation [Tri,Tri,Tri] → 27 scaled-width venturis ──
            DesignTopology::TripleTrifurcationVenturi => {
                let trunk_len = TREATMENT_HEIGHT_MM * 0.15e-3;
                let branch1_len = TREATMENT_HEIGHT_MM * 0.12e-3;
                let branch2_len = TREATMENT_HEIGHT_MM * 0.09e-3;
                let frac = self.trifurcation_center_frac;
                triple_trifurcation_venturi_rect(
                    &self.id,
                    trunk_len,
                    branch1_len,
                    branch2_len,
                    w,
                    frac,
                    dt,
                    h,
                    tl,
                )
            }

            // ── TrifurcationBifurcationBifurcation [Tri,Bi,Bi] → 12 scaled-width venturis ──
            DesignTopology::TrifurcationBifurcationBifurcationVenturi => {
                let trunk_len = TREATMENT_HEIGHT_MM * 0.18e-3;
                let branch1_len = TREATMENT_HEIGHT_MM * 0.15e-3;
                let branch2_len = TREATMENT_HEIGHT_MM * 0.11e-3;
                let frac = self.trifurcation_center_frac;
                trifurcation_bifurcation_bifurcation_venturi_rect(
                    &self.id,
                    trunk_len,
                    branch1_len,
                    branch2_len,
                    w,
                    frac,
                    dt,
                    h,
                    tl,
                )
            }

            // ── QuadTrifurcation [Tri^4] → 81 scaled-width venturis ──
            DesignTopology::QuadTrifurcationVenturi => {
                let trunk_len = TREATMENT_HEIGHT_MM * 0.12e-3;
                let branch1_len = TREATMENT_HEIGHT_MM * 0.10e-3;
                let branch2_len = TREATMENT_HEIGHT_MM * 0.08e-3;
                let branch3_len = TREATMENT_HEIGHT_MM * 0.06e-3;
                let frac = self.trifurcation_center_frac;
                quad_trifurcation_venturi_rect(
                    &self.id,
                    trunk_len,
                    branch1_len,
                    branch2_len,
                    branch3_len,
                    w,
                    frac,
                    dt,
                    h,
                    tl,
                )
            }

            // ── CascadeCenterTrifurcationSeparator → single venturi on center arm ──
            DesignTopology::CascadeCenterTrifurcationSeparator { n_levels } => {
                let trunk_len = TREATMENT_HEIGHT_MM * 0.20e-3;
                let branch_len = TREATMENT_HEIGHT_MM * 0.15e-3;
                let frac = self.trifurcation_center_frac;
                cascade_center_trifurcation_rect(
                    &self.id, trunk_len, branch_len, n_levels, w, frac, dt, tl, h,
                )
            }

            // ── CIF tri-first then tri/bi skimming → single treatment venturi ──
            DesignTopology::IncrementalFiltrationTriBiSeparator { n_pretri } => {
                use cfd_1d::cell_separation::{cif_pretri_stage_q_fracs, tri_center_q_frac};

                let trunk_len = TREATMENT_HEIGHT_MM * 0.20e-3;
                let pretri_len = TREATMENT_HEIGHT_MM * 0.15e-3;
                let hybrid_len = TREATMENT_HEIGHT_MM * 0.12e-3;
                // Keep CIF streams separated through most of the device and
                // remerge close to outlet to preserve selective venturi exposure.
                // Higher selective venturi-flow fractions tolerate a slightly longer
                // post-remerge tail; strongly selective layouts keep the tail short.
                let q_pretri_product = cif_pretri_stage_q_fracs(
                    n_pretri,
                    self.cif_pretri_center_frac(),
                    self.cif_terminal_tri_center_frac(),
                )
                .into_iter()
                .product::<f64>();
                let q_terminal_tri = tri_center_q_frac(self.cif_terminal_tri_center_frac());
                let q_bi_treat = self.cif_terminal_bi_treat_frac();
                let selective_qfrac =
                    (q_pretri_product * q_terminal_tri * q_bi_treat).clamp(0.0, 1.0);
                let outlet_tail_factor = (0.08 + 0.25 * selective_qfrac).clamp(0.08, 0.30);
                let outlet_tail_len = trunk_len * outlet_tail_factor;
                incremental_filtration_tri_bi_rect_staged_remerge(
                    &self.id,
                    trunk_len,
                    pretri_len,
                    hybrid_len,
                    n_pretri,
                    w,
                    self.cif_pretri_center_frac(),
                    self.cif_terminal_tri_center_frac(),
                    self.cif_terminal_bi_treat_frac(),
                    dt,
                    tl,
                    outlet_tail_len,
                    h,
                )
            }

            // ── Two venturis in series on one path (closed-loop) ──
            DesignTopology::SerialDoubleVenturi => {
                let inter_len = TREATMENT_HEIGHT_MM * 0.40e-3; // 18 mm inter-venturi
                serial_double_venturi_rect(&self.id, w, dt, h, tl, inter_len)
            }

            // ── Bifurcation + serpentine in each arm (closed-loop) ──
            DesignTopology::BifurcationSerpentine => {
                let trunk_len = TREATMENT_HEIGHT_MM * 0.25e-3;
                bifurcation_serpentine_rect(&self.id, trunk_len, n, sl, w, h)
            }

            // ── Trifurcation + serpentine in each arm (closed-loop) ──
            DesignTopology::TrifurcationSerpentine => {
                let trunk_len = TREATMENT_HEIGHT_MM * 0.33e-3;
                trifurcation_serpentine_rect(&self.id, trunk_len, n, sl, w, h)
            }

            // ── Asymmetric bifurcation → wide serpentine (cancer/WBC) + narrow bypass (RBC) ──
            DesignTopology::AsymmetricBifurcationSerpentine => {
                let trunk_len = TREATMENT_HEIGHT_MM * 1e-3 / 6.0;
                asymmetric_bifurcation_serpentine_rect(
                    &self.id,
                    trunk_len,
                    n,
                    sl,
                    w,
                    self.asymmetric_narrow_frac,
                    h,
                )
            }

            // ── Asymmetric 3-stream trifurcation with center-only venturi ──
            DesignTopology::AsymmetricTrifurcationVenturi => {
                let trunk_len = TREATMENT_HEIGHT_MM * 1e-3 / 6.0;
                let branch_len = TREATMENT_HEIGHT_MM * 1e-3 * 0.25;
                asymmetric_trifurcation_venturi_rect(
                    &self.id,
                    trunk_len,
                    branch_len,
                    w,
                    self.trifurcation_center_frac,
                    self.trifurcation_left_frac,
                    dt,
                    tl,
                    h,
                )
            }

            // ── Tri→Bi→Tri selective center venturi ──
            DesignTopology::TriBiTriSelectiveVenturi => {
                let trunk_len = TREATMENT_HEIGHT_MM * 1e-3 / 6.0;
                let branch_len = TREATMENT_HEIGHT_MM * 1e-3 * 0.15;
                cascade_tri_bi_tri_selective_rect(
                    &self.id,
                    trunk_len,
                    branch_len,
                    w,
                    self.trifurcation_center_frac,
                    self.cif_terminal_bi_treat_frac,
                    self.trifurcation_center_frac,
                    dt,
                    tl,
                    h,
                )
            }

            // ── Constriction-expansion array (wide→narrow cycles, WBC margination) ──
            DesignTopology::ConstrictionExpansionArray { n_cycles } => {
                let wide_w = w;
                let narrow_w = w * 0.40;
                // Each cycle has a wide + narrow segment; split treatment length evenly
                let seg_len = TREATMENT_HEIGHT_MM * 1e-3 / (n_cycles as f64 * 2.0);
                constriction_expansion_array_rect(
                    &self.id,
                    n_cycles,
                    seg_len,
                    seg_len * 0.5,
                    wide_w,
                    narrow_w,
                    h,
                )
            }

            // ── Spiral serpentine (N turns, high Dean number, WBC/RBC separation) ──
            DesignTopology::SpiralSerpentine { n_turns } => {
                let turn_len = TREATMENT_HEIGHT_MM * 1e-3 / n_turns as f64;
                spiral_channel_rect(&self.id, n_turns, turn_len, w, h)
            }

            // ── Parallel microchannel array (N identical channels, clinical throughput) ──
            DesignTopology::ParallelMicrochannelArray { n_channels } => {
                let ch_len = TREATMENT_HEIGHT_MM * 1e-3;
                parallel_microchannel_array_rect(&self.id, n_channels, ch_len, w, h)
            }

            // ── AdaptiveTree: variable-depth split tree → venturi at every leaf ──
            DesignTopology::AdaptiveTree {
                levels,
                split_types,
            } => {
                let trunk_len = TREATMENT_HEIGHT_MM * 1e-3 * 0.5 / (f64::from(levels) + 1.0);
                let branch_len = trunk_len * 0.8;
                match levels {
                    0 => venturi_rect(&self.id, w, dt, h, tl),
                    1 => {
                        let level_len = TREATMENT_HEIGHT_MM * 1e-3 / 4.0;
                        if split_types & 1 == 0 {
                            bifurcation_venturi_rect(&self.id, level_len, w, dt, h, tl)
                        } else {
                            trifurcation_venturi_rect(&self.id, level_len, w, dt, h, tl)
                        }
                    }
                    2 => match split_types & 0b11 {
                        // Bi→Bi = 4 leaves
                        0b00 => double_bifurcation_venturi_rect(
                            &self.id, trunk_len, branch_len, w, dt, h, tl,
                        ),
                        // Bi→Tri = 6 leaves
                        0b10 => bifurcation_trifurcation_venturi_rect(
                            &self.id, trunk_len, branch_len, w, dt, h, tl,
                        ),
                        // Tri→Bi = 6 leaves (new preset)
                        0b01 => trifurcation_bifurcation_venturi_rect(
                            &self.id, trunk_len, branch_len, w, dt, h, tl,
                        ),
                        // Tri→Tri = 9 leaves
                        _ => double_trifurcation_venturi_rect(
                            &self.id, trunk_len, branch_len, w, dt, h, tl,
                        ),
                    },
                    // levels ≥3: use CCT-style selective center-only venturi
                    // (protects RBCs while cancer-enriched center stream gets treatment)
                    _ => {
                        let lv = levels.min(5);
                        let cf = self.trifurcation_center_frac;
                        cascade_center_trifurcation_rect(
                            &self.id, trunk_len, branch_len, lv, w, cf, dt, tl, h,
                        )
                    }
                }
            }
        }
    }

    /// Build a [`cfd_schematics::geometry::ChannelSystem`] for 2D schematic
    /// visualisation using `cfd-schematics`.
    ///
    /// All physical units are converted from metres (candidate fields) to
    /// millimetres (the schematics API).  The generated system mirrors the
    /// x/y symmetry shown in existing `cfd-schematics` examples.
    ///
    /// Returns a fully-populated `ChannelSystem` that can be passed directly
    /// to [`cfd_schematics::plot_geometry`].
    #[must_use]
    pub fn to_channel_system(&self) -> cfd_schematics::geometry::ChannelSystem {
        use cfd_schematics::{
            config::{
                ArcConfig, ChannelTypeConfig, FrustumConfig, GeometryConfig, SerpentineConfig,
                TaperProfile,
            },
            geometry::{create_geometry, SplitType},
        };

        let w_mm = self.channel_width_m * 1e3; // e.g. 2.0 mm
        let h_mm = self.channel_height_m * 1e3; // e.g. 0.5 mm
                                                // Clamp throat to valid FrustumConfig range: strictly inside (0, w_mm)
        let dt_mm = (self.throat_diameter_m * 1e3).max(0.06).min(w_mm * 0.9);

        let gc = GeometryConfig {
            wall_clearance: 4.0,
            channel_width: w_mm,
            channel_height: h_mm,
            ..Default::default()
        };

        // Full-plate layout: use the entire 96-well footprint for schematic
        // generation rather than constraining topology to the 6x6 center zone.
        let box_dims = (PLATE_WIDTH_MM, PLATE_HEIGHT_MM);

        // ── Wave-shape configs ───────────────────────────────────────────
        // Sine wave — smooth Dean-flow visualisation for tree / venturi topologies.
        let sine_wave = SerpentineConfig {
            fill_factor: 0.75,
            wave_density_factor: 2.5,
            gaussian_width_factor: 4.0,
            ..SerpentineConfig::default()
        };

        // Square wave — high-density sharp wave for cell-size separation designs.
        let square_wave = SerpentineConfig {
            fill_factor: 0.80,
            wave_density_factor: 4.0,
            gaussian_width_factor: 4.0,
            ..SerpentineConfig::default()
        }
        .with_square_wave();

        let system = match self.topology {
            // ── Single venturi: sine wave, no splits ──
            DesignTopology::SingleVenturi | DesignTopology::SerialDoubleVenturi => create_geometry(
                box_dims,
                &[],
                &gc,
                &ChannelTypeConfig::AllSerpentine(sine_wave),
            ),

            // ── Cell separation topologies: square wave + asymmetric bifurcation ──
            DesignTopology::CellSeparationVenturi | DesignTopology::WbcCancerSeparationVenturi => {
                create_geometry(
                    box_dims,
                    &[SplitType::AsymmetricBifurcation { ratio: 0.67 }],
                    &gc,
                    &ChannelTypeConfig::AllSerpentine(square_wave),
                )
            }

            // ── Pure serpentine grid ──
            DesignTopology::SerpentineGrid => create_geometry(
                box_dims,
                &[],
                &gc,
                &ChannelTypeConfig::AllSerpentine(sine_wave),
            ),

            // ── Venturi + serpentine (series): adaptive (sine + frustum) ──
            DesignTopology::VenturiSerpentine => {
                let frustum_cfg = FrustumConfig {
                    inlet_width: w_mm,
                    throat_width: dt_mm,
                    outlet_width: w_mm,
                    taper_profile: TaperProfile::Smooth,
                    smoothness: 50,
                    throat_position: 0.5,
                };
                create_geometry(
                    box_dims,
                    &[],
                    &gc,
                    &ChannelTypeConfig::Adaptive {
                        serpentine_config: sine_wave,
                        arc_config: ArcConfig::default(),
                        frustum_config: frustum_cfg,
                    },
                )
            }

            // ── Binary venturi trees → sine wave ──
            DesignTopology::BifurcationVenturi => create_geometry(
                box_dims,
                &[SplitType::Bifurcation],
                &gc,
                &ChannelTypeConfig::AllSerpentine(sine_wave),
            ),
            DesignTopology::TrifurcationVenturi => create_geometry(
                box_dims,
                &[SplitType::Trifurcation],
                &gc,
                &ChannelTypeConfig::AllSerpentine(sine_wave),
            ),
            DesignTopology::DoubleBifurcationVenturi => create_geometry(
                box_dims,
                &[SplitType::Bifurcation, SplitType::Bifurcation],
                &gc,
                &ChannelTypeConfig::AllSerpentine(sine_wave),
            ),
            DesignTopology::TripleBifurcationVenturi => create_geometry(
                box_dims,
                &[
                    SplitType::Bifurcation,
                    SplitType::Bifurcation,
                    SplitType::Bifurcation,
                ],
                &gc,
                &ChannelTypeConfig::AllSerpentine(sine_wave),
            ),
            DesignTopology::DoubleTrifurcationVenturi => create_geometry(
                box_dims,
                &[SplitType::Trifurcation, SplitType::Trifurcation],
                &gc,
                &ChannelTypeConfig::AllSerpentine(sine_wave),
            ),
            DesignTopology::BifurcationTrifurcationVenturi => create_geometry(
                box_dims,
                &[SplitType::Bifurcation, SplitType::Trifurcation],
                &gc,
                &ChannelTypeConfig::AllSerpentine(sine_wave),
            ),

            // ── Symmetric serpentine arms → sine wave ──
            DesignTopology::BifurcationSerpentine => create_geometry(
                box_dims,
                &[SplitType::Bifurcation],
                &gc,
                &ChannelTypeConfig::AllSerpentine(sine_wave),
            ),
            DesignTopology::TrifurcationSerpentine => create_geometry(
                box_dims,
                &[SplitType::Trifurcation],
                &gc,
                &ChannelTypeConfig::AllSerpentine(sine_wave),
            ),

            // ── Asymmetric bifurcation → square wave (cell-size separation) ──
            DesignTopology::AsymmetricBifurcationSerpentine => create_geometry(
                box_dims,
                &[SplitType::AsymmetricBifurcation { ratio: 0.67 }],
                &gc,
                &ChannelTypeConfig::AllSerpentine(square_wave),
            ),

            // ── Constriction-expansion array: square wave (WBC margination) ──
            DesignTopology::ConstrictionExpansionArray { .. } => create_geometry(
                box_dims,
                &[],
                &gc,
                &ChannelTypeConfig::AllSerpentine(square_wave),
            ),

            // ── Spiral serpentine: sine wave (Dean-flow separation) ──
            DesignTopology::SpiralSerpentine { .. } => create_geometry(
                box_dims,
                &[],
                &gc,
                &ChannelTypeConfig::AllSerpentine(sine_wave),
            ),

            // ── Parallel microchannel array: sine wave ──
            DesignTopology::ParallelMicrochannelArray { .. } => create_geometry(
                box_dims,
                &[],
                &gc,
                &ChannelTypeConfig::AllSerpentine(sine_wave),
            ),

            // ── New multi-level trifurcation / cascade topologies → sine wave ──
            DesignTopology::TrifurcationBifurcationVenturi => create_geometry(
                box_dims,
                &[SplitType::Trifurcation, SplitType::Bifurcation],
                &gc,
                &ChannelTypeConfig::AllSerpentine(sine_wave),
            ),
            DesignTopology::TripleTrifurcationVenturi => create_geometry(
                box_dims,
                &[
                    SplitType::Trifurcation,
                    SplitType::Trifurcation,
                    SplitType::Trifurcation,
                ],
                &gc,
                &ChannelTypeConfig::AllSerpentine(sine_wave),
            ),
            DesignTopology::TrifurcationBifurcationBifurcationVenturi => create_geometry(
                box_dims,
                &[
                    SplitType::Trifurcation,
                    SplitType::Bifurcation,
                    SplitType::Bifurcation,
                ],
                &gc,
                &ChannelTypeConfig::AllSerpentine(sine_wave),
            ),
            DesignTopology::QuadTrifurcationVenturi => create_geometry(
                box_dims,
                &[
                    SplitType::Trifurcation,
                    SplitType::Trifurcation,
                    SplitType::Trifurcation,
                    SplitType::Trifurcation,
                ],
                &gc,
                &ChannelTypeConfig::AllSerpentine(sine_wave),
            ),
            DesignTopology::CascadeCenterTrifurcationSeparator { n_levels } => {
                // Cascade topology: cap at 3 schematic splits to keep ≤ 27
                // parallel channels and prevent visual overlap.
                let schematic_levels = (n_levels as usize).min(3);
                let splits: Vec<SplitType> = (0..schematic_levels)
                    .map(|_| SplitType::SymmetricTrifurcation {
                        center_ratio: self.trifurcation_center_frac,
                    })
                    .collect();
                create_geometry(
                    box_dims,
                    &splits,
                    &gc,
                    &ChannelTypeConfig::AllSerpentine(sine_wave),
                )
            }
            DesignTopology::IncrementalFiltrationTriBiSeparator { n_pretri } => {
                // CIF schematic: cap pre-tri splits so parallel channels remain
                // visually distinct.  Full n_pretri can produce 3^(n+1)*2 leaf
                // channels; at n_pretri=3 that is 162 — far too many for
                // 85 mm plate height.  Show at most 1 representative pre-tri
                // to keep the total ≤ 18 parallel paths.
                let schematic_pretri = n_pretri.min(1);
                let mut splits: Vec<SplitType> = (0..schematic_pretri as usize)
                    .map(|_| SplitType::SymmetricTrifurcation {
                        center_ratio: self.cif_pretri_center_frac(),
                    })
                    .collect();
                splits.push(SplitType::SymmetricTrifurcation {
                    center_ratio: self.cif_terminal_tri_center_frac(),
                });
                splits.push(SplitType::AsymmetricBifurcation {
                    ratio: self.cif_terminal_bi_treat_frac(),
                });
                create_geometry(
                    box_dims,
                    &splits,
                    &gc,
                    &ChannelTypeConfig::AllSerpentine(sine_wave),
                )
            }

            // ── Asymmetric 3-stream trifurcation → sine wave ──
            DesignTopology::AsymmetricTrifurcationVenturi => create_geometry(
                box_dims,
                &[SplitType::Trifurcation],
                &gc,
                &ChannelTypeConfig::AllSerpentine(sine_wave),
            ),

            // ── Tri→Bi→Tri selective center venturi → sine wave ──
            DesignTopology::TriBiTriSelectiveVenturi => create_geometry(
                box_dims,
                &[
                    SplitType::Trifurcation,
                    SplitType::Bifurcation,
                    SplitType::Trifurcation,
                ],
                &gc,
                &ChannelTypeConfig::AllSerpentine(sine_wave),
            ),

            // ── Adaptive tree: decode split sequence from packed bits → sine wave ──
            DesignTopology::AdaptiveTree {
                levels,
                split_types,
            } => {
                let splits: Vec<SplitType> = (0..levels as usize)
                    .map(|i| {
                        if (split_types >> i) & 1 == 0 {
                            SplitType::Bifurcation
                        } else {
                            SplitType::Trifurcation
                        }
                    })
                    .collect();
                create_geometry(
                    box_dims,
                    &splits,
                    &gc,
                    &ChannelTypeConfig::AllSerpentine(sine_wave),
                )
            }
        };
        map_to_plate_coords(system, box_dims, w_mm, h_mm)
    }

    /// Total approximate channel path length [mm].
    #[must_use]
    pub fn total_path_length_mm(&self) -> f64 {
        let mut len = 0.0_f64;
        // Venturi section (inlet cone + throat + diffuser, approximated as 10× d_inlet)
        if self.topology.has_venturi() {
            len += self.inlet_diameter_m * 10.0 * self.topology.venturi_count() as f64 * 1000.0;
        }
        // Serpentine section
        if self.topology.has_serpentine() {
            let n = self.serpentine_segments as f64;
            let bends = (self.serpentine_segments.saturating_sub(1)) as f64;
            len += (n * self.segment_length_m + bends * std::f64::consts::PI * self.bend_radius_m)
                * 1000.0;
        }
        len
    }
}

/// Remap a treatment-zone–local `ChannelSystem` into full 96-well plate
/// coordinates (ANSI/SLAS 1-2004, 127.76 × 85.47 mm).
///
/// The geometry generator produces channels in a coordinate frame matching
/// `local_dims`.  When `local_dims` already equals the plate envelope the
/// remap is an identity; only the inlet/outlet trunk insertion is needed.
///
/// Previous implementation computed the bounding box from all path points
/// (including serpentine oscillation extremes) and rescaled the entire
/// system to fill the plate, which undid the edge-padding that keeps
/// peripheral channels away from the plate walls.  Now we use `local_dims`
/// as the authoritative source-of-truth for the coordinate frame.
fn map_to_plate_coords(
    mut sys: cfd_schematics::geometry::ChannelSystem,
    local_dims: (f64, f64),
    default_channel_width_mm: f64,
    default_channel_height_mm: f64,
) -> cfd_schematics::geometry::ChannelSystem {
    use cfd_schematics::geometry::{Channel, ChannelType, Node};

    let (box_w, box_h) = local_dims;

    // Use the generator's own coordinate frame rather than the point
    // bounding-box so that edge-padding is preserved.
    let (x_min, x_max) = (0.0, box_w);
    let (y_min, y_max) = (0.0, box_h);

    let span_x = (x_max - x_min).abs();
    let span_y = (y_max - y_min).abs();
    let sx = if span_x > 1e-9 {
        PLATE_WIDTH_MM / span_x
    } else {
        1.0
    };
    let sy = if span_y > 1e-9 {
        PLATE_HEIGHT_MM / span_y
    } else {
        1.0
    };

    let remap_point = |pt: &mut (f64, f64)| {
        pt.0 = (pt.0 - x_min) * sx;
        pt.1 = (pt.1 - y_min) * sy;
    };

    for node in &mut sys.nodes {
        remap_point(&mut node.point);
    }

    for channel in &mut sys.channels {
        match &mut channel.channel_type {
            ChannelType::Straight => {} // coordinates resolved from nodes
            ChannelType::SmoothStraight { path }
            | ChannelType::Serpentine { path }
            | ChannelType::Arc { path }
            | ChannelType::Frustum { path, .. } => {
                for pt in path.iter_mut() {
                    remap_point(pt);
                }
            }
        }
    }

    sys.box_dims = (PLATE_WIDTH_MM, PLATE_HEIGHT_MM);
    sys.box_outline = vec![
        ((0.0, 0.0), (PLATE_WIDTH_MM, 0.0)),
        ((PLATE_WIDTH_MM, 0.0), (PLATE_WIDTH_MM, PLATE_HEIGHT_MM)),
        ((PLATE_WIDTH_MM, PLATE_HEIGHT_MM), (0.0, PLATE_HEIGHT_MM)),
        ((0.0, PLATE_HEIGHT_MM), (0.0, 0.0)),
    ];

    // Ensure the rendered schematic has explicit plate-edge inlet/outlet
    // trunks so the flow path connects into and out of the full cuboid.
    if !sys.channels.is_empty() && !sys.nodes.is_empty() {
        let mut degree = vec![0usize; sys.nodes.len()];
        for ch in &sys.channels {
            if ch.from_node < degree.len() {
                degree[ch.from_node] += 1;
            }
            if ch.to_node < degree.len() {
                degree[ch.to_node] += 1;
            }
        }

        let mut terminals: Vec<usize> = degree
            .iter()
            .enumerate()
            .filter_map(|(idx, d)| if *d == 1 { Some(idx) } else { None })
            .collect();

        if terminals.len() >= 2 {
            terminals.sort_by(|a, b| sys.nodes[*a].point.0.total_cmp(&sys.nodes[*b].point.0));
            let inlet_node = terminals[0];
            let outlet_node = terminals[terminals.len() - 1];
            let next_node_id_start = sys.nodes.iter().map(|n| n.id).max().unwrap_or(0) + 1;
            let next_channel_id_start = sys.channels.iter().map(|c| c.id).max().unwrap_or(0) + 1;
            let mut next_node_id = next_node_id_start;
            let mut next_channel_id = next_channel_id_start;

            const EDGE_EPS: f64 = 1e-6;

            let inlet_x = sys.nodes[inlet_node].point.0;
            let inlet_y = sys.nodes[inlet_node].point.1;
            if inlet_x > EDGE_EPS {
                let port_id = next_node_id;
                next_node_id += 1;
                sys.nodes.push(Node {
                    id: port_id,
                    point: (0.0, inlet_y),
                    metadata: None,
                });
                let (w, h) = sys
                    .channels
                    .iter()
                    .find(|c| c.from_node == inlet_node || c.to_node == inlet_node)
                    .map_or((default_channel_width_mm, default_channel_height_mm), |c| {
                        (c.width, c.height)
                    });
                sys.channels.push(Channel {
                    id: next_channel_id,
                    from_node: port_id,
                    to_node: inlet_node,
                    width: w,
                    height: h,
                    channel_type: ChannelType::Straight,
                    metadata: None,
                });
                next_channel_id += 1;
            }

            let outlet_x = sys.nodes[outlet_node].point.0;
            let outlet_y = sys.nodes[outlet_node].point.1;
            if outlet_x < PLATE_WIDTH_MM - EDGE_EPS {
                let port_id = next_node_id;
                sys.nodes.push(Node {
                    id: port_id,
                    point: (PLATE_WIDTH_MM, outlet_y),
                    metadata: None,
                });
                let (w, h) = sys
                    .channels
                    .iter()
                    .find(|c| c.from_node == outlet_node || c.to_node == outlet_node)
                    .map_or((default_channel_width_mm, default_channel_height_mm), |c| {
                        (c.width, c.height)
                    });
                sys.channels.push(Channel {
                    id: next_channel_id,
                    from_node: outlet_node,
                    to_node: port_id,
                    width: w,
                    height: h,
                    channel_type: ChannelType::Straight,
                    metadata: None,
                });
            }
        }
    }

    sys
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::design::{CrossSectionShape, DesignCandidate, DesignTopology};

    fn sample_candidate(topology: DesignTopology) -> DesignCandidate {
        DesignCandidate {
            id: "sample".to_string(),
            topology,
            flow_rate_m3_s: 1.667e-6,
            inlet_gauge_pa: 100_000.0,
            throat_diameter_m: 50e-6,
            inlet_diameter_m: 4e-3,
            throat_length_m: 250e-6,
            channel_width_m: 6e-3,
            channel_height_m: 1e-3,
            serpentine_segments: 4,
            segment_length_m: 0.01,
            bend_radius_m: 0.004,
            feed_hematocrit: 0.45,
            trifurcation_center_frac: 1.0 / 3.0,
            cif_pretri_center_frac: 1.0 / 3.0,
            cif_terminal_tri_center_frac: 1.0 / 3.0,
            cif_terminal_bi_treat_frac: 0.68,
            asymmetric_narrow_frac: 0.5,
            trifurcation_left_frac: 1.0 / 3.0,
            cross_section_shape: CrossSectionShape::Rectangular,
        }
    }

    #[test]
    fn channel_system_has_plate_edge_inlet_and_outlet_ports() {
        let candidate = sample_candidate(DesignTopology::TrifurcationSerpentine);
        let system = candidate.to_channel_system();

        let has_left_edge_port = system.nodes.iter().any(|n| n.point.0.abs() < 1e-9);
        let has_right_edge_port = system
            .nodes
            .iter()
            .any(|n| (n.point.0 - PLATE_WIDTH_MM).abs() < 1e-9);

        assert!(has_left_edge_port, "missing left plate-edge inlet trunk");
        assert!(has_right_edge_port, "missing right plate-edge outlet trunk");
    }

    #[test]
    fn channel_system_plate_box_matches_full_plate_dimensions() {
        let candidate = sample_candidate(DesignTopology::BifurcationVenturi);
        let system = candidate.to_channel_system();
        assert!((system.box_dims.0 - PLATE_WIDTH_MM).abs() < 1e-9);
        assert!((system.box_dims.1 - PLATE_HEIGHT_MM).abs() < 1e-9);
        assert_eq!(system.box_outline.len(), 4);
    }

    #[test]
    fn channel_geometry_is_remapped_to_full_plate_bounds() {
        use cfd_schematics::geometry::ChannelType;

        let candidate =
            sample_candidate(DesignTopology::IncrementalFiltrationTriBiSeparator { n_pretri: 2 });
        let system = candidate.to_channel_system();

        let mut x_min = f64::INFINITY;
        let mut x_max = f64::NEG_INFINITY;
        let mut y_min = f64::INFINITY;
        let mut y_max = f64::NEG_INFINITY;

        let mut include = |x: f64, y: f64| {
            if x < x_min {
                x_min = x;
            }
            if x > x_max {
                x_max = x;
            }
            if y < y_min {
                y_min = y;
            }
            if y > y_max {
                y_max = y;
            }
        };

        for n in &system.nodes {
            include(n.point.0, n.point.1);
        }
        for ch in &system.channels {
            match &ch.channel_type {
                ChannelType::Straight => {}
                ChannelType::SmoothStraight { path }
                | ChannelType::Serpentine { path }
                | ChannelType::Arc { path }
                | ChannelType::Frustum { path, .. } => {
                    for &(x, y) in path {
                        include(x, y);
                    }
                }
            }
        }

        // X should span full plate (inlet trunk at x=0, outlet trunk at x=PLATE_WIDTH_MM)
        assert!((x_min - 0.0).abs() < 1e-6, "x_min={x_min}");
        assert!((x_max - PLATE_WIDTH_MM).abs() < 1e-6, "x_max={x_max}");
        // Y channels should be within plate bounds but padded away from walls
        // (edge_padding keeps peripheral channels inward)
        assert!(y_min >= 0.0, "y_min={y_min} should be >= 0");
        assert!(y_max <= PLATE_HEIGHT_MM, "y_max={y_max} should be <= {PLATE_HEIGHT_MM}");
        assert!(y_min < PLATE_HEIGHT_MM / 2.0, "y_min={y_min} should be below center");
        assert!(y_max > PLATE_HEIGHT_MM / 2.0, "y_max={y_max} should be above center");
    }
}
