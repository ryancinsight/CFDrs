//! Resistance analysis for network components

use super::traits::NetworkAnalyzer;
use crate::domain::network::{Network, NetworkGraphExt};
use crate::solver::analysis::ResistanceAnalysis;
use aequitas::systems::si::quantities::HydraulicResistance;
use cfd_core::CfdScalar;
use cfd_core::conversion::{SafeFromF64, SafeFromUsize};
use cfd_core::error::ResistanceCalculationErrorKind as ResistanceCalculationError;
use cfd_core::error::Result;
use cfd_core::physics::constants::physics::thermo::{P_ATM, T_STANDARD};
use eunomia::NumericElement;
use petgraph::Direction::Outgoing;
use petgraph::graph::EdgeIndex;
use std::iter::Sum;

/// Resistance analyzer for network components
pub struct ResistanceAnalyzer<T: CfdScalar + Copy> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: CfdScalar + Copy> Default for ResistanceAnalyzer<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: CfdScalar + Copy> ResistanceAnalyzer<T> {
    /// Create new resistance analyzer
    #[must_use]
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T: CfdScalar + Copy + SafeFromF64 + SafeFromUsize + Sum> NetworkAnalyzer<T>
    for ResistanceAnalyzer<T>
{
    type Result = ResistanceAnalysis<T>;

    fn analyze(&mut self, network: &Network<T>) -> Result<ResistanceAnalysis<T>> {
        let mut analysis = ResistanceAnalysis::new();
        self.populate_edge_resistances(network, &mut analysis)?;

        let critical = self.critical_paths(network, &analysis);
        for path in critical.paths() {
            analysis.add_critical_path(
                path.iter()
                    .map(|&edge_idx| network.graph[edge_idx].id.clone())
                    .collect(),
            );
        }

        Ok(analysis)
    }

    fn name(&self) -> &'static str {
        "ResistanceAnalyzer"
    }
}

impl<T: CfdScalar + Copy + SafeFromF64 + Sum> ResistanceAnalyzer<T> {
    fn populate_edge_resistances(
        &self,
        network: &Network<T>,
        analysis: &mut ResistanceAnalysis<T>,
    ) -> Result<()> {
        let fluid = network.fluid();

        for edge in network.edges_with_properties() {
            let flow_rate = if edge.flow_rate.into_base() == T::ZERO {
                None
            } else {
                Some(edge.flow_rate.into_base())
            };

            let resistance = self
                .calculate_resistance(edge.properties, fluid, flow_rate)
                .map_err(|e| {
                    cfd_core::error::Error::InvalidInput(format!(
                        "Failed to analyze resistance for edge '{}': {}",
                        edge.id, e
                    ))
                })?;

            analysis.add_resistance(edge.id.clone(), HydraulicResistance::from_base(resistance));

            let component_type = edge.properties.component_type;
            analysis.add_resistance_by_type(
                component_type.as_str().to_string(),
                HydraulicResistance::from_base(resistance),
            );
        }

        Ok(())
    }

    /// Enumerate every simple inlet-to-outlet edge path by depth-first search
    /// over outgoing edges and keep those tied at the maximal series
    /// resistance.
    ///
    /// Parallel edges are distinct branches, so each edge path is visited
    /// exactly once. The running sum `prefix` adds resistances from the inlet
    /// forward, the same left-to-right order as summing each path afresh.
    fn critical_paths(
        &self,
        network: &Network<T>,
        analysis: &ResistanceAnalysis<T>,
    ) -> CriticalEdgePaths<T> {
        let mut critical = CriticalEdgePaths::default();
        let inlet_nodes = network.graph.inlet_nodes();
        let outlet_nodes = network.graph.outlet_nodes();
        if inlet_nodes.is_empty() || outlet_nodes.is_empty() {
            return critical;
        }

        let graph = &network.graph;
        let edge_resistance: Vec<T> = graph
            .edge_indices()
            .map(|edge_idx| {
                let edge = &graph[edge_idx];
                analysis
                    .resistances
                    .get(&edge.id)
                    .copied()
                    .unwrap_or(edge.resistance)
                    .into_base()
            })
            .collect();
        let edge_target = |edge_idx: EdgeIndex| graph.raw_edges()[edge_idx.index()].target();

        let mut on_path = vec![false; graph.node_count()];
        let mut path: Vec<EdgeIndex> = Vec::new();
        let mut prefix: Vec<T> = Vec::new();
        let mut cursors: Vec<Option<EdgeIndex>> = Vec::new();

        for &inlet in &inlet_nodes {
            for &outlet in &outlet_nodes {
                on_path.fill(false);
                on_path[inlet.index()] = true;
                path.clear();
                prefix.clear();
                prefix.push(T::ZERO);
                cursors.clear();
                cursors.push(graph.first_edge(inlet, Outgoing));

                while let Some(cursor) = cursors.last_mut() {
                    let Some(edge_idx) = *cursor else {
                        cursors.pop();
                        if let Some(entered) = path.pop() {
                            prefix.pop();
                            on_path[edge_target(entered).index()] = false;
                        }
                        continue;
                    };
                    *cursor = graph.next_edge(edge_idx, Outgoing);

                    let child = edge_target(edge_idx);
                    let resistance_sum = prefix[path.len()] + edge_resistance[edge_idx.index()];
                    if child == outlet {
                        critical.offer(resistance_sum, &path, edge_idx);
                    } else if !on_path[child.index()] {
                        on_path[child.index()] = true;
                        path.push(edge_idx);
                        prefix.push(resistance_sum);
                        cursors.push(graph.first_edge(child, Outgoing));
                    }
                }
            }
        }

        critical
    }
}

/// Edge paths tied at the maximal series resistance, stored contiguously:
/// path `k` is `edges[offsets[k]..offsets[k + 1]]`.
struct CriticalEdgePaths<T> {
    best: Option<T>,
    edges: Vec<EdgeIndex>,
    offsets: Vec<usize>,
}

impl<T> Default for CriticalEdgePaths<T> {
    fn default() -> Self {
        Self {
            best: None,
            edges: Vec::new(),
            offsets: vec![0],
        }
    }
}

impl<T: CfdScalar + Copy> CriticalEdgePaths<T> {
    /// Offer the path `prefix` followed by `last`, of series resistance
    /// `resistance_sum`: a larger sum replaces the tied set, and a sum within
    /// `epsilon * (|best| + 1)` of the best joins it.
    fn offer(&mut self, resistance_sum: T, prefix: &[EdgeIndex], last: EdgeIndex) {
        let replaces = match self.best {
            None => true,
            Some(best) if resistance_sum > best => true,
            Some(best) => {
                let diff = <T as NumericElement>::abs(resistance_sum - best);
                let scale = <T as NumericElement>::abs(best) + T::ONE;
                if diff <= T::default_epsilon() * scale {
                    false
                } else {
                    return;
                }
            }
        };
        if replaces {
            self.best = Some(resistance_sum);
            self.edges.clear();
            self.offsets.truncate(1);
        }
        self.edges.extend_from_slice(prefix);
        self.edges.push(last);
        self.offsets.push(self.edges.len());
    }

    fn paths(&self) -> impl Iterator<Item = &[EdgeIndex]> {
        self.offsets
            .windows(2)
            .map(|bounds| &self.edges[bounds[0]..bounds[1]])
    }
}

impl<T: CfdScalar + Copy + SafeFromF64> ResistanceAnalyzer<T> {
    fn calculate_resistance(
        &self,
        properties: &crate::domain::network::EdgeProperties<T>,
        fluid: &cfd_core::physics::fluid::ConstantPropertyFluid<T>,
        flow_rate: Option<T>,
    ) -> std::result::Result<T, ResistanceCalculationError> {
        use crate::physics::resistance::{FlowConditions, HagenPoiseuilleModel, ResistanceModel};

        // Require hydraulic diameter - no silent fallbacks
        let hydraulic_diameter = properties
            .hydraulic_diameter
            .ok_or(ResistanceCalculationError::MissingHydraulicDiameter)?
            .into_base();

        // Create resistance model with validated parameters
        let model = HagenPoiseuilleModel::new(hydraulic_diameter, properties.length.into_base());

        const REF_TEMPERATURE_KEY: &str = "reference_temperature";
        const REF_PRESSURE_KEY: &str = "reference_pressure";

        let temperature = properties
            .properties
            .get(REF_TEMPERATURE_KEY)
            .copied()
            .unwrap_or_else(|| T::from_f64_or_one(T_STANDARD));
        let pressure = properties
            .properties
            .get(REF_PRESSURE_KEY)
            .copied()
            .unwrap_or_else(|| T::from_f64_or_one(P_ATM));

        let conditions = FlowConditions {
            reynolds_number: flow_rate.map(|q| {
                let velocity = q / properties.area.into_base();
                fluid.density.into_base() * velocity * hydraulic_diameter
                    / fluid.viscosity.into_base()
            }),
            velocity: flow_rate.map(|q| q / properties.area.into_base()),
            flow_rate,
            shear_rate: None,
            temperature,
            pressure,
        };

        // Calculate resistance and propagate any errors
        model
            .calculate_resistance(fluid, &conditions)
            .map_err(|e| ResistanceCalculationError::ModelError(e.to_string()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::network::{
        ComponentType, EdgeProperties, Network, NetworkBuilder, ResistanceUpdatePolicy,
    };
    use aequitas::systems::si::quantities::{Area, Length};
    use cfd_core::physics::fluid::database::water_20c;
    use std::collections::HashMap;

    fn pipe_properties(id: &str, diameter: f64) -> EdgeProperties<f64> {
        EdgeProperties {
            id: id.to_string(),
            component_type: ComponentType::Pipe,
            length: Length::from_base(0.01),
            area: Area::from_base(std::f64::consts::FRAC_PI_4 * diameter * diameter),
            hydraulic_diameter: Some(Length::from_base(diameter)),
            resistance: HydraulicResistance::from_base(1.0),
            geometry: None,
            resistance_update_policy: ResistanceUpdatePolicy::FlowInvariant,
            properties: HashMap::new(),
        }
    }

    /// inlet =(a1|a2)=> junction =(b1|b2)=> outlet, plus a direct pipe `d`.
    /// `a1` and `a2` are identical, `b2` is narrower than `b1`, and `d` is
    /// wide, so the two maximal paths tie through `b2`. Outgoing edges are
    /// visited newest first, so `a2` precedes `a1`; each path appears once.
    fn parallel_hop_network() -> Network<f64> {
        let mut builder = NetworkBuilder::new();
        let inlet = builder.add_inlet("inlet".to_string());
        let junction = builder.add_junction("junction".to_string());
        let outlet = builder.add_outlet("outlet".to_string());
        let pipes = [
            (inlet, junction, "a1", 1.0e-3),
            (inlet, junction, "a2", 1.0e-3),
            (junction, outlet, "b1", 1.0e-3),
            (junction, outlet, "b2", 0.5e-3),
            (inlet, outlet, "d", 2.0e-3),
        ];
        let edges: Vec<_> = pipes
            .iter()
            .map(|&(from, to, id, _)| builder.connect_with_pipe(from, to, id.to_string()))
            .collect();
        let graph = builder.build().expect("parallel-hop network must validate");
        let mut network = Network::new(graph, water_20c::<f64>().expect("water properties exist"));
        for (edge, &(_, _, id, diameter)) in edges.into_iter().zip(&pipes) {
            network.add_edge_properties(edge, pipe_properties(id, diameter));
        }
        network
    }

    #[test]
    fn critical_paths_are_the_tied_maximal_series_paths_once_each() -> Result<()> {
        let analysis = ResistanceAnalyzer::<f64>::new().analyze(&parallel_hop_network())?;

        assert_eq!(
            analysis.critical_paths,
            vec![
                vec!["a2".to_string(), "b2".to_string()],
                vec!["a1".to_string(), "b2".to_string()],
            ]
        );
        Ok(())
    }

    #[test]
    fn unreachable_outlet_has_no_critical_path() -> Result<()> {
        let mut builder = NetworkBuilder::new();
        let inlet = builder.add_inlet("inlet".to_string());
        let outlet = builder.add_outlet("outlet".to_string());
        let edge = builder.connect_with_pipe(outlet, inlet, "a".to_string());
        let mut network = Network::new(
            builder
                .build()
                .expect("reversed two-node network must validate"),
            water_20c::<f64>().expect("water properties exist"),
        );
        network.add_edge_properties(edge, pipe_properties("a", 1.0e-3));

        let analysis = ResistanceAnalyzer::<f64>::new().analyze(&network)?;

        assert!(analysis.critical_paths.is_empty());
        assert_eq!(analysis.resistances.len(), 1);
        Ok(())
    }
}
