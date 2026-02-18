use super::{EdgeId, NodeId};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum NodeKind {
    Inlet,
    Outlet,
    Reservoir,
    Junction,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NodeSpec {
    pub id: NodeId,
    pub kind: NodeKind,
}

impl NodeSpec {
    #[must_use]
    pub fn new(id: impl Into<String>, kind: NodeKind) -> Self {
        Self {
            id: NodeId::new(id),
            kind,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum EdgeKind {
    Pipe,
    Valve,
    Pump,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChannelSpec {
    pub id: EdgeId,
    pub kind: EdgeKind,
    pub from: NodeId,
    pub to: NodeId,
    pub length_m: f64,
    pub diameter_m: f64,
    pub resistance: f64,
    pub quad_coeff: f64,
    // Component properties
    pub valve_cv: Option<f64>,
    pub pump_max_flow: Option<f64>,
    pub pump_max_pressure: Option<f64>,
}

impl ChannelSpec {
    #[must_use]
    pub fn new_pipe(
        id: impl Into<String>,
        from: impl Into<String>,
        to: impl Into<String>,
        length_m: f64,
        diameter_m: f64,
        resistance: f64,
        quad_coeff: f64,
    ) -> Self {
        Self {
            id: EdgeId::new(id),
            kind: EdgeKind::Pipe,
            from: NodeId::new(from),
            to: NodeId::new(to),
            length_m,
            diameter_m,
            resistance,
            quad_coeff,
            valve_cv: None,
            pump_max_flow: None,
            pump_max_pressure: None,
        }
    }

    #[must_use]
    pub fn new_valve(
        id: impl Into<String>,
        from: impl Into<String>,
        to: impl Into<String>,
        cv: f64,
    ) -> Self {
        Self {
            id: EdgeId::new(id),
            kind: EdgeKind::Valve,
            from: NodeId::new(from),
            to: NodeId::new(to),
            length_m: 0.0, // Valves are 0-length in 1D usually, or negligible
            diameter_m: 0.0,
            resistance: 0.0, // Calculated from Cv
            quad_coeff: 0.0, // Calculated from Cv
            valve_cv: Some(cv),
            pump_max_flow: None,
            pump_max_pressure: None,
        }
    }

    #[must_use]
    pub fn new_pump(
        id: impl Into<String>,
        from: impl Into<String>,
        to: impl Into<String>,
        max_flow: f64,
        max_pressure: f64,
    ) -> Self {
        Self {
            id: EdgeId::new(id),
            kind: EdgeKind::Pump,
            from: NodeId::new(from),
            to: NodeId::new(to),
            length_m: 0.0,
            diameter_m: 0.0,
            resistance: 0.0,
            quad_coeff: 0.0,
            valve_cv: None,
            pump_max_flow: Some(max_flow),
            pump_max_pressure: Some(max_pressure),
        }
    }
}
