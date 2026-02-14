use super::{EdgeId, NodeId};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NodeKind {
    Inlet,
    Outlet,
    Junction,
}

#[derive(Debug, Clone)]
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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EdgeKind {
    Pipe,
}

#[derive(Debug, Clone)]
pub struct ChannelSpec {
    pub id: EdgeId,
    pub kind: EdgeKind,
    pub from: NodeId,
    pub to: NodeId,
    pub length_m: f64,
    pub diameter_m: f64,
    pub resistance: f64,
    pub quad_coeff: f64,
}

impl ChannelSpec {
    #[must_use]
    pub fn new(
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
        }
    }
}
