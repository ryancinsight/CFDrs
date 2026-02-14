//! Error types for the scheme↔cfd-1d bridge.

use std::fmt;

/// Errors that can occur during scheme→cfd-1d conversion.
#[derive(Debug)]
pub enum BridgeError {
    /// The scheme `ChannelSystem` has no channels.
    EmptyNetwork,
    /// A channel references a node id that does not exist.
    InvalidNodeReference {
        /// Channel id in the scheme system
        channel_id: usize,
        /// The missing node id
        node_id: usize,
    },
    /// The converter could not infer any inlet nodes.
    NoInlets,
    /// The converter could not infer any outlet nodes.
    NoOutlets,
    /// A channel has zero or negative computed path length.
    ZeroLengthChannel {
        /// Channel id in the scheme system
        channel_id: usize,
    },
    /// An error from the cfd-core layer.
    CfdError(cfd_core::error::Error),
}

impl fmt::Display for BridgeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyNetwork => write!(f, "Scheme ChannelSystem contains no channels"),
            Self::InvalidNodeReference {
                channel_id,
                node_id,
            } => write!(
                f,
                "Channel {} references non-existent node {}",
                channel_id, node_id
            ),
            Self::NoInlets => write!(
                f,
                "Could not infer any inlet nodes from the scheme topology"
            ),
            Self::NoOutlets => write!(
                f,
                "Could not infer any outlet nodes from the scheme topology"
            ),
            Self::ZeroLengthChannel { channel_id } => {
                write!(f, "Channel {} has zero or negative path length", channel_id)
            }
            Self::CfdError(e) => write!(f, "CFD error during conversion: {}", e),
        }
    }
}

impl std::error::Error for BridgeError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::CfdError(e) => Some(e),
            _ => None,
        }
    }
}

impl From<cfd_core::error::Error> for BridgeError {
    fn from(e: cfd_core::error::Error) -> Self {
        Self::CfdError(e)
    }
}
