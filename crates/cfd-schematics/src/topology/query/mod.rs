//! Derived topology queries computed from [`BlueprintTopologySpec`] structure.
//!
//! These methods replace the per-variant `match` dispatch formerly in
//! `cfd-optim::DesignTopology`. Every query is derived from the declarative
//! spec, not from an enum variant name, enabling the GA to compose arbitrary
//! topologies without extending an enum.

mod classification;
mod distribution;
mod naming;
mod venturi;

#[cfg(test)]
mod tests;
