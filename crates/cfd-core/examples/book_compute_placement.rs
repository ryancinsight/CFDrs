//! Book example: selecting compute backends and carrying placement hints.
//!
//! This keeps backend selection and NUMA placement at the compute seam, while
//! downstream domain layers remain backend-agnostic.

use cfd_core::compute::{ComputeBackend, ComputeCapability};

fn main() {
    let capability = ComputeCapability::detect();

    let workloads = [256usize, 8_192, 250_000];
    let selected_backends: Vec<ComputeBackend> = workloads
        .into_iter()
        .map(|problem_size| capability.select_backend(problem_size))
        .collect();

    let placement_hint = capability.placement_hint();

    let _ = (selected_backends, placement_hint);
}
