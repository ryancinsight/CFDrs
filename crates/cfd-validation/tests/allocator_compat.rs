//! Proves a consumer can select its own process allocator.

use std::alloc::System;

#[global_allocator]
static CONSUMER_ALLOCATOR: System = System;

#[test]
fn consumer_allocator_coexists_with_validation_library() {
    let stats = cfd_validation::benchmarking::MemoryStats::new();
    assert_eq!(stats.snapshot().allocation_count, 0);

    let value = Box::new(42_u64);
    assert_eq!(*value, 42);
}
