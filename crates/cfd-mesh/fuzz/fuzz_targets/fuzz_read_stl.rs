// Fuzz target: STL parser (IndexedMesh path).
//
// Run with:
//   cargo +nightly fuzz run fuzz_read_stl -- -max_len=65536
//
// The parser must never panic on arbitrary input; it may return `Err`.
#![no_main]

use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    // All errors are acceptable; panics are not.
    let _ = cfd_mesh::io::stl::fuzz_read_stl(data);
});
