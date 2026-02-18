// Fuzz target: STL parser into HalfEdgeMesh.
//
// Run with:
//   cargo +nightly fuzz run fuzz_read_stl_he -- -max_len=65536
//
// The parser must never panic on arbitrary input; it may return `Err`.
#![no_main]

use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    use cfd_mesh::mesh::with_mesh;
    use cfd_mesh::io::stl::read_stl_he;

    let _ = with_mesh(|mut mesh, mut token| {
        let _ = read_stl_he(data, &mut mesh, &mut token);
        mesh.face_count()
    });
});
