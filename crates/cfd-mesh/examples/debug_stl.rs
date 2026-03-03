use cfd_mesh::infrastructure::io::stl::read_stl;
use std::fs::File;
use std::path::Path;

fn main() {
    let solid_path = Path::new("crates/cfd-mesh/outputs/schematic_to_3d/mirrored_bifurcation/mirrored_bifurcation_solid.stl");
    let chan_path = Path::new("crates/cfd-mesh/outputs/schematic_to_3d/mirrored_bifurcation/mirrored_bifurcation_channels.stl");

    let solid = read_stl(File::open(solid_path).unwrap()).unwrap();
    let chans = read_stl(File::open(chan_path).unwrap()).unwrap();

    let s_aabb = solid.bounding_box();
    let c_aabb = chans.bounding_box();

    println!("Solid AABB: {:?}", s_aabb);
    println!("Chans AABB: {:?}", c_aabb);

    // Look for vertices near X=0
    let mut zero_verts = 0;
    for pos in solid.vertices.positions() {
        if pos.x.abs() < 1e-3 {
            zero_verts += 1;
        }
    }
    println!("Solid vertices exactly at X=0: {}", zero_verts);

    // Faces touching X=0
    let mut faces_at_0 = 0;
    for face in solid.faces.iter() {
        let p0 = solid.vertices.position(face.vertices[0]);
        let p1 = solid.vertices.position(face.vertices[1]);
        let p2 = solid.vertices.position(face.vertices[2]);
        if p0.x.abs() < 1e-3 || p1.x.abs() < 1e-3 || p2.x.abs() < 1e-3 {
            faces_at_0 += 1;
        }
    }
    println!("Solid faces touching X=0: {}", faces_at_0);
}
