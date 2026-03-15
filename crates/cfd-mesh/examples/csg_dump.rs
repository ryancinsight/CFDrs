use cfd_mesh::application::csg::boolean::{csg_boolean, BooleanOp};
use cfd_mesh::domain::core::scalar::Point3r;
use cfd_mesh::domain::geometry::primitives::PrimitiveMesh;
use cfd_mesh::Cube;

fn main() {
    let cube_a = Cube {
        origin: Point3r::new(0.0, 0.0, 0.0),
        width: 2.0,
        height: 2.0,
        depth: 2.0,
    }
    .build()
    .unwrap();
    let cube_b = Cube {
        origin: Point3r::new(1.0, 0.0, 0.0),
        width: 2.0,
        height: 2.0,
        depth: 2.0,
    }
    .build()
    .unwrap();

    let mesh = csg_boolean(BooleanOp::Union, &cube_a, &cube_b).unwrap();

    let mut fractional = 0;
    for face in mesh.faces.iter() {
        for &vid in &face.vertices {
            let p = mesh.vertices.position(vid);
            if p.x.fract() != 0.0 || p.y.fract() != 0.0 || p.z.fract() != 0.0 {
                println!("Fractional Coordinate found: {:?}", p);
                fractional += 1;
            }
        }
    }
    println!("Total fractional vertices: {}", fractional);
}
