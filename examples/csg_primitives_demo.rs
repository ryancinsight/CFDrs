//! CSG demonstration using csgrs
//!
//! This example demonstrates basic usage of the csgrs library
//! for creating and manipulating 3D geometries.

use std::fs;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🔧 CSG Demonstration");
    println!("===================");
    
    // Create output directory
    fs::create_dir_all("output/csg_primitives")?;
    
    // Create basic shapes using csgrs
    println!("Creating basic shapes...");
    
    // Create a cube
    let cube = csgrs::Cube::new(2.0, 2.0, 2.0);
    let cube_stl = cube.to_stl_ascii("cube");
    fs::write("output/csg_primitives/cube.stl", cube_stl)?;
    println!("✓ Cube created and exported");
    
    // Create a sphere
    let sphere = csgrs::Sphere::new(1.0, 32, 16);
    let sphere_stl = sphere.to_stl_ascii("sphere");
    fs::write("output/csg_primitives/sphere.stl", sphere_stl)?;
    println!("✓ Sphere created and exported");
    
    // Create a cylinder
    let cylinder = csgrs::Cylinder::new(0.8, 3.0, 24);
    let cylinder_stl = cylinder.to_stl_ascii("cylinder");
    fs::write("output/simple_csg/cylinder.stl", cylinder_stl)?;
    println!("✓ Cylinder created and exported");
    
    // Perform boolean operations
    println!("\nPerforming boolean operations...");
    
    // Create two cubes for boolean operations
    let cube1 = csgrs::Cube::new(2.0, 2.0, 2.0);
    let mut cube2 = csgrs::Cube::new(1.5, 1.5, 1.5);
    cube2.translate(0.5, 0.5, 0.5);
    
    // Union
    let union_result = cube1.union(&cube2);
    let union_stl = union_result.to_stl_ascii("union");
    fs::write("output/simple_csg/union.stl", union_stl)?;
    println!("✓ Union operation completed");
    
    // Difference
    let cube3 = csgrs::Cube::new(2.0, 2.0, 2.0);
    let sphere_hole = csgrs::Sphere::new(1.0, 24, 12);
    let difference_result = cube3.difference(&sphere_hole);
    let diff_stl = difference_result.to_stl_ascii("cube_with_hole");
    fs::write("output/simple_csg/cube_with_hole.stl", diff_stl)?;
    println!("✓ Difference operation completed");
    
    // Intersection
    let cube4 = csgrs::Cube::new(2.0, 2.0, 2.0);
    let sphere2 = csgrs::Sphere::new(1.5, 24, 12);
    let intersection_result = cube4.intersection(&sphere2);
    let intersect_stl = intersection_result.to_stl_ascii("intersection");
    fs::write("output/simple_csg/intersection.stl", intersect_stl)?;
    println!("✓ Intersection operation completed");
    
    // Transformations
    println!("\nApplying transformations...");
    
    let mut transformed_cube = csgrs::Cube::new(1.0, 1.0, 1.0);
    
    // Translation
    transformed_cube.translate(2.0, 0.0, 0.0);
    let translated_stl = transformed_cube.to_stl_ascii("translated_cube");
    fs::write("output/simple_csg/translated_cube.stl", translated_stl)?;
    println!("✓ Translation applied");
    
    // Rotation
    transformed_cube.rotate(0.0, 0.0, 1.0, std::f64::consts::PI / 4.0);
    let rotated_stl = transformed_cube.to_stl_ascii("rotated_cube");
    fs::write("output/simple_csg/rotated_cube.stl", rotated_stl)?;
    println!("✓ Rotation applied");
    
    // Scaling
    transformed_cube.scale(1.5, 1.0, 0.8);
    let scaled_stl = transformed_cube.to_stl_ascii("scaled_cube");
    fs::write("output/simple_csg/scaled_cube.stl", scaled_stl)?;
    println!("✓ Scaling applied");
    
    // Create a complex CFD-relevant geometry
    println!("\nCreating CFD-relevant geometry...");
    
    // Flow domain with obstacle
    let flow_domain = csgrs::Cube::new(10.0, 4.0, 1.0);
    let obstacle = csgrs::Cylinder::new(0.5, 1.2, 32);
    let flow_geometry = flow_domain.difference(&obstacle);
    let flow_stl = flow_geometry.to_stl_ascii("flow_around_cylinder");
    fs::write("output/simple_csg/flow_around_cylinder.stl", flow_stl)?;
    println!("✓ Flow domain with cylindrical obstacle created");
    
    // Pipe with elbow
    let straight_pipe = csgrs::Cylinder::new(0.5, 4.0, 24);
    let mut elbow_pipe = csgrs::Cylinder::new(0.5, 3.0, 24);
    elbow_pipe.rotate(0.0, 1.0, 0.0, std::f64::consts::PI / 2.0);
    elbow_pipe.translate(1.5, 0.0, 2.0);
    
    let pipe_system = straight_pipe.union(&elbow_pipe);
    let pipe_stl = pipe_system.to_stl_ascii("pipe_elbow");
    fs::write("output/simple_csg/pipe_elbow.stl", pipe_stl)?;
    println!("✓ Pipe with elbow created");
    
    // Heat exchanger fin array
    let base_plate = csgrs::Cube::new(3.0, 0.1, 2.0);
    let mut fin_assembly = base_plate;
    
    for i in 0..5 {
        let mut fin = csgrs::Cube::new(0.05, 1.0, 1.8);
        fin.translate(-1.0 + i as f64 * 0.5, 0.5, 0.0);
        fin_assembly = fin_assembly.union(&fin);
    }
    
    let fin_stl = fin_assembly.to_stl_ascii("heat_exchanger");
    fs::write("output/simple_csg/heat_exchanger.stl", fin_stl)?;
    println!("✓ Heat exchanger fin array created");
    
    println!("\n✅ All CSG operations completed successfully!");
    println!("📁 STL files saved to: output/simple_csg/");
    println!("\nFiles created:");
    println!("  - Primitive shapes: cube.stl, sphere.stl, cylinder.stl");
    println!("  - Boolean operations: union.stl, cube_with_hole.stl, intersection.stl");
    println!("  - Transformations: translated_cube.stl, rotated_cube.stl, scaled_cube.stl");
    println!("  - CFD geometries: flow_around_cylinder.stl, pipe_elbow.stl, heat_exchanger.stl");
    
    Ok(())
}