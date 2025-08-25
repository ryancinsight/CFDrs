//! CSG demonstration using cfd_mesh CSG wrapper
//!
//! This example demonstrates basic usage of the CFD mesh CSG wrapper
//! for creating and manipulating 3D geometries.

#![allow(missing_docs)]

use cfd_mesh::csg::{CsgBuilder, CsgError, CsgOperator};
use nalgebra::Vector3;
use std::fs;

fn main() -> Result<(), CsgError> {
    println!("🔧 CSG Primitives Demonstration");
    println!("==============================");

    // Create output directory
    fs::create_dir_all("output/csg_primitives")
        .map_err(|e| CsgError::ExportError(format!("Failed to create output directory: {}", e)))?;

    // Create CSG operator
    let operator = CsgOperator::<f64>::new();

    // Create basic shapes
    println!("\nCreating basic shapes...");

    // Create a cube
    let cube = operator.create_cube(2.0, 2.0, 2.0)?;
    let cube_stl = cube.to_stl("cube")?;
    fs::write("output/csg_primitives/cube.stl", cube_stl)
        .map_err(|e| CsgError::ExportError(format!("Failed to write cube STL: {}", e)))?;
    println!("✓ Cube created and exported");

    // Create a sphere
    let sphere = operator.create_sphere(1.0, 32, 16)?;
    let sphere_stl = sphere.to_stl("sphere")?;
    fs::write("output/csg_primitives/sphere.stl", sphere_stl)
        .map_err(|e| CsgError::ExportError(format!("Failed to write sphere STL: {}", e)))?;
    println!("✓ Sphere created and exported");

    // Create a cylinder
    let cylinder = operator.create_cylinder(0.8, 3.0, 24)?;
    let cylinder_stl = cylinder.to_stl("cylinder")?;
    fs::write("output/csg_primitives/cylinder.stl", cylinder_stl)
        .map_err(|e| CsgError::ExportError(format!("Failed to write cylinder STL: {}", e)))?;
    println!("✓ Cylinder created and exported");

    // Create a frustum (cone)
    let cone = operator.create_frustum(1.0, 0.2, 2.0, 16)?;
    let cone_stl = cone.to_stl("cone")?;
    fs::write("output/csg_primitives/cone.stl", cone_stl)
        .map_err(|e| CsgError::ExportError(format!("Failed to write cone STL: {}", e)))?;
    println!("✓ Cone created and exported");

    // Perform boolean operations
    println!("\nPerforming boolean operations...");

    // Create two cubes for boolean operations
    let cube1 = operator.create_cube(2.0, 2.0, 2.0)?;
    let mut cube2 = operator.create_cube(1.5, 1.5, 1.5)?;
    cube2.translate(&Vector3::new(0.5, 0.5, 0.5))?;

    // Union
    let union_result = cube1.union(&cube2);
    let union_stl = union_result.to_stl("union")?;
    fs::write("output/csg_primitives/union.stl", union_stl)
        .map_err(|e| CsgError::ExportError(format!("Failed to write union STL: {}", e)))?;
    println!("✓ Union operation completed");

    // Create fresh cubes for difference
    let cube3 = operator.create_cube(2.0, 2.0, 2.0)?;
    let mut cube4 = operator.create_cube(1.5, 1.5, 1.5)?;
    cube4.translate(&Vector3::new(0.5, 0.5, 0.5))?;

    // Difference
    let difference_result = cube3.difference(&cube4);
    let difference_stl = difference_result.to_stl("difference")?;
    fs::write("output/csg_primitives/difference.stl", difference_stl)
        .map_err(|e| CsgError::ExportError(format!("Failed to write difference STL: {}", e)))?;
    println!("✓ Difference operation completed");

    // Create fresh cubes for intersection
    let cube5 = operator.create_cube(2.0, 2.0, 2.0)?;
    let mut cube6 = operator.create_cube(1.5, 1.5, 1.5)?;
    cube6.translate(&Vector3::new(0.5, 0.5, 0.5))?;

    // Intersection
    let intersection_result = cube5.intersection(&cube6);
    let intersection_stl = intersection_result.to_stl("intersection")?;
    fs::write("output/csg_primitives/intersection.stl", intersection_stl)
        .map_err(|e| CsgError::ExportError(format!("Failed to write intersection STL: {}", e)))?;
    println!("✓ Intersection operation completed");

    // Complex example: Create a hollowed sphere
    println!("\nCreating complex geometry...");

    let outer_sphere = operator.create_sphere(1.5, 32, 16)?;
    let inner_sphere = operator.create_sphere(1.2, 32, 16)?;
    let hollow_sphere = outer_sphere.difference(&inner_sphere);

    // Add a window by cutting with a cube
    let mut window_cube = operator.create_cube(1.0, 2.0, 1.0)?;
    window_cube.translate(&Vector3::new(0.8, 0.0, 0.0))?;

    let windowed_sphere = hollow_sphere.difference(&window_cube);
    let complex_stl = windowed_sphere.to_stl("windowed_hollow_sphere")?;
    fs::write(
        "output/csg_primitives/windowed_hollow_sphere.stl",
        complex_stl,
    )
    .map_err(|e| CsgError::ExportError(format!("Failed to write complex STL: {}", e)))?;
    println!("✓ Complex geometry created");

    // Using the builder pattern
    println!("\nUsing builder pattern...");

    let built_geometry = CsgBuilder::<f64>::new()
        .cube(1.0, 1.0, 1.0)?
        .translate(Vector3::new(0.0, 0.0, 1.0))?
        .build()?;

    let builder_stl = built_geometry.to_stl("builder_result")?;
    fs::write("output/csg_primitives/builder_result.stl", builder_stl)
        .map_err(|e| CsgError::ExportError(format!("Failed to write builder STL: {}", e)))?;
    println!("✓ Builder pattern geometry created");

    println!("\n✅ All CSG primitives demonstrations completed successfully!");
    println!("   Check the output/csg_primitives directory for STL files.");

    Ok(())
}
