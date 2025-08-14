//! Test the csgrs API to understand its structure

fn main() {
    println!("Testing csgrs API...");
    
    // Test basic shape creation
    let cube = csgrs::Cube::new(1.0, 1.0, 1.0);
    println!("Created cube: {:?}", std::any::type_name_of_val(&cube));
    
    // Test STL export
    let stl = cube.to_stl_ascii("test");
    println!("STL length: {}", stl.len());
    
    // Test sphere
    let sphere = csgrs::Sphere::new(1.0, 16, 8);
    println!("Created sphere: {:?}", std::any::type_name_of_val(&sphere));
    
    // Test boolean operation
    let union_result = cube.union(&sphere);
    println!("Union result: {:?}", std::any::type_name_of_val(&union_result));
    
    println!("csgrs API test completed successfully!");
}