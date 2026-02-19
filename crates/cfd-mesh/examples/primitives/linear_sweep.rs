//! Export: Linear Sweep primitives → `outputs/primitives/linear_sweep_*.stl`
//!
//! Demonstrates three profiles swept along +Y using [`cfd_mesh::LinearSweep`]:
//!
//! 1. **Square prism** — 2×2 mm profile, 3 mm tall (V = 12 mm³, exact)
//! 2. **Hexagonal prism** — r=1 mm regular hexagon, 2 mm tall (V ≈ 5.196 mm³)
//! 3. **Triangular prism** — equilateral-ish triangle, 2 mm tall
//!
//! All meshes are validated for watertightness (closed + oriented, χ=2)
//! and positive signed volume.
//!
//! Run with:
//! ```sh
//! cargo run -p cfd-mesh --example primitives_linear_sweep
//! ```

use std::f64::consts::PI;
use std::fs;
use std::io::BufWriter;

use cfd_mesh::LinearSweep;
use cfd_mesh::geometry::primitives::{PrimitiveMesh, linear_sweep::Point2};
use cfd_mesh::io::stl;
use cfd_mesh::storage::edge_store::EdgeStore;
use cfd_mesh::watertight::check::check_watertight;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  Primitive: Linear Sweep  (via LinearSweep builder)");
    println!("=================================================================");

    let crate_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let out_dir = crate_dir.join("outputs").join("primitives");
    fs::create_dir_all(&out_dir)?;

    // ── 1. Square prism ────────────────────────────────────────────────────
    {
        let profile = vec![
            Point2::new(-1.0, -1.0),
            Point2::new( 1.0, -1.0),
            Point2::new( 1.0,  1.0),
            Point2::new(-1.0,  1.0),
        ];
        let sweep = LinearSweep { profile, height: 3.0 };
        let mesh = sweep.build()?;

        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        let expected_vol = 4.0 * 3.0; // 2×2 × 3

        println!("  [Square Prism 2×2×3]");
        println!("    Vertices  : {}", mesh.vertices.len());
        println!("    Faces     : {}", mesh.faces.len());
        println!("    Closed    : {}", report.is_closed);
        println!("    Oriented  : {}", report.orientation_consistent);
        println!("    Watertight: {}", report.is_watertight);
        println!("    Volume    : {:.6} mm³  (expected {:.6})", report.signed_volume, expected_vol);
        println!("    Euler χ   : {:?}  (expected 2)", report.euler_characteristic);

        let stl_path = out_dir.join("linear_sweep_square.stl");
        let file = fs::File::create(&stl_path)?;
        let mut writer = BufWriter::new(file);
        stl::write_binary_stl(&mut writer, &mesh.vertices, &mesh.faces)?;
        println!("    STL       : {}", stl_path.display());
        println!();
    }

    // ── 2. Hexagonal prism ────────────────────────────────────────────────
    {
        let profile = LinearSweep::regular_polygon(6, 1.0);
        let height = 2.0_f64;
        let sweep = LinearSweep { profile, height };
        let mesh = sweep.build()?;

        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        // Hexagon area: 3√3/2 × r²
        let hex_area = 3.0 * (3.0_f64).sqrt() / 2.0;
        let expected_vol = hex_area * height;

        println!("  [Hexagonal Prism r=1, h=2]");
        println!("    Vertices  : {}", mesh.vertices.len());
        println!("    Faces     : {}", mesh.faces.len());
        println!("    Closed    : {}", report.is_closed);
        println!("    Oriented  : {}", report.orientation_consistent);
        println!("    Watertight: {}", report.is_watertight);
        println!("    Volume    : {:.6} mm³  (expected {:.6})", report.signed_volume, expected_vol);
        println!("    Euler χ   : {:?}  (expected 2)", report.euler_characteristic);

        let stl_path = out_dir.join("linear_sweep_hexagon.stl");
        let file = fs::File::create(&stl_path)?;
        let mut writer = BufWriter::new(file);
        stl::write_binary_stl(&mut writer, &mesh.vertices, &mesh.faces)?;
        println!("    STL       : {}", stl_path.display());
        println!();
    }

    // ── 3. Triangular prism ───────────────────────────────────────────────
    {
        // Equilateral triangle: side ≈ 2, so r (circumradius) ≈ 1.155
        // Place vertices at angles 90°, 210°, 330° (= 90 + 0, 120, 240)
        let r = 1.0_f64;
        let profile = vec![
            Point2::new(r * (PI / 2.0).cos(),       r * (PI / 2.0).sin()),
            Point2::new(r * (PI / 2.0 + 2.0 * PI / 3.0).cos(), r * (PI / 2.0 + 2.0 * PI / 3.0).sin()),
            Point2::new(r * (PI / 2.0 + 4.0 * PI / 3.0).cos(), r * (PI / 2.0 + 4.0 * PI / 3.0).sin()),
        ];
        let height = 2.0_f64;
        // Area of equilateral triangle with circumradius r: A = 3√3/4 × side²
        // side = r√3, so A = 3√3/4 × 3r² = 9√3/4 × r² ... or use signed_area
        let area = LinearSweep::signed_area(&profile);
        let expected_vol = area * height;
        let sweep = LinearSweep { profile, height };
        let mesh = sweep.build()?;

        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);

        println!("  [Triangular Prism (equilateral r=1, h=2)]");
        println!("    Vertices  : {}", mesh.vertices.len());
        println!("    Faces     : {}", mesh.faces.len());
        println!("    Closed    : {}", report.is_closed);
        println!("    Oriented  : {}", report.orientation_consistent);
        println!("    Watertight: {}", report.is_watertight);
        println!("    Volume    : {:.6} mm³  (expected {:.6})", report.signed_volume, expected_vol);
        println!("    Euler χ   : {:?}  (expected 2)", report.euler_characteristic);

        let stl_path = out_dir.join("linear_sweep_triangle.stl");
        let file = fs::File::create(&stl_path)?;
        let mut writer = BufWriter::new(file);
        stl::write_binary_stl(&mut writer, &mesh.vertices, &mesh.faces)?;
        println!("    STL       : {}", stl_path.display());
        println!();
    }

    println!("=================================================================");
    Ok(())
}
