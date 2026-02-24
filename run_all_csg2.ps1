$examples = @(
    "csg_cube_cube",
    "csg_cube_sphere",
    "csg_cube_cylinder",
    "csg_sphere_sphere",
    "csg_sphere_cylinder",
    "csg_cylinder_cylinder",
    "csg_cylinder_cylinder_asymmetric",
    "csg_cylinder_cylinder_l_shape",
    "csg_cylinder_cylinder_perpendicular",
    "csg_cylinder_cylinder_symmetric",
    "csg_cylinder_cylinder_t_junction",
    "csg_cylinder_cylinder_v_shape",
    "csg_coplanar_disk_disk"
)

# These are the Cargo.toml example names (without csg_ prefix in file)
# Check if union/intersection/difference/compound exist as targets
$extra = @("union", "intersection", "difference", "compound")
# Probably not registered targets -- the available names use csg_ prefix

$results = @()
foreach ($ex in $examples) {
    $results += "===== $ex ====="
    $out = & cargo run --features csg -p cfd-mesh --example $ex 2>&1
    foreach ($line in $out) {
        $s = $line.ToString()
        if ($s -match "Watertight|PASS|FAIL|Vol error|Volume\s+:|Area error|CSG|panick|thread.*panicked|error\[") {
            $results += $s
        }
    }
    $results += ""
}
$results | Set-Content .\all_csg_results2.txt -Encoding UTF8
