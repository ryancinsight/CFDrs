$examples = @(
    "csg_cylinder_cylinder",
    "csg_cylinder_cylinder_asymmetric",
    "csg_cylinder_cylinder_l_shape",
    "csg_cylinder_cylinder_perpendicular",
    "csg_cylinder_cylinder_symmetric",
    "csg_cylinder_cylinder_t_junction",
    "csg_cylinder_cylinder_v_shape",
    "csg_coplanar_disk_disk"
)

$results = @()
foreach ($ex in $examples) {
    $results += "===== $ex ====="
    $out = & cargo run --features csg -p cfd-mesh --example $ex 2>&1
    foreach ($line in $out) {
        $s = $line.ToString()
        if ($s -match "Watertight|PASS|FAIL|Vol error|Volume\s+:|CSG|panick|thread.*panicked") {
            $results += $s
        }
    }
    $results += ""
}
$results | Set-Content .\all_csg_results.txt -Encoding UTF8
