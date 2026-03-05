Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$workspaceRoot = Split-Path -Parent $PSScriptRoot

$files = @(
    "report/milestone12_results.md",
    "report/ARPA-H_SonALAsense_Milestone 12 Report.md",
    "report/milestone12_scenario_audit.md",
    "report/milestone12_leukapheresis_infeasibility.md",
    "report/low_flow_band_analysis.md",
    "report/wall_shear_diffuser_analysis.md",
    "report/figures/selected_ga_schematic.svg",
    "report/figures/selected_cifx_combined_schematic.svg",
    "report/figures/selected_cif_schematic.svg",
    "report/figures/top_hydrosdt_schematic.svg",
    "report/figures/m12_cross_mode_scoring.svg",
    "report/figures/m12_head_to_head_scores.svg",
    "report/figures/m12_cavitation_distribution.svg",
    "report/figures/m12_pareto_oncology.svg",
    "report/figures/m12_multifidelity_dp.svg",
    "report/figures/m12_ga_convergence.svg"
)

$dirs = @(
    "report/milestone12",
    "crates/cfd-optim/outputs/sdt_therapy",
    "crates/cfd-optim/outputs/sdt_scenario_audit",
    "crates/cfd-optim/outputs/milestone12_analysis"
)

Write-Host "Cleaning Milestone 12 generated artifacts..."

foreach ($relative in $files) {
    $path = Join-Path $workspaceRoot $relative
    if (Test-Path $path) {
        Remove-Item -LiteralPath $path -Force
        Write-Host "  removed file: $relative"
    }
}

foreach ($relative in $dirs) {
    $path = Join-Path $workspaceRoot $relative
    if (Test-Path $path) {
        Remove-Item -LiteralPath $path -Recurse -Force
        Write-Host "  removed dir : $relative"
    }
}

$canonicalTargets = @(
    "report/milestone12_results.md",
    "report/ARPA-H_SonALAsense_Milestone 12 Report.md",
    "report/milestone12/figure_manifest.json",
    "report/milestone12",
    "report/figures/selected_ga_schematic.svg",
    "report/figures/selected_cifx_combined_schematic.svg",
    "report/figures/selected_cif_schematic.svg",
    "report/figures/top_hydrosdt_schematic.svg",
    "report/figures/m12_cross_mode_scoring.svg",
    "report/figures/m12_head_to_head_scores.svg",
    "report/figures/m12_cavitation_distribution.svg",
    "report/figures/m12_pareto_oncology.svg",
    "report/figures/m12_multifidelity_dp.svg",
    "report/figures/m12_ga_convergence.svg"
)

$remaining = @()
foreach ($relative in $canonicalTargets) {
    $path = Join-Path $workspaceRoot $relative
    if (Test-Path $path) {
        $remaining += $relative
    }
}

if ($remaining.Count -eq 0) {
    Write-Host "Post-clean verification: 0 canonical Milestone 12 artifacts remain."
} else {
    Write-Host "Post-clean verification: found remaining canonical artifacts:"
    foreach ($relative in $remaining) {
        Write-Host "  remains    : $relative"
    }
    throw "Milestone 12 cleanup verification failed."
}

Write-Host "Milestone 12 cleanup complete. Regeneration can proceed."
