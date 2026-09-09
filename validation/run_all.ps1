
$tests = @(
    "validation/validate_bifurcation.py",
    "validation/validate_trifurcation.py",
    "validation/validate_venturi.py",
    "validation/validate_serpentine.py"
)

foreach ($test in $tests) {
    Write-Host "Running $test..."
    python $test
    if ($LASTEXITCODE -ne 0) {
        Write-Error "$test failed!"
        exit 1
    }
}
Write-Host "All validation tests passed!"
