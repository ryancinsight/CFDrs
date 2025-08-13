#!/bin/bash
# Comprehensive verification script for CFD Suite

echo "========================================="
echo "CFD Suite - Complete Verification Script"
echo "========================================="
echo ""

# Set up Rust environment
source "/usr/local/cargo/env" 2>/dev/null || true

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Function to check command result
check_result() {
    if [ $1 -eq 0 ]; then
        echo -e "${GREEN}✓ $2 PASSED${NC}"
    else
        echo -e "${RED}✗ $2 FAILED${NC}"
        exit 1
    fi
}

# 1. Check Rust installation
echo "1. Checking Rust installation..."
rustc --version > /dev/null 2>&1
check_result $? "Rust installation"

# 2. Build all crates
echo ""
echo "2. Building all crates..."
cargo build --all --release 2>&1 | tail -3
check_result ${PIPESTATUS[0]} "Build all crates"

# 3. Run all tests
echo ""
echo "3. Running all tests..."
TEST_OUTPUT=$(cargo test --all --release 2>&1)
TEST_RESULT=$?
echo "$TEST_OUTPUT" | grep "test result:" | head -5
check_result $TEST_RESULT "All tests"

# 4. Count test statistics
TOTAL_TESTS=$(echo "$TEST_OUTPUT" | grep "passed" | awk '{sum+=$1} END {print sum}')
echo "   Total tests passed: $TOTAL_TESTS"

# 5. Build examples
echo ""
echo "4. Building all examples..."
cargo build --examples --release 2>&1 | tail -2
check_result ${PIPESTATUS[0]} "Build examples"

# 6. Run sample examples
echo ""
echo "5. Running sample examples..."
echo "   - Running pipe_flow_1d..."
cargo run --release --example pipe_flow_1d 2>&1 | grep "Converged:" > /dev/null
check_result $? "pipe_flow_1d example"

echo "   - Running 2d_heat_diffusion..."
cargo run --release --example 2d_heat_diffusion 2>&1 | grep "Example Complete" > /dev/null
check_result $? "2d_heat_diffusion example"

# 7. Check for documentation
echo ""
echo "6. Checking documentation generation..."
cargo doc --no-deps --quiet 2>&1
check_result $? "Documentation generation"

# 8. Summary
echo ""
echo "========================================="
echo "         VERIFICATION COMPLETE"
echo "========================================="
echo ""
echo "Summary:"
echo "  - Rust: $(rustc --version | cut -d' ' -f2)"
echo "  - Cargo: $(cargo --version | cut -d' ' -f2)"
echo "  - Total tests passed: $TOTAL_TESTS"
echo "  - All crates: BUILD SUCCESS"
echo "  - All tests: PASS"
echo "  - All examples: WORKING"
echo "  - Documentation: GENERATED"
echo ""
echo -e "${GREEN}✓ CFD Suite is fully operational!${NC}"