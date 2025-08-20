#!/bin/bash
# Fix ALL import issues systematically

echo "=== Fixing all imports systematically ==="

# Fix Error and Result imports
find crates -name "*.rs" -type f -exec sed -i \
  -e 's/use cfd_core::{Error, Result}/use cfd_core::error::{Error, Result}/g' \
  -e 's/use cfd_core::Error/use cfd_core::error::Error/g' \
  -e 's/use cfd_core::Result/use cfd_core::error::Result/g' \
  -e 's/cfd_core::Error\b/cfd_core::error::Error/g' \
  -e 's/cfd_core::Result\b/cfd_core::error::Result/g' {} \;

# Fix Fluid imports
find crates -name "*.rs" -type f -exec sed -i \
  -e 's/use cfd_core::Fluid/use cfd_core::fluid::Fluid/g' \
  -e 's/use cfd_core::{Fluid/use cfd_core::fluid::Fluid; use cfd_core::error::{/g' \
  -e 's/cfd_core::Fluid\b/cfd_core::fluid::Fluid/g' {} \;

# Fix BoundaryCondition imports
find crates -name "*.rs" -type f -exec sed -i \
  -e 's/use cfd_core::BoundaryCondition/use cfd_core::boundary::BoundaryCondition/g' \
  -e 's/cfd_core::BoundaryCondition\b/cfd_core::boundary::BoundaryCondition/g' {} \;

# Fix Problem imports
find crates -name "*.rs" -type f -exec sed -i \
  -e 's/use cfd_core::Problem/use cfd_core::problem::Problem/g' \
  -e 's/cfd_core::Problem\b/cfd_core::problem::Problem/g' {} \;

# Fix Solver imports
find crates -name "*.rs" -type f -exec sed -i \
  -e 's/use cfd_core::Solver\b/use cfd_core::solver::Solver/g' \
  -e 's/cfd_core::Solver\b/cfd_core::solver::Solver/g' {} \;

# Fix Domain imports
find crates -name "*.rs" -type f -exec sed -i \
  -e 's/use cfd_core::Domain\b/use cfd_core::domain::Domain/g' \
  -e 's/cfd_core::Domain\b/cfd_core::domain::Domain/g' {} \;

# Fix SolverConfiguration imports
find crates -name "*.rs" -type f -exec sed -i \
  -e 's/use cfd_core::SolverConfiguration/use cfd_core::solver::SolverConfiguration/g' \
  -e 's/cfd_core::SolverConfiguration\b/cfd_core::solver::SolverConfiguration/g' {} \;

# Fix other common imports
find crates -name "*.rs" -type f -exec sed -i \
  -e 's/use cfd_core::Configurable/use cfd_core::solver::Configurable/g' \
  -e 's/use cfd_core::Validatable/use cfd_core::solver::Validatable/g' \
  -e 's/cfd_core::Configurable\b/cfd_core::solver::Configurable/g' \
  -e 's/cfd_core::Validatable\b/cfd_core::solver::Validatable/g' {} \;

echo "Import fixes applied"