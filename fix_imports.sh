#!/bin/bash

# Fix combined imports with Result
sed -i 's/use cfd_core::{Error, Result};/use cfd_core::error::{Error, Result};/g' $(find crates -name "*.rs")
sed -i 's/use cfd_core::{Result, Error};/use cfd_core::error::{Error, Result};/g' $(find crates -name "*.rs")

# Fix Fluid imports in combinations
sed -i 's/use cfd_core::{Fluid, Result};/use cfd_core::error::Result;\nuse cfd_core::fluid::Fluid;/g' $(find crates -name "*.rs")
sed -i 's/use cfd_core::{Result, Fluid};/use cfd_core::error::Result;\nuse cfd_core::fluid::Fluid;/g' $(find crates -name "*.rs")

# Fix Problem imports
sed -i 's/use cfd_core::{Problem, Result};/use cfd_core::error::Result;\nuse cfd_core::problem::Problem;/g' $(find crates -name "*.rs")

# Fix BoundaryCondition imports
sed -i 's/use cfd_core::{BoundaryCondition, Result};/use cfd_core::error::Result;\nuse cfd_core::boundary::BoundaryCondition;/g' $(find crates -name "*.rs")
sed -i 's/use cfd_core::{BoundaryCondition, Fluid, Result};/use cfd_core::error::Result;\nuse cfd_core::boundary::BoundaryCondition;\nuse cfd_core::fluid::Fluid;/g' $(find crates -name "*.rs")

# Fix other Problem combinations
sed -i 's/use cfd_core::{BoundaryCondition, Fluid, Problem};/use cfd_core::boundary::BoundaryCondition;\nuse cfd_core::fluid::Fluid;\nuse cfd_core::problem::Problem;/g' $(find crates -name "*.rs")

# Fix Solver trait imports
sed -i 's/use cfd_core::Solver;/use cfd_core::solver::Solver;/g' $(find crates -name "*.rs")
sed -i 's/use cfd_core::Configurable;/use cfd_core::solver::Configurable;/g' $(find crates -name "*.rs")
sed -i 's/use cfd_core::Validatable;/use cfd_core::solver::Validatable;/g' $(find crates -name "*.rs")

echo "Import fixes applied"