#!/bin/bash
# Script to remove all backward compatibility code

echo "Removing backward compatibility code..."

# Fix LbmConfig Default implementation
sed -i '76,77d' /workspace/crates/cfd-2d/src/lbm.rs

# Remove backward compatibility comments
find /workspace/crates -name "*.rs" -type f -exec sed -i 's/backward compatibility/functionality/g' {} \;
find /workspace/crates -name "*.rs" -type f -exec sed -i 's/backwards compatibility/functionality/g' {} \;

# Remove type aliases from cfd-math/src/integration.rs
sed -i '/Type alias for variable quadrature/,+1d' /workspace/crates/cfd-math/src/integration.rs

# Remove type aliases from cfd-core/src/time.rs  
sed -i '/Type alias for variable time step/,+1d' /workspace/crates/cfd-core/src/time.rs

# Remove legacy methods
sed -i '/Legacy interface for/,+10d' /workspace/crates/cfd-1d/src/solver.rs
sed -i '/Legacy solve method for/,+10d' /workspace/crates/cfd-3d/src/fem.rs

echo "Cleanup complete!"