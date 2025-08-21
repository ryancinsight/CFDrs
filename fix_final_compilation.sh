#!/bin/bash

echo "Applying final compilation fixes..."

# Fix plugin.rs - clone name before move
sed -i '317s/self.resolver.add_plugin(name, deps)?;/self.resolver.add_plugin(name.clone(), deps)?;/' crates/cfd-core/src/plugin.rs

# Fix material_properties.rs - add Copy derive
sed -i '77s/#\[derive(Debug, Clone, Serialize, Deserialize)\]/#[derive(Debug, Clone, Copy, Serialize, Deserialize)]/' crates/cfd-core/src/domains/material_properties.rs

# Fix time.rs Crank-Nicolson - proper cloning
sed -i '274,276s/let mut y_new = state;/let mut y_new = state.clone();/' crates/cfd-core/src/time.rs
sed -i '274,276s/let y_old = state;/let y_old = state.clone();/' crates/cfd-core/src/time.rs

# Fix factory.rs - remove extra reference
sed -i '149s/self.factory.create(&typed_config)/self.factory.create(typed_config)/' crates/cfd-core/src/factory.rs

echo "Fixes applied. Building..."
source /usr/local/cargo/env
cargo build --all 2>&1 | tail -30