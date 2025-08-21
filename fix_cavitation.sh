#!/bin/bash
# Fix the cavitation.rs file properly

cat > /tmp/cavitation_fix.sed << 'EOF'
# Fix line 246
s/let r_b = bubble_radius;/let r_b = *bubble_radius;/

# Fix line 250
s/let rate = f_vap \*/let rate = *f_vap */

# Fix line 257
s/let rate = f_cond \*/let rate = *f_cond */

# Fix line 251
s/nucleation_fraction \* (T::one()/*nucleation_fraction * (T::one()/
EOF

sed -i -f /tmp/cavitation_fix.sed /workspace/crates/cfd-core/src/cavitation.rs