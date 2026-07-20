# Example: venturi_cavitation

**Crate**: `cfd-suite` (workspace root)  
**Run**: `cargo run --example venturi_cavitation`  
**Source**: [`examples/venturi_cavitation.rs`](../../examples/venturi_cavitation.rs)

## What This Example Demonstrates

Selective cavitation screening for a venturi flow with biologically-realistic
cell populations — healthy RBCs, CTCs (circulating tumour cells), and WBCs —
showing differential inception thresholds driven by membrane mechanics.

| Population | Radius | Stiffness | Nucleation weight |
|---|---|---|---|
| Healthy RBC | 3.3 μm | 72 kPa | 1.0 |
| CTC | 8.0 μm | 28 kPa | 1.0 |
| Healthy WBC | varies | varies | 1.0 |

## Key Code Snippet

```rust
use cfd_1d::{evaluate_venturi_screening, VenturiScreeningInput};
use cfd_core::physics::cavitation::{
    SelectiveCavitationInput, SelectiveCavitationPopulation, CellMechanicalState,
};

let selective = SelectiveCavitationInput {
    base_vapor_pressure_pa: 3_170.0,
    density_kg_m3:          1_025.0,
    populations: vec![
        SelectiveCavitationPopulation { identity: CellPopulationIdentity::HealthyRbc, .. },
        SelectiveCavitationPopulation { identity: CellPopulationIdentity::CirculatingTumorCell, .. },
    ],
};
let result = evaluate_venturi_screening(VenturiScreeningInput {
    selective: Some(selective), ..
})?;
```

## Physics Background

**Selective cavitation** exploits the difference in membrane stiffness and size
between healthy and malignant cells. A CTC (stiffer, larger nucleus) has a lower
cavitation inception number than a healthy RBC. A well-tuned venturi geometry can
therefore selectively lyse CTCs while leaving erythrocytes intact.

## Book Chapter

[← Turbulence Models and Cavitation](../turbulence_multiphase.md)

