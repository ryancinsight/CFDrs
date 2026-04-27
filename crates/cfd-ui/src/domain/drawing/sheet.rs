//! Drawing sheet — container for an engineering drawing with title block and views.

use crate::domain::drawing::annotation::Annotation;
use crate::domain::drawing::dimension::DimensionSpec;
use crate::domain::drawing::title_block::TitleBlock;
use crate::domain::drawing::view::ProjectedView;
use serde::{Deserialize, Serialize};

/// Standard sheet sizes for engineering drawings.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub enum SheetSize {
    // ISO sizes
    A0,
    A1,
    A2,
    A3,
    A4,
    // ANSI sizes
    AnsiA,
    AnsiB,
    AnsiC,
    AnsiD,
    AnsiE,
    /// Custom size with explicit dimensions.
    Custom {
        /// Width in millimeters.
        width_mm: f64,
        /// Height in millimeters.
        height_mm: f64,
    },
}

impl SheetSize {
    /// Return the sheet dimensions in millimeters as (width, height).
    #[must_use]
    pub fn dimensions_mm(self) -> (f64, f64) {
        match self {
            Self::A0 => (1189.0, 841.0),
            Self::A1 => (841.0, 594.0),
            Self::A2 => (594.0, 420.0),
            Self::A3 => (420.0, 297.0),
            Self::A4 => (297.0, 210.0),
            Self::AnsiA => (279.4, 215.9),
            Self::AnsiB => (431.8, 279.4),
            Self::AnsiC => (558.8, 431.8),
            Self::AnsiD => (863.6, 558.8),
            Self::AnsiE => (1117.6, 863.6),
            Self::Custom {
                width_mm,
                height_mm,
            } => (width_mm, height_mm),
        }
    }
}

/// A complete engineering drawing sheet.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct DrawingSheet {
    /// The physical sheet size.
    pub size: SheetSize,
    /// Title block information.
    pub title_block: TitleBlock,
    /// The projected views placed on this sheet.
    pub views: Vec<ProjectedView>,
    /// Dimensions placed on the drawing.
    pub dimensions: Vec<DimensionSpec>,
    /// Text annotations and notes.
    pub annotations: Vec<Annotation>,
    /// Border margin from sheet edge in millimeters.
    pub border_margin_mm: f64,
}

impl DrawingSheet {
    /// Create a new drawing sheet with the given size and default title block.
    #[must_use]
    pub fn new(size: SheetSize) -> Self {
        Self {
            size,
            title_block: TitleBlock::default(),
            views: Vec::new(),
            dimensions: Vec::new(),
            annotations: Vec::new(),
            border_margin_mm: 10.0,
        }
    }

    /// The drawable area dimensions (sheet minus margins) in mm.
    #[must_use]
    pub fn drawable_area_mm(&self) -> (f64, f64) {
        let (w, h) = self.size.dimensions_mm();
        (
            w - 2.0 * self.border_margin_mm,
            h - 2.0 * self.border_margin_mm,
        )
    }
}
