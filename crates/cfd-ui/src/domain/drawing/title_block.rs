//! Title block — engineering drawing identification and revision data (ISO 7200).

use serde::{Deserialize, Serialize};

/// Revision history entry.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RevisionEntry {
    /// Revision letter or number (e.g. "A", "B", "1").
    pub revision: String,
    /// Date of the revision.
    pub date: String,
    /// Description of changes.
    pub description: String,
    /// Who approved the revision.
    pub approved_by: String,
}

/// Title block data for an engineering drawing.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct TitleBlock {
    /// Organization or company name.
    pub company_name: String,
    /// Drawing title / part name.
    pub drawing_title: String,
    /// Drawing number / part number.
    pub drawing_number: String,
    /// Current revision letter.
    pub revision: String,
    /// Author / drawn by.
    pub author: String,
    /// Date created.
    pub date: String,
    /// Drawing scale (e.g. "1:1", "2:1").
    pub scale: String,
    /// Sheet number and total (e.g. 1 of 3).
    pub sheet_number: u32,
    /// Total number of sheets.
    pub sheet_total: u32,
    /// Material specification.
    pub material: Option<String>,
    /// Surface finish specification.
    pub finish: Option<String>,
    /// General tolerances note.
    pub tolerances: Option<String>,
    /// Approval signatures (role, name) pairs.
    pub approvals: Vec<(String, String)>,
    /// Revision history.
    pub revisions: Vec<RevisionEntry>,
}

impl Default for TitleBlock {
    fn default() -> Self {
        Self {
            company_name: String::new(),
            drawing_title: "Untitled Drawing".to_owned(),
            drawing_number: "DWG-001".to_owned(),
            revision: "A".to_owned(),
            author: String::new(),
            date: String::new(),
            scale: "1:1".to_owned(),
            sheet_number: 1,
            sheet_total: 1,
            material: None,
            finish: None,
            tolerances: None,
            approvals: Vec::new(),
            revisions: Vec::new(),
        }
    }
}
