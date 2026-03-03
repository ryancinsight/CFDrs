//! Property editor panel — context-sensitive property display.

/// Property value types for the property editor.
#[derive(Clone, Debug)]
pub enum PropertyValue {
    Text(String),
    Float(f64),
    Integer(i64),
    Boolean(bool),
}

/// A single property entry in the editor.
#[derive(Clone, Debug)]
pub struct PropertyEntry {
    /// Property name/label.
    pub name: String,
    /// Current value.
    pub value: PropertyValue,
    /// Whether this property is read-only.
    pub read_only: bool,
}

/// Collected properties for the currently selected object.
#[derive(Clone, Debug, Default)]
pub struct PropertySheet {
    /// The name of the selected object.
    pub object_name: String,
    /// Properties of the selected object.
    pub properties: Vec<PropertyEntry>,
}

impl PropertySheet {
    /// Create an empty property sheet (nothing selected).
    #[must_use]
    pub fn empty() -> Self {
        Self::default()
    }

    /// Create a property sheet for a mesh object.
    #[must_use]
    pub fn for_mesh(name: &str, vertex_count: usize, face_count: usize) -> Self {
        Self {
            object_name: name.to_owned(),
            properties: vec![
                PropertyEntry {
                    name: "Vertices".to_owned(),
                    value: PropertyValue::Integer(vertex_count as i64),
                    read_only: true,
                },
                PropertyEntry {
                    name: "Faces".to_owned(),
                    value: PropertyValue::Integer(face_count as i64),
                    read_only: true,
                },
            ],
        }
    }
}
