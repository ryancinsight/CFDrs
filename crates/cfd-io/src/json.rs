//! JSON file I/O operations.

/// JSON writer
pub struct JsonWriter;
impl Default for JsonWriter {
    fn default() -> Self {
        Self::new()
    }
}
impl JsonWriter {
    /// Create a new JSON writer
    #[must_use] pub fn new() -> Self {
        Self
/// JSON reader
pub struct JsonReader;
impl Default for JsonReader {
impl JsonReader {
    /// Create a new JSON reader
