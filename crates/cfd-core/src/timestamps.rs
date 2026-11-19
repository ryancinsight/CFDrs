//! UTC timestamp rendering for the suite's metadata and report surfaces.
//!
//! Every surface renders the same two shapes: an RFC 3339 instant for the
//! artifacts the suite writes, and a `YYYY-MM-DD HH:MM:SS UTC` line for the
//! report headers a reader sees. The calendar fields are written directly
//! rather than through a formatter, because the formatter's result is
//! fallible and every call site would carry that error path for an instant
//! that is already a valid date and time.

use time::OffsetDateTime;

/// The calendar fields a UTC instant renders, in the order the formats use.
fn calendar_fields(now: OffsetDateTime) -> (i32, u8, u8, u8, u8, u8) {
    (
        now.year(),
        u8::from(now.month()),
        now.day(),
        now.hour(),
        now.minute(),
        now.second(),
    )
}

/// Renders a UTC instant as an RFC 3339 timestamp, such as
/// `2026-11-19T12:34:56Z`.
///
/// # Example
///
/// ```rust
/// use cfd_core::timestamps::rfc3339_utc;
/// use time::OffsetDateTime;
///
/// let instant = OffsetDateTime::from_unix_timestamp(1_700_000_000).unwrap();
/// assert_eq!(rfc3339_utc(instant), "2023-11-14T22:13:20Z");
/// ```
#[must_use]
pub fn rfc3339_utc(now: OffsetDateTime) -> String {
    let (year, month, day, hour, minute, second) = calendar_fields(now);
    format!("{year:04}-{month:02}-{day:02}T{hour:02}:{minute:02}:{second:02}Z")
}

/// Renders a UTC instant for a report header, such as
/// `2026-11-19 12:34:56 UTC`.
///
/// # Example
///
/// ```rust
/// use cfd_core::timestamps::display_utc;
/// use time::OffsetDateTime;
///
/// let instant = OffsetDateTime::from_unix_timestamp(1_700_000_000).unwrap();
/// assert_eq!(display_utc(instant), "2023-11-14 22:13:20 UTC");
/// ```
#[must_use]
pub fn display_utc(now: OffsetDateTime) -> String {
    let (year, month, day, hour, minute, second) = calendar_fields(now);
    format!("{year:04}-{month:02}-{day:02} {hour:02}:{minute:02}:{second:02} UTC")
}

#[cfg(test)]
mod tests {
    use super::{display_utc, rfc3339_utc};
    use time::OffsetDateTime;

    #[test]
    fn rfc3339_renders_the_utc_calendar_fields() {
        let instant = OffsetDateTime::from_unix_timestamp(1_700_000_000)
            .expect("a unix timestamp inside the calendar range");
        assert_eq!(rfc3339_utc(instant), "2023-11-14T22:13:20Z");
    }

    #[test]
    fn display_renders_the_utc_calendar_fields() {
        let instant = OffsetDateTime::from_unix_timestamp(1_700_000_000)
            .expect("a unix timestamp inside the calendar range");
        assert_eq!(display_utc(instant), "2023-11-14 22:13:20 UTC");
    }
}
