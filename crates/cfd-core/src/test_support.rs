//! Assertions that name the constraint a call violated.
//!
//! An assertion that only checks a result is an error passes whenever the call
//! fails, including for a reason unrelated to the one the test exists to
//! catch. Asserting the rendered error instead makes a test fail when the code
//! rejects its input for the wrong reason.
//!
//! [`Error`](crate::error::Error) renders as `category: detail`, so a fragment
//! spanning the colon pins the variant as well as the cause.

use core::fmt::Display;

/// Assert that `result` is an error whose rendered form contains `fragment`.
///
/// Generic over the error type so it serves the crates that carry their own,
/// and taking the result by reference so a test that inspects the value
/// afterwards keeps its binding.
///
/// # Panics
///
/// Panics when `result` is `Ok`, or when its error does not contain
/// `fragment`.
#[track_caller]
pub fn assert_rejects<T, E: Display>(result: &Result<T, E>, fragment: &str) {
    match result {
        Err(error) => {
            let rendered = format!("{error}");
            assert!(
                rendered.contains(fragment),
                "rejection must name the violated constraint: expected a \
                 message containing {fragment:?}, got {rendered:?}"
            );
        }
        Ok(_) => panic!("expected a rejection containing {fragment:?}, got Ok"),
    }
}
