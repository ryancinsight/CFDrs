//! Runtime diagnostics control for CSG pipelines.
//!
//! Diagnostics are disabled by default to keep CSG operations and tests
//! deterministic and I/O-light. Enable with:
//! `CFD_MESH_CSG_TRACE=1` (also accepts `true`, `yes`, `on`).

use std::sync::OnceLock;

/// Returns `true` when CSG trace diagnostics are enabled.
#[inline]
pub(crate) fn trace_enabled() -> bool {
    static TRACE_ENABLED: OnceLock<bool> = OnceLock::new();
    *TRACE_ENABLED.get_or_init(|| {
        std::env::var("CFD_MESH_CSG_TRACE")
            .is_ok_and(|v| {
                let s = v.trim();
                s == "1"
                    || s.eq_ignore_ascii_case("true")
                    || s.eq_ignore_ascii_case("yes")
                    || s.eq_ignore_ascii_case("on")
            })
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn trace_disabled_by_default() {
        // Test process env for this key is unset in normal CI/test runs.
        // If it is set externally, this assertion is still coherent with
        // runtime behavior by checking the key first.
        if std::env::var("CFD_MESH_CSG_TRACE").is_err() {
            assert!(!trace_enabled());
        }
    }
}
