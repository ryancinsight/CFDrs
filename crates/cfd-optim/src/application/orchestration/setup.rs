//! Environment setup and configuration for Milestone 12 example pipelines.

use std::path::PathBuf;

use crate::scoring::OptimMode;

/// Initialise the `tracing_subscriber` with INFO default and env-filter.
///
/// Safe to call multiple times — subsequent calls are no-ops.
pub fn init_tracing() {
    let _ = tracing_subscriber::fmt()
        .with_env_filter(
            tracing_subscriber::EnvFilter::builder()
                .with_default_directive(tracing::Level::INFO.into())
                .from_env_lossy(),
        )
        .try_init();
}

/// Resolve workspace root and create `report/milestone12/` + `report/figures/`.
///
/// Returns `(workspace_root, out_dir, figures_dir)`.
pub fn resolve_output_directories(
) -> Result<(PathBuf, PathBuf, PathBuf), Box<dyn std::error::Error>> {
    let workspace_root = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("cfd-optim crate has a parent")
        .parent()
        .expect("crates/ has a workspace root")
        .to_path_buf();

    let report_root = workspace_root.join("report");
    let out_dir = report_root.join("milestone12");
    let figures_dir = report_root.join("figures");
    std::fs::create_dir_all(&out_dir)?;
    std::fs::create_dir_all(&figures_dir)?;

    Ok((workspace_root, out_dir, figures_dir))
}

/// Canonical Milestone 12 reporting must run in release mode because the full
/// audit/report pipeline materializes long-running evidence products.
///
/// # Errors
/// Returns an error when invoked from a debug build.
pub fn ensure_release_reports() -> Result<(), Box<dyn std::error::Error>> {
    if cfg!(debug_assertions) {
        return Err(
            "Milestone 12 report orchestration is release-only; rerun with --release".into(),
        );
    }
    Ok(())
}

/// Check the `M12_FAST` environment variable.
///
/// Fast mode is the default; set `M12_FAST=0|false|no` to force a full run.
pub fn fast_mode() -> bool {
    std::env::var("M12_FAST")
        .ok()
        .map(|v| {
            let s = v.trim().to_ascii_lowercase();
            !(s == "0" || s == "false" || s == "no")
        })
        .unwrap_or(true)
}

/// Parse an environment variable as `usize`, falling back to `default` (min 1).
pub fn fast_env(name: &str, default: usize) -> usize {
    std::env::var(name)
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(default)
        .max(1)
}

/// Number of ranked full-design artifacts to retain per Milestone 12 stage.
///
/// Defaults to 1500 so report generation and downstream analysis can compare a
/// materially larger ranked pool than the human-facing top-5 tables.
pub fn milestone12_ranked_pool_size() -> usize {
    std::env::var("M12_RANKED_POOL_SIZE")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(1500)
        .max(1)
}

/// Standard Option 2 scoring mode: adult oncology patient (70 kg), 50/50
/// CombinedSdtLeukapheresis.
pub fn option2_mode() -> OptimMode {
    OptimMode::CombinedSdtLeukapheresis {
        leuka_weight: 0.5,
        sdt_weight: 0.5,
        patient_weight_kg: 70.0,
    }
}

#[cfg(test)]
mod tests {
    use super::ensure_release_reports;

    #[test]
    fn milestone12_reports_are_debug_guarded() {
        if cfg!(debug_assertions) {
            assert!(ensure_release_reports().is_err());
        } else {
            assert!(ensure_release_reports().is_ok());
        }
    }
}
