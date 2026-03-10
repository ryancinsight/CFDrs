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

    let out_dir = workspace_root.join("report").join("milestone12");
    let figures_dir = workspace_root.join("report").join("figures");
    std::fs::create_dir_all(&out_dir)?;
    std::fs::create_dir_all(&figures_dir)?;

    Ok((workspace_root, out_dir, figures_dir))
}

/// Check the `M12_FAST` environment variable.
pub fn fast_mode() -> bool {
    std::env::var("M12_FAST")
        .ok()
        .map(|v| {
            let s = v.trim().to_ascii_lowercase();
            s == "1" || s == "true" || s == "yes"
        })
        .unwrap_or(false)
}

/// Parse an environment variable as `usize`, falling back to `default` (min 1).
pub fn fast_env(name: &str, default: usize) -> usize {
    std::env::var(name)
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(default)
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
