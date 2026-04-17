//! Generated asset-review manifest support for Milestone 12 report releases.
//!
//! # Theorem
//! The emitted review manifest is a deterministic scaffold over the current
//! generated asset set: every current figure and the current narrative report
//! appear exactly once in the review manifest, and prior review findings are
//! preserved only for assets with the same canonical path, asset type, and
//! content digest.
//!
//! **Proof sketch**
//! The builder performs a single pass over the current figure list and the
//! narrative report to create one review entry per current asset path. Existing
//! review entries are keyed by their canonical path string and reused only on
//! exact path/type/digest matches; unmatched prior entries are dropped. Because
//! the current asset list is the only source for output entries and each current
//! path is visited once, the result cannot omit, duplicate, or alias assets
//! across runs, and any content change invalidates stale review state.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use chrono::Utc;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

use super::figures::NarrativeFigureSpec;

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "lowercase")]
pub(crate) enum AssetType {
    Svg,
    Png,
    Jpg,
    Pdf,
    Md,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "lowercase")]
pub(crate) enum ReviewStatus {
    Pass,
    Fail,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub(crate) struct AssetReviewRecord {
    pub opened: bool,
    pub status: ReviewStatus,
    pub actual_finding: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reviewed_at: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub(crate) struct AssetReviewEntry {
    pub path: String,
    pub asset_type: AssetType,
    #[serde(default)]
    pub content_sha256: String,
    pub expected_story: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub expected_geometry_or_data: Option<String>,
    pub review: AssetReviewRecord,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub(crate) struct AssetReviewGenerator {
    pub path: String,
    pub kind: String,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub inputs: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub(crate) struct AssetReviewManifest {
    pub generated_at: String,
    pub generator: AssetReviewGenerator,
    pub assets: Vec<AssetReviewEntry>,
}

#[derive(Debug, Clone)]
pub(crate) struct AssetReviewOutcome {
    pub manifest_path: PathBuf,
    pub complete: bool,
}

#[derive(Debug, Clone)]
struct AssetReviewSeed {
    path: String,
    asset_type: AssetType,
    content_sha256: String,
    expected_story: String,
    expected_geometry_or_data: Option<String>,
}

pub(crate) fn write_asset_review_manifest(
    report_dir: &Path,
    narrative_path: &Path,
    figure_specs: &[NarrativeFigureSpec],
) -> Result<AssetReviewOutcome, Box<dyn std::error::Error>> {
    let manifest_path = report_dir
        .join("milestone12")
        .join("asset_review_manifest.json");
    let existing: Option<AssetReviewManifest> = if manifest_path.exists() {
        Some(serde_json::from_str(&std::fs::read_to_string(
            &manifest_path,
        )?)?)
    } else {
        None
    };
    let mut existing_reviews = existing
        .map(|manifest| {
            manifest
                .assets
                .into_iter()
                .map(|mut asset| {
                    asset.path = normalize_manifest_path_string(&asset.path);
                    (asset.path.clone(), asset)
                })
                .collect::<HashMap<_, _>>()
        })
        .unwrap_or_default();

    let mut seeds = figure_specs
        .iter()
        .map(|spec| {
            let path = asset_path_string(report_dir, &spec.path);
            Ok(AssetReviewSeed {
                content_sha256: sha256_file(Path::new(&path))?,
                path,
                asset_type: AssetType::Svg,
                expected_story: spec.caption.clone(),
                expected_geometry_or_data: Some(expected_geometry_or_data(spec)),
            })
        })
        .collect::<Result<Vec<_>, Box<dyn std::error::Error>>>()?;
    seeds.push(AssetReviewSeed {
        content_sha256: sha256_file(narrative_path)?,
        path: normalize_manifest_path(narrative_path),
        asset_type: AssetType::Md,
        expected_story: "Milestone 12 narrative report should accurately summarize the selected designs, cite the generated figures, and avoid placeholder or unavailable-asset language.".to_string(),
        expected_geometry_or_data: Some(
            "Markdown report should contain the canonical Milestone 12 sections, align with the current figure manifest, and match the current report-manifest authority semantics.".to_string(),
        ),
    });

    let mut assets = Vec::with_capacity(seeds.len());
    for seed in seeds {
        let entry = existing_reviews
            .remove(&seed.path)
            .filter(|existing| {
                existing.asset_type == seed.asset_type
                    && existing.content_sha256 == seed.content_sha256
            })
            .map_or_else(
                || scaffold_entry(&seed),
                |mut existing| {
                    existing.content_sha256.clone_from(&seed.content_sha256);
                    existing.expected_story.clone_from(&seed.expected_story);
                    existing
                        .expected_geometry_or_data
                        .clone_from(&seed.expected_geometry_or_data);
                    existing
                },
            );
        assets.push(entry);
    }

    let manifest = AssetReviewManifest {
        generated_at: Utc::now().to_rfc3339(),
        generator: AssetReviewGenerator {
            path: "crates/cfd-optim/examples/milestone12_report.rs".to_string(),
            kind: "report".to_string(),
            inputs: vec![
                narrative_path.display().to_string(),
                report_dir
                    .join("milestone12")
                    .join("figure_manifest.json")
                    .display()
                    .to_string(),
            ],
        },
        assets,
    };

    if let Some(parent) = manifest_path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    std::fs::write(&manifest_path, serde_json::to_string_pretty(&manifest)?)?;

    Ok(AssetReviewOutcome {
        complete: manifest
            .assets
            .iter()
            .all(|asset| asset.review.opened && asset.review.status == ReviewStatus::Pass),
        manifest_path,
    })
}

fn asset_path_string(report_dir: &Path, manifest_path: &str) -> String {
    let normalized = manifest_path.replace('\\', "/");
    if let Some(stripped) = normalized.strip_prefix("../report/") {
        normalize_manifest_path(&report_dir.join(stripped))
    } else {
        normalize_manifest_path(&report_dir.join(normalized))
    }
}

fn scaffold_entry(seed: &AssetReviewSeed) -> AssetReviewEntry {
    AssetReviewEntry {
        path: seed.path.clone(),
        asset_type: seed.asset_type.clone(),
        content_sha256: seed.content_sha256.clone(),
        expected_story: seed.expected_story.clone(),
        expected_geometry_or_data: seed.expected_geometry_or_data.clone(),
        review: AssetReviewRecord {
            opened: false,
            status: ReviewStatus::Fail,
            actual_finding: "Pending manual visual review of the generated asset.".to_string(),
            reviewed_at: None,
        },
    }
}

fn sha256_file(path: &Path) -> Result<String, Box<dyn std::error::Error>> {
    let bytes = std::fs::read(path)?;
    let mut hasher = Sha256::new();
    hasher.update(bytes);
    Ok(format!("{:x}", hasher.finalize()))
}

fn normalize_manifest_path(path: &Path) -> String {
    path.display().to_string().replace('\\', "/")
}

fn normalize_manifest_path_string(path: &str) -> String {
    path.replace('\\', "/")
}

fn expected_geometry_or_data(spec: &NarrativeFigureSpec) -> String {
    match spec.number {
        1 => "Process figure should show the canonical pipeline boxes and arrows from topology generation through report export without missing stages.".to_string(),
        2 => "Bifurcation concept figure should depict a geometry-authored bifurcation-root treatment layout within the report figure frame.".to_string(),
        3 => "Trifurcation concept figure should depict a geometry-authored trifurcation-root treatment layout with selective-routing structure clearly visible.".to_string(),
        4..=6 => "Selected-design schematic should be geometry-authored, include inlet and outlet markers, and include venturi throat markers where the caption claims active throats.".to_string(),
        7 | 8 | 11 => "Bar heights, labels, and ordering should match the ranked design metrics described in the caption.".to_string(),
        9 => "Scatter plot should include Option 2 and GA points with visible σ=0 and σ=1 reference lines and legible point labels.".to_string(),
        10 => "Pareto plot should show the highlighted shortlist inside the broader eligible background cloud with a visible frontier.".to_string(),
        12 => "Convergence plot should show a generation-indexed fitness trajectory with a trailing-plateau or trailing-improvement annotation.".to_string(),
        13 => "Dean-versus-cavitation figure should show paired per-placement physics with Dean magnitude, cavitation labels, and upstream pressure trend.".to_string(),
        14 => "Treatment-lane zoom should show side-by-side local treatment geometry with venturi throat callouts and a delta summary box.".to_string(),
        _ => "Generated asset should visually match its caption and report role.".to_string(),
    }
}

#[cfg(test)]
mod tests {
    use std::path::{Path, PathBuf};

    use super::{write_asset_review_manifest, ReviewStatus};
    use crate::reporting::figures::NarrativeFigureSpec;

    #[test]
    fn asset_review_manifest_scaffolds_all_current_figures() {
        let unique = format!(
            "m12-asset-review-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .expect("clock")
                .as_nanos()
        );
        let report_dir = std::env::temp_dir().join(unique);
        std::fs::create_dir_all(report_dir.join("milestone12")).expect("create review dir");
        write_svg(
            &report_dir,
            "milestone12_creation_optimization_process.svg",
            "<svg>process</svg>",
        );
        write_svg(
            &report_dir,
            "selected_option1_schematic.svg",
            "<svg>option1</svg>",
        );
        let narrative_path = write_narrative(&report_dir, "# Narrative");

        let figure_specs = vec![
            NarrativeFigureSpec {
                number: 1,
                title: "Process".to_string(),
                path: "../report/figures/milestone12_creation_optimization_process.svg".to_string(),
                caption: "Process caption".to_string(),
                alt: "Process".to_string(),
            },
            NarrativeFigureSpec {
                number: 4,
                title: "Option 1".to_string(),
                path: "../report/figures/selected_option1_schematic.svg".to_string(),
                caption: "Option 1 caption".to_string(),
                alt: "Option 1".to_string(),
            },
        ];

        let outcome = write_asset_review_manifest(&report_dir, &narrative_path, &figure_specs)
            .expect("asset review manifest should write");
        let manifest = std::fs::read_to_string(outcome.manifest_path).expect("manifest readable");

        assert!(manifest.contains("selected_option1_schematic.svg"));
        assert!(manifest.contains("milestone12_creation_optimization_process.svg"));
        assert!(manifest.contains("ARPA-H_SonALAsense_Milestone 12 Report.md"));
        assert!(manifest.contains("\"content_sha256\""));
        assert!(manifest.contains("\"opened\": false"));
        assert!(manifest.contains("\"status\": \"fail\""));
        assert!(!outcome.complete);
    }

    #[test]
    fn asset_review_manifest_preserves_matching_passed_review_entries() {
        let unique = format!(
            "m12-asset-review-preserve-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .expect("clock")
                .as_nanos()
        );
        let report_dir = std::env::temp_dir().join(unique);
        let milestone_dir = report_dir.join("milestone12");
        std::fs::create_dir_all(&milestone_dir).expect("create review dir");
        let manifest_path = milestone_dir.join("asset_review_manifest.json");
        let figure_path = write_svg(
            &report_dir,
            "selected_option1_schematic.svg",
            "<svg>option1</svg>",
        );
        let narrative_path = write_narrative(&report_dir, "# Narrative");
        let reviewed_asset_path = figure_path.display().to_string();
        let reviewed_asset_sha = super::sha256_file(&figure_path).expect("hash");
        let narrative_asset_path = narrative_path.display().to_string();
        let narrative_sha = super::sha256_file(&narrative_path).expect("hash");
        let seeded_manifest = serde_json::json!({
            "generated_at": "2026-04-13T00:00:00Z",
            "generator": { "path": "x", "kind": "report" },
            "assets": [
                {
                    "path": reviewed_asset_path,
                    "asset_type": "svg",
                    "content_sha256": reviewed_asset_sha,
                    "expected_story": "old",
                    "review": {
                        "opened": true,
                        "status": "pass",
                        "actual_finding": "Representative.",
                        "reviewed_at": "2026-04-13T01:00:00Z"
                    }
                },
                {
                    "path": narrative_asset_path,
                    "asset_type": "md",
                    "content_sha256": narrative_sha,
                    "expected_story": "Narrative",
                    "review": {
                        "opened": true,
                        "status": "pass",
                        "actual_finding": "Narrative consistent.",
                        "reviewed_at": "2026-04-13T01:00:00Z"
                    }
                }
            ]
        });
        std::fs::write(
            &manifest_path,
            serde_json::to_string_pretty(&seeded_manifest).expect("serialize"),
        )
        .expect("seed manifest");

        let figure_specs = vec![NarrativeFigureSpec {
            number: 4,
            title: "Option 1".to_string(),
            path: "../report/figures/selected_option1_schematic.svg".to_string(),
            caption: "Option 1 caption".to_string(),
            alt: "Option 1".to_string(),
        }];
        let outcome = write_asset_review_manifest(&report_dir, &narrative_path, &figure_specs)
            .expect("asset review manifest should rewrite");
        let manifest: super::AssetReviewManifest =
            serde_json::from_str(&std::fs::read_to_string(outcome.manifest_path).expect("read"))
                .expect("deserialize");

        assert_eq!(manifest.assets.len(), 2);
        let figure_entry = manifest
            .assets
            .iter()
            .find(|entry| entry.asset_type == super::AssetType::Svg)
            .expect("figure entry");
        assert!(figure_entry.review.opened);
        assert_eq!(figure_entry.review.status, ReviewStatus::Pass);
        assert_eq!(figure_entry.review.actual_finding, "Representative.");
        let narrative_entry = manifest
            .assets
            .iter()
            .find(|entry| entry.asset_type == super::AssetType::Md)
            .expect("narrative entry");
        assert!(narrative_entry.review.opened);
        assert_eq!(narrative_entry.review.status, ReviewStatus::Pass);
    }

    #[test]
    fn asset_review_manifest_invalidates_reviews_when_content_changes() {
        let unique = format!(
            "m12-asset-review-invalidate-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .expect("clock")
                .as_nanos()
        );
        let report_dir = std::env::temp_dir().join(unique);
        let milestone_dir = report_dir.join("milestone12");
        std::fs::create_dir_all(&milestone_dir).expect("create review dir");
        let figure_path = write_svg(
            &report_dir,
            "selected_option1_schematic.svg",
            "<svg>option1-v1</svg>",
        );
        let narrative_path = write_narrative(&report_dir, "# Narrative");
        let manifest_path = milestone_dir.join("asset_review_manifest.json");
        let seeded_manifest = serde_json::json!({
            "generated_at": "2026-04-13T00:00:00Z",
            "generator": { "path": "x", "kind": "report" },
            "assets": [
                {
                    "path": figure_path.display().to_string(),
                    "asset_type": "svg",
                    "content_sha256": super::sha256_file(&figure_path).expect("hash"),
                    "expected_story": "old",
                    "review": {
                        "opened": true,
                        "status": "pass",
                        "actual_finding": "Representative.",
                        "reviewed_at": "2026-04-13T01:00:00Z"
                    }
                }
            ]
        });
        std::fs::write(
            &manifest_path,
            serde_json::to_string_pretty(&seeded_manifest).expect("serialize"),
        )
        .expect("seed manifest");
        std::fs::write(&figure_path, "<svg>option1-v2</svg>").expect("mutate figure");

        let figure_specs = vec![NarrativeFigureSpec {
            number: 4,
            title: "Option 1".to_string(),
            path: "../report/figures/selected_option1_schematic.svg".to_string(),
            caption: "Option 1 caption".to_string(),
            alt: "Option 1".to_string(),
        }];
        let outcome = write_asset_review_manifest(&report_dir, &narrative_path, &figure_specs)
            .expect("asset review manifest should rewrite");
        let manifest: super::AssetReviewManifest =
            serde_json::from_str(&std::fs::read_to_string(outcome.manifest_path).expect("read"))
                .expect("deserialize");

        let figure_entry = manifest
            .assets
            .iter()
            .find(|entry| entry.asset_type == super::AssetType::Svg)
            .expect("figure entry");
        assert!(!figure_entry.review.opened);
        assert_eq!(figure_entry.review.status, ReviewStatus::Fail);
    }

    fn write_svg(report_dir: &Path, name: &str, body: &str) -> PathBuf {
        let path = report_dir.join("figures").join(name);
        std::fs::create_dir_all(path.parent().expect("figure parent")).expect("mkdir");
        std::fs::write(&path, body).expect("write figure");
        path
    }

    fn write_narrative(report_dir: &Path, body: &str) -> PathBuf {
        let path = report_dir.join("ARPA-H_SonALAsense_Milestone 12 Report.md");
        std::fs::write(&path, body).expect("write narrative");
        path
    }
}
