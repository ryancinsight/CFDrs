//! Persistent on-disk cache for expensive `SdtMetrics` evaluations.

use crate::{DesignCandidate, SdtMetrics};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

/// Bump this when the metric model or blueprint semantics change.
pub const METRICS_CACHE_VERSION: &str = "m12_metrics_v1";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MetricsCacheStats {
    pub hits: usize,
    pub misses: usize,
}

#[derive(Debug, Clone)]
pub struct MetricsCache {
    root: Arc<PathBuf>,
    version: &'static str,
    hits: Arc<AtomicUsize>,
    misses: Arc<AtomicUsize>,
}

#[derive(Debug, Serialize, Deserialize)]
struct MetricsCacheEntry {
    cache_version: String,
    candidate_json: String,
    metrics: SdtMetrics,
}

impl MetricsCache {
    /// Create a cache rooted at `root`.
    ///
    /// # Errors
    /// Returns an error if the cache directory cannot be created.
    pub fn new(root: impl AsRef<Path>) -> Result<Self, Box<dyn std::error::Error>> {
        let root = root.as_ref().to_path_buf();
        std::fs::create_dir_all(&root)?;
        Ok(Self {
            root: Arc::new(root),
            version: METRICS_CACHE_VERSION,
            hits: Arc::new(AtomicUsize::new(0)),
            misses: Arc::new(AtomicUsize::new(0)),
        })
    }

    /// Fetch cached metrics or compute and persist them on a miss.
    ///
    /// # Errors
    /// Returns an error if metric computation fails.
    pub fn get_or_compute<E, F>(
        &self,
        candidate: &DesignCandidate,
        compute: F,
    ) -> Result<SdtMetrics, Box<dyn std::error::Error>>
    where
        E: std::error::Error + Send + Sync + 'static,
        F: FnOnce() -> Result<SdtMetrics, E>,
    {
        let candidate_json = serde_json::to_string(candidate)?;
        let cache_path = self.cache_path(&candidate_json);

        if let Ok(cached) = self.try_load(&cache_path, &candidate_json) {
            self.hits.fetch_add(1, Ordering::Relaxed);
            return Ok(cached);
        }

        self.misses.fetch_add(1, Ordering::Relaxed);
        let metrics = compute().map_err(|e| -> Box<dyn std::error::Error> { Box::new(e) })?;
        let entry = MetricsCacheEntry {
            cache_version: self.version.to_string(),
            candidate_json,
            metrics: metrics.clone(),
        };
        self.persist(&cache_path, &entry)?;
        Ok(metrics)
    }

    #[must_use]
    pub fn stats(&self) -> MetricsCacheStats {
        MetricsCacheStats {
            hits: self.hits.load(Ordering::Relaxed),
            misses: self.misses.load(Ordering::Relaxed),
        }
    }

    fn try_load(
        &self,
        path: &Path,
        candidate_json: &str,
    ) -> Result<SdtMetrics, Box<dyn std::error::Error>> {
        let raw = std::fs::read_to_string(path)?;
        let entry: MetricsCacheEntry = serde_json::from_str(&raw)?;
        if entry.cache_version == self.version && entry.candidate_json == candidate_json {
            Ok(entry.metrics)
        } else {
            Err("stale cache entry".into())
        }
    }

    fn cache_path(&self, candidate_json: &str) -> PathBuf {
        let mut hasher = Sha256::new();
        hasher.update(candidate_json.as_bytes());
        hasher.update(self.version.as_bytes());
        let digest = hasher.finalize();
        let mut hex = String::with_capacity(digest.len() * 2);
        for byte in digest {
            let _ = write!(hex, "{byte:02x}");
        }
        self.root.join(format!("{hex}.json"))
    }

    fn persist(
        &self,
        path: &Path,
        entry: &MetricsCacheEntry,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let payload = serde_json::to_string(entry)?;
        let stamp = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)?
            .as_nanos();
        let tmp_path = path.with_extension(format!("tmp-{}-{stamp}", std::process::id()));
        std::fs::write(&tmp_path, payload)?;
        match std::fs::rename(&tmp_path, path) {
            Ok(()) => Ok(()),
            Err(_) => {
                if path.exists() {
                    let _ = std::fs::remove_file(&tmp_path);
                    Ok(())
                } else {
                    std::fs::rename(&tmp_path, path)?;
                    Ok(())
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::design::{
        CrossSectionShape, DesignTopology, PrimitiveSplitSequence, TreatmentZoneMode,
    };

    fn sample_candidate() -> DesignCandidate {
        DesignCandidate {
            id: "cache-test".to_string(),
            topology: DesignTopology::PrimitiveSelectiveTree {
                sequence: PrimitiveSplitSequence::TriTri,
            },
            flow_rate_m3_s: 1.667e-6,
            inlet_gauge_pa: 100_000.0,
            throat_diameter_m: 45e-6,
            inlet_diameter_m: 4e-3,
            throat_length_m: 135e-6,
            channel_width_m: 6e-3,
            channel_height_m: 1e-3,
            serpentine_segments: 5,
            segment_length_m: 8e-3,
            bend_radius_m: 3.5e-3,
            feed_hematocrit: 0.45,
            trifurcation_center_frac: 0.45,
            cif_pretri_center_frac: 0.45,
            cif_terminal_tri_center_frac: 0.45,
            cif_terminal_bi_treat_frac: 0.76,
            asymmetric_narrow_frac: 0.5,
            trifurcation_left_frac: 1.0 / 3.0,
            cross_section_shape: CrossSectionShape::Rectangular,
            treatment_zone_mode: TreatmentZoneMode::UltrasoundOnly,
            centerline_venturi_throat_count: 1,
        }
    }

    #[test]
    fn cache_reuses_existing_metrics() {
        let cache_root = std::env::temp_dir().join(format!(
            "cfd-optim-cache-test-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .expect("system time must be after unix epoch")
                .as_nanos()
        ));
        let cache = MetricsCache::new(&cache_root).expect("cache must initialize");
        let candidate = sample_candidate();
        let expected = SdtMetrics {
            flow_rate_ml_min: 100.0,
            therapy_channel_fraction: 0.25,
            ..SdtMetrics::default()
        };

        let first = cache
            .get_or_compute(&candidate, || Ok::<_, std::io::Error>(expected.clone()))
            .expect("first compute should succeed");
        let second = cache
            .get_or_compute(&candidate, || -> Result<SdtMetrics, std::io::Error> {
                panic!("cached lookup should not recompute");
            })
            .expect("cached compute should succeed");

        assert_eq!(first.flow_rate_ml_min, expected.flow_rate_ml_min);
        assert_eq!(
            second.therapy_channel_fraction,
            expected.therapy_channel_fraction
        );
        assert_eq!(cache.stats(), MetricsCacheStats { hits: 1, misses: 1 });

        let _ = std::fs::remove_dir_all(cache_root);
    }
}
