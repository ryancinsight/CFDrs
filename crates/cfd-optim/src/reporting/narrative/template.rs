//! Lightweight strict template renderer for Milestone 12 narrative markdown.
//!
//! Placeholders use `{{TOKEN}}` form. Rendering fails if any placeholder
//! remains unresolved to prevent silent drift in milestone reports.

use std::collections::{BTreeMap, BTreeSet};

/// Render a template using strict placeholder replacement.
///
/// Uses a single-pass scan over the template to resolve `{{KEY}}` placeholders
/// in one allocation, instead of N successive `String::replace` calls which
/// each clone the entire template.
///
/// # Errors
/// Returns an error if unresolved placeholders remain in the rendered output.
pub fn render_template_strict(
    template: &str,
    values: &BTreeMap<String, String>,
) -> Result<String, Box<dyn std::error::Error>> {
    // Single-pass: scan for `{{` ... `}}`, look up the key, and emit either
    // the replacement value or the raw token.  This replaces the previous
    // N × String::replace pattern which allocated N intermediate Strings.
    let mut rendered = String::with_capacity(template.len());
    let mut unresolved = BTreeSet::new();
    let mut idx = 0;
    while idx < template.len() {
        if let Some(open_rel) = template[idx..].find("{{") {
            let open = idx + open_rel;
            rendered.push_str(&template[idx..open]);
            let search_start = open + 2;
            if let Some(close_rel) = template[search_start..].find("}}") {
                let close = search_start + close_rel;
                let key = template[search_start..close].trim();
                if let Some(value) = values.get(key) {
                    rendered.push_str(value);
                } else {
                    // Preserve unresolved token verbatim for error reporting.
                    rendered.push_str(&template[open..close + 2]);
                    if !key.is_empty() {
                        unresolved.insert(key.to_string());
                    }
                }
                idx = close + 2;
            } else {
                // No closing `}}` — emit the rest verbatim.
                rendered.push_str(&template[open..]);
                idx = template.len();
            }
        } else {
            rendered.push_str(&template[idx..]);
            break;
        }
    }

    if unresolved.is_empty() {
        Ok(rendered)
    } else {
        Err(format!(
            "unresolved template placeholders: {}",
            unresolved
                .into_iter()
                .enumerate()
                .fold(String::new(), |mut acc, (i, token)| {
                    if i > 0 {
                        acc.push_str(", ");
                    }
                    acc.push_str(&token);
                    acc
                })
        )
        .into())
    }
}

/// Extract all `{{TOKEN}}` placeholders from a template string.
///
/// Retained for testing; production rendering uses the single-pass scanner
/// in [`render_template_strict`] which detects unresolved tokens inline.
#[cfg(test)]
#[must_use]
pub fn extract_placeholders(template: &str) -> BTreeSet<String> {
    let mut out = BTreeSet::new();
    let mut idx = 0usize;
    while let Some(open_rel) = template[idx..].find("{{") {
        let open = idx + open_rel;
        let search_start = open + 2;
        if let Some(close_rel) = template[search_start..].find("}}") {
            let close = search_start + close_rel;
            let token = template[search_start..close].trim();
            if !token.is_empty() {
                out.insert(token.to_string());
            }
            idx = close + 2;
        } else {
            break;
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn strict_renderer_replaces_all_tokens() {
        let tpl = "A {{ONE}} B {{TWO}}";
        let mut map = BTreeMap::new();
        map.insert("ONE".to_string(), "1".to_string());
        map.insert("TWO".to_string(), "2".to_string());

        let rendered = render_template_strict(tpl, &map).expect("template should render");
        assert_eq!(rendered, "A 1 B 2");
    }

    #[test]
    fn strict_renderer_fails_on_unresolved() {
        let tpl = "A {{ONE}} B {{TWO}}";
        let mut map = BTreeMap::new();
        map.insert("ONE".to_string(), "1".to_string());
        let err = render_template_strict(tpl, &map).expect_err("must fail");
        assert!(err.to_string().contains("TWO"));
    }
}
