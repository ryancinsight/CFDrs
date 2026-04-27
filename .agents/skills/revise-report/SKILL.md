---
name: revise-report
description: Process all pending report review annotations — find the generating Rust code, apply each requested change, rebuild, and mark annotations resolved. Use when the user wants to apply accumulated report review comments.
allowed-tools: Read, Edit, Write, Grep, Glob, Bash, Agent
---

Process every `status = "pending"` annotation in `report/review/pending.toml` and apply the requested changes to the **Rust generation code** (not the generated `.md` files).

## Step-by-step

1. **Read** `report/review/pending.toml` and collect all entries where `status = "pending"`.
   - If none are pending, tell the user and stop.

2. **For each pending annotation**, in order:
   a. **Locate the source** — grep for the `text` value (or a distinctive substring) in the narrative generation code under `crates/cfd-optim/src/reporting/narrative/` and `report/templates/`. Also check `crates/cfd-optim/src/reporting/canonical/writer.rs` and `crates/cfd-optim/src/reporting/figures/`.
   b. **Apply the change** described in `comment` to the Rust source (format strings, static text, conditional logic, etc.). Preserve all existing `{:.Nf}` / `{:.Ne}` format specifiers unless the comment explicitly asks to change them.
   c. **Update the annotation** in `pending.toml`: set `status = "resolved"` and fill in `source` with the file path that was changed.
   d. **Report** what you changed (file, line, before → after summary).

3. **Build** — run `cargo build -p cfd-optim --release --no-default-features` and fix any errors.

4. **Run tests** — run `cargo test -p cfd-optim --no-default-features --lib` and confirm all pass.

5. **Summarize** — print a table of all annotations processed: text snippet, comment, source file, status.

## Rules

- Never edit the generated report `.md` files directly — only edit Rust source and templates.
- If you cannot locate the generating code for an annotation, leave it as `status = "pending"` and note the issue in your summary.
- If a comment is ambiguous, apply your best interpretation but flag it in the summary for the user to verify.
- When the template file `report/templates/m12_narrative_template.md` contains the text, check whether the text is a `{{PLACEHOLDER}}` — if so, trace it to the writer code that fills the placeholder value.
- Keep changes minimal — only modify what the annotation asks for.
