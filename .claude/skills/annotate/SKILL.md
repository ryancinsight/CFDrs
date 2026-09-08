---
name: annotate
description: Add a review annotation to the Milestone 12 report. Select text in the report, then run /annotate with your comment. Use when the user selects text in a report .md file and wants to mark it for revision.
argument-hint: "<comment describing the desired change>"
allowed-tools: Read, Edit, Write, Grep, Glob
---

The user has selected text in a generated Milestone 12 report and wants to annotate it for a future revision pass.

## Your task

1. **Identify the selected text** from the IDE selection context in this conversation. If no text is selected, ask the user to select the text they want to annotate.

2. **Identify the section** by looking at the nearest markdown heading above the selected text in the report file.

3. **Append a new `[[annotation]]` entry** to `report/review/pending.toml` with:
   - `text` — the selected text (trimmed, max ~200 chars; use a distinctive substring if longer)
   - `comment` — the user's comment from `$ARGUMENTS`
   - `status` — `"pending"`
   - `section` — the report section heading
   - `source` — `""` (filled when resolved)

4. **Confirm** by printing the annotation you added (text snippet + comment) so the user can verify.

## Rules

- Do NOT modify the report `.md` file — it is auto-generated and will be overwritten.
- Do NOT modify any Rust source code — that happens during `/revise-report`.
- Only append to `report/review/pending.toml`. Never remove or reorder existing entries.
- If the user provides no argument text, ask what change they want.
- Keep `text` values short but unique enough to grep for in the generating code.
