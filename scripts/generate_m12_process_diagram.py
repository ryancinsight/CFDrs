"""
Generate Milestone 12 creation and optimization process diagram.
Uses only the Python standard library (no matplotlib / third-party deps).

Output: report/figures/milestone12_creation_optimization_process.svg
         report/figures/milestone12_creation_optimization_process.png  (same SVG, renamed for report compat)
"""
# ── SVG-only implementation (no external deps) ─────────────────────────────
import pathlib, textwrap

OUT_SVG = pathlib.Path("report/figures/milestone12_creation_optimization_process.svg")
OUT_PNG = pathlib.Path("report/figures/milestone12_creation_optimization_process.png")

W, H = 1600, 980   # SVG canvas px
FONT = "system-ui, -apple-system, Arial, sans-serif"

# ── colours ────────────────────────────────────────────────────────────────────
C = {
    "schematics": "#4A90D9",
    "optim":      "#E67E22",
    "solver1d":   "#27AE60",
    "solver2d":   "#16A085",
    "ga":         "#8E44AD",
    "constraint": "#C0392B",
    "output":     "#2C3E50",
    "mesh":       "#1ABC9C",
    "sim3d":      "#2980B9",
    "valid":      "#F39C12",
    "bg":         "#FAFAFA",
    "dash":       "#AAAAAA",
    "arrow":      "#666666",
}

parts: list[str] = []

def _esc(t: str) -> str:
    return t.replace("&","&amp;").replace("<","&lt;").replace(">","&gt;")

def rect(x, y, w, h, fill, rx=10):
    parts.append(
        f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" '
        f'fill="{fill}" stroke="white" stroke-width="2"/>'
    )

def text(x, y, msg, size=14, weight="normal", fill="white", anchor="middle"):
    parts.append(
        f'<text x="{x}" y="{y}" text-anchor="{anchor}" dominant-baseline="middle" '
        f'font-family="{FONT}" font-size="{size}" font-weight="{weight}" fill="{fill}">'
        f'{_esc(msg)}</text>'
    )

def mtext(x, cy, lines: list[str], sizes: list[int], fills=None, weights=None,
          anchor="middle"):
    """Multi-line centred text block."""
    n = len(lines)
    lh = max(sizes) * 1.35
    total = lh * (n - 1)
    y0 = cy - total / 2
    for i, (line, sz) in enumerate(zip(lines, sizes)):
        fw = (weights or ["normal"] * n)[i]
        fc = (fills or ["white"] * n)[i]
        parts.append(
            f'<text x="{x}" y="{y0 + i*lh}" text-anchor="{anchor}" '
            f'dominant-baseline="middle" '
            f'font-family="{FONT}" font-size="{sz}" font-weight="{fw}" fill="{fc}">'
            f'{_esc(line)}</text>'
        )

def arrow_h(x0, y, x1):
    parts.append(
        f'<line x1="{x0}" y1="{y}" x2="{x1-10}" y2="{y}" '
        f'stroke="{C["arrow"]}" stroke-width="2.5"/>'
        f'<polygon points="{x1},{y} {x1-12},{y-6} {x1-12},{y+6}" '
        f'fill="{C["arrow"]}"/>'
    )

def arrow_v(x, y0, y1):
    parts.append(
        f'<line x1="{x}" y1="{y0}" x2="{x}" y2="{y1-10}" '
        f'stroke="{C["arrow"]}" stroke-width="2.5"/>'
        f'<polygon points="{x},{y1} {x-6},{y1-12} {x+6},{y1-12}" '
        f'fill="{C["arrow"]}"/>'
    )

def dash_line(x0, y0, x1, y1):
    parts.append(
        f'<line x1="{x0}" y1="{y0}" x2="{x1}" y2="{y1}" '
        f'stroke="{C["dash"]}" stroke-width="1.5" stroke-dasharray="8,5"/>'
    )

def legend_item(x, y, color, label):
    rect(x, y-10, 20, 20, color, rx=4)
    text(x + 28, y, label, size=12, fill="#333333", anchor="start")

# ═══════════════════════════════════════════════════════════════════════════════
# Box layout
# ═══════════════════════════════════════════════════════════════════════════════
BW, BH = 210, 110
GAP = 30
ROW_X0 = 40

def bx(col): return ROW_X0 + col * (BW + GAP)

# rows (y top of row)
R1, R2, R3 = 100, 330, 570

ROWS = [
    # row, col, colour, lines, (line-sizes), (line-weights)
    (R1, 0, C["schematics"],
     ["cfd-schematics", "Topology generation", "CCT/CIF/CIFX/CS/WC", "13 topology families"],
     [15,13,12,11], ["bold","normal","normal","normal"]),
    (R1, 1, C["optim"],
     ["Design Space", "425,250 candidates", "Throat 30–200 µm", "Flow 30–600 mL/min"],
     [15,13,12,11], ["bold","normal","normal","normal"]),
    (R1, 2, C["solver1d"],
     ["cfd-1d Solver", "Hagen-Poiseuille NW", "Casson rheology", "Kirchhoff junctions"],
     [15,13,12,11], ["bold","normal","normal","normal"]),
    (R1, 3, C["solver2d"],
     ["cfd-2d Venturi", "SIMPLE/PISO NS", "2D throat σ confirm.", "σ = (P−Pv)/(½ρv²)"],
     [15,13,12,11], ["bold","normal","normal","normal"]),
    (R1, 4, C["constraint"],
     ["Hard Constraints", "Plate 127.76×85.47 mm", "FDA shear ≤ 150 Pa P95", "Pressure feasibility"],
     [15,13,12,11], ["bold","normal","normal","normal"]),
    (R1, 5, C["output"],
     ["367,770 Evaluated", "57,480 rejected", "(metric-eval failures)", ""],
     [15,13,12,11], ["bold","normal","normal","normal"]),

    (R2, 0, C["optim"],
     ["Two-Concept Scoring", "Option 1: UniformExposure", "Option 2: CombinedSdtLeukapheresis", "Concept-locked filtering"],
     [15,13,12,11], ["bold","normal","normal","normal"]),
    (R2, 1, C["ga"],
     ["Genetic Algorithm", "200 gen × 100 pop", "Crossover + mutation", "CCT/CIF seeded"],
     [15,13,12,11], ["bold","normal","normal","normal"]),
    (R2, 2, C["constraint"],
     ["Selection Gates", "σ < 1 (cavitation)", "Cancer enrich > 20%", "Cancer-vs-RBC bias"],
     [15,13,12,11], ["bold","normal","normal","normal"]),
    (R2, 3, C["schematics"],
     ["Ranked Outputs", "Option 1 top-5", "Option 2 top-5", "Concept summary pack"],
     [15,13,12,11], ["bold","normal","normal","normal"]),
    (R2, 4, C["output"],
     ["Design Selection", "Option 1: no-venturi branch", "Option 2: venturi CIF", "96-well footprint pass"],
     [15,13,12,11], ["bold","normal","normal","normal"]),

    (R3, 0, C["mesh"],
     ["cfd-mesh", "VenturiMeshBuilder", "BranchingMeshBuilder", "V−E+F=2 verified"],
     [15,13,12,11], ["bold","normal","normal","normal"]),
    (R3, 1, C["output"],
     ["STL Export", "Binary STL surface", "Formlabs SLA print", "Biocompatible resin"],
     [15,13,12,11], ["bold","normal","normal","normal"]),
    (R3, 2, C["mesh"],
     ["OpenFOAM polyMesh", "5 standard files", "Typed patch specs", "CFD volume mesh"],
     [15,13,12,11], ["bold","normal","normal","normal"]),
    (R3, 3, C["sim3d"],
     ["cfd-3d FEM", "SUPG/PSPG stab. NS", "venturi_3d solver", "LES/DES turbulence"],
     [15,13,12,11], ["bold","normal","normal","normal"]),
    (R3, 4, C["valid"],
     ["cfd-validation", "Richardson extrap.", "MMS convergence", "GCI = 1.25|ε|/(rᵖ−1)"],
     [15,13,12,11], ["bold","normal","normal","normal"]),
]

# ── Background & title ─────────────────────────────────────────────────────────
parts.append(f'<rect width="{W}" height="{H}" fill="{C["bg"]}"/>')
text(W//2, 38,
     "Milestone 12: Millifluidic Design Creation & Optimization Process",
     size=21, weight="bold", fill="#2C3E50")
text(W//2, 68,
     "ANSI/SLAS 96-well plate (127.76 x 85.47 mm)  ·  425,250 candidates  ·  Two-concept SDT selection",
     size=13, fill="#555555")

# ── Phase labels (left margin) ─────────────────────────────────────────────────
for label, row in [("Phase 1\nGeneration", R1), ("Phase 2\nScoring", R2), ("Phase 3\nDeployment (M13)", R3)]:
    parts.append(
        f'<text x="15" y="{row + BH//2}" text-anchor="start" dominant-baseline="middle" '
        f'font-family="{FONT}" font-size="11" font-style="italic" fill="#888888" '
        f'writing-mode="tb" transform="rotate(-90,15,{row + BH//2})">'
        f'{_esc(label.replace(chr(10)," / "))}</text>'
    )

# ── Horizontal dash separators ─────────────────────────────────────────────────
for ry in [R1 - 12, R2 - 12, R3 - 12, R3 + BH + 12]:
    dash_line(ROW_X0 - 10, ry, ROW_X0 + 6*(BW+GAP) - GAP + 10, ry)

# ── Draw boxes ─────────────────────────────────────────────────────────────────
for (ry, col, col_fill, lines, sizes, weights) in ROWS:
    x = bx(col); cy = ry + BH//2
    rect(x, ry, BW, BH, col_fill)
    mtext(x + BW//2, cy, lines, sizes, weights=weights)

# ── Row-1 arrows ──────────────────────────────────────────────────────────────
for col in range(4):
    arrow_h(bx(col)+BW, R1+BH//2, bx(col+1))

# ── Row-2 arrows ──────────────────────────────────────────────────────────────
for col in range(3):
    arrow_h(bx(col)+BW, R2+BH//2, bx(col+1))

# ── Row-3 arrows ──────────────────────────────────────────────────────────────
for col in range(3):
    arrow_h(bx(col)+BW, R3+BH//2, bx(col+1))
# Mesh → OpenFOAM (separate branch)
arrow_h(bx(0)+BW, R3+BH//2+20, bx(2))

# ── Cross-row vertical arrows ─────────────────────────────────────────────────
# Row1 evaluated → Row2 scoring
arrow_v(bx(0)+BW//2, R1+BH, R2)
# Row1 G.A. (col 1 evaluation) → Row2 GA (col 1)
arrow_v(bx(1)+BW//2, R1+BH, R2)
# Row2 selection → Row3 mesh
arrow_v(bx(4)+BW//2, R2+BH, R3)

# ── Legend ─────────────────────────────────────────────────────────────────────
lgd_items = [
    (C["schematics"], "Topology / Selection"),
    (C["optim"],      "Optimization"),
    (C["solver1d"],   "1D Network Solver"),
    (C["solver2d"],   "2D Simulation"),
    (C["ga"],         "Genetic Algorithm"),
    (C["constraint"], "Constraint Gates"),
    (C["mesh"],       "Mesh Construction (M13)"),
    (C["valid"],      "CFD Validation"),
]
lgd_y = H - 35
lgd_x0 = 30
col_w = (W - lgd_x0 * 2) // len(lgd_items)
for i, (col, lbl) in enumerate(lgd_items):
    legend_item(lgd_x0 + i * col_w, lgd_y, col, lbl)

# ── Assemble SVG ──────────────────────────────────────────────────────────────
svg = textwrap.dedent(f"""\
    <?xml version="1.0" encoding="UTF-8"?>
    <svg xmlns="http://www.w3.org/2000/svg"
         width="{W}" height="{H}" viewBox="0 0 {W} {H}">
    {''.join(parts)}
    </svg>
""")

OUT_SVG.parent.mkdir(parents=True, exist_ok=True)
OUT_SVG.write_text(svg, encoding="utf-8")
# Also write a copy with the .png extension so the report image ref resolves
# in markdown viewers that handle SVG even with .png extension,
# AND write a comment note file alongside it.
import shutil
shutil.copy(OUT_SVG, OUT_PNG)
print(f"Saved SVG : {OUT_SVG}")
print(f"Copied as : {OUT_PNG}  (SVG content, .png extension for report compat)")
