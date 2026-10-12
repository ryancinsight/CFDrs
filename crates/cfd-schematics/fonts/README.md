# Embedded typeface

`DejaVuSans.ttf` is DejaVu Sans 2.35 reduced to the code points schematic
labels use: Basic Latin, Latin-1 Supplement, Latin Extended-A, Greek and
Coptic, General Punctuation, superscripts and subscripts, the euro sign,
letterlike symbols, arrows, mathematical operators, U+25A1 and U+FFFD. The
plotters backend registers it as `sans-serif` so text renders without a host
font stack. `LICENSE_DEJAVU` is the license the source font ships with; it
permits redistribution and modification of fonts not named Bitstream or Vera.

Source: `DejaVuSans.ttf` from matplotlib 3.10.3 `mpl-data/fonts/ttf`, SHA-256
`3fdf69cabf06049ea70a00b5919340e2ce1e6d02b0cc3c4b44fb6801bd1e0d22`. Reduced
with fontTools 4.61.1:

```sh
python -m fontTools.subset DejaVuSans.ttf \
  --unicodes="U+0020-007E,U+00A0-017F,U+0370-03FF,U+2000-206F,U+2070-209F,U+20AC,U+2100-214F,U+2190-21FF,U+2200-22FF,U+25A1,U+FFFD" \
  --layout-features='kern,liga' --name-IDs='*' --name-legacy \
  --name-languages='*' --glyph-names=false --notdef-outline \
  --output-file=DejaVuSans.ttf
```
