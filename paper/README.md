# Manuscript 1: figures and tables

`make_figures.py` builds every display item from the repository's own model outputs, so the figures always match the code.

Format follows Elsevier artwork guidelines (Field Crops Research):
- widths 90, 140 or 190 mm;
- Arial-metric sans-serif, 6.5–8 pt;
- Okabe–Ito colour-blind-safe palette;
- vector PDF plus 600 dpi PNG.

Rebuild after changing model parameters:

```bash
python rice-crf-map/build.py && python wheat-crf-map/build.py && python maize-crf-map/build.py
python analysis/mixed_shapes.py && python analysis/beta_scan.py
python paper/make_figures.py
```

| File | Content |
|---|---|
| `figures/Fig1_worked_example` | Static-water curves (β 1.3 vs 2.5), field supply vs demand and soil pool vs floor for Wuhan mid-season rice |
| `figures/Fig2_global_recipes` | Recommended main control period and total N for every rice, wheat and maize season |
| `figures/Fig3_grades_vs_shape` | N saved vs split urea; further saving from a 2nd linear grade, an S-shaped grade, or linear + S |
| `figures/Fig4_beta_scan` | β scan: yield-loss proxy at 85% N (mean, worst) and N needed at equal yield, by crop group |
| `figures/Fig5_robustness` | Pool trajectories in cool, normal and warm years for β 1.3, 2.0, 4.0 (NE China maize) |
| `figures/Fig6_environment` | Change in N, NH3, N2O, leaching, GHG and damage cost at the optimal β vs β 1.3 |
| `figures/FigS1_validation_demo` | Validation workflow on simulated burial-bag data (clearly labelled as simulated) |
| `tables/Table2_summary_by_crop.csv` | Per-crop summary: grades, urea share, N saved, optimal β |
