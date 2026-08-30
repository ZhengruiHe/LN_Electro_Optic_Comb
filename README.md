# LN Electro-Optic Comb

Design workspace for a compact thin-film lithium-niobate (TFLN)
electro-optic comb and poling-free second-harmonic platform.

The first tapeout is organized around three structures:

1. a PDK-aligned 1550-nm travelling-wave phase modulator for 25-GHz
   electro-optic-comb characterization;
2. a compact x-cut LN spontaneous-quasi-phase-matched (SQPM) micro-racetrack
   for poling-free 1550-to-775-nm SHG;
3. a monolithic `CW SHG -> dual-wavelength travelling-wave PM` sequence that
   creates 1550- and 775-nm EO combs with the same RF line spacing.

The EO modulation remains non-resonant; only the compact SHG block is
resonant. SHG occurs before modulation, so the SHG racetrack FSR does not need
to equal the 25-GHz RF drive. The resonator still requires a 1550/775-nm
double-resonant mode pair and thermal control.

## Key documents

- [`docs/three_structure_layout_plan.md`](docs/three_structure_layout_plan.md):
  literature-backed architecture, dimensions, area budget, DOE, simulation,
  and measurement plan
- [`docs/three_structure_floorplan.svg`](docs/three_structure_floorplan.svg):
  compact 21.8 mm x 3.8 mm PDK-area floorplan
- [`docs/simulation_roadmap.md`](docs/simulation_roadmap.md): staged execution
  gates from stack definition to tapeout
- [`docs/pdk_compatibility_review.md`](docs/pdk_compatibility_review.md):
  supported components, custom 775-nm gaps, and foundry questions
- [`models/system/design_targets.json`](models/system/design_targets.json):
  machine-readable preliminary design targets
- [`scripts/design_budget.py`](scripts/design_budget.py): reproducible comb
  span, QPM bandwidth, conversion, and area estimates

## Repository layout

- `docs/`: design decisions, evidence tables, and simulation roadmap
- `models/optical/`: optical modes, QPM, couplers, and nonlinear propagation
- `models/rf/`: HFSS or COMSOL RF electrode models
- `models/system/`: EO-comb and SHG/SFG system models
- `scripts/`: reproducible parameter sweeps and post-processing
- `results/`: generated figures and exported data; ignored by default

The foundry PDK archive is confidential local input and is intentionally
excluded from Git. The current S2/S3 layouts are simulation concepts, not
fabrication-ready cells, because the supplied PDK qualifies 1550-nm devices
but does not qualify SHG or 775-nm components.
