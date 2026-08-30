# LN Electro-Optic Comb

Design workspace for a compact thin-film lithium-niobate (TFLN)
electro-optic comb and poling-free second-harmonic platform.

The first design study is organized around three fully custom structures:

1. a 1550-nm four-pass optical-recycling travelling-wave phase modulator for
   low-`Vpi` 25-GHz electro-optic-comb generation;
2. a compact x-cut LN spontaneous-quasi-phase-matched (SQPM) micro-racetrack
   for poling-free 1550-to-775-nm SHG;
3. a monolithic `CW SHG -> dual-rail travelling-wave PM` sequence in which
   separately optimized 1550- and 775-nm waveguides share one RF electrode.

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
excluded from Git. It is used only for the process stack, layer rules and die
boundary; none of the three structures relies on a supplied device BlackBox.
All dimensions are simulation starting points rather than fabrication-ready
mask dimensions.
