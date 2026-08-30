# LN Electro-Optic Comb

Design workspace for a non-resonant thin-film lithium-niobate (TFLN)
electro-optic comb and periodically poled lithium-niobate (PPLN)
second-harmonic platform.

The first tapeout is organized around three structures:

1. a low-`Vpi` travelling-wave, optical-recycling phase modulator for
   25-GHz electro-optic comb generation near 1550 nm;
2. a stand-alone PPLN quasi-phase-matched 1550-to-775-nm SHG test array;
3. a monolithic travelling-wave EO-comb-to-PPLN-SHG cascade that generates
   fundamental and second-harmonic combs with the same RF line spacing.

The current baseline is deliberately non-resonant. It does not require an
optical-ring FSR condition at either wavelength. A resonant implementation is
kept as a later comparison, after the travelling-wave and QPM building blocks
have been measured independently.

## Key documents

- [`docs/three_structure_layout_plan.md`](docs/three_structure_layout_plan.md):
  literature-backed architecture, dimensions, area budget, DOE, simulation,
  and measurement plan
- [`docs/three_structure_floorplan.svg`](docs/three_structure_floorplan.svg):
  preliminary 20 mm x 6 mm die floorplan
- [`docs/simulation_roadmap.md`](docs/simulation_roadmap.md): staged execution
  gates from stack definition to tapeout
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

