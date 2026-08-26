# LN Electro-Optic Comb

Design workspace for a lithium-niobate optical resonator combining:

- resonant electro-optic comb generation near the fundamental wavelength;
- PPLN quasi-phase-matched second-harmonic generation;
- simultaneous fundamental and second-harmonic optical resonances;
- RF electrode and electro-optic overlap design;
- coupled-mode simulation of the two spectral bands.

The current baseline is a single optical racetrack resonator with separate
functional sections for PPLN frequency conversion and electro-optic modulation.

## Repository layout

- `docs/`: design decisions, parameter tables, and simulation roadmap
- `models/optical/`: optical mode, resonator, coupler, and QPM models
- `models/rf/`: HFSS or COMSOL RF electrode models
- `models/system/`: coupled-mode and comb-dynamics models
- `scripts/`: reproducible parameter sweeps and post-processing
- `results/`: generated figures and exported data; ignored by default

See `docs/simulation_roadmap.md` for the staged implementation plan.

