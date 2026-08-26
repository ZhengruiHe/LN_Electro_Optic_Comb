# Simulation Roadmap

## Baseline concept

Use one optical racetrack resonator with two spatially separated functions:

1. an unpoled electro-optic modulation section with RF electrodes;
2. a periodically poled section for quasi-phase-matched SHG.

The first design target is deliberately limited: generate a resonant EO comb
around the fundamental wavelength and convert it into a second-harmonic comb.
Only the central fundamental and second-harmonic modes are initially required
to be doubly resonant. Equal FSRs across the two bands are not assumed.

## Stage 0: freeze the fabrication stack

Record the foundry and design-rule inputs before optimizing geometry:

- LN cut and film thickness;
- etch depth, sidewall angle, and minimum feature size;
- BOX and top-cladding materials and thicknesses;
- permitted metal stack, minimum metal-to-waveguide gap, and heater rules;
- target fundamental wavelength and available RF band.

Deliverable: `docs/design_inputs.md` with values, sources, and uncertainties.

Gate: no final waveguide or electrode dimensions are selected until the stack
and fabrication constraints are known.

## Stage 1: optical cross-section and mode-family scan

Use Lumerical MODE or COMSOL Wave Optics at the fundamental and second harmonic.
Sweep waveguide width and etch depth, including all plausible mode families.

Extract:

- effective index and group index;
- mode profiles and polarization;
- propagation and bend loss estimates;
- group-velocity dispersion;
- nonlinear mode-overlap integral;
- electro-optic field component relevant to the selected LN tensor element.

Deliverables:

- mode maps at both wavelengths;
- `n_eff`, `n_g`, and dispersion versus geometry;
- ranked fundamental/second-harmonic mode pairs.

Gate: select a mode pair that has acceptable confinement, overlap, loss, and
fabrication tolerance. Do not require equal group indices unless later analysis
shows that line-by-line dual-band resonance is essential.

## Stage 2: PPLN quasi-phase matching

For each selected mode pair, calculate

`Delta k = beta_2w - 2 beta_w`

and the first-order QPM period

`Lambda = 2*pi/abs(Delta k)`.

Include duty-cycle, domain-wall placement, temperature, width, and etch-depth
sensitivity. Compute the nonlinear overlap using the full LN tensor rather than
a scalar effective index approximation.

Deliverables:

- nominal poling period and permitted fabrication window;
- normalized SHG coupling coefficient;
- predicted phase-matching bandwidth and temperature sensitivity.

Gate: the required poling period and duty cycle must comply with the foundry
rules and remain tolerant to expected fabrication variation.

## Stage 3: racetrack resonance and double-resonance search

Choose a target fundamental FSR from the available RF source. Estimate the
racetrack length from the fundamental group index, then refine it using the
frequency-dependent propagation constant.

Search for longitudinal mode pairs satisfying the central double-resonance
condition near `omega_2w = 2*omega_w`. Keep separate records of:

- fundamental FSR;
- second-harmonic FSR;
- central double-resonance mismatch;
- integrated dispersion in each band;
- loaded and intrinsic optical Q targets.

Deliverables:

- racetrack radius and straight-section lengths;
- resonance-frequency tables for both mode families;
- thermal and electro-optic trimming ranges required after fabrication.

Gate: at least one central mode pair must be alignable within the loaded cavity
linewidth using realistic tuning. Equal FSRs are not a first-version gate.

## Stage 4: bus-to-resonator couplers

Simulate local coupling regions with FDTD or EME instead of a full three-
dimensional millimetre-scale resonator. Prefer separate wavelength-selective
couplers for the fundamental input and second-harmonic extraction.

Extract coupling coefficients versus wavelength, gap, coupling length, and
fabrication error. Set coupling targets separately for the two bands.

Gate: both bands must have usable external coupling without forcing excessive
loss or an impractical shared coupler geometry.

## Stage 5: RF electrode simulation

Use HFSS or COMSOL RF for the unpoled modulation section. Start with a 2D or
quasi-static cross-section, then validate pads, bends, feeds, and discontinuities
with a 3D model.

Extract:

- characteristic impedance and propagation constant;
- RF loss and electric-field distribution;
- S11 and S21 over the target band;
- microwave resonance and Q if a resonant electrode is used;
- optical-RF overlap at both optical wavelengths.

Gate: the electrode must provide adequate EO coupling without unacceptable
metal absorption, heating, or impedance mismatch.

## Stage 6: electro-optic and nonlinear coupled-mode model

Implement a reproducible Python or MATLAB model containing:

- the fundamental optical mode family;
- the second-harmonic mode family;
- pump detuning and external coupling;
- optical loss and integrated dispersion;
- RF-driven EO coupling;
- SHG and sum-frequency coupling between comb lines;
- thermal or photorefractive detuning where needed.

Begin with the carrier and first sidebands, verify energy conservation and
limiting cases, then expand the number of modes.

Gate: reproduce single-band resonant EO comb behaviour and doubly resonant SHG
separately before enabling all coupling terms simultaneously.

## Stage 7: tolerance and control-loop design

Run Monte Carlo or corner sweeps for waveguide width, film thickness, etch depth,
poling period, coupling gap, and temperature. Determine the required tuning and
locking architecture:

- laser-to-fundamental-resonance lock;
- central SH double-resonance alignment;
- RF-to-fundamental-FSR adjustment;
- heater and fast EO actuator ranges;
- lock acquisition and recovery sequence.

Gate: the tuning range must cover expected fabrication offsets with margin, and
the actuator bandwidths must be separated enough to avoid control-loop conflict.

## Stage 8: layout and verification

Only after the previous gates pass:

- generate the optical, PPLN, RF, heater, and pad layout;
- run project-level geometric and connectivity checks;
- run foundry DRC when the official rule deck is available;
- archive simulator versions, material models, and all parameter sources.

## Immediate next task

Create `docs/design_inputs.md` from the actual wafer stack and laboratory RF
capabilities. Then run the Stage 1 optical mode scan; HFSS electrode optimization
starts only after the fundamental mode and target RF frequency are selected.

