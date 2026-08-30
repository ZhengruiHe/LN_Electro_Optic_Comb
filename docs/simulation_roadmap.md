# Simulation Roadmap

## Baseline concept

The compact baseline contains three independently measurable structures:

1. a standard-PDK 1550-nm travelling-wave phase modulator for EO-comb
   characterization;
2. a poling-free x-cut LN spontaneous-quasi-phase-matched (SQPM)
   micro-racetrack for 1550-to-775-nm SHG;
3. an integrated sequence in which CW SHG occurs first and a custom
   dual-wavelength travelling-wave phase modulator creates the two combs
   afterward.

This order deliberately removes an RF-to-SHG-ring-FSR constraint.  The ring
still requires fundamental/SH double resonance, but the 25-GHz EO sidebands are
generated downstream and need not coincide with the ring modes.

All dimensions in `docs/three_structure_layout_plan.md` are planning values.
Foundry confirmation, the optical/RF simulations below, and official DRC take
precedence.

## Stage 0: reconstruct the real PDK stack

Obtain and record:

- exact LN1 and LN2 etch depths and their allowed combinations;
- LN sidewall angle and thickness tolerances;
- final top-cladding thickness and refractive-index model;
- residual SiN geometry beneath every LN region;
- M1 material, thickness, conductivity/sheet resistance, and RF limits;
- permitted 775-nm power, propagation and custom-coupler rules.

Recreate the layer booleans in a small KLayout test cell and verify that the
required SiN isolation trench follows every LN guide, bend and coupling region.

Gate: no final ring width, SQPM order, straight length, bus gap, or electrode
cross-section is released until these values are known.

## Stage 1: 1550/775-nm anisotropic optical modes

Use COMSOL Wave Optics or Lumerical MODE with a full anisotropic LN tensor.
Sweep both LN etches, top width, residual slab, sidewall angle, oxide, and
underlying SiN at 1550 and 775 nm.

Extract:

- all relevant quasi-TE mode families and any polarization hybridization;
- `n_eff`, `n_g`, confinement, bend loss and radiation leakage;
- nonlinear overlap using the full `chi(2)` tensor;
- EO overlap candidates at both wavelengths;
- sensitivities to width, thickness, etch and temperature.

Preferred interaction: an x-cut, all-quasi-TE process that can access a strong
`d33` contribution and remain compatible with downstream modulation.  A
type-I CPM solution using different polarizations is a fallback, not the first
choice, because the nonlinear coefficient and two-colour EO overlap may be
worse.

Gate: select one manufacturable fundamental/SH mode pair with adequate
nonlinear overlap and a 775-nm mode that can be routed and coupled.

## Stage 2: S1 travelling-wave EO characterization model

Treat the supplied 9.1-mm PM as a fixed BlackBox for layout.  Build an external
RF fixture model from its published ports/bounding box and the actual probe,
launch, termination and package geometry.

Use HFSS 2D/Q3D for the line cross-section when sufficient geometry is
available, then 3D HFSS for pads, transitions and termination.  Extract:

- `Z0(f)`, microwave effective index and conductor/dielectric loss;
- `S11`, `S21`, launch discontinuity and expected voltage along the line;
- estimated `Vpi(f)` only after combining the RF field with the optical mode;
- metal-induced optical loss.

The supplied PM's width, gap and length are locked and cannot form an electrode
DOE.  Add separate custom CPW coupons if low-`Vpi` optimization is required.

Gate at 25 GHz: `|Z0-50 ohm| <= 5 ohm`, `S11 < -10 dB`, and a predicted drive
that is safe for the probes and amplifier.  The PM `Vpi(f)` remains a measured
quantity until sufficient internal geometry is available.

## Stage 3: SQPM racetrack geometry

For each candidate mode pair, calculate the angle-dependent mismatch

`Delta_k(theta) = beta_SH(theta) - 2 beta_FF(theta)`

around an x-cut bend.  Find PDK-compliant integer-phase solutions satisfying

`Delta_phi_straight = Delta_k_y L0 = m pi`

and

`Delta_phi_arc = integral_0^pi Delta_k(theta) R dtheta = 2 N pi`,

with `R >= 80 um`.  Do not mechanically copy the 129.03-um radius and
81.8/245.4-um straight lengths of the 600-nm-film paper.

For every integer solution, calculate:

- round-trip phase and azimuthal/momentum selection;
- nonlinear overlap integrated around arcs and straights;
- expected FSR and mode density at both wavelengths;
- intrinsic-Q bounds from straight, bend, transition and sidewall loss;
- process-corner drift of the SQPM wavelength.

Gate: keep the three smallest integer-phase solutions whose phase error remains
inside the selected tolerance window across fabrication corners.

## Stage 4: double resonance and thermal control

Solve the cavity spectra near 1550 and 775 nm and search for mode pairs with

`delta_2 = omega_SH_mode - 2 omega_FF_mode`

inside the combined loaded-linewidth/tuning range.  Use a coupled-mode model
that includes pump depletion, intrinsic/external Q, nonlinear coupling,
photothermal shift, photorefraction, and the heater response.

Include a heater on every racetrack.  Simulate heater resistance, optical
absorption risk and thermal crosstalk in COMSOL Heat Transfer or an equivalent
model.

Gate: at least one process-corner solution reaches double resonance without
exceeding the heater-power or temperature limit, and has a practical detuning
capture range.

## Stage 5: dual-band ring coupling and extraction

Use 3D FDTD or EME to compare:

1. one pulley bus carrying both 1550 and 775 nm;
2. separate pump and SH buses; and
3. a pump bus plus first-mask visible-scatter collection.

Sweep ring/bus widths, gap, pulley angle and taper.  Enforce the 0.30-um LN
spacing rule.  Extract external Q, insertion loss, unwanted higher-order modes,
and process tolerance at both wavelengths.

Gate: choose the smallest topology that gives measurable SH output without
destroying pump buildup.  A common bus is preferred only if both external-Q
targets can be met.

## Stage 6: dual-wavelength travelling-wave PM

The standard 1550-nm PM is not a validated 775-nm device.  Design a custom
waveguide/electrode section for S3 and compute separately

- 1550- and 775-nm optical-metal loss;
- EO tensor projection and RF-field overlap;
- `VpiL_1550(f)` and `VpiL_775(f)`;
- RF phase velocity relative to both optical group indices;
- phase mismatch and modulation roll-off along the 9-mm-class line.

The two modulation indices are

`beta_lambda = pi V_RF / Vpi_lambda`

only after the appropriate frequency-dependent travelling-wave voltage is
defined.  Do not assume `beta_775 = 2 beta_1550`.

Gate: both wavelengths obtain a useful modulation index at 25 GHz with
acceptable loss and RF power.

## Stage 7: system-level spectrum model

Implement the sequence

`CW cavity SHG -> two-colour phase modulation`.

For each carrier,

`E_lambda(t) = A_lambda exp[i omega_lambda t + i beta_lambda sin(Omega t)]`.

Validate single-tone cavity SHG, Bessel-function PM spectra, energy
conservation, RF voltage calibration and the limiting case of zero SHG before
combining the blocks.  Add measured S1/S2 transfer functions when available.

Gate: reproduce the observed S1 `Vpi(f)` and S2 SHG/detuning curves before
claiming an S3 prediction.

## Stage 8: layout and signoff

Generate the three lanes plus:

- CPW thru/reflect/line calibration;
- passive PM-length optical reference;
- unheated and uncoupled racetrack controls;
- 1550/775 straight/bend/coupler tests;
- SiN and LN cutbacks;
- heater and metal sheet-resistance structures.

The supplied IPKISS/Latitudeda wrappers require their corresponding commercial
EDA environments; the current general Python environment does not provide
`ipkiss3` or `fnpcell`.  The PDK also contains no complete automatic DRC deck.
Final assembly and signoff must therefore run in the supported PDK environment
and through the foundry's official DRC service.

Gate: hierarchy, BlackBox keep-outs, connectivity, die boundary, probe/fibre
clearance and official DRC all pass.

## Measurement order

1. Measure passive 1550-nm loss and all RF calibration coupons.
2. Measure S1 `S11/S21`, `Vpi(f)` and EO comb spectra.
3. Identify S2 resonances, Q, heater tuning and visible SH scatter at low pump.
4. Calibrate 775-nm collection and measure SH power, detuning, stability and
   process-variant dependence.
5. Measure S3 SHG before applying RF.
6. Apply 25 GHz and verify equal tooth spacing, but separately fitted
   modulation indices, at both optical carriers.

Stop/go criterion: S3 is interpreted only after S1 and S2 independently agree
with their models.  A missing 775-nm comb must first be separated into SHG,
coupling, PM overlap, and RF-delivery failure modes.
