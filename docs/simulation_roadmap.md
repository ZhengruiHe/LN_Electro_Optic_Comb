# Simulation Roadmap

## Baseline concept

The first version is a non-resonant, all-LN platform with three independently
measurable structures:

1. a travelling-wave optical-recycling phase modulator for a 1550-nm EO comb;
2. a straight PPLN waveguide array for 1550-to-775-nm SHG;
3. a monolithic cascade of the first two functions.

The integrated device first produces the fundamental EO comb and then performs
SHG/SFG in PPLN. It therefore does **not** require
`FSR(omega) = FSR(2 omega)` or any optical-ring resonance. The second-harmonic
comb inherits the RF line spacing from the fundamental comb.

The numerical values in `docs/three_structure_layout_plan.md` are preliminary
targets, not released mask dimensions. Foundry rules and simulation gates below
take precedence.

## Stage 0: freeze the common fabrication stack

Record and confirm:

- x-cut MgO:LNOI availability and crystal-axis marker;
- LN thickness, etch depth, sidewall angle, BOX, and top cladding;
- minimum waveguide, gap, domain, and alignment tolerances;
- poling sequence and maximum high-voltage field;
- permitted metal thickness and minimum metal-to-waveguide gap;
- RF probe pitch and whether a deposited 50-ohm termination is available.

Deliverable: a signed-off stack table and foundry-rule checklist.

Gate: do not release a final QPM period or CPW geometry before this step.

## Stage 1: dual-wavelength optical modes

Use Lumerical MODE or COMSOL Wave Optics. Sweep waveguide top width, etch
depth, film thickness, and sidewall angle at 1550 nm and 775 nm.

Extract:

- `n_eff`, `n_g`, field profiles, polarization, and higher-order modes;
- propagation and bend-loss estimates;
- group-velocity mismatch and dispersion;
- full-tensor EO and nonlinear overlap integrals;
- width/thickness sensitivities of the QPM wavelength.

Gate: select a common cross-section that supports low-loss 1550-nm routing and
a fabricable, high-overlap SHG mode pair.

## Stage 2: travelling-wave electrode

Use HFSS 2D/quasi-static or COMSOL Electrostatics first, then a 3D HFSS model
of the complete 10-mm electrode, pads, and termination transition.

Optimize:

- 50-ohm characteristic impedance;
- microwave effective index matched to optical group index;
- conductor and dielectric RF loss;
- `S11`, `S21`, and EO overlap from 5 to 40 GHz;
- metal-induced optical loss and RF heating;
- phase accumulation through the four-pass optical recycling path.

Gate targets at 25 GHz: `|Z0 - 50 ohm| <= 5 ohm`, `S11 < -10 dB`, effective
`Vpi <= 2.5 V`, and usable EO response through at least 35 GHz. The last two
targets must be validated with the complete recycling path, not just a CPW
cross-section.

## Stage 3: PPLN QPM and bandwidth

Calculate

`Delta k = beta_2w - 2 beta_w`

and

`Lambda = 2*pi/abs(Delta k)`.

Then sweep poling period, duty cycle, waveguide width, film thickness,
temperature, and device length. Use a coupled-wave or split-step solver to
include pump depletion and multiple comb lines.

Two different optimizations are required:

- the stand-alone SHG array uses a nominal 4-mm uniform PPLN section for high
  conversion efficiency;
- the integrated comb converter uses a nominal 2-mm uniform PPLN section for
  bandwidth, with a 4-mm weakly chirped PPLN variant kept in the DOE.

Gate: the simulated QPM FWHM of the integrated path must cover the chosen EO
comb span with at least 20% spectral margin.

## Stage 4: system-level comb model

For a phase-modulated pump,

`E_w(t) = E0*exp[i*w0*t + i*beta*sin(Omega*t)]`.

In the ideal instantaneous, phase-matched limit,

`P_2w(t) proportional to E_w(t)^2`

and the second harmonic has modulation index `2*beta`, while its line spacing
remains `Omega/(2*pi)`. The numerical model must add finite QPM bandwidth,
loss, dispersion, RF phase error, and pump depletion.

Gate: reproduce the single-tone SHG limit, the Bessel-function EO-comb limit,
and power conservation before running the full multi-line cascade.

## Stage 5: couplers, tapers, and layout verification

Use EME/FDTD for:

- 1550-nm edge couplers and tapers;
- the 1550-nm EO-waveguide to dual-wavelength PPLN transition;
- optional 1550/775-nm output demultiplexer;
- TE0/TE1 mode multiplexers, crossing, and loop-back structures.

Generate the three floorplan blocks plus controls: single-pass PM, unpoled
waveguides, CPW thru/open/short, and waveguide cutbacks. Run connectivity,
minimum-spacing, density, and official foundry DRC.

Gate: all critical structures must have independent controls and accessible
optical/RF ports.

## Stage 6: experimental validation order

1. Measure passive loss, mode-multiplexer loss, and crossing loss.
2. Measure CPW `S11/S21`, microwave index, and loss.
3. Measure single-pass and four-pass `Vpi(f)` and EO comb spectra.
4. Map SHG versus wavelength, temperature, period, width, length, and power.
5. Measure the integrated fundamental and second-harmonic combs.
6. Compare line spacing, line-to-line power, conversion, and phase coherence.

Stop/go criterion for the integrated structure: independently measured EO and
PPLN blocks must agree with their respective models before attributing any
integrated-device failure to the cascade physics.

## Later resonant branch

Only after the non-resonant blocks work should a doubly resonant optical-ring
version be reconsidered. Its central `omega/2omega` double resonance, two-band
dispersion, coupling, thermal tuning, and RF-to-fundamental-FSR condition form
a separate risk set and are not part of the first tapeout baseline.

