# Compact Three-Structure Layout Plan

## 1. Scope and design decision

All three devices are custom. The supplied PDK is used only for the material
stack, drawing layers, minimum geometry, die boundary and eventual foundry
DRC. No phase-modulator or resonator BlackBox is used.

The selected integrated order is

`1550-nm CW -> poling-free SQPM racetrack -> separate 1550/775 rails -> common travelling-wave RF electrode -> two EO combs`.

Only the SHG block is resonant. EO modulation occurs after the ring, so the
ring FSR does not need to equal 25 GHz. The ring must instead provide one
fundamental/second-harmonic double-resonant mode pair.

| ID | Custom structure | Rough reserved box | Rough area | Main purpose |
|---|---|---:|---:|---|
| S1 | four-pass 1550-nm travelling-wave phase modulator | 12.0 mm x 0.9 mm | 10.80 mm2 | low-`Vpi` 25-GHz EO comb |
| S2 | three-cell poling-free SQPM/MPM racetrack DOE | 2.2 mm x 0.7 mm | 1.54 mm2 | compact 1550-to-775-nm CW SHG |
| S3 | one SHG racetrack plus dual-rail common-electrode PM | 12.5 mm x 1.0 mm | 12.50 mm2 | two combs driven by the same RF tone |

The functional reservations total about **24.84 mm2**, or **30.0%** of the
21.8 mm x 3.8 mm effective area. The three lane widths total 2.6 mm, leaving
about 1.2 mm of die width for inter-lane clearance and controls. These are
engineering envelopes, not final GDS dimensions.

## 2. PDK constraints applied to every custom cell

The currently extracted planning constraints are:

- x-cut TFLN on SiN, nominally 400-nm LN, 300-nm SiN and about 400-nm
  interlayer oxide;
- 21.8 mm x 3.8 mm effective design region;
- LN minimum line/space 0.30/0.30 um and minimum bend radius 80 um;
- M1 minimum line/space 2/3 um and RF electrode gap no smaller than 3 um.

Exact LN etch depths, sidewall angle, final top cladding and M1
material/thickness/conductivity are still missing. Consequently, no waveguide
width, coupling gap or CPW cross-section in this document is frozen.

## 3. S1: custom low-Vpi travelling-wave EO comb

### 3.1 Physical structure

Use a 10-mm, nominally 50-ohm GSG travelling-wave electrode with a 50-ohm
termination. A custom 1550-nm waveguide traverses the same active RF region
four times. TE0/TE1 adiabatic mode multiplexers, crossings and loopback delay
lines recycle the optical wave while the microwave propagates only once.

The important point is coherent accumulation: each optical pass must meet the
correct RF phase. The loop lengths are therefore simulation variables tied to
the 25-GHz RF period and the measured optical group delay, not simply the
shortest possible connections.

A useful first RF sweep is signal width about 43 um and electrode-waveguide
gap about 5.5 um, following the literature device below. These values are not
copied to mask because the present M1 stack and cladding differ. Sweep signal
width, ground width, gap, metal thickness and waveguide offset subject to the
PDK limits.

### 3.2 Why four passes

For phase modulation, the comb modulation index is approximately

`beta = pi V_peak / Vpi(f)`,

and the ideal line powers follow `P_n = P_0 J_n(beta)^2`. Four coherent passes
increase interaction length without requiring a 40-mm RF line, whose
microwave loss would be excessive. The price is additional optical loss,
mode-multiplexer crosstalk and strict loop-delay control.

### 3.3 Rough size and first design gate

Reserve 12.0 mm x 0.9 mm:

- 10.0 mm for the active electrode;
- about 1.0 mm at each end in total for RF taper/termination and optical
  transitions;
- 0.9 mm width for GSG pads/line, four optical tracks, 80-um-radius loopbacks
  and routing clearance.

The first gate at 25 GHz is `45-55 ohm`, `S11 < -10 dB`, acceptable `S21`,
coherent four-pass addition, and predicted effective `Vpi <= 2.5 V`. The
`Vpi` value is a target; it cannot be claimed before RF/EO overlap simulation
and measurement.

## 4. S2: custom poling-free SQPM racetrack SHG

### 4.1 Physical mechanism and structure

In an x-cut LN bend, the propagation direction rotates relative to the crystal
axis. Both the effective index and effective nonlinear coefficient can vary
periodically around the path. A racetrack can therefore reset or compensate
the generated SH phase without periodic domain inversion.

For each candidate FF/SH mode pair, solve

`Delta_k(theta) = beta_2w(theta) - 2 beta_w(theta)`,

then choose a straight length `L0` and bend radius `R` satisfying approximate
integer-phase conditions

`Delta_k_y L0 = m pi`

and

`integral_0^pi Delta_k(theta) R dtheta = 2 N pi`.

The racetrack contains the LN ring, one pump bus and preferably a separate
775-nm pickup bus, a metal heater and optical reference ports. A one-bus
variant is retained only if it gives useful external Q at both wavelengths.

### 4.2 DOE and rough size

Use three cells rather than blindly copying three radii:

1. the smallest robust SQPM integer solution from the actual 400-nm stack;
2. the next integer solution with lower phase-error sensitivity;
3. a modal-phase-matching fallback using a different SH mode family.

Start numerical searches around `R = 120-160 um`, while enforcing `R >= 80
um`; choose the final radius and straight length from the phase integral and
double-resonance calculation. A single ring core with buses and heater is
expected to fit roughly 0.8 mm x 0.5 mm. Reserve 2.2 mm x 0.7 mm for three
compact cells and shared fanout, about 1.54 mm2 total.

The first-mask success criterion is measurable 775-nm output with quadratic
low-power scaling and heater-accessible double resonance, not record
conversion efficiency.

## 5. S3: SHG plus a dual-rail common-electrode modulator

### 5.1 Physical structure

Use one S2-derived racetrack to generate two CW carriers. Separate the
residual 1550-nm pump and generated 775-nm light with two ring buses or a custom
dual-band extraction section. Route them into two independent optical rails:

- a 1550-nm rail along one side of the GSG signal electrode;
- a 775-nm rail along the other side.

Both rails see the same 25-GHz travelling RF field, but each waveguide width,
etch choice and electrode offset is optimized independently. The transverse
field reverses sign on opposite sides of the signal electrode; this reverses
the relative phase modulation sign, but not the comb spacing or ideal line
power envelope.

For each colour,

`E_lambda(t) = A_lambda exp[i omega_lambda t + i beta_lambda sin(Omega t)]`.

Both combs have spacing `Omega/(2 pi) = 25 GHz`, while generally
`beta_775 != 2 beta_1550`. No equality between the ring FSR and RF frequency
is required.

### 5.2 Rough size and risk

Reserve 12.5 mm x 1.0 mm, about 12.50 mm2:

- up to 0.8 mm x 0.6 mm for the racetrack, buses and heater;
- 10.0 mm for the common travelling-wave electrode;
- remaining length for WDM/tapers, RF pads and two optical outputs.

This is the highest-risk cell. It requires low-loss 775-nm routing, usable EO
overlap at both wavelengths, acceptable 775-nm metal absorption and one RF
phase velocity that is useful relative to both optical group indices. A
single-pass dual-rail PM is selected for the first design; a four-pass
dual-colour recycler would add unverified 775-nm multiplexers and too much
system risk.

## 6. Simulation workflow

Do not simulate the entire chip in one full-3D model.

| Level | Recommended tool | Model and required outputs |
|---|---|---|
| optical cross-section | Ansys Lumerical MODE FDE or COMSOL Wave Optics | anisotropic modes at 1550/775 nm, `n_eff`, `n_g`, confinement, bend loss, dispersion, EO/nonlinear overlaps and process corners |
| RF cross-section | Ansys Q3D or HFSS 2D Extractor | `Z0(f)`, `n_RF(f)`, conductor/dielectric loss and field distribution versus GSG geometry |
| RF 3D details | HFSS | probe pad, taper, 10-mm line, termination, `S11/S21`, current density and launch discontinuity |
| EO overlap | COMSOL Electrostatics plus optical modes, or field export/integration | `d(n_eff)/dV`, `VpiL`, metal loss and travelling-wave `Vpi(f)` for each rail |
| optical transitions | EME/FDTD | TE0/TE1 multiplexer, crossing, loopback, ring/bus coupler, WDM/taper insertion loss and crosstalk |
| SQPM cavity | mode solver plus Python | angle-dependent `Delta_k`, nonlinear phase integral, cavity modes, double resonance, loaded/intrinsic/external Q and SHG coupled-mode response |
| thermal tuning | COMSOL Heat Transfer | heater efficiency, temperature rise, crosstalk and resonance capture range |
| system spectrum | Python | Bessel combs, measured RF voltage, SHG detuning, losses and separate `beta_1550/beta_775` |
| layout/signoff | PDK-native EDA plus official DRC | custom hierarchy, layer booleans, connectivity, die/probe/facet clearances and official DRC |

HFSS is therefore the correct tool for the microwave electrode, pads and
termination, but it is not the tool for anisotropic optical modes or nonlinear
SHG. A practical minimum toolchain is MODE + Q3D/HFSS + Python; COMSOL is most
useful for EO overlap and thermal tuning.

For the travelling-wave response use

`H(f) proportional to [1-exp(-(alpha_RF+j Delta_beta)L)] /(alpha_RF+j Delta_beta)`,

with `Delta_beta = 2 pi f (n_RF-n_g)/c`. For S3 evaluate the expression once
with `n_g_1550` and once with `n_g_775`.

For the nonlinear ring, numerically integrate

`A_2w proportional to integral kappa(s) exp[i integral_0^s Delta_k(s') ds'] ds`

around the complete racetrack, then place the result in a temporal
coupled-mode model including pump/SH detuning, intrinsic/external Q, heater,
photothermal shift and pump depletion.

## 7. Execution order and measurable gates

1. Obtain the missing foundry cross-section and M1 parameters.
2. Solve 1550/775-nm optical modes and select manufacturable waveguide families.
3. Close the S1 RF cross-section and 3D launch, then calculate `Vpi(f)`.
4. Close the four-pass optical multiplexer, loss and delay budget.
5. Search S2 integer-phase and double-resonance solutions; simulate buses and
   heater only for the best three cells.
6. Combine one selected ring with the S3 dual-rail electrode model.
7. Generate custom GDS only after the optical, RF, thermal and testability
   gates pass; submit it to official foundry DRC.

Mandatory controls are RF thru/reflect/line standards, single-pass and passive
S1 routes, optical cutbacks, unheated/uncoupled rings, 1550/775 straight/bend
tests, heater resistance cells and separate access to the S3 SHG output before
RF modulation.

## 8. Evidence, engineering inference and unknowns

### Source-backed results

- [Zhang et al., Communications Physics 6, 17
  (2023)](https://www.nature.com/articles/s42005-023-01137-9) demonstrated a
  four-pass, 1-cm travelling-wave TFLN EO comb: 47 lines at 25 GHz and 28 dBm,
  effective RF `Vpi = 1.90 V` near 24.95 GHz, versus 15 lines for its
  single-pass control. It supports the S1 architecture, not the exact area
  estimate or transfer of its electrode dimensions to this process.
- [Yuan et al., Science China Physics, Mechanics & Astronomy 66, 284211
  (2023)](https://doi.org/10.1007/s11433-023-2145-6) demonstrated poling-free
  SQPM in x-cut 600-nm-LN racetracks with 129.03-um outer radius and
  81.8/245.4-um straight sections. The measured normalized on-chip
  efficiencies were only `1.01e-4/W` and `0.43e-4/W`, so it is a geometry and
  mechanism reference rather than an efficiency promise.
- [Zhu et al., Chinese Optics Letters 22, 031903
  (2024)](https://doi.org/10.3788/COL202422.031903) reported efficient CQPM in
  a roughly 100-um x-cut microdisk, but its CMP/suspended-disk process and
  tapered-fibre coupling are not transferable to this PDK.
- [Luo et al., Physical Review Applied 11, 034026
  (2019)](https://doi.org/10.1103/PhysRevApplied.11.034026) demonstrated a
  modal-phase-matched LN ring and motivates the third S2 fallback cell; its
  Z-cut 600-nm geometry is not copied.

The solver roles above follow the official [Ansys MODE FDE
description](https://optics.ansys.com/hc/en-us/articles/360034917233-MODE-Finite-Difference-Eigenmode-FDE-solver-introduction),
[Ansys Q3D impedance documentation](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v252/en/Subsystems/Q3DExtractor/Content/Q3D/CharacteristicImpedance.htm),
and [COMSOL anisotropic-waveguide
example](https://www.comsol.com/model/optically-anisotropic-waveguide-57481).

### Engineering estimates

- The 12.0 x 0.9, 2.2 x 0.7 and 12.5 x 1.0 mm boxes are routing reserves based
  on a 10-mm active electrode, PDK bend radius and local controls. They are not
  measured paper footprints.
- `Vpi <= 2.5 V`, `45-55 ohm` and `S11 < -10 dB` are first-pass design gates,
  not demonstrated device results.
- A dual-rail common electrode is physically plausible and reduces optical
  compromise, but the exact S3 combination has not been demonstrated by the
  cited papers.

### Unknowns before fabrication

- exact LN etches, sidewall, top oxide and M1 RF material parameters;
- foundry permission and propagation/coupling loss at 775 nm;
- the final SQPM mode pair, radius, straight length, double-resonance tolerance
  and heater capture range;
- four-pass loop delay, accumulated loss and fabrication sensitivity;
- `VpiL_1550`, `VpiL_775`, RF power handling and optical-metal loss.

## 9. Current feasibility conclusion

S1 is medium-to-high feasibility because a closely related four-pass device
has been demonstrated, although every waveguide and RF dimension must be
re-solved for this stack. S2 is medium feasibility for observing SHG but low
confidence for high conversion efficiency on a first mask. S3 is physically
consistent and compact, but remains the highest-risk structure because it
combines 775-nm extraction, two EO overlaps, microwave velocity matching and
thermal resonance control.

The next real calculation should be the anisotropic 1550/775-nm mode sweep and
the custom GSG cross-section sweep, not a whole-chip HFSS model and not a final
GDS layout.
