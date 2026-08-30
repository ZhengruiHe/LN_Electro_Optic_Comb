# Compact Three-Structure Layout Plan

## 1. Design decision

The baseline has been changed from a long PPLN waveguide cascade to a compact,
poling-free SHG architecture compatible with the geometric envelope of the
supplied x-cut TFLN-on-SiN PDK.

| ID | Structure | Reserved box | Area | Purpose |
|---|---|---:|---:|---|
| S1 | PDK 1550-nm travelling-wave PM EO comb | 10.8 mm x 0.65 mm | 7.02 mm2 | establish `Vpi(f)`, RF loss and EO-comb span using the lowest-risk PDK block |
| S2 | poling-free SQPM SHG racetrack DOE | 2.2 mm x 0.65 mm | 1.43 mm2 | find an integer-phase, 1550/775-nm double-resonant solution without domain inversion |
| S3 | SQPM-SHG first, then dual-band travelling-wave PM | 11.8 mm x 0.75 mm | 8.85 mm2 | generate 1550- and 775-nm CW carriers first, then create two EO combs with one RF tone |

The functional reservations total **17.30 mm2**, only **20.9%** of the PDK's
21.8 mm x 3.8 mm effective design area.  The remaining area is intentionally
kept for optical fanout, probe clearance, heaters, cutbacks, RF calibration,
alignment marks, and tolerance variants.  These are planning boxes, not
foundry-approved mask dimensions.

The key system ordering is:

`1550-nm CW -> poling-free SQPM racetrack -> 1550/775-nm CW -> travelling-wave phase modulation -> two EO combs`

This ordering is smaller and physically cleaner than sending an already broad
EO comb into a high-Q SHG ring.  Only the CW pump and one SH cavity mode must be
aligned before modulation.  The ring FSR is not required to equal the 25-GHz
RF frequency.

## 2. What the supplied PDK changes

The earlier all-LN 500-nm-film assumptions are obsolete.  The supplied process
uses a nominal 400-nm x-cut TFLN film over a nominal 300-nm SiN layer, separated
by about 400 nm of oxide, and provides a 21.8 mm x 3.8 mm effective design
area.  The PDK is aimed at 1550 nm.

The standard extraordinary-polarized Y-propagating phase-modulator BlackBox
has a 9.1-mm optical length and a protected bounding box of approximately
9.58 mm x 0.37 mm.  It is fixed geometry: length, signal width and electrode
gap are not exposed as parameters.  The manual quotes greater than 67 GHz and
less than 3 V.cm `VpiL` for the supplied MZI, but gives no corresponding
`Vpi(f)` metric for the PM.  The MZI number is therefore not used as a claimed
PM result.

No periodic-poling, SHG, 775-nm, visible coupler, nonlinear ring, or dual-band
component is supplied.  S2 and S3 are custom exploratory cells and require
written foundry permission.  The detailed compatibility gate is recorded in
`docs/pdk_compatibility_review.md`.

## 3. Why spontaneous/cyclic phase matching is plausible

For an x-cut LN ring, the propagation direction rotates relative to the
crystal optical axis.  An in-plane TE field therefore experiences an
azimuth-dependent effective index and effective nonlinear coefficient.  One
form used in the literature is

`d_eff(theta) = -d22 cos^3(theta) + 3 d31 cos^2(theta) sin(theta) + d33 sin^3(theta)`.

The sign and magnitude variation acts like a natural nonlinear grating.  A
racetrack adds straight sections whose length can be chosen so that the phase
accumulated in the arcs and straights repeats constructively.  In the notation
used by the racetrack literature, candidate geometries are selected from

`Delta_phi_straight = Delta_k L0 = m pi`

and

`Delta_phi_arc = integral_0^pi Delta_k(theta) R dtheta = 2 N pi`.

This is spontaneous quasi-phase matching (SQPM), closely related to cyclic
phase matching.  It removes the fabricated domain grating, not the optical
resonance requirements.

For useful cavity-enhanced SHG, the selected modes still need

`omega_SH ~= 2 omega_FF`

within their loaded linewidths, adequate nonlinear overlap, the appropriate
azimuthal selection rule, and usable coupling at both wavelengths.  Geometry
and a heater provide coarse and fine tuning, respectively.

## 4. Why the ring comes before the modulator

If an EO comb enters a high-Q SHG ring, only the teeth that coincide with ring
modes couple efficiently.  A 2024 CQPM microdisk experiment did demonstrate
broadband SHG, but it explicitly reported that wideband-light coupling was
inefficient and that the source repetition rate was not fully compatible with
the cavity resonances over a wide range.  Its high CW efficiency and its
wide wavelength tuning range must not be interpreted as simultaneous,
uniform conversion of every tooth of a 25-GHz EO comb.

Putting SHG first avoids that trap.  The ring produces two continuous-wave
carriers.  A downstream travelling-wave PM driven at angular frequency
`Omega` gives

`E_w(t) = A_w exp[i omega t + i beta_w sin(Omega t)]`

`E_2w(t) = A_2w exp[i 2 omega t + i beta_2w sin(Omega t)]`.

Both spectra have line spacing `Omega/(2 pi)`.  Their modulation indices need
not be equal; approximately, the shorter wavelength tends to accumulate more
phase for the same index perturbation, but the actual ratio must be computed
from the two optical modes and RF overlap.

## 5. S1: compact 1550-nm travelling-wave EO comb

### 5.1 First-tapeout baseline

Use the standard PDK `PM_e_Y_1550_LN` BlackBox, the provided 1550-nm LN edge
couplers, and a 50-ohm termination.  Reserve 10.8 mm x 0.65 mm, including the
9.58 mm x 0.37 mm protected PM box, two approximately 0.52-mm edge-coupler
allowances, and routing tolerance.

This is smaller and lower risk than immediately recreating a custom four-pass
device.  It is also diagnostic: the first measurements establish the actual
PM `Vpi(f)`, RF loss, optical loss and phase-modulation index available from
this foundry stack.

### 5.2 Expected, not guaranteed, comb scale

A published 1-cm single-pass TFLN phase modulator produced 15 lines at
24.95 GHz and 28 dBm.  The same paper's custom four-pass optical-recycling
device produced 47 lines, reduced RF power by about 15 times for an equivalent
span, and achieved an effective RF `Vpi` of 1.90 V at 24.95 GHz.  Those results
set a useful scale, but the PDK BlackBox is not that four-pass device.

Use the following first-pass goals:

| Quantity | Initial target | Status |
|---|---:|---|
| RF frequency | 25 GHz | chosen system point; sweep 10-40 GHz |
| impedance | 50 ohm +/- 5 ohm | simulation/measurement gate |
| return loss | `S11 < -10 dB` around 25 GHz | simulation/measurement gate |
| EO lines | at least 15 observable lines at the available amplifier power | engineering target, not PDK guarantee |
| `Vpi(f)` | measure, do not infer from MZI | mandatory result |

### 5.3 Low-`Vpi` upgrade

If S1 does not provide enough modulation index, the next cell is a custom
multi-pass travelling-wave PM.  It requires TE0/TE1 mode multiplexers and
precise RF-delay matching that are not supplied by this PDK.  A custom CPW is
also required because the PDK exposes no parameterized GSG travelling-wave
electrode PCell.  This upgrade follows measured S1 data rather than being
placed on the critical path of the first compact mask.

## 6. S2: poling-free SQPM racetrack DOE

### 6.1 Selected structure

Use an x-cut LN micro-racetrack rather than a released microdisk.  The
foundry process supports etched LN guides but does not advertise the pedestal
release and chemo-mechanical-polishing process used by the highest-Q microdisk
papers.  A ring is therefore the closest manufacturable implementation of the
CPM physics, although its Q and efficiency may be substantially lower.

The closest published process-compatible geometry used a 129.03-um outer
radius, 1-um top width and 0.8-um coupling gap, all above this PDK's geometric
minimums.  It used 600-nm LN with a 380-nm etch, so none of its phase-compensation
lengths may be copied to the present 400-nm stack.  The actual DOE radii and
straight lengths must be integer-phase solutions of the two equations above,
with `R >= 80 um`.  The 80, 129 and 160-um values used in the budget script are
only scale points.  Their approximate FSRs also do not constrain the downstream
25-GHz modulator.

### 6.2 Minimal three-cell DOE

After the optical-mode sweep, place three integer-phase racetrack solutions
inside 2.2 mm x 0.65 mm:

- the smallest valid solution above the 80-um bend rule;
- the solution nearest the literature-scale 129-um radius;
- one larger-radius solution that improves mode density or tolerance;
- a compact heater on every racetrack;
- a simulated pulley-bus gap, never below the PDK's 0.30-um LN spacing rule;
- one unheated passive reference and one straight dual-band guide in adjacent
  control space.

Do not freeze `w0` or the coupling geometry from a paper.  Solve the actual
400-nm-film, two-etch, oxide-clad, SiN-underlayer cross-section at both
1550 and 775 nm first.  The 775-nm mode will likely be higher order and is
especially sensitive to etch depth, sidewall angle and residual SiN.

### 6.3 Success criterion

S2 succeeds when at least one cell shows all of the following:

1. a reproducible 1550-nm resonance and measurable SH near 775 nm;
2. quadratic low-power SH scaling;
3. a heater-accessible double-resonance point;
4. separately calibrated pump coupling and SH collection; and
5. a measured thermal/photorefractive stability window.

Maximum paper efficiency is not the first-mask criterion.

## 7. S3: SHG first, dual-band EO modulation second

### 7.1 Optical path

`1550 input -> dual-band ring coupler -> CQPM/MPM ring -> residual 1550 + generated 775 -> dual-band travelling-wave PM -> common output -> external dichroic separation`

A common output and external dichroic filter minimize the first-chip area.  A
two-bus ring may be substituted if a single bus cannot simultaneously extract
the selected FF and SH modes.

### 7.2 Reserved area

Reserve 11.8 mm x 0.75 mm:

- up to 0.8 mm in length and 0.65 mm in width for the racetrack, heater and dual-band transition;
- about 9.6 mm for the PDK-scale travelling-wave electrode envelope;
- the remaining length for input/output tapering and isolation.

The standard PM BlackBox is qualified only at 1550 nm and cannot simply be
declared dual-band.  S3 therefore needs a custom waveguide under a custom M1
electrode with 1550- and 775-nm EO overlap, optical-metal loss, microwave
index, and RF attenuation all simulated.  The PDK layer rules permit a custom
drawing in principle; foundry approval is still required.

### 7.3 Why this produces two combs

The same 25-GHz RF tone modulates both optical carriers after SHG.  The two
comb centers are near 193 THz and 387 THz, while both tooth spacings are
25 GHz.  There is no requirement that either comb spacing equal the SHG ring
FSR because the sidebands are created outside the ring.

## 8. Compact die floorplan

Use three horizontal lanes inside the 21.8 mm x 3.8 mm effective region:

- top lane: S1 and RF calibration structures;
- middle lane: S3 with clear GSG probe access;
- bottom-left: S2 ring DOE;
- bottom-right: 1550/775 straight guides, bend/coupler DOE, heater references,
  CPW thru/reflect/line, and waveguide cutbacks.

Place optical interfaces on the left/right facets and RF probe pads so that a
probe body does not collide with lensed fibers.  The final pad pitch must be
set from the actual probe, not from the optical floorplan.

## 9. Simulation route

| Stage | Tool | Required outputs |
|---|---|---|
| PDK geometry reconstruction | KLayout + foundry values | LN1/LN2 cross-sections, layer booleans, hierarchy |
| two-colour modes | COMSOL Wave Optics or Lumerical MODE | anisotropic modes, `n_eff`, `n_g`, loss proxy, overlap, width/thickness corners |
| ring double resonance | MODE/FEM + Python coupled-mode model | mode families, azimuthal selection, `2 omega` mismatch, FSRs, heater range |
| dual-band coupler | 3D FDTD or EME | external Q at both wavelengths, parasitic-mode content, tolerance |
| travelling-wave electrode | HFSS 2D/Q3D then 3D HFSS | `Z0`, RF index/loss, `S11/S21`, pad transition, EO overlap, optical-metal loss |
| two-colour PM | optical modes + electrostatic/RF field integration | `VpiL` and modulation index at both wavelengths |
| system spectrum | Python coupled-mode/Bessel model | CW SHG, two EO combs, power conservation, thermal detuning |
| layout signoff | PDK-native EDA + official DRC | hierarchy, BlackBox keep-outs, connectivity and foundry DRC |

HFSS is appropriate for the electrical line, pads and microwave-optical
velocity matching.  It cannot solve the anisotropic optical modes or nonlinear
SHG cavity by itself.

## 10. Mandatory measurement controls

### S1

- RF thru/reflect/line coupons and de-embedded `S11/S21`;
- optical cutbacks and a passive route matching the PM path;
- `Vpi(f)` from 10-40 GHz;
- EO spectrum versus RF power and optical wavelength.

### S2

- ring transmission, loaded/intrinsic Q and heater tuning coefficient;
- visible-scatter camera check followed by calibrated 775-nm extraction;
- SH power versus detuning, pump power and heater power;
- long-term drift, thermal bistability and photorefraction check.

### S3

- spectra before and after the travelling-wave PM at both wavelengths;
- exact 25-GHz tooth spacing at 1550 and 775 nm;
- RF power-to-modulation-index curves for both colours;
- a comparison against independently measured S1 and S2 transfer functions.

## 11. Mask-release gates

1. **Foundry gate:** written approval for custom 775-nm LN routing and visible
   collection; exact etch/cladding/metal stack received.
2. **Mode gate:** one FF/SH mode pair has adequate nonlinear overlap and a
   double-resonance mismatch reachable by heater tuning across process corners.
3. **Coupling gate:** one/two-bus design gives usable external Q at both
   wavelengths without violating 0.30-um LN spacing.
4. **RF gate:** S1 or the custom S3 line meets the 50-ohm, reflection, loss and
   velocity-match targets at 25 GHz.
5. **Dual-band PM gate:** computed `VpiL` and metal loss are acceptable at both
   1550 and 775 nm.
6. **Testability gate:** every critical block has an independent passive or
   electrical control.
7. **DRC gate:** official foundry DRC passes.  The PDK package does not contain
   a complete automatic DRC deck, so local geometry checks are not signoff.

## 12. Evidence and limits

### Source-backed results

1. Ke Zhang et al., *A power-efficient integrated lithium niobate
   electro-optic comb generator*, Communications Physics 6, 17 (2023),
   https://www.nature.com/articles/s42005-023-01137-9 .  A 1-cm, four-pass
   travelling-wave device produced 47 lines at 25 GHz and 28 dBm with an
   effective RF `Vpi` of 1.90 V at 24.95 GHz; its single-pass control produced
   15 lines under the same drive.
2. Tingge Yuan et al., *Chip-scale spontaneous quasi-phase matched second
   harmonic generation in a micro-racetrack resonator*, Science China Physics,
   Mechanics & Astronomy 66, 284211 (2023),
   https://doi.org/10.1007/s11433-023-2145-6 .  The x-cut 600-nm-LN device used
   a 129.03-um outer radius and 81.8/245.4-um straight sections for 37th/111th
   order SQPM.  Its core boxes were about 0.34 x 0.258 mm and
   0.503 x 0.258 mm, while measured normalized on-chip efficiencies were only
   `1.01e-4/W` and `0.43e-4/W`.  It is the closest geometric reference and also
   shows the efficiency cost of high-order SQPM.
3. Jintian Lin et al., *Broadband Quasi-Phase-Matched Harmonic Generation in
   an On-Chip Monocrystalline Lithium Niobate Microdisk Resonator*, Physical
   Review Letters 122, 173903 (2019),
   https://doi.org/10.1103/PhysRevLett.122.173903 .  An x-cut approximately
   30-um-diameter microdisk used natural cyclic QPM and reported normalized SHG
   efficiency up to 9.9%/mW without domain engineering.
4. Jiefu Zhu et al., *Broadband second-harmonic generation in thin-film
   lithium niobate microdisk via cyclic quasi-phase matching*, Chinese Optics
   Letters 22, 031903 (2024),
   https://doi.org/10.3788/COL202422.031903 .  A roughly 100-um-diameter,
   550-nm-thick x-cut microdisk reported loaded Q values of `3.32e7` near the
   pump and `2.83e6` near the SH, 15.2%/mW CW normalized efficiency, and SHG
   over a pump scan wider than 100 nm.  It also documents the much lower
   conversion of broadband sources and the comb-spacing/cavity-resonance
   mismatch limitation.
5. Rui Luo et al., *Optical Parametric Generation in a Lithium Niobate
   Microring with Modal Phase Matching*, Physical Review Applied 11, 034026
   (2019), https://doi.org/10.1103/PhysRevApplied.11.034026 .  A 50-um-radius
   Z-cut ring used TM00/TM20 modal phase matching and reported 1500%/W SHG;
   it is an alternative mechanism, not a geometry to copy into this x-cut
   400-nm PDK.
6. Xingze Song et al., *Broadband birefringence phase-matched second-harmonic
   generation in a slightly curved lithium niobate-on-insulator waveguide*,
   Applied Optics 65, 1511-1515 (2026),
   https://doi.org/10.1364/AO.586578 .  A poling-free x-cut curved-waveguide
   method exceeded 100-nm bandwidth but reported only 1.38%/(W.cm2)
   normalized efficiency; it is a broadband pulsed-light fallback rather than
   the preferred CW resonant converter.

### Engineering inferences

- The final racetrack radii are integer-phase solutions from the actual PDK
  stack; the 80/129/160-um values are only scale points.
- The 17.30-mm2 total reservation is an area budget derived from the supplied
  BlackBox bounds and routing margins, not a measured paper footprint.
- S3's common-bus dual-band modulation is physically plausible, but no cited
  paper or supplied PDK component demonstrates this exact monolithic sequence.

### Unknowns

- Exact LN etch depths, sidewall angle, top-cladding thickness and M1 RF
  properties are missing.
- The foundry has not qualified Q, loss, photorefraction, couplers, or power
  handling at 775 nm.
- The final double-resonant mode pair, ring width, bus gap, heater range and
  dual-band `VpiL` cannot be specified before simulation.

## 13. Current conclusion

The user's proposed non-PPLN direction is physically sound and materially
reduces the layout area.  The most defensible compact architecture is not to
modulate inside the SHG ring and not to send a broad EO comb into it.  It is to
generate 1550/775-nm CW light in a small poling-free resonator and then use a
separate travelling-wave phase modulator.  This removes the RF-to-ring-FSR
constraint while preserving identical comb spacing at the two optical
carriers.

The design is ready for cross-section reconstruction and optical/RF
simulation, but S2 and S3 remain outside the qualified PDK scope and are not
ready for fabrication.
