# Three-Structure Layout Plan

## 1. Decision summary

The recommended first tapeout contains three non-resonant structures on a
common x-cut MgO:TFLN platform:

| ID | Structure | Main purpose | Reserved box | Reserved area |
|---|---|---|---:|---:|
| S1 | Four-pass travelling-wave EO phase modulator | Low-`Vpi`, 25-GHz fundamental EO comb | 13.0 mm x 1.4 mm | 18.2 mm^2 |
| S2 | Uniform/chirped PPLN-SHG DOE array | Independently verify QPM, efficiency, bandwidth, and tolerance | 7.0 mm x 1.6 mm | 11.2 mm^2 |
| S3 | Travelling-wave EO comb + broadband PPLN cascade | Simultaneous 1550-nm and 775-nm comb output | 18.0 mm x 1.6 mm | 28.8 mm^2 |

The three functional reservations total 58.2 mm^2. A preliminary
20 mm x 6 mm die is 120 mm^2, leaving 61.8 mm^2 for facet margins, dicing
streets, RF calibration, passive controls, alignment marks, heaters, and metal
density fill. These are **floorplan reservations**, not measured paper
footprints or final foundry-approved mask dimensions.

The baseline die drawing is in `docs/three_structure_floorplan.svg`, and the
machine-readable targets are in `models/system/design_targets.json`.

## 2. Why these three structures

The order separates the two major technical risks:

- S1 answers whether the travelling-wave electrode and optical recycling
  actually achieve low RF `Vpi` and a broad EO comb.
- S2 answers whether domain inversion, dual-wavelength modes, and QPM work on
  the selected wafer stack.
- S3 tests only the final cascade after S1 and S2 can be diagnosed separately.

This is lower risk than starting with a doubly resonant optical ring. In S3,
there is no optical FSR to match. A single RF drive at `f_RF` produces the
fundamental comb and the PPLN maps it to a second-harmonic comb with the same
line spacing.

If the fundamental field is

`E_w(t) = E0 exp[i*w0*t + i*beta sin(Omega*t)]`,

then ideal SHG gives

`E_2w(t) proportional to E_w(t)^2`

`= E0^2 exp[i*2*w0*t + i*2*beta sin(Omega*t)]`.

Therefore the SH comb is centered at `2*w0`, has line spacing `Omega`, and has
an ideal phase-modulation index of `2*beta`. Equivalently, all input tooth pairs
`(n,m)` contribute through SHG/SFG at
`2*w0 + (n+m)*Omega`; the allowed output grid is still spaced by `Omega`, not
`2*Omega`.

## 3. Common fabrication platform

### 3.1 Provisional stack

| Parameter | Preliminary value | Reason/status |
|---|---:|---|
| Substrate | x-cut MgO:LNOI | x-cut supports strong in-plane `r33`; MgO is preferred for visible-power robustness |
| LN film | 500 nm | Common value demonstrated in the EO-comb and PPLN literature |
| Ridge etch | 300 nm | Leaves a 200-nm slab and matches a published 1550/775-nm PPLN geometry |
| BOX | 2.0 um SiO2 | Literature baseline |
| Top cladding | 0.7-1.0 um SiO2 | Optical isolation from metal; final value is an RF/optical co-design variable |
| Waveguide axis | Along LN y axis | Makes the dominant TE optical field and lateral RF field use the large EO tensor element |
| EO waveguide width | Start at 1.0-1.3 um | Must be selected by the dual-wavelength optical-mode sweep |
| PPLN waveguide width | 1.20 um nominal | Published 500-nm-film PPLN starting point; sweep required |
| First-order QPM period | `Lambda0 = 3.70 um` nominal | Published starting point only; recalculate from the actual stack |
| Domain duty cycle | 50% target | Maximum first-order Fourier component; include fabrication error |
| Electrode metal | Au, 1.0-1.6 um target | Thicker than the 520-nm literature demonstration to reduce conductor loss if the process permits |

The final common stack cannot be frozen from literature alone. Film thickness,
etch depth, sidewall angle, cladding, MgO content, and temperature all shift
the effective indices and QPM period.

### 3.2 Mandatory fabrication sequence

1. Pattern temporary Cr/Au poling electrodes and invert the domains.
2. Inspect domain fidelity using SH microscopy or a qualified equivalent.
3. Remove the poling metal.
4. Align and etch the optical waveguides inside the poled stripes.
5. Deposit the top cladding.
6. Form the permanent travelling-wave electrodes, probe pads, and optional
   heaters.

The PPLN interaction region must not retain the temporary poling fingers. The
permanent EO electrode belongs in the unpoled modulation region. Keeping the
functions spatially separated avoids RF-metal loss in the PPLN section and
avoids EO sign cancellation across inverted domains.

## 4. S1: low-Vpi travelling-wave EO comb

### 4.1 Chosen architecture

Use a non-resonant 10-mm, 50-ohm GSG travelling-wave electrode with four-pass
optical recycling. TE0/TE1 mode multiplexers, one crossing, and loop-back
waveguides make the optical field coherently traverse the same RF interaction
region four times. The RF remains a normal two-port travelling wave and is
terminated in 50 ohm.

This architecture is selected because a published TFLN device demonstrated:

- a 1-cm travelling-wave electrode;
- 5.5-um electrode gap and 43-um signal width;
- RF index about 2.30 and optical group index about 2.27;
- effective RF `Vpi = 1.90 V` at 24.95 GHz;
- 47 comb lines at 25 GHz using 28 dBm RF drive;
- approximately 4 dB total on-chip loss for its 12-cm recycled optical path.

### 4.2 Preliminary layout values

| Item | Baseline | Sweep/control |
|---|---:|---|
| Physical electrode length | 10.0 mm | 5, 7.5, and 10 mm CPW coupons |
| Optical passes | 4 | single-pass and double-pass controls |
| RF center frequency | 25 GHz | characterize 5-40 GHz |
| GSG signal width | 43 um starting point | RF optimizer changes this |
| Waveguide-to-electrode gap | 5.5 um starting point | 4.5, 5.5, 6.5 um coupons if metal loss allows |
| Target impedance | 50 ohm | +/-5 ohm acceptance at 25 GHz |
| Target effective `Vpi` | <=2.5 V at 25 GHz | stretch goal <=2.0 V |
| Target EO bandwidth | useful through >=35 GHz | characterize full 5-40 GHz band |
| Core footprint | about 12 mm x 0.8 mm | literature-scale reference |
| Reserved mask box | 13.0 mm x 1.4 mm | includes RF pads, optical bends, and keep-outs |

The 5.5-um gap is not automatically optimal. Reducing it lowers `Vpi` but can
increase metal absorption and change impedance. Increasing electrode length
also lowers `Vpi` but increases RF attenuation and the velocity-mismatch
penalty. HFSS and the measured passive controls decide the final point.

### 4.3 Expected comb behavior

For a pure phase modulator, the line powers follow Bessel functions rather than
forming a flat top. S1 is intended to prove low-power broad comb generation,
not flatness. If a flat-top spectrum becomes a requirement, add a travelling-
wave MZM before two phase modulators in a later mask; a 2025 TFLN experiment
using one intensity modulator and two phase modulators reported over 70 lines
within a 10-dB bandwidth and 5-25-GHz repetition-rate tuning.

## 5. S2: PPLN quasi-phase-matched SHG array

### 5.1 Chosen architecture

Use straight, non-resonant ridge waveguides with removable poling electrodes.
The nominal interaction is type-0-like, using the strongest practical tensor
combination allowed by the selected x-cut geometry. Start from a 1.20-um top
width, 300-nm ridge, 200-nm slab, 700-nm cladding, and 3.70-um period, then
recalculate with the full anisotropic LN model.

The literature gives two important scale references:

- a 4-mm-long, 75-um-wide poled region accommodated three ridge waveguides,
  used an approximately 4.1-um period, reached 2600% W^-1 cm^-2 normalized
  efficiency, and produced 117 mW at 775 nm from 220 mW on-chip pump;
- another 500-nm-film geometry used a 1.20-um-wide waveguide, 300-nm ridge,
  200-nm slab, 700-nm cladding, and a 3.70-um period.

The different periods are not contradictory; they arise from different film,
waveguide, cladding, and modal dispersion.

### 5.2 DOE inside the 7.0 mm x 1.6 mm block

Use 11 lanes at a nominal 120-um pitch:

- 7 period lanes:
  `Lambda = Lambda0 + {-0.15,-0.10,-0.05,0,+0.05,+0.10,+0.15} um`,
  with nominal width and 4-mm length;
- 2 width lanes: `w = w0 - 0.10 um` and `w0 + 0.10 um`, with nominal period
  and 4-mm length;
- 2 length lanes: 2 mm and 6 mm, with nominal width and period;
- place an unpoled 4-mm reference in the adjacent unused floorplan area.

This is a center-plus-one-factor sweep, not a full factorial. Before mask
release, set `Lambda0` from the measured wafer stack and eigenmode simulation;
do not blindly retain 3.70 um.

### 5.3 Length and bandwidth choice

For uniform QPM, an approximate intensity FWHM is

`Delta f_FWHM ~= 0.442/(GVM*L)`.

Using the 150 fs/mm group-velocity mismatch extracted in the high-efficiency
PPLN experiment gives:

| PPLN length | Estimated FWHM at 1550 nm | Low-depletion ideal conversion at 100 mW, using 2600% W^-1 cm^-2 |
|---:|---:|---:|
| 2 mm | 1.47 THz, about 11.8 nm | 10.4% |
| 4 mm | 0.74 THz, about 5.9 nm | 41.6% |
| 6 mm | 0.49 THz, about 3.9 nm | 93.6%, outside the safe low-depletion regime |

These are planning calculations, not guaranteed performance. Loss, duty-cycle
error, photorefraction, coupling, pump depletion, nonuniformity, and actual GVM
must be included in the final model. The standalone baseline is 4 mm because
its purpose is maximum CW SHG efficiency; the 2-mm and 6-mm lanes identify the
efficiency-bandwidth and uniformity tradeoff experimentally.

## 6. S3: monolithic EO-comb-to-PPLN-SHG cascade

No source identified in this survey demonstrates this exact monolithic
combination of a four-pass travelling-wave TFLN EO comb and a 1550-to-775-nm
PPLN section on the same die. The two component technologies are independently
demonstrated, and system-level EO-comb-to-PPLN harmonic transfer has been
demonstrated, but S3 itself is an engineering synthesis. Its main unresolved
risks are the common-stack compromise, the EO-to-PPLN mode transition, and
fabrication compatibility between domain poling and permanent RF metal.

### 6.1 Optical and RF path

`1550-nm CW input -> four-pass TW phase modulator -> adiabatic width/mode`

`transition -> 2-mm PPLN -> co-propagating 1550/775-nm output`

Use end-fire collection for the first mask. An on-chip 1550/775-nm
demultiplexer is optional and should only be inserted after a dual-wavelength
EME/FDTD design proves low loss and sufficient isolation. A simple external
dichroic filter is less risky for the first measurement.

### 6.2 Why the integrated PPLN is 2 mm

A 47-line, 25-GHz comb occupies approximately

`(47 - 1)*25 GHz = 1.15 THz`,

which corresponds to about 9.2 nm near 1550 nm. A 4-mm uniform PPLN section is
estimated to provide only about 5.9 nm FWHM for a 150-fs/mm GVM, so it would
preferentially convert the center of the comb. A 2-mm section gives about
11.8 nm FWHM and covers the planned comb with roughly 28% bandwidth margin.

Therefore:

- **S3 baseline:** 2-mm uniform PPLN, optimized for comb coverage;
- **S2 efficiency reference:** 4-mm uniform PPLN;
- **S2/S3 advanced coupon:** 4-mm weakly chirped PPLN, whose `Lambda(z)` is
  generated from the simulated propagation constants rather than chosen by
  an arbitrary linear chirp.

### 6.3 Preliminary layout values

| Item | Baseline |
|---|---:|
| EO block | 10-mm electrode, four optical passes |
| EO-to-PPLN transition | 0.3-0.5 mm adiabatic starting length |
| Integrated PPLN | 2.0 mm uniform first-order QPM |
| Optical input/output fanout | about 1.5-2.0 mm total allowance |
| Reserved box | 18.0 mm x 1.6 mm = 28.8 mm^2 |
| Fundamental comb target | >=40 observable teeth near 1550 nm at 25 GHz |
| SH comb target | same 25-GHz spacing, QPM-limited span near 775 nm |
| On-chip pump for first test | 20-100 mW sweep, then increase under thermal/photorefractive monitoring |

The integrated 2-mm PPLN has an ideal low-depletion conversion estimate of
about 10.4% at 100 mW using the best published normalized efficiency. A first
device achieving lower conversion can still validate the two-comb mechanism;
the measured S2 efficiency should be used to predict S3 before comparison.

## 7. Preliminary die floorplan

Use a 20 mm x 6 mm die with all end-fire optical ports on the left and right
facets. Put RF probe pads along the top edge to keep the probes clear of lensed
fibers.

- top row: S1, 13.0 mm x 1.4 mm;
- middle row: S3, 18.0 mm x 1.6 mm;
- bottom-left row: S2, 7.0 mm x 1.6 mm;
- bottom-right free region: CPW thru/open/short, unpoled reference, mode-
  multiplexer/crossing controls, and waveguide cutbacks.

Keep at least 0.5 mm longitudinal facet margin before the first taper and use
the foundry dicing-street rule in the transverse direction. The 20 mm x 6 mm
die is a planning envelope; a foundry reticle or maximum die size can force a
different arrangement.

## 8. Required simulations and tools

| Problem | Recommended tool | Model/output |
|---|---|---|
| 1550/775 modes and dispersion | Lumerical MODE or COMSOL Wave Optics | `n_eff`, `n_g`, modes, GVM, overlap, sensitivity |
| QPM period and SHG | MODE/COMSOL + Python coupled-wave model or pyChi | `Lambda`, efficiency, bandwidth, chirp, pump depletion |
| CPW cross-section | HFSS 2D/Q3D or COMSOL | `Z0`, microwave index, loss, RF field overlap |
| Full electrode/pads | 3D HFSS | `S11`, `S21`, discontinuities, RF heating proxy |
| Mode multiplexer/crossing/taper | EME/FDTD | insertion loss, crosstalk, fabrication tolerance |
| Full EO-SHG comb | Python/MATLAB | multi-line EO + SHG/SFG spectrum and power conservation |
| GDS generation/check | KLayout or gdsfactory | hierarchy, connectivity, DRC-ready layout |

HFSS is appropriate for the electrical structure, but it cannot replace the
optical eigenmode, QPM, and nonlinear comb simulations.

## 9. Verification and measurement matrix

### S1

- RF `S11/S21` for the complete line and calibration coupons;
- extracted RF loss and microwave index;
- `Vpi(f)` from 5-40 GHz for single-, double-, and four-pass devices;
- EO comb spectra versus RF frequency and power;
- insertion loss separated into straight waveguide, mode multiplexers,
  crossing, and loops.

### S2

- SH power versus pump wavelength and power;
- QPM center/FWHM versus period, width, length, and temperature;
- normalized and absolute on-chip efficiency after coupling calibration;
- domain duty cycle and longitudinal uniformity;
- photorefractive drift or thermal bistability at 775-nm output.

### S3

- simultaneous 1550-nm and 775-nm spectra;
- exact RF spacing at both wavelengths;
- SH line-to-line response relative to the S2 QPM transfer function;
- power scaling, conversion, long-term drift, and thermal sensitivity;
- coherence/phase-noise test only after the optical mechanism is verified.

### Minimum laboratory equipment

- 1520-1580-nm tunable CW laser, polarization controller, and calibrated
  1550-nm power meter;
- EDFA and variable attenuator for PPLN power sweeps, with an interlock or
  conservative software limit during initial tests;
- microwave synthesizer, amplifier, bias tee if required, 40-GHz-class VNA,
  two compatible GSG probes, and a broadband 50-ohm termination;
- telecom OSA and a visible OSA or scanning spectrometer whose resolution is
  better than 25 GHz (about 0.05 nm at 775 nm) if individual SH teeth are to be
  resolved;
- 775-nm power meter/detector, dichroic filters, visible camera, and calibrated
  775-nm coupling reference;
- temperature-controlled stage and a top-view microscope suitable for
  alignment and monitoring visible scatter.

Absolute on-chip SHG efficiency is not credible until input/output coupling at
both wavelengths has been calibrated separately.

## 10. Mask-release gates

1. **Foundry gate:** stack, poling, metal, probe, and dicing rules confirmed.
2. **Optical gate:** a dual-wavelength mode pair with fabricable QPM period and
   acceptable overlap is selected.
3. **RF gate:** complete 3D electrode meets impedance, reflection, loss,
   velocity-match, and `Vpi` targets.
4. **Bandwidth gate:** S3 PPLN covers the simulated S1 comb span with at least
   20% FWHM margin.
5. **Tolerance gate:** period/width/thickness/temperature corners remain inside
   the DOE and available heater tuning.
6. **Testability gate:** every critical building block has an independent
   control and an accessible optical or RF port.
7. **DRC gate:** official foundry DRC passes; project-level geometry checks are
   not a substitute.

## 11. Main risks and mitigations

| Risk | Consequence | Mitigation in this plan |
|---|---|---|
| PPLN bandwidth narrower than the EO comb | Only central SH teeth survive | 2-mm S3 baseline; simulated chirped coupon; explicit bandwidth gate |
| Domain nonuniformity | Distorted QPM curve and reduced efficiency | 2/4/6-mm lanes, period sweep, SH microscopy |
| Four-pass phase error or mode crosstalk | Low effective modulation index | single/double/four-pass controls and isolated mux/crossing coupons |
| RF attenuation or mismatch | `Vpi` rises at 25 GHz | thick metal if allowed, HFSS pads/line/termination co-design, TRL coupons |
| Metal absorption | High insertion loss | gap sweep and optical-metal loss simulation |
| Photorefraction at 775 nm | drift or instability | MgO:LNOI preference, start at low power, temperature monitor/control |
| Facet coupling at two wavelengths | uncertain absolute efficiency | calibrated 1550/775 coupling references; external dichroic first |
| Literature dimensions copied to a different stack | QPM/RF miss | use all values only as starting points and enforce simulation gates |

## 12. Evidence table

The following are source-backed facts used to set the starting scale; all
reserved boxes and performance gates above are engineering proposals.

1. Ke Zhang et al., *A power-efficient integrated lithium niobate
   electro-optic comb generator*, Communications Physics 6, 17 (2023),
   https://www.nature.com/articles/s42005-023-01137-9 . Four-pass optical
   recycling, 1-cm travelling-wave electrode, 1.90-V effective RF `Vpi` at
   24.95 GHz, and 47 lines at 25 GHz and 28 dBm.
2. Zixuan Huang et al., *High-Efficiency and Ultra-Compact Thin-Film Lithium
   Niobate Phase Modulator Based on a Dual-Noncentered Waveguides Recycling
   Structure*, JLT 44, 4120-4128 (2026),
   https://opg.optica.org/jlt/abstract.cfm?uri=jlt-44-10-4120 . Demonstrates a
   compact travelling-wave recycling alternative and tabulates the
   approximately 0.8 mm x 12 mm footprint scale of the 1.9-Vcm multi-pass
   predecessor.
3. Cheng Wang et al., *Ultrahigh-efficiency second-harmonic generation in
   nanophotonic PPLN waveguides*, Optica 5, 1438-1441 (2018),
   https://arxiv.org/abs/1810.09235 . Reports 4-mm PPLN, 2600% W^-1 cm^-2,
   approximately 4.1-um period, 150-fs/mm GVM, and 53% absolute conversion.
4. Hubert S. Stokowski et al., *Integrated frequency-modulated optical parametric
   oscillator*, Nature 627, 95-100 (2024),
   https://www.nature.com/articles/s41586-024-07071-2 . Extended data gives the
   1.20-um width, 300-nm ridge, 200-nm slab, 700-nm cladding, and 3.70-um PPLN
   starting geometry.
5. Markus Ludwig et al., *Ultraviolet astronomical spectrograph calibration with
   laser frequency combs from nanophotonic lithium niobate waveguides*, Nature
   Communications 15, 7614 (2024),
   https://www.nature.com/articles/s41467-024-51560-x . Demonstrates that an
   18-GHz EO comb can be transferred through tailored PPLN harmonic generation
   and explicitly uses chirped poling for broadband QPM.
6. Hao Wen et al., *Broadband and flat-top integrated electro-optic frequency
   combs on a thin-film lithium niobate platform*, Optics Letters 50,
   3692-3695 (2025),
   https://opg.optica.org/ol/abstract.cfm?URI=ol-50-11-3692 . One intensity
   modulator plus two phase modulators produced more than 70 lines within a
   10-dB bandwidth, supporting the later flat-top upgrade path.

## 13. Current conclusion

The three-structure tapeout is physically plausible and experimentally
diagnosable. The most defensible first integrated device is not a doubly
resonant ring; it is a four-pass, low-`Vpi` travelling-wave phase modulator
followed by a short broadband PPLN section. The approximate 20 mm x 6 mm die
envelope is sufficient for all three functional blocks plus essential controls.

The dimensions are now detailed enough to start optical and RF simulation, but
not yet to fabricate. The next irreversible step is to obtain the actual
foundry stack/design rules and run the Stage 1 dual-wavelength mode sweep.
