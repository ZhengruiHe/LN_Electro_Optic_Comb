# HFSS Travelling-Wave Electrode Simulation Plan

## 1. Simulation objective

The HFSS task is to design the custom GSG travelling-wave electrode shared by
the S1 and S3 devices. It must establish, rather than assume:

- characteristic impedance `Z0(f)`;
- microwave propagation constant `gamma(f) = alpha(f) + j beta(f)`;
- microwave effective index `n_RF(f)` and conductor/dielectric loss;
- two-port `S11/S21` including GSG launches;
- normalized RF electric field around the 1550- and 775-nm waveguides;
- the RF part of the frequency-dependent EO response.

HFSS does not solve the optical mode or the Pockels overlap by itself. Optical
fields and `n_g` come from MODE/COMSOL and are combined with the HFSS field in
Python after both models have converged.

## 2. Why the final model is HFSS 3D Driven Terminal

Use **HFSS 3D, solution type Driven Terminal**, with one signal terminal and
the two GSG ground conductors as references. This directly represents the
multi-conductor CPW and provides terminal impedance and S parameters.

Q3D may be used only for a fast isotropic pre-screen. Ansys documents that
anisotropic materials are not available in Q3D, whereas x-cut LN has an
orientation-dependent RF permittivity tensor. The final `Z0`, `n_RF`, loss and
field distribution must therefore come from the anisotropic HFSS model.

## 3. Model hierarchy

Do not begin with the complete 10-mm device.

| Model | Geometry | Purpose | Required output |
|---|---|---|---|
| M0a | uniform CPW, 0.5 mm | short-line calibration | unrenormalized S2P, terminal `Zpi`, 25-GHz field |
| M0b | same cross-section, 1.0 mm | two-length extraction | unrenormalized S2P |
| M1 | uniform CPW, 10 mm | verify accumulated RF loss and phase | S2P from 1-50 GHz |
| M2 | one GSG probe pad and taper | optimize launch independently | launch return/insertion loss and current crowding |
| M3 | pad + 10-mm line + pad | final electrical model | 50-ohm-renormalized and native S2P |
| M4 | selected M0 cross-section with both optical rail positions | EO-field export | complex `E_RF(x,z,f)/V_line` at 5, 25 and 50 GHz |

M0a/M0b are the main parameter-sweep models. Only the best two or three
cross-sections proceed to M1. Only one cross-section proceeds to M3.

## 4. Coordinate system and stack

Use the following global coordinates consistently in HFSS and the optical mode
solver:

- global `X`: lateral direction across the GSG gaps; crystal `z` axis;
- global `Y`: RF/optical propagation direction; crystal `y` axis;
- global `Z`: wafer normal; crystal `x` axis for x-cut LN.

The LN RF tensor must be entered after this coordinate rotation. The current
configuration uses the high-frequency literature values `[28, 43, 43]` in
global `X/Y/Z` only as placeholders. Replace them with foundry-approved or
measured values before accepting any result.

The explicit stack is, from bottom to top: handle substrate, bottom oxide,
SiN, LN-SiN interlayer oxide, LN, top cladding and M1. The active-region case
with SiN removed below the LN guide must be compared against the continuous
SiN case if the foundry requires that trench.

## 5. Starting CPW geometry and parameter sweep

The conventional CPW is selected first. Periodic capacitively loaded or
T-rail electrodes are a backup only if the conventional line cannot meet both
impedance and velocity targets.

Starting point:

- signal width `Ws = 43 um`;
- signal-ground gap `G = 5.5 um`;
- each ground width `Wg = 100 um`;
- active length `L = 10 mm`;
- optical rail placed approximately at the centre of each signal-ground gap.

The 43/5.5-um pair is taken only as a sweep seed from Zhang et al.; it was
fabricated on a different 500-nm-LN/2-um-oxide stack and is not a mask value
for this PDK.

Coarse screening grid for M0a/M0b:

| Variable | Values |
|---|---|
| `Ws` | 35, 43, 50, 60 um |
| `G` | 4, 5.5, 7, 9 um |
| `Wg` | 80, 100, 120 um |

This is 48 geometric cases. Use the 0.5-mm model first and retain only cases
with plausible terminal impedance and field concentration. Build the paired
1.0-mm model only for approximately the best 8-12 cases. Then run a local
fine sweep around the best point. Do not sweep M1 thickness or oxide values as
free optimization variables after the foundry has fixed them; use them only as
process corners.

## 6. Ports and boundaries

### Wave ports

- Assign one wave port at each Y end of the uniform CPW.
- Use Driven Terminal with the signal as the excited conductor and both grounds
  as reference conductors.
- Use one CPW mode. Export the native terminal impedance (`Zpi`) and also a
  separately 50-ohm-renormalized S2P for comparison with VNA measurements.
- The port rectangle must include substrate and air, and its left/right edges
  must touch the two outer ground conductors. Start with the Ansys CPW sizing
  guidance, then increase port width and height until `Z0` and `n_RF` change by
  less than 0.5%.
- Set de-embedding to zero for the two-length M0 extraction. De-embed only a
  known uniform section in the launch study.

### Open boundaries

Follow the official Ansys CPW Driven Terminal example:

- radiation boundary on the top/bottom thickness faces of the air region;
- Perfect H/open treatment on the lateral and end faces behind the wave ports;
- port sheet faces touch the air region but the port edges do not touch the
  outer air-box edges.

Run an air-box convergence test using at least 200, 300 and 450 um lateral/top
padding. The final result is accepted only when `Z0`, loss and `n_RF` change by
less than 0.5% between the last two cases.

### Backside condition

Simulate both physical measurement cases if the chip mounting is not fixed:

1. floating/open backside;
2. metal chuck or package ground at the measured chip-to-metal spacing.

Silicon substrate leakage and microwave index can change substantially between
these cases. The final mask decision must use the actual probe/chuck/package
condition.

## 7. Materials and conductor model

- Enter LN as an anisotropic RF dielectric using the rotated tensor.
- Enter Si, SiO2, SiN and top cladding with frequency-appropriate permittivity
  and loss tangent.
- Use the actual M1 conductivity, thickness and roughness. Do not silently use
  ideal PEC for the loss model.
- Model M1 as a volumetric conductor when its thickness is comparable with a
  few skin depths. A finite-conductivity surface approximation is acceptable
  only after checking the HFSS validity condition that conductor thickness is
  sufficiently larger than skin depth.
- Run material corners for metal conductivity/roughness and dielectric loss,
  but do not mix geometry and material uncertainty in the first sweep.

## 8. Mesh and convergence

Use multi-frequency adaptive refinement at 5, 25 and 50 GHz:

- maximum Delta-S `0.01`;
- maximum 12 adaptive passes;
- minimum two converged passes;
- interpolating sweep 1-50 GHz for S parameters;
- separate discrete 25-GHz sweep with fields saved.

Apply a skin-depth mesh to M1 with at least three layers and a surface triangle
length no larger than `min(G/2, 5 um)`. Add local refinement in both CPW gaps
and the optical-field export boxes. Do not force a micrometre mesh through the
entire 500-um handle substrate or the complete 10-mm length.

Convergence is not established merely because HFSS stops. Record:

- final Delta-S and whether the pass limit was reached;
- tetrahedron count and peak memory;
- change in `Z0`, `n_RF`, loss and 25-GHz gap field over the last two passes;
- port-size and air-box convergence.

## 9. Extracting transmission-line quantities

Export **unrenormalized** M0a/M0b Touchstone data on identical frequency grids.
For each frequency, convert both networks to ABCD matrices and form the
incremental line matrix. Its eigenvalues are `exp(+/- gamma Delta L)`, giving

`gamma(f) = alpha(f) + j beta(f)`

and

`n_RF(f) = c beta(f)/(2 pi f)`.

The repository script uses the 0.5-mm length difference to suppress common
port error. Report RF loss in dB/cm, terminal impedance, and `n_RF`. Compare the
extracted impedance with HFSS `Zpi`; disagreement greater than 2% triggers a
port/de-embedding review.

Commands after HFSS is installed and the stack is confirmed:

```powershell
.\.venv\Scripts\python.exe models\rf\build_hfss_cpw.py build --length-um 500
.\.venv\Scripts\python.exe models\rf\build_hfss_cpw.py solve --length-um 500
.\.venv\Scripts\python.exe models\rf\build_hfss_cpw.py solve --length-um 1000
.\.venv\Scripts\python.exe models\rf\postprocess_tline.py `
  --short results\hfss\CPW_L500um.s2p `
  --long results\hfss\CPW_L1000um.s2p `
  --delta-length-mm 0.5 `
  --output results\hfss\cpw_extracted.csv
```

## 10. Combining HFSS with the EO calculation

For each optical rail, export the complex RF electric field over a small X-Z
box around the optical waveguide. Normalize it by the corresponding travelling
wave terminal voltage, not simply by the default 1-W port excitation.

Combine that field with the normalized optical field and Pockels tensor using
first-order perturbation:

`Delta beta_per_V = (omega/(4 P_opt)) integral E_opt* dot Delta-epsilon(E_RF_per_V) dot E_opt dA`.

Then

`VpiL = pi / Delta beta_per_V`.

The frequency-dependent travelling-wave response for a matched line is

`H(f) proportional to [1-exp(-(alpha+j Delta_beta)L)]/(alpha+j Delta_beta)`,

where

`Delta_beta = 2 pi f (n_RF-n_g)/c`.

For S1, calculate this with `n_g_1550`, then apply the independently simulated
four-pass delay/loss factors. For S3, evaluate it twice using `n_g_1550` and
`n_g_775`; the two modulation indices are not forced to be equal.

## 11. S1 and S3 interpretation

### S1

Both CPW gaps contain 1550-nm optical modulation regions. The four passes use
the two gap field signs together with TE0/TE1 conversion and designed delay
lines. HFSS supplies the local RF amplitude and phase; optical EME/FDTD supplies
mode conversion, crossing loss and loop delay. Do not multiply efficiency by
four until the optical phase condition is satisfied.

### S3

Place the 1550-nm rail in one CPW gap and the 775-nm rail in the opposite gap.
Their RF field signs are opposite, which changes their relative comb phase but
not the common 25-GHz spacing. Use separate optical modes and separate overlap
integrals. If the calculated velocity-mismatch penalty exceeds 1 dB at 25 GHz
for either colour, compare a split-electrode fallback.

## 12. Electrical acceptance gates

The first design advances only if all are met:

- `45 ohm <= Re(Z0) <= 55 ohm` around 20-30 GHz;
- `S11 < -10 dB` over 20-30 GHz for the final pad-line-pad model;
- RF attenuation no greater than 4 dB/cm at 25 GHz as an initial screening
  target;
- S1 `abs(n_RF-n_g_1550) <= 0.1` as a preferred velocity target;
- S3 calculated velocity-mismatch penalty no greater than 1 dB at 25 GHz for
  both rails;
- no severe current crowding or local field hot spot at the selected RF drive;
- mesh, port, air-box and backside convergence documented.

These are engineering gates, not source-backed device results. Revise them
after amplifier power, probe rating and the optical `n_g` values are known.

## 13. Current environment status

PyAEDT 1.4.0 and its optional analysis dependencies are installed in the
project `.venv`. This computer currently has Ansys Optics 2024 R1 but no
detected `ANSYSEM_ROOTxxx` Electronics Desktop/HFSS installation. The scripts
can therefore be imported, checked and post-processed now, but cannot create
or solve an AEDT project on this host yet.

## 14. Evidence and limitations

### Source-backed

- [Ansys CPW Driven Terminal example](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v251/en/Subsystems/HFSS/Content/GettingStarted/CoPlanarWaveguideDrivenTerminal.htm)
  supports the terminal solution, port and open-boundary strategy.
- [Ansys wave-port sizing guidance](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v252/en/Subsystems/HFSS/Content/HFSS/WaveportSize.htm)
  requires deliberate CPW port sizing and ground contact.
- [Ansys anisotropic-material documentation](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v251/en/Subsystems/HFSS/Content/Materials/DefiningAnisotropicTensors.htm)
  states that anisotropic materials are unavailable in Q3D.
- [Ansys adaptive-meshing documentation](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v242/en/Subsystems/HFSS3DLayout/Content/HFSS/AdaptiveMeshingMethodsInHFSS.htm)
  supports multi-frequency adaptation and Delta-S convergence.
- [Ansys finite-conductivity documentation](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v252/en/Subsystems/HFSS/Content/HFSS/AssigningFiniteConductivityBoundaries.htm)
  defines the thickness-versus-skin-depth validity condition.
- [Zhang et al., Communications Physics 6, 17 (2023)](https://www.nature.com/articles/s42005-023-01137-9)
  used a 1-cm GSG line with 43-um signal width and 5.5-um signal-ground gap,
  simulated `n_RF` about 2.30 versus optical `n_g` about 2.27, and demonstrated
  four-pass 25-GHz EO-comb generation. Its stack differs from this PDK.
- [PyAEDT installation documentation](https://aedt.docs.pyansys.com/version/stable/Getting_started/Installation.html)
  supports the project-local CPython installation used here.

### Unknown before final simulation

The exact PDK etches, upper/bottom oxide, substrate, M1 stack, RF dielectric
tensor/loss, SiN trench, GSG probe and mounting boundary are still unconfirmed.
The checked-in JSON values for those fields are placeholders and the automation
script refuses a final solve unless this condition is explicitly overridden.
