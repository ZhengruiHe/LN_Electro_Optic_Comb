# PDK Compatibility Review

Review date: 2026-08-30

## Outcome

The supplied process is used only as the physical envelope: material stack,
layer definitions, minimum geometry and die boundary.  The three proposed
devices are all custom.  In particular, no supplied phase-modulator BlackBox
is used.  The process does **not** presently qualify periodic poling,
second-harmonic generation, 775-nm routing, dual-band couplers or the proposed
four-pass modulator, so this remains a custom-device study rather than a
foundry-supported PCell flow.

The PDK archive and extracted files are NDA/confidential source material.  They
remain local and are excluded from Git.  This repository records only the
derived design constraints needed for project planning.

## Confirmed process envelope

| Item | Supplied value | Design consequence |
|---|---:|---|
| Platform | x-cut TFLN on LPCVD SiN | suitable crystal cut for in-plane EO use and cyclic phase matching |
| Nominal wavelength | 1550 nm | every 775-nm component is custom and unqualified |
| LN film | 400 nm typical | replaces the former 500-nm literature assumption |
| SiN film | 300 nm typical | useful for passive routing, but must be removed below active LN guides where required |
| LN-SiN oxide spacer | about 400 nm typical | must be included in optical and RF simulations |
| Effective design area | 21.8 mm x 3.8 mm | replaces the former 20 mm x 6 mm floorplan |
| LN minimum line/space | 0.30/0.30 um | applies to custom rings, buses, and tapers |
| LN minimum bend radius | 80 um | hard lower bound for custom rings and loopbacks |
| M1 minimum line/space | 2/3 um | applies to travelling-wave electrodes and heaters/fanout checks |
| RF minimum electrode gap | 3 um | lower geometric bound, not the optical-loss optimum |
| GDS grid | 0.001 um | use hierarchical custom cells and preserve the grid |

The manual does not disclose the two LN etch depths, sidewall angle, final top
cladding thickness, or M1 material/thickness.  These are mandatory simulation
inputs and must be obtained before a mask can be frozen.

## How the PDK is and is not used

- Reuse the declared LN/SiN/oxide/metal layers, die boundary, minimum
  line/space, minimum bend radius and final foundry DRC.
- Draw the S1 waveguide, optical-recycling loops, mode multiplexers, GSG
  electrode and termination as custom geometry.
- Draw the S2/S3 racetracks, buses, heaters, 775-nm routes and dual-rail
  modulation section as custom geometry.
- Supplied couplers or passive cells may be used only as measurement-interface
  references after checking their wavelength and polarization qualification;
  they are not part of the device physics assumed here.
- Do not transfer a bandwidth or `VpiL` number from a supplied MZI or PM into
  these custom structures.

## Structure-by-structure compatibility

| Structure | Status | What is reusable | What remains custom |
|---|---|---|---|
| S1: four-pass 1550-nm travelling-wave EO comb | custom | x-cut stack, M1 rules, bend-rule envelope | CPW geometry, optical/RF velocity match, EO overlap, mode multiplexers, delay loops and `Vpi(f)` |
| S2: poling-free SQPM SHG racetrack DOE | custom, not qualified | x-cut stack, heater/M1 rules, bend-rule envelope | 775-nm modes/loss, integer-phase geometry, double resonance, ring/bus coupling and visible extraction |
| S3: SHG then dual-rail EO comb | custom, highest risk | the process envelope only | two wavelength-specific optical rails, dual-band extraction, common RF electrode and thermal control |

## Why poling-free does not mean phase-matching-free

In an x-cut LN ring, the propagation direction rotates relative to the optical
axis.  The TE effective index and effective nonlinear coefficient therefore
vary periodically around the ring.  This cyclic quasi-phase matching (CQPM)
can make successive generated SH fields add constructively without domain
inversion.

It still requires:

1. energy matching between a fundamental cavity mode and an SH cavity mode;
2. adequate azimuthal/momentum selection and nonlinear mode overlap;
3. coupling at both wavelengths; and
4. resonance control against fabrication and thermal detuning.

The selected integrated order is consequently:

`1550-nm CW -> SQPM-SHG racetrack -> separate 1550/775 rails -> common-electrode travelling-wave EO phase modulation`

The ring converts a single optical carrier before comb generation.  Its FSR
does not need to equal the 25-GHz RF drive, and the EO comb teeth do not have to
land on a complete two-colour cavity grid.

## Foundry questions that block release

Request written answers for all of the following:

1. Which archive/manual version is current?  The outer filename and the manual
   revision label do not agree.
2. What are the exact LN1/LN2 etch depths, sidewall-angle distribution, final
   top-cladding thickness, and wafer thickness tolerances?
3. What are the M1 metal stack, thickness, conductivity/sheet resistance, and
   maximum recommended RF frequency?
4. May users draw custom LN waveguides, mode multiplexers, travelling-wave
   electrodes, micro-rings/micro-disks and dual-band buses?
5. Is propagation and optical-power handling near 775 nm permitted, and are
   there measured LN/SiN losses or photorefractive limits?
6. Can the foundry support a 1550/775-nm edge coupler or permit a custom one at
   the die facet?
7. Does the required SiN-removal trench below LN guides apply identically to
   the ring, bus, coupling gap, and dual-band phase-modulator region?
8. What is the current MPW schedule?  The dates printed in the supplied manual
   are already past on the review date.

## Release decision

The compact layout may proceed through simulation and test-cell planning.  It
must remain marked **not for fabrication** until the foundry answers the eight
questions above, the dual-band optical model closes, and the official DRC
passes.

The supplied IPKISS and Latitudeda wrappers also require their corresponding
commercial runtimes.  The current generic Python environment has neither
`ipkiss3` nor `fnpcell`, and the package does not contain a complete automatic
DRC deck.  Final custom-cell assembly and signoff therefore need the supported
PDK environment rather than a generic GDS-only script.
