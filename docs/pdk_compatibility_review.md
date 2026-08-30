# PDK Compatibility Review

Review date: 2026-08-30

## Outcome

The supplied process is a useful baseline for the 1550-nm travelling-wave
electro-optic section, but it does **not** presently qualify periodic poling,
second-harmonic generation, 775-nm routing, or dual-band couplers.  The compact
poling-free SHG branch is therefore a custom-device proposal inside the layer
rules, not a foundry-supported PCell flow.

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
| LN minimum bend radius | 80 um | sets the first CQPM ring radius in the DOE |
| M1 minimum line/space | 2/3 um | applies to travelling-wave electrodes and heaters/fanout checks |
| RF minimum electrode gap | 3 um | lower geometric bound, not the optical-loss optimum |
| GDS grid | 0.001 um | use one hierarchical top cell; do not flatten foundry BlackBoxes |

The manual does not disclose the two LN etch depths, sidewall angle, final top
cladding thickness, or M1 material/thickness.  These are mandatory simulation
inputs and must be obtained before a mask can be frozen.

## Supplied building blocks relevant to this project

- The extraordinary-polarized, Y-propagating 1550-nm phase-modulator
  BlackBox has a 9.1-mm optical interaction length and an approximately
  9.58 mm x 0.37 mm protected bounding box.
- The supplied 1550-nm MZI is specified in the manual at greater than 67 GHz
  bandwidth and less than 3 V.cm `VpiL`.  This number must **not** be silently
  assigned to the phase-modulator BlackBox; its `Vpi(f)` is not stated.
- 1550-nm LN/SiN edge and grating couplers, MMIs, crossings, heaters, and an
  LN-SiN interlayer coupler are provided.
- No PPLN, poling, SHG, 775-nm, visible coupler, dual-band WDM, or nonlinear
  resonator component appears in the supplied layer list or BlackBox library.

## Structure-by-structure compatibility

| Structure | Status | What is reusable | What remains custom |
|---|---|---|---|
| S1: 1550-nm travelling-wave EO comb | PDK-aligned concept | standard LN PM, M1, 1550-nm couplers | measured `Vpi(f)`, RF loss/index and actual comb span |
| S2: poling-free SQPM SHG racetrack DOE | custom, not qualified | x-cut LN layers, heater/M1, bend-rule envelope | 775-nm modes/loss, integer-phase geometry, double resonance, ring/bus coupling and visible extraction |
| S3: SHG then dual-band EO comb | custom, highest risk | S1 electrode concept and S2 racetrack concept | one bus or two buses at 1550/775 nm, dual-band PM, couplers and thermal control |

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

`1550-nm CW -> SQPM-SHG racetrack -> co-propagating 1550/775 CW -> travelling-wave EO phase modulation`

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
4. May users draw custom LN micro-rings/micro-disks and dual-band buses rather
   than only connecting qualified BlackBoxes?
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
DRC deck.  BlackBox placement and final signoff therefore need the supported
PDK environment rather than a generic GDS-only script.
