# Worked examples

One folder per molecule, each carried through the full pipeline
(`.gbw` -> `orca_2mkl` -> `molden2molden` -> `janpa` -> viewer/analysis).
Every folder is **self-contained**: the viewer and analysis modes run in
place using only the files in the folder:

```powershell
cd ethene
# re-create the viewer file from the substrate:
python ../orca_to_janpa.py --to-cart ethene_CLPO_spherical.molden --sort-energy --avogadro
# re-run the interaction table:
python ../orca_to_janpa.py --e2 ethene_CLPO_spherical.molden
```

(The originals of these files live in the ORCA calculation folders; these
copies are the pipeline artifacts as produced.)

## Files in each folder

| file | produced by | what it is |
| --- | --- | --- |
| `<mol>.molden.input` | `orca_2mkl` | SCF MOs + basis, ORCA flavour |
| `<mol>.PURE` | `molden2molden` | canonical Molden; janpa input, and the canonical reference for `--sort-energy` / `--e2` |
| `<mol>_CLPO_spherical.molden` | `janpa` + label fixes | the **substrate**: spherical d/f, `Ene=` = sequential numbers, markers + `Spin=` corrected; the analysis input |
| `<mol>.JANPA` | janpa stdout | NPA charges, Wiberg bond indices, full CLPO labels, JANPA's own charge-transfer table |
| `<mol>.S.txt`, `<mol>.fock_ao.txt` | `janpa -doFock` | overlap / Fock in the AO basis (spherical, [GTO] order) |
| `<mol>.fock_nao.txt`, `<mol>.clpo2lho.txt`, `<mol>.lho2nao.txt` | `janpa -doFock` | NAO Fock and the LHO transformation chain (route-B cross-check inside `--e2`) |
| `<mol>_CLPO_E2.txt` | `--e2` | pair-interaction table: E2 (kcal/mol) + charge transfer q (e) |
| `<mol>_CLPO.molden` | `--clpo` | **the viewer file**: cartesian d/f, markers clean, real Fock energies, occupied-first order; generated with `--avogadro` in these copies (integer `Occup`) — **open this one in Avogadro** |
| `water_CLPO_Alpha.molden` (water only) | `--clpo` without `--avogadro` | cartesian + fractional `Occup`: the repro file for the Avogadro electron-counting bug |
| `<mol>.xyz` | ORCA | geometry |

Every folder carries the complete set: the substrate, the viewer file, the
dumps, the E2 table and the labels. To regenerate the viewer file with the
true fractional `Occup` (instead of the Avogadro `2/0`), rerun it from the
substrate:

```powershell
python ../orca_to_janpa.py --to-cart ethene_CLPO_spherical.molden --sort-energy
```

## The molecules

| example | level of theory | what it demonstrates |
| --- | --- | --- |
| `ethene/` | wB97X-D3/def2-TZVP | the reference example: the pi CLPO **is** the canonical HOMO; sigma(C-H) -> sigma*(C-H) hyperconjugation ~5.6 kcal/mol; substrate + viewer pair |
| `water/` | **HF**/def2-SVP | the Hartree-Fock case (where E(2) is defensible); weak delocalization, table tops out at ~2 kcal/mol |
| `formaldehyde/` | wB97X-D3/def2-TZVP | lone-pair donor: O n -> sigma*(C-H) ~29 kcal/mol x2 (classic negative hyperconjugation) |
| `isobutene/` | wB97X-D3/def2-TZVP | sigma(C-H) -> pi*(C=C) hyperconjugation ~5.4 kcal/mol per methyl C-H; vinylic C-H -> sigma*(C-C) ~8.7 |
| `ethylium/` | wB97X-D3/def2-TZVP | bridged 3c-2e ethyl cation: the E2 table marks the strong interaction with `*` instead of pretending it is a hyperconjugation energy |
| `tbutyl/` | wB97X-D3/def2-TZVP | carbocation hyperconjugation: 3x C-H -> empty p on C+ ~34 kcal/mol |

Two sanity signals worth knowing: each `_E2.txt` run re-checks JANPA's
printed CT pairs live (formaldehyde 3/3, isobutene 6/6, tbutyl 9/9 in the
copies here), and the `sum q` column reproduces the per-molecule CT totals
quoted in `<mol>.JANPA` (e.g. ethene 0.05933 e).

Note on `ethylium/`: this is the **bridged** (nonclassical) structure — H
bridging both carbons at ~1.30 A, C-C = 1.37 A. The Freq run in the
original calculation shows no imaginary frequencies, which is the known
DFT-functional behaviour (bridged = minimum, e.g. JPCA 2002,
doi 10.1021/jp0215264); wavefunction methods make the classical structure
the minimum and this one a transition state.
