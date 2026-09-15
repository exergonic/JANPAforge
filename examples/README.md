# Worked examples

One folder per molecule, each carried through the full pipeline
(`.gbw` -> `orca_2mkl` -> `molden2molden` -> `janpa` -> viewer/analysis).
Every folder is **self-contained**: the viewer and analysis modes run in
place using only the files in the folder:

```powershell
cd ethene
python ../orca_to_janpa.py --to-cart ethene_CLPO.molden --avogadro --sort-energy
python ../orca_to_janpa.py --e2 ethene_CLPO.molden
```

(The originals of these files live in the ORCA calculation folders; these
copies are the pipeline artifacts as produced.)

## Files in each folder

| file | produced by | what it is |
| --- | --- | --- |
| `<mol>.molden.input` | `orca_2mkl` | SCF MOs + basis, ORCA flavour |
| `<mol>.PURE` | `molden2molden` | canonical Molden; janpa input, and the canonical reference for `--sort-energy` / `--e2` |
| `<mol>_CLPO.molden` | `janpa` | the CLPO export (spherical, `Ene=` = sequential numbers) |
| `<mol>.JANPA` | janpa stdout | NPA charges, Wiberg bond indices, full CLPO labels, JANPA's own charge-transfer table |
| `<mol>.S.txt`, `<mol>.fock_ao.txt` | `janpa -doFock` | overlap / Fock in the AO basis (spherical, [GTO] order) |
| `<mol>.fock_nao.txt`, `<mol>.clpo2lho.txt`, `<mol>.lho2nao.txt` | `janpa -doFock` | NAO Fock and the LHO transformation chain (route-B cross-check inside `--e2`) |
| `<mol>_CLPO_E2.txt` | `--e2` | pair-interaction table: E2 (kcal/mol) + charge transfer q (e) |
| `<mol>_CLPO_cart.molden` | `--to-cart --sort-energy` | cartesian d/f, real energies, fractional Occup (viewers needing cartesian only) |
| `<mol>_CLPO_Avogadro.molden` | `--to-cart --avogadro --sort-energy` | **open this one in Avogadro**: cartesian + integer `Occup` + `Spin= Alpha` + energy order |
| `<mol>_CLPO_7F.molden` | `--fix-markers --sort-energy` | spherical with corrected `[5D]`/`[7F]` markers (MOrbVis) |
| `<mol>.xyz` | ORCA | geometry |

Not every folder carries every view file yet — the two earliest examples
(ethene, formaldehyde) have them all; the later analysis-only ones
(water, isobutene, ethylium, tbutyl) have `_CLPO_E2.txt` but no
`_cart`/`_Avogadro`/`_7F` files yet.

## The molecules

| example | level of theory | what it demonstrates |
| --- | --- | --- |
| `ethene/` | wB97X-D3/def2-TZVP | the reference example: the pi CLPO **is** the canonical HOMO; sigma(C-H) -> sigma*(C-H) hyperconjugation ~5.6 kcal/mol; all view files |
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
