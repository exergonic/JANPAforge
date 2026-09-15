# Worked examples

One folder per molecule, each carried through the full pipeline
(`.gbw` -> `orca_2mkl` -> `molden2molden` -> `janpa` -> viewer/analysis).
Every folder is **self-contained**: the viewer and analysis modes run in
place using only the files in the folder:

```powershell
cd ethene
# re-create the viewer file from the substrate:
python ../../orca_to_janpa.py --to-cart ethene_CLPO_spherical.molden --sort-energy --avogadro
# re-run the interaction table:
python ../../orca_to_janpa.py --e2 ethene_CLPO_spherical.molden
```

(The originals of these files live in the ORCA calculation folders; these
copies are the pipeline artifacts as produced. Every set works the same
way — swap `CLPO` for `LHO`, `AHO`, `LPO`, `NAO` or `PNAO`.)

## Files in each folder

| file | produced by | what it is |
| --- | --- | --- |
| `<mol>.molden.input` | `orca_2mkl` | SCF MOs + basis, ORCA flavour |
| `<mol>.PURE` | `molden2molden` | canonical Molden; janpa input, and the canonical reference for `--sort-energy` / `--e2` |
| `<mol>_<SET>_spherical.molden` | `janpa` + label fixes | the **substrate** for set `<SET>` (CLPO, LHO, AHO, LPO, NAO, PNAO): spherical d/f, `Ene=` = sequential numbers, markers + `Spin=` corrected; the analysis input |
| `<mol>.JANPA` | janpa stdout | NPA charges, Wiberg bond indices, full CLPO labels, JANPA's own charge-transfer table |
| `<mol>.S.txt`, `<mol>.fock_ao.txt` | `janpa -doFock` | overlap / Fock in the AO basis (spherical, [GTO] order) |
| `<mol>.fock_nao.txt`, `<mol>.clpo2lho.txt`, `<mol>.lho2nao.txt`, `<mol>.aho2nao.txt`, `<mol>.lpo2aho.txt` | `janpa -doFock` | NAO Fock and the transformation chains (per-set route-B cross-checks inside `--e2`) |
| `<mol>_CLPO_E2.txt` | `--e2` | pair-interaction table: E2 (kcal/mol) + charge transfer q (e); any other set runs the same way |
| `<mol>_<SET>.molden` | `<SET>` flag | **the viewer file** for set `<SET>`: cartesian d/f, markers clean, real Fock energies, occupied-first order; generated with `--avogadro` in these copies (integer `Occup`) — **open these in Avogadro**. `PNAO` is energy-sorted like the rest; it is the pre-orthogonalization set (normalized but mutually non-orthogonal, `Occup` not summing to the electron count), so its run report states that and `--e2` is undefined for it |
| `water_CLPO_Alpha.molden.txt` (water only) | `--clpo` without `--avogadro` | cartesian + fractional `Occup`: the repro file for the Avogadro electron-counting bug ([OpenChemistry/avogadrolibs#3005](https://github.com/OpenChemistry/avogadrolibs/issues/3005)); the `.txt` suffix is the renamed form GitHub accepts as an attachment — bytes identical |
| `<mol>.xyz` | ORCA | geometry |

Every folder carries the complete set for all six JANPA sets (viewer +
substrate each), the dumps and transformation chains, the CLPO E2 table
and the labels. To regenerate a viewer file with the true fractional
`Occup` (instead of the Avogadro `2/0`), rerun it from the substrate:

```powershell
python ../../orca_to_janpa.py --to-cart ethene_CLPO_spherical.molden --sort-energy
```

## The molecules

| example | level of theory | what it demonstrates |
| --- | --- | --- |
| `ethene/` | wB97X-D3/def2-TZVP | the reference example: the pi CLPO **is** the canonical HOMO; sigma(C-H) -> sigma*(C-H) hyperconjugation, q = 0.0080 e (E2 ~5.6 kcal/mol, indicative); substrate + viewer pair |
| `water/` | **HF**/def2-SVP | the Hartree-Fock case -- HF, the regime `VALIDATION.md` uses for its NBO 3.1 comparison (cc-pVDZ in that run, def2-SVP here); weak delocalization, table tops out at ~2 kcal/mol (q = 0.0012 e) |
| `formaldehyde/` | wB97X-D3/def2-TZVP | lone-pair donor: O n -> sigma*(C-H) x2, **q = 0.059 e** each (E2 ~29 kcal/mol, indicative -- and the one channel where CLPO and NBO disagree most; see `../VALIDATION.md`) |
| `isobutene/` | wB97X-D3/def2-TZVP | sigma(C-H) -> pi*(C=C) hyperconjugation, **q = 0.015 e** per methyl C-H (E2 ~5.4, indicative); vinylic C-H -> sigma*(C-C) E2 ~8.7 |
| `ethylium/` | wB97X-D3/def2-TZVP | bridged 3c-2e ethyl cation: the E2 table marks the strong interaction with `*` instead of pretending it is a hyperconjugation energy |
| `tbutyl/` | wB97X-D3/def2-TZVP | carbocation hyperconjugation: 3x C-H -> empty p on C+, **q = 0.084 e** each (E2 ~34 kcal/mol, indicative) |

**Reading the E(2) numbers.** For the wB97X-D3 examples, E2 is
**indicative**: the Kohn-Sham "Fock" matrix is not the HF Fock (each
`_E2.txt` header says so), and the robust per-pair quantity there is the
charge q = D_ij^2/D_ii -- the number JANPA's own CT analysis prints. The
one channel where the two localization schemes genuinely part ways is
formaldehyde's O n -> sigma*(C-H): at HF, CLPO gives ~26.5 kcal/mol where
NBO 3.1 gives 0.85 (different acceptors; the full story is in
`../VALIDATION.md`). Read the DFT ~29 as CLPO-side and scheme-dependent,
not as an NBO-grade number. The HF row (water) is the one with a direct
cross-program check.

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
