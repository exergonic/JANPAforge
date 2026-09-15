<div align="center">

# JANPAforge

**Localized orbitals out of ORCA that actually open.**

One stdlib-only script: `orca_to_janpa.py`.

`.gbw` → Molden → JANPA → Molden files your viewer can trust:
real energies, clean labels, the right electron count.

[![python](https://img.shields.io/badge/python-3.10%2B-3776AB?style=flat-square&logo=python&logoColor=white)](https://www.python.org/)
[![dependencies](https://img.shields.io/badge/dependencies-0-success?style=flat-square)](#project-rules)
![output: Molden](https://img.shields.io/badge/output-Molden-e36209?style=flat-square)
[![tested with ORCA 6.1.1](https://img.shields.io/badge/tested_with-ORCA_6.1.1-0b5f8a?style=flat-square)](https://www.faccts.de/orca/)
[![tested with JANPA 2.02](https://img.shields.io/badge/tested_with-JANPA_2.02-6f42c1?style=flat-square)](http://janpa.sourceforge.net/)
[![NPA + Wiberg validated vs NBO 3.1](https://img.shields.io/badge/NPA_%2B_Wiberg_validated_vs-NBO_3.1-brightgreen?style=flat-square)](VALIDATION.md)
[![no .47 file needed](https://img.shields.io/badge/.47_file-not_needed-blueviolet?style=flat-square)](#quick-start)
[![license: MIT](https://img.shields.io/badge/license-MIT-yellow?style=flat-square)](LICENSE)

**A plain SCF `.gbw` is all it takes — no NPA/NBO keywords, no extra ORCA output.**

</div>

## What it does

[JANPA](http://janpa.sourceforge.net/) gives you the free, open-source route to
Lewis-like localized orbitals — NBO-style bonding analysis without the license.
But its raw export fights every viewer: spherical d/f coefficients, broken
`[5D]/[7F]/[9G]` marker lines, a sometimes-wrong `Spin=`, `Ene=` values that are
really just 0, 1, 2…, orbitals in JANPA's internal order, and fractional
`Occup` values that Avogadro reads as integers.

**JANPAforge** is the single-file script `orca_to_janpa.py`: it takes a plain
ORCA SCF `.gbw` the whole way — conversion, JANPA run, and cleanup — and
writes, for each JANPA orbital set, **exactly two Molden files**: one
*viewer* file that opens correctly anywhere, and the *spherical substrate*
kept for further analysis. It also computes the
**NBO-style E(2) and charge-transfer table** for any set, and it is
**validated against Gaussian 09's NBO 3.1** ([details](VALIDATION.md)).

```text
ORCA (.gbw)
   │  orca_2mkl
   ▼
.molden.input
   │  molden2molden
   ▼
.PURE
   │  janpa
   ▼
JANPA localized orbitals   (CLPO · LHO · AHO · LPO · NAO · PNAO)
   │  orca_to_janpa.py --clpo --avogadro --e2
   ▼
<base>_CLPO.molden              ← open this one in your viewer
<base>_CLPO_spherical.molden    ← the analysis input (--e2, --sort-energy, JANPA)
<base>_CLPO_E2.txt              ← pair interactions: E2 (kcal/mol) + charge transfer (e)
```

## Quick start

Everything the script itself needs is already in your Python: no packages,
no venv, nothing to install. You bring **ORCA** (for `orca_2mkl`) and the
**JANPA package** (`janpa.jar` + `molden2molden.jar` side by side).

```bash
python orca_to_janpa.py path/to/molecule --clpo --avogadro --e2
```

`--clpo` picks JANPA's NBO-analog localized-orbital set — the one for
bonding analysis; `--avogadro` applies the viewer fix for Avogadro's
integer-`Occup` parsing; `--e2` adds the pair-interaction table. (All six
sets, and what each one is for, are glossed
[below](#which-janpa-orbital-set).)

Per `<base>`, that is three steps:

- `[1/3]` `orca_2mkl <base> -molden` → `<base>.molden.input`
- `[2/3]` `molden2molden -fromorca3bf -orca3signs` → `<base>.PURE`
- `[3/3]` `janpa -i <base>.PURE` → saved as `<base>.JANPA`

…and with the set flags, the files you came for:

```text
<base>_<SET>.molden            the viewer file: cartesian d/f, markers clean,
                               Spin= corrected, real Fock energies, occupied-
                               first order; integer Occup with --avogadro
<base>_<SET>_spherical.molden  the substrate: JANPA's export with the two
                               corrected labels; the analysis input
<base>.S.txt, <base>.fock_ao.txt, <base>.fock_nao.txt
                               plus the transformation chains of the
                               requested sets — the janpa dumps
<base>_<SET>_E2.txt            with --e2: pair-interaction table
                               and full CLPO labels in <base>.JANPA
```

One set flag per set — `--clpo`, `--lho`, `--aho`, `--lpo`, `--nao`,
`--pnao` — or `--all-sets` for all six at once. Extra flags after `--` are
passed straight to `janpa.jar`:

```bash
python orca_to_janpa.py molecule -- --npacharges charges.txt
```

Point the tool at your installs once: `--orca-dir` (default
`C:/ORCA_6.1.1`, with a `PATH` fallback for `orca_2mkl`) and `--janpa-dir`
(the folder holding `janpa.jar` and `molden2molden.jar` — **not** searched
on `PATH`; default: the current directory) — or set the `ORCA_DIR` /
`JANPA_DIR` environment variables. Checking an ORCA output by hand:
`python orca_to_janpa.py --diagnose molecule.out`.

## Verified

End-to-end on `ethene_NBO.gbw` (wB97X-D3/def2-TZVP, ORCA 6.1.1):
`orca_2mkl` → `molden2molden` → `janpa` converges, 16.000000 electrons,
charge sum 0.00000, CLPO C=C / C–H bonding graph correct. No `NPA`/`NBO`
keyword needed in the ORCA input for this route — a plain SCF `.gbw`
suffices.

## Examples

`examples/` carries six molecules through the whole pipeline — every file
in place and self-contained, so the viewer and analysis modes rerun
in-folder **without ORCA**:

```bash
cd examples/ethene
python ../../orca_to_janpa.py --to-cart ethene_CLPO_spherical.molden --sort-energy --avogadro
python ../../orca_to_janpa.py --e2 ethene_CLPO_spherical.molden
```

See [examples/README.md](examples/README.md) for what each molecule shows —
hyperconjugation in isobutene, the bridged ethyl cation's `*` flag,
carbocation σ(C–H)→empty-p, and more.

## Validation

Cross-checked against Gaussian 09's NBO 3.1 on HF/cc-pVDZ for water,
formaldehyde, and isobutene: populations agree within 0.004 e, Wiberg bond
indices within 0.004, orbital occupancies within 6e-4, and isobutene's
hyperconjugation channels within 3 % — with one honest divergence
documented. The full record, plus `validation/compare.py` and the raw
evidence, is in [VALIDATION.md](VALIDATION.md).

## Viewing orbitals (both fixes are on by default)

Raw JANPA exports have two universal defects — spherical d/f coefficients
(Avogadro renders cartesian only) and wrong marker lines (JANPA omits
`[7F]` and writes a spurious `[9G]`, which mistracks f shells in readers
that honor markers). Both fixes are applied **by default** by the pipeline.
Every set flag writes exactly two Molden files for its set:

- `<base>_<SET>.molden` — **the viewer file**: cartesian d/f, markers
  clean, `Spin=` corrected from the canonical source, real Fock energies,
  occupied-first order (see "Orbital order and energies"). One file per
  set, one name each, works in every Molden reader. **Avogadro users:
  always pass `--avogadro`** — the automatic fixes are markers + cartesian;
  the occupancy workaround is opt-in. It keeps the same filename but writes
  integer `Occup= 2/0`, because Avogadro parses `Occup` as int (1.986 → 1,
  miscounting 16 electrons as 8) and fills orbitals positionally, so
  fractional occupations mislabel occupied/virtual
  ([upstream issue #3005](https://github.com/OpenChemistry/avogadrolibs/issues/3005)).
  Drop the flag when the upstream bug is fixed.
- `<base>_<SET>_spherical.molden` — the **analysis substrate**: JANPA's
  export with only the two corrected labels (markers, `Spin=`); JANPA's
  order, sequential `Ene=` values and fractional `Occup` stay exactly as
  written. It is the input for `--e2` / `--sort-energy` (same spherical
  basis as the `-doFock` dumps) and for further analysis with JANPA or
  other tools; a rerun regenerates it byte-identically from `<base>.PURE`.

`--pnao` is the exception, and only in one respect: PNAO is the
pre-orthogonalization intermediate of the NAO construction — its orbitals
are normalized but **mutually non-orthogonal** (water: max |CᵀS C − I| =
4.8e-01) and its occupancies do not sum to the electron count (14.54 e
for water's 10). The viewer file gets the full treatment including the
energy sort — E = ⟨φ|F|φ⟩ is a well-defined single-orbital expectation
value, gate-verified against the canonical expansion (5.2e-09 on water) —
but the sort report states the non-orthogonality instead of gating on it,
and the pair-interaction analysis (`--e2`) refuses PNAO files with that
explanation.

The standalone converters remain for re-processing existing files:

- `python orca_to_janpa.py --to-cart FILE.molden [--avogadro]` — the
  cartesian conversion (d 5→6, f 7→10 in Molden order, markers dropped).
- `python orca_to_janpa.py --fix-markers FILE.molden` — the spherical
  convention with corrected markers (`[5D]` if d/f shells exist, `[7F]`
  if f shells exist, `[9G]` only for real g shells), for readers that
  prefer spherical files.

Both also correct `Spin=` by default: a uniform closed-shell set takes the
label of the sibling `<base>.PURE` (`--spin Alpha|Beta` overrides).

The spherical↔cartesian maps are derived numerically from Gaussian moments
for normalized functions on both sides and cross-checked against
`molden2molden -cart2pure` unit responses; every run is gated and the
script refuses to write on mismatch.

### Orbital order and energies

A JANPA export carries no orbital energies -- every `Ene=` is a sequential
number (0, 1, 2, ...) and the MOs sit in JANPA's internal hybrid-pairing
order, so the viewer's orbital list is unreadable. The viewer file fixes
both, always (the pipeline applies this by default; `--sort-energy` remains
as the modifier for standalone `--to-cart` / `--fix-markers`
re-processing):

1. it computes the Fock expectation energy of every orbital, E = <phi|F|phi>
   -- the diagonal Fock element, i.e. the physically meaningful energy of a
   localized orbital -- and writes it into `Ene=`;
2. it orders the MOs occupied-first (Occup > 1.0), then by ascending E. The
   occupied block stays at the top, so Avogadro's positional electron
   assignment still counts 2n electrons.

Inputs are the `<base>.S.txt` / `<base>.fock_ao.txt` dumps (JANPA `-doFock`
Fock + overlap in the spherical AO basis) and the `<base>.PURE` canonical
set -- all three written by `--clpo`. Gates, all before anything is
written: orthonormality of both orbital sets under S, the canonical
residual `F C = S C diag(eps)`, and an independent cross-check of E via the
canonical-MO expansion weights. A wrong or mismatched dump fails loudly
instead of writing numbers.

When the energy inputs are absent, it falls back to class order from the
labels in `<base>.JANPA` (needs `-RyOccPrintThreshold -1`, which `--clpo`
sets): `[BD + LP] [NB] [RY]`, i.e. bonding and lone pairs, then
antibonding, then Rydberg -- with the stdout-index→file-position mapping
verified against the occupancies. `Ene=` then holds sequential numbers, not
energies.

Example -- ethene (wB97X-D3/def2-TZVP), `--sort-energy` occupied block:
core x2 (-10.13 Ha), C=C sigma (-0.89), 4x C-H (-0.62), pi (-0.37); then
the C-C / C-H antibonding partners, then 72 Rydberg orbitals. Works on any
JANPA export (CLPO, LHO, NAO) because it only uses the exported AO
coefficients.

Caveats inherited from JANPA's export: the substrate's `[MO]` energies are
still sequential numbers (only the viewer file carries real energies), and
molden2molden's reader crashes on blank lines inside `[Atoms]` / a missing
blank at the end of `[GTO]` — outputs of this script already use the safe
layout with LF endings.

### Orbital-interaction analysis (`--e2`)

JANPA has no E(2) feature — its wiki page *"Can I compute the 'E(2) energy'
using the output of JANPA?"* explains why: the Fock matrix is only fully
meaningful for Hartree-Fock, under DFT it belongs to the auxiliary
Kohn-Sham system, and E(2) is not an observable. What JANPA does print is
an experimental charge-transfer table (donor -> acceptor pairs above a
fixed 0.01 e threshold). `--e2` computes both quantities for **every** pair
of any JANPA export and writes `<stem>_E2.txt`:

    python orca_to_janpa.py --e2 ethene_CLPO_spherical.molden
    # or in one shot: python orca_to_janpa.py ethene --clpo --e2

- `E2 = n_i F_ij^2/(F_jj - F_ii)` [kcal/mol] — the NBO-style second-order
  perturbation estimate in the localized basis. Read it with the caveat
  printed in the report: well-defined for HF, indicative for DFT (prefer
  `q` there).
- `q = D_ij^2/D_ii` [e] — the charge transferred between the two orbitals.
  This is exactly the number JANPA prints in its "Approximate charge
  transfer analysis"; `--e2` reproduces every printed pair (checked at run
  time whenever the log describes the export) and shows the pairs below
  JANPA's print threshold.
- Rows marked `*` have |F_ij|/(F_jj-F_ii) >= 0.25: the two orbitals are
  strongly mixed and the second-order estimate is not meaningful for them
  (often a sign the Lewis-like reference itself is inadequate — e.g. the
  3c-2e bridged ethyl cation surfaces as one strongly mixed pair instead of
  a hyperconjugation energy).

`--e2` gates everything the way `--sort-energy` does, plus an independent
cross-check of the whole Fock matrix against JANPA's own NAO data (the
`E2_pert` recipe: `<base>.fock_nao.txt` transformed with the set's chain —
`clpo2lho` @ `lho2nao` for CLPO, `lho2nao` for LHO, `aho2nao` for AHO,
`lpo2aho` @ `aho2nao` for LPO, the F_NAO dump itself for NAO), and it
refuses to write on mismatch. PNAO has no NAO-space path *and* is not
orthonormal: `--e2` refuses it with that explanation. JANPA's own CT
table only describes the CLPO set, so the printed-pair cross-check and
the labels apply there; other sets verify through the route-B chain and
report "labels/CT check skipped".

Numbers from the examples folder (wB97X-D3/def2-TZVP, so E2 is
**indicative** there — prefer `q`): isobutene sigma(C-H) -> pi*(C=C) ~5
kcal/mol per methyl C-H (q = 0.015 e); formaldehyde O lone pair ->
sigma*(C-H) ~29 kcal/mol (q = 0.059 e; this is the one channel where
CLPO and NBO genuinely diverge -- see [VALIDATION.md](VALIDATION.md));
tert-butyl cation sigma(C-H) -> empty-p ~34 kcal/mol x3 (q = 0.084 e);
water (HF, where E2 is directly meaningful) tops out at ~2 kcal/mol
(q = 0.0012 e).

### Which JANPA orbital set (all work with the viewer/analysis modes)

One flag per set — `--nao`, `--pnao`, `--lho`, `--aho`, `--lpo`, `--clpo`
(or `--all-sets`) — each producing the viewer + substrate pair described
above:

- **NAO** (natural atomic orbitals) — the orthonormal atomic set behind NPA
  charges, Wiberg bond indices, angular-momentum populations
  (`-NAO_Molden_File`, Fock: `<base>.fock_nao.txt`; its `--e2` route-B
  check compares the F_NAO dump directly).
- **PNAO** — pre-orthogonalization NAOs; NAO-construction intermediate.
  Normalized but mutually non-orthogonal by construction: viewer +
  substrate with the energy sort, but no pair analysis (see "Viewing
  orbitals").
- **LHO** (localized hybrid orbitals) — the atom-centred hybrids the CLPOs
  are built from (`-LHO_Molden_File`; `<base>.lho2nao.txt` = LHOs in the
  NAO basis).
- **AHO / LPO** — the LPO family (Int J Quantum Chem 2019, e25798);
  property-optimized localized orbitals and their atomic hybrids
  (`-AHO2NAO_File`, `-LPO2AHO_File` = their transformation chains).
- **CLPO** ("chemist's LPO") — the NBO-analog Lewis-like set (BD/NB/LP/RY);
  the one to use for bonding analysis, orbital visualization and the full
  `--e2` treatment (CLPO labels and the printed CT table live in
  `<base>.JANPA`).

The naive approach of just running `janpa` and reading the files also
needs: the janpa stdout (`<base>.JANPA`) for labels, and a `-doFock` run
for the matrices (`-Fock_AO_File`, `-Fock_NAO_File`, `-S_Matrix_File` plus
the transformation chain of whichever set: `-CLPO2LHO_File`,
`-LHO2NAO_File`, `-AHO2NAO_File`, `-LPO2AHO_File`; `-doFock` is a bare
flag).

## Limits

- **Closed-shell molecules, in practice.** The occupied/virtual split is an
  occupancy cut (Occup > 1.0), `--avogadro` writes integer `Occup= 2/0`,
  and `Spin=` is retagged from one uniform label — all closed-shell
  assumptions; open-shell references are out of scope.
- **No g shells.** `--to-cart` stops with that message (instead of a
  traceback) and the spherical file from `--fix-markers` is the artifact
  for such molecules.
- PNAO gets the viewer treatment but no pair analysis (see "Viewing
  orbitals").

## Project rules

- Single file, stdlib only: `orca_to_janpa.py` runs with plain `python`,
  nothing to install, no venv, global Python stays clean.
- Correctness fixes are on by default — not flags you have to remember.
- Every number the script writes is gated first; on a failed check it
  refuses to write rather than emit plausible-looking values.

## License

MIT — see [LICENSE](LICENSE).
