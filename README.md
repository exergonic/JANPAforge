# orca-to-janpa

ORCA `.gbw` → JANPA → NBO Visualization Pipeline

**No `.47` file needed.**

## Usage (single script, stdlib only — plain `python`, nothing to install)

```powershell
cd C:/Users/mccan/Code/orca-to-janpa
python orca_to_janpa.py C:/Users/mccan/orca_calcs/NBO/ethene_NBO
python orca_to_janpa.py --diagnose C:/Users/mccan/orca_calcs/NBO/ethene_NBO.out
```

What it does per `<base>`:

- `[1/3]` runs `C:/ORCA_6.1.1/orca_2mkl.exe <base> -molden`
  (uses `<base>.mp2nos` via `-anyorbs` when present)
- `[2/3]` runs `molden2molden -fromorca3bf -orca3signs` → `<base>.PURE`
  (add `--dot47 file.47` for the correlated `-ds47` route)
- `[3/3]` runs `janpa -i <base>.PURE`, saves `<base>.JANPA`

Add `--clpo` to also export the CLPOs and the data the analysis modes need
(see "Orbital order" and "Orbital-interaction analysis" below):

```powershell
python orca_to_janpa.py ethene --clpo
# -> ethene_CLPO.molden, ethene.S.txt, ethene.fock_ao.txt,
#    ethene.fock_nao.txt, ethene.clpo2lho.txt, ethene.lho2nao.txt,
#    full CLPO labels in ethene.JANPA
python orca_to_janpa.py --to-cart ethene_CLPO.molden --avogadro --sort-energy
python orca_to_janpa.py --e2 ethene_CLPO.molden
```

Extra flags after `--` are passed to `janpa.jar`, e.g.:

```powershell
python orca_to_janpa.py ethene_NBO -- --npacharges charges.txt
```

## Verified

End-to-end on `ethene_NBO.gbw` (wB97X-D3/def2-TZVP, ORCA 6.1.1):
`orca_2mkl` → `molden2molden` → `janpa` converges, 16.000000 electrons,
charge sum 0.00000, CLPO C=C / C–H bonding graph correct. No `NPA`/`NBO`
keyword needed in the ORCA input for this route — a plain SCF `.gbw`
suffices.

## Viewing orbitals (Avogadro needs cartesian, MOrbVis needs markers)

JANPA writes spherical MOs. Two viewer quirks, both fixed by this script:

- **MOrbVis** honors `[5D]`/`[7F]`/`[9G]` (checked against its WebGPU
  evaluator source), but JANPA omits `[7F]` — f shells then evaluate as
  10-component cartesian and every orbital mistracks. Fix:
  `python orca_to_janpa.py --fix-markers ethene_CLPO.molden`
  (also drops JANPA's spurious `[9G]` when no g shells exist).
- **Avogadro** renders cartesian only (same reason `puream=0` was needed
  for avo_ibo). Fix:
  `python orca_to_janpa.py --to-cart ethene_CLPO.molden --avogadro`
  (d 5→6, f 7→10 in Molden order, markers dropped, `--spin Alpha`
  implied). `--avogadro` additionally reorders MOs occupied-first and
  writes integer `Occup= 2/0`: Avogadro parses `Occup` as int
  (1.986 → 1, miscounting 16 electrons as 8) and fills orbitals
  positionally by energy order, so fractional occupations with
  all-zero energies mislabel virtuals as occupied. The maps are derived numerically
  from Gaussian moments for normalized functions on both sides and
  cross-checked against `molden2molden -cart2pure` unit responses;
  every run is gated and the script refuses to write on mismatch.

### Orbital order and energies (`--sort-energy`)

A JANPA export carries no orbital energies -- every `Ene=` is a sequential
number (0, 1, 2, ...) and the MOs sit in JANPA's internal hybrid-pairing
order, so the viewer's orbital list is unreadable. `--sort-energy` (a
`--to-cart` / `--fix-markers` modifier) fixes both:

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

Caveats inherited from JANPA's export: `[MO]` energies are sequential
numbers unless `--sort-energy` rewrote them, and molden2molden's reader
crashes on blank lines inside `[Atoms]` / a missing blank at the end of
`[GTO]` — outputs of this script already use the safe layout with LF
endings.

### Orbital-interaction analysis (`--e2`)

JANPA has no E(2) feature — its wiki page *"Can I compute the 'E(2) energy'
using the output of JANPA?"* explains why: the Fock matrix is only fully
meaningful for Hartree-Fock, under DFT it belongs to the auxiliary
Kohn-Sham system, and E(2) is not an observable. What JANPA does print is
an experimental charge-transfer table (donor -> acceptor pairs above a
fixed 0.01 e threshold). `--e2` computes both quantities for **every** pair
of any JANPA export and writes `<stem>_E2.txt`:

    python orca_to_janpa.py --e2 ethene_CLPO.molden

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
`E2_pert` recipe: `<base>.fock_nao.txt` transformed with
`<base>.clpo2lho.txt` @ `<base>.lho2nao.txt`), and it refuses to write on
mismatch.

Numbers from this folder (wB97X-D3/def2-TZVP): isobutene sigma(C-H) ->
pi*(C=C) ~5 kcal/mol per methyl C-H; formaldehyde O lone pair ->
sigma*(C-H) ~29 kcal/mol; tert-butyl cation sigma(C-H) -> empty-p
~34 kcal/mol x3; water (HF) tops out at ~2 kcal/mol.

### Which JANPA orbital set (all work with the viewer/analysis modes)

- **NAO** (natural atomic orbitals) — the orthonormal atomic set behind NPA
  charges, Wiberg bond indices, angular-momentum populations
  (`-NAO_Molden_File`, Fock: `<base>.fock_nao.txt`).
- **PNAO** — pre-orthogonalization NAOs; NAO-construction intermediate.
- **LHO** (localized hybrid orbitals) — the atom-centred hybrids the CLPOs
  are built from (`-LHO_Molden_File`; `<base>.lho2nao.txt` = LHOs in the
  NAO basis).
- **AHO / LPO** — the LPO family (Int J Quantum Chem 2019, e25798);
  property-optimized localized orbitals and their atomic hybrids.
- **CLPO** ("chemist's LPO") — the NBO-analog Lewis-like set (BD/NB/LP/RY);
  the one to use for bonding analysis, orbital visualization and `--e2`.

The naive approach of just running `janpa` and reading the files also
needs: the janpa stdout (`<base>.JANPA`) for labels, and a `-doFock` run
for the matrices (`-Fock_AO_File`, `-Fock_NAO_File`, `-S_Matrix_File`,
`-CLPO2LHO_File`, `-LHO2NAO_File`; `-doFock` is a bare flag).

## Project rules

- Single file, stdlib only: `orca_to_janpa.py` runs with plain `python`,
  nothing to install, no venv, global Python stays clean.
