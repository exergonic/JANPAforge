# orca-to-janpa

ORCA `.gbw` → JANPA pipeline for HF/DFT. **No `.47` file needed.**

## Why you get no `.47` file (the short version)

Nothing is wrong with your inputs. Your `.out` files show the tell-tale
signature:

```
Now starting NBO....
/c/Users/mccan/orca_calcs/NBO
Stopping NBO...-------
```

That middle line is the stdout of `pwd` — your `NBOEXE`/`GENEXE` dummy.
The `NBOEXE=pwd` trick dates to the **ORCA 3.0.x era**, when ORCA wrote
the `.47` file itself (via its internal `gennbo`) and only needed *some*
zero-exit program in `NBOEXE` to proceed. The JANPA wiki page you followed
says exactly that: tested on 3.0.0 / 3.0.2, with the trick documented for
"post-3.0.0 (but not v4.x)".

Since ORCA 3.1 the `!NPA` / `!NBO` keywords drive the **licensed Weinhold
NBO6/NBO7 binary**. The ORCA 6.1 manual (§5.2) states `NBOEXE` must be the
real `nbo7` executable (i4 integer build). A placeholder exits 0, so ORCA
takes it as success and moves on — writing **no `.47` file and no NBO
analysis**. That is what you see.

To get a true `.47` from ORCA 6.1.1 you need the NBO license from the
University of Wisconsin. There is no free knob that re-enables the old
behavior.

## Good news: for HF/DFT you don't need the `.47` at all

The JANPA wiki's HF/DFT route never touches a `.47`:

1. `orca_2mkl <base> -molden` → `<base>.molden.input` (SCF MOs + basis)
2. `molden2molden -fromorca3bf -orca3signs` → `<base>.PURE` (conventional Molden)
3. `janpa -i <base>.PURE` → NPA / CLPO analysis

The `-ds47` option (density from a `.47`) is only for **correlated**
densities: MP2-relaxed, QCISD/CCSD `orbopt`. And MP2 has its own
`.47`-free path: run with `%mp2 Density relaxed / NatOrbs true`, then feed
`orca_2mkl` the resulting `.mp2nos` instead of the `.gbw` (this tool picks
`.mp2nos` automatically when present). Only CI/CC-orbopt densities truly
strand you without the NBO license.

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

Add `--clpo` to also export the CLPOs and the data the viewer-order step
needs (see "Orbital order" below):

```powershell
python orca_to_janpa.py ethene --clpo
# -> ethene_CLPO.molden, ethene.S.txt, ethene.fock_ao.txt,
#    full CLPO labels in ethene.JANPA
python orca_to_janpa.py --to-cart ethene_CLPO.molden --avogadro --sort-energy
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

## Project rules

- Single file, stdlib only: `orca_to_janpa.py` runs with plain `python`,
  nothing to install, no venv, global Python stays clean.
