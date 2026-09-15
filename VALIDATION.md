# Cross-program validation: JANPA/CLPO analysis vs NBO 3.1

Purpose: establish that the orbital-interaction numbers this toolchain
produces (ORCA → JANPA → CLPO analysis) are comparable to an independent,
canonical implementation of the same physics (Gaussian's NBO 3.1). This
document is a validation record: the questions it answers are whether the
numbers are right in kind, and by how much the two implementations differ.
It is not a regression test suite.

## What is compared, and why in these tiers

NBO 3.1 and JANPA are not clones of each other: JANPA is an open
implementation of the **NPA** layer (NAO construction, populations,
Wiberg indices — the same well-defined algorithm), while for *localization*
NBO builds NHOs/NBOs with its own search and JANPA builds its own LPO
family (LHO/AHO/LPO/CLPO). So the comparison has two kinds of rows:

| tier | NBO program | JANPA side | expectation |
| --- | --- | --- | --- |
| populations | NPA algorithm | same algorithm reimplemented | **tight** (≤ ~0.005 e) |
| per-atom hybrids | NHO | LHO | role-analogous, qualitative |
| bonds / lone pairs | NBO | CLPO | same job, different machinery |
| E(2) table | Fock matrix in the NBO basis | this tool's E2 in the CLPO basis | **interpretive** (same partners/ranking, magnitudes within ~20–30 %) |

Agreement in the tight rows (populations, bond indices) validates the whole
ORCA → Molden → JANPA chain.
The E(2) row documents the cross-algorithm delta rather than asserting equality.

## Method

- **Level:** HF/cc-pVDZ, single point in both programs, identical fixed
  geometry in each case — the exact coordinates are in the input files
  under `validation/` (`.gjf` for G09, `.inp` for ORCA). HF because it is
  unambiguous across programs and E(2) has its cleanest meaning there.
- **Gaussian:** G09W Rev. B.01 (ships NBO 3.1), route
  `#p hf/cc-pvdz pop=nboread` plus a `$NBO BNDIDX $END` block.
  Two notes for reproducing these runs on G09W: the NBO key list is read
  only with `pop=nboread` (with plain `pop=nbo` the `$NBO` block is
  ignored), and `%Mem` is kept at 1000 MB, the reliable maximum for this
  32-bit build.
- **ORCA 6.1.1:** `! HF cc-pVDZ SP`, then this repo's pipeline:
  `python orca_to_janpa.py <base> --clpo --e2`.
- **Basis-convention notes.** 6-31G\* cannot be compared directly: Gaussian
  uses 6D cartesian d functions for Pople sets, ORCA 5D, and ORCA 6 rejects
  every simple cartesian keyword tried. cc-pVDZ is 5D/7F by default in both.
  ORCA's header prints a *decontracted* basis count for these cc-pVDZ
  runs (water: 30 listed vs the true SCF dimension of 24 — confirmed via
  the Molden conversion and the 1e-8 energy agreement).
- **Evidence layout.** Each `validation/<molecule>/` folder holds the
  exact inputs and outputs: the G09 `.gjf`/`.out`, the ORCA `.inp`/`.out`,
  the JANPA log, and the CLPO E2 table. `validation/compare.py`
  recomputes every table in this document from those files
  (`validation/compare_all.txt` is its saved output).

## Cross-molecule summary

| molecule | atoms / SCF functions | ΔE(RHF) | max Δq | max ΔW | max Δoccup | E(2), main comparison |
| --- | --- | --- | --- | --- | --- | --- |
| water | 3 / 24 | 1.1e-8 Ha | 0.0011 e | 0.0019 | 5e-4 | LP→H-Rydberg 2.33 vs 2.00 (14 %) |
| formaldehyde | 4 / 38 | 1.1e-8 Ha | 0.0040 e | 0.0032 | 1.6e-4 | LP→C-Rydberg 14.22 vs 15.90 (12 %); one interaction diverges (see below) |
| isobutene | 12 / 96 | 2.4e-9 Ha | 0.0010 e | 0.0008 | 5.6e-4 | four hyperconjugation interactions within 3 % |

## Water, HF/cc-pVDZ

Geometry (Å): O `0 0 0.0626`; H `−0.792 0 −0.4973`; H `0.792 0 −0.4973`.

### Energy

| quantity | G09 | this pipeline | delta | tolerance | verdict |
| --- | --- | --- | --- | --- | --- |
| E(RHF) / Ha | −76.0250404980 | −76.0250404875 | **1.1e-8** | ≤ 1e-5 | pass |

### Populations

| quantity | G09 NBO 3.1 | JANPA | delta | tolerance | verdict |
| --- | --- | --- | --- | --- | --- |
| NPA charge O | −0.92791 | −0.92902 | 0.0011 e | ≤ 0.005 | pass |
| NPA charge H | +0.46396 | +0.46451 | 0.0006 e | ≤ 0.005 | pass |
| Wiberg index O–H | 0.7871 | 0.7862 | 0.0009 | ≤ 0.01 | pass |
| Wiberg index H–H | 0.0020 | 0.0014 | 0.0006 | – | – |
| Wiberg total, O | 1.5742 | 1.5723 | 0.0019 | ≤ 0.01 | pass |

### Orbital occupancies

| orbital | G09 NBO | JANPA CLPO | delta |
| --- | --- | --- | --- |
| core (O 1s) | 1.99994 (CR) | 1.99997 | 0.00003 |
| O–H bond (×2) | 1.99834 | 1.99887 | 0.0005 |
| lone pairs (O) | 1.99756 / 1.99695 | 1.99789 / 1.99742 | ≤ 0.0005 |

### E(2)

| interaction | G09 | this pipeline | delta |
| --- | --- | --- | --- |
| LP → H-Rydberg | 2.33 (×2) | 2.00 (×2) | −14 % |
| second LP → H-Rydberg | 0.87 (×2) | 1.32 (×2) | same class |
| BD(O–H) → H-Rydberg | 0.64–0.67 (×4) | 0.42 / 0.76 | same class |

**Verdict:** same dominant interaction, same ranking, top pair 14 % apart —
inside the interpretive tolerance.

Internal consistency checks: orthonormality 2.9e-09, canonical
residual 9.9e-09, route A vs B 1.1e-08.

## Formaldehyde, HF/cc-pVDZ

Geometry: fixed single-point coordinates, identical in both input files
(`validation/formaldehyde/g09_formaldehyde_hf_ccpvdz.gjf` and
`orca_formaldehyde_hf_ccpvdz.inp`; 4 atoms: C, O, H, H).

### Energy and populations

| quantity | G09 | this pipeline | delta | tolerance | verdict |
| --- | --- | --- | --- | --- | --- |
| E(RHF) / Ha | −113.8768187500 | −113.8768187387 | 1.1e-8 | ≤ 1e-5 | pass |
| NPA charge C | +0.36698 | +0.37042 | 0.0034 e | ≤ 0.005 | pass |
| NPA charge O | −0.57961 | −0.58363 | 0.0040 e | ≤ 0.005 | pass |
| NPA charge H (×2) | +0.10631 | +0.10660 | 0.0003 e | ≤ 0.005 | pass |
| Wiberg C=O | 1.8806 | 1.8774 | 0.0032 | ≤ 0.01 | pass |
| Wiberg C–H (×2) | 0.9356 | 0.9357 | 0.0001 | ≤ 0.01 | pass |
| Wiberg totals C / O / H | 3.7518 / 1.9768 / 0.9951 | 3.7489 / 1.9735 / 0.9951 | ≤ 0.0033 | – | – |

### Orbital occupancies

| orbital | G09 NBO | JANPA CLPO | delta |
| --- | --- | --- | --- |
| π(C=O) | 1.99998 | 1.99993 | 0.00005 |
| σ(C=O) | 1.99961 | 1.99966 | 0.00005 |
| C–H (×2) | 1.99437 | 1.99441 | 0.00004 |
| core C / core O | 1.99971 / 1.99981 | 1.99974 / 1.99997 | ≤ 0.0002 |
| O lone pairs (σ / π) | 1.98636 / 1.90947 | 1.98627 / 1.90972 | ≤ 0.0003 |

### E(2)

| interaction | G09 | this pipeline | delta |
| --- | --- | --- | --- |
| in-plane LP(O) → C-Rydberg | 14.22 | 15.90 | +12 % |
| π LP(O) → C-Rydberg | 4.90 | 5.23 | +7 % |
| BD(C–H) → O-Rydberg (×2) | 2.28 | 2.19 | −4 % |
| BD(C–H) ↔ BD\*(C–H) (×2) | 1.21 | 1.20 | −1 % |
| in-plane LP(O) → σ\*(C–H) (×2) | 0.85 | **26.49** | see below |

**Documented divergence.** The LP(O)→σ\*(C–H) interaction is 0.85 kcal/mol
in the NBO basis and 26.5 kcal/mol in the CLPO basis. The cause is in the
coupling, not the energies: F_ij = 0.034 Ha (NBO) vs 0.159 Ha (CLPO)
— 4.7× — with similar ΔE (1.64 vs 1.15 Ha), and E2 ∝ F². The two sets'
"σ\*(C–H)" acceptors are *different functions*: NBO's BD\* is a
tightly-localized antibond, while the CLPO NB is more diffuse and overlaps
the oxygen lone pair far more strongly. This is a property of the two
localization schemes, not a numerical error — every other interaction
above, and the donor identification, agree. This is documented explicitly so the
divergence is not mistaken for a defect.

Internal consistency checks: orthonormality 8.3e-09, canonical
residual 1.0e-08, E cross-check 2.4e-08.

## Isobutene, HF/cc-pVDZ

Geometry: fixed single-point coordinates, identical in both input files
(`validation/isobutene/g09_isobutene_hf_ccpvdz.gjf` and
`orca_isobutene_hf_ccpvdz.inp`; 12 atoms).

### Energy and populations

| quantity | G09 | this pipeline | delta | tolerance | verdict |
| --- | --- | --- | --- | --- | --- |
| E(RHF) / Ha | −156.1241243460 | −156.1241243484 | 2.4e-9 | ≤ 1e-5 | pass |
| NPA charges, all 12 atoms | — | — | max 0.0010 e | ≤ 0.005 | pass |
| Wiberg C=C | 1.9256 | 1.9264 | 0.0008 | ≤ 0.01 | pass |
| Wiberg C–C (×2) | 1.0253 | 1.0253 | 0.0000 | ≤ 0.01 | pass |
| Wiberg C–H (×8) | 0.9354–0.9458 | 0.9353–0.9456 | ≤ 0.0002 | ≤ 0.01 | pass |
| Wiberg totals, all atoms | — | — | ≤ 0.002 | – | – |

Representative NPA charges (G09 / this pipeline): C1 `+0.00034 / −0.00021`;
methyl C `−0.59562 / −0.59636`; =CH₂ C `−0.41847 / −0.41743`;
H `+0.19 … +0.21` within 0.0003.

### Orbital occupancies

| orbital | G09 NBO | JANPA CLPO | delta |
| --- | --- | --- | --- |
| C1–C2 / C1–C3 σ | 1.98070 | 1.98061 | 0.00009 |
| σ(C1=C4) | 1.98498 | 1.98503 | 0.00005 |
| π(C1=C4) | 1.96415 | 1.96467 | 0.0005 |
| methyl C–H | 1.98294 / 1.99018 | 1.98300 / 1.99029 | ≤ 0.0001 |
| vinylic C–H (×2) | 1.98558 | 1.98614 | 0.0006 |
| cores | 1.99900–1.99937 | 1.99912–1.99938 | ≤ 0.0001 |

### E(2) — hyperconjugation interactions

| interaction | G09 | this pipeline | delta |
| --- | --- | --- | --- |
| vinylic C–H → σ\*(C–C) (×2) | 8.94 | 8.67 | −3 % |
| methyl C–H → σ\*(C–C) (×2) | 5.72 | 5.65 | −1.2 % |
| methyl C–H → π\*(C=C) (×4) | 5.47 | 5.51 | +0.7 % |
| C–C σ → σ\*(C–H) (×2) | 4.00 | 4.03 | +0.8 % |

F_ij and ΔE match to three digits on the methyl C–H → π\*(C=C)
interaction (0.064 / 0.93 Ha in both bases). The count of small
interactions differs (76 rows ≥ 0.5 in NBO's basis, 48 in this
pipeline's) because the Rydberg manifolds are constructed differently —
expected, and not part of the claim.

Internal consistency checks: orthonormality 6.3e-08, canonical residual
2.0e-08, E cross-check 1.3e-07.

## Limits

- Validates numbers of the kind compared here; **does not** validate the
  Molden file-format layer, the viewer fixes, or DFT (all runs here are HF).
- The E(2) rows compare two *different* localized bases. Interactions
  whose acceptors are built differently can differ by several-fold
  (formaldehyde LP→σ\*(C–H)); that is a property of the constructions, not
  an error signal. The robust claims are: same donors, same interaction
  ranking, and — when ΔE and F_ij agree — numerically matching energies.
- NBO 3.1 is a baseline, not a modern NBO-grade reference; and the JANPA
  side has no E(2) of its own — the E2 table is this tool's computation
  on CLPOs.
- Pinned versions: ORCA 6.1.1 / G09W B.01 / JANPA 2.02.

## Coverage

- water — tight tiers and E(2) compared
- formaldehyde — plus one documented E(2) divergence
- isobutene — its four hyperconjugation interactions within 3 %

Evidence: `validation/water/`, `validation/formaldehyde/`,
`validation/isobutene/` — inputs and outputs exactly as produced, plus
`compare.py` (regenerates every table) and `compare_all.txt`.
