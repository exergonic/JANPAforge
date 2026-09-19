# JANPA orbital sets

JANPAforge is an ORCA → JANPA pipeline. It exports orbital sets JANPA
already builds. It does not implement a new localization algorithm.

## Why localize

Canonical molecular orbitals are eigenfunctions of a symmetry-adapted
Fock or Kohn–Sham operator. They spread over equivalent atoms by
construction. Bonding analysis wants objects centred on atoms and on
bonds: charges, two-centre occupancies, donor–acceptor pairs. That is
the reason to leave the canonical basis.

## Where JANPA sits

JANPA is a Weinhold-style natural / Lewis stack: atomic orbitals, then
hybrids, then a Lewis-like pairing. It is not a Foster–Boys /
Pipek–Mezey family of unitary rotations of the occupied MO block. IBO,
IAO, Boys, Edmiston–Ruedenberg, and Pipek–Mezey orbitals are out of
scope; this software does not produce them.

JANPA reimplements the Natural Population Analysis (NPA) layer —
natural atomic orbital (NAO) construction, populations, Wiberg
indices — the same well-defined algorithm as Natural Bond Orbitals
(NBO) (Nikolaienko, Bulavin, and Hovorun, *Comput. Theor. Chem.*
**1050** (2014) 15–22). On HF/cc-pVDZ water, formaldehyde, and
isobutene, NPA charges agree with the NBO 3.1 module of Gaussian 09
within 0.004 e and Wiberg indices within 0.004 (that module is not
authentic NBO — see the [caveat](VALIDATION.md)). Localization is
JANPA’s own family
of localized property-optimized orbitals (LPO) (Nikolaienko and
Bulavin, *Int. J. Quantum Chem.* **119** (2019) e25798). The
chemist’s LPO (CLPO) does the same job as NBO — bonding (BD),
antibonding (NB), lone pair (LP), and Rydberg (RY) — with different
machinery. Treat the bonding graph as comparable, the charge-transfer
numbers as the primary comparison, and the E(2) magnitudes as
interpretive ([VALIDATION.md](VALIDATION.md)).

## PNAO — `--pnao`

**What it is.** Pre-orthogonalization natural atomic orbitals. This is
the NAO-construction intermediate: each function is normalized, and
the set is not mutually orthogonal. On water, max |CᵀS C − I| =
4.8×10⁻¹. Occupancies do not sum to the electron count (14.54 e
against water’s 10).

**What it is for.** The construction step before orthogonalization.
The viewer still sorts by E = ⟨φ|F|φ⟩, a well-defined single-orbital
expectation (cross-check 5.2×10⁻⁹ on water).

**What it is not for.** Pair analysis. `--e2` refuses PNAO files.
Do not read `Occup` as an electron count.

**Files.** `--pnao` writes `<base>_PNAO.molden` (viewer: cartesian,
real Fock energies, occupied-first) and
`<base>_PNAO_spherical.molden` (substrate: JANPA’s export with
markers and `Spin=` corrected).

## NAO — `--nao`

**What it is.** Orthonormal natural atomic orbitals. This is the
atomic set behind NPA charges, Wiberg bond indices, and
angular-momentum populations.

**What it is for.** Atomic charges and bond indices. That is the tight
validation tier against NBO 3.1: populations ≤ 0.005 e, Wiberg ≤
0.01.

**What it is not for.** A Lewis bonding picture. NAOs are
atom-centred. They are not BD / LP objects. Bond pictures and the
labelled E(2) table belong on CLPO.

**Files.** `--nao` writes `<base>_NAO.molden` and
`<base>_NAO_spherical.molden`. `--e2` is defined on the substrate;
its route-B check compares `<base>.fock_nao.txt` directly.

## LHO — `--lho`

**What it is.** Localized hybrid orbitals. These are the atom-centred
hybrids that CLPOs are built from. Against NBO they are the natural
hybrid orbital (NHO) analogue: same role, qualitative agreement only.

**What it is for.** The hybrid picture before Lewis pairing.
`<base>.lho2nao.txt` is that set in the NAO basis.

**What it is not for.** The BD / NB / LP / RY graph. Pairing is a
CLPO step.

**Files.** `--lho` writes `<base>_LHO.molden` and
`<base>_LHO_spherical.molden`. `--e2` runs; JANPA’s printed
charge-transfer table and type labels describe CLPO, not this export.

## CLPO — `--clpo`

**What it is.** Chemist’s localized property-optimized orbitals: a
Lewis-like set labelled BD / NB / LP / RY (bonding, antibonding,
lone pair, Rydberg). This is the NBO analogue in the stack.

**What it is for.** The bonding graph, orbital pictures, and the full
`--e2` treatment (charge transfer first, E(2) second). CLPO labels and
JANPA’s printed charge-transfer table live in `<base>.JANPA`.
Occupied-block energies on ethene
(wB97X-D3/def2-TZVP): cores −10.13 Ha, C=C σ −0.89, four C–H −0.62,
π −0.37.

**What it is not for.** A drop-in NBO E(2) table. The robust
cross-program claims are the same donors and the same ranking.
Magnitudes can differ when the two schemes build different acceptors.
At HF/cc-pVDZ, formaldehyde LP(O)→σ\*(C–H) is 26.49 kcal/mol in the
CLPO basis and 0.85 kcal/mol in Gaussian's NBO 3.1 module, because
F_ij is 0.159 Ha
against 0.034 Ha (E2 ∝ F²). Isobutene’s four hyperconjugation
channels stay within 3% on the same geometry. Documented as a
property of the constructions, not as a numerical error
([VALIDATION.md](VALIDATION.md)).

**Files.** `--clpo` writes `<base>_CLPO.molden` and
`<base>_CLPO_spherical.molden`. `--e2` reads the substrate and writes
`<base>_CLPO_E2.txt`.

## AHO / LPO — `--aho`, `--lpo`

**What they are.** Localized property-optimized orbitals (Nikolaienko
and Bulavin, *Int. J. Quantum Chem.* **119** (2019) e25798). These
are the localized orbitals that best approximate the first-order
reduced density matrix (1-RDM) as a sum of one- and two-centre
terms (Frobenius norm). Atomic hybrid orbitals (AHO) are their
atomic hybrids. CLPO is the same construction with extra chemical
constraints, built from localized hybrid orbitals (LHOs) rather
than AHOs. The NAO-space chains are `<base>.aho2nao.txt` and
`<base>.lpo2aho.txt`.

**What they are for.** One-electron property partitioning into atom
and pair contributions. That is the question `--lpo` answers and
`--clpo` does not. LPO has only a localization constraint.
Unconstrained pairing of a lone-pair AHO with a weakly occupied
hybrid on another centre produces high-ionicity “pseudobonds” —
formally better for the 1-RDM, not a Lewis structure (Nikolaienko,
*Phys. Chem. Chem. Phys.* **21** (2019) 5285–5294).

**What they are not for.** The BD / NB / LP / RY Lewis graph. Those
labels and JANPA’s printed CT table describe CLPO. CLPO caps BD/NB
ionicity at 0.90 and keeps BD/LP occupancy above 1.0 so that one
orbital remains per electron pair.

**Files.** `--aho` writes `<base>_AHO.molden` and
`<base>_AHO_spherical.molden`. `--lpo` writes `<base>_LPO.molden` and
`<base>_LPO_spherical.molden`. `--e2` runs on either substrate;
labels and charge-transfer checks are skipped.

## Which flag

If the question is NPA charges or Wiberg indices, `--nao`. Export
PNAO only to inspect the construction intermediate; it is not a
second charge set. Hybrids before Lewis pairing: `--lho`. A
one-electron property partitioned into atom and pair terms, without
Lewis constraints: `--lpo` (hybrids: `--aho`). Bonding graph,
pictures, and q / E(2) pairs: `--clpo`. Every set flag writes the
same pair — cartesian viewer, spherical substrate — and `--e2`
takes the substrate. `--all-sets` writes all six.

## Limits

The pipeline is closed-shell SCF in practice. The occupied/virtual
split is Occup > 1.0. `--avogadro` writes integer `Occup = 2/0`.
`Spin=` is retagged from one uniform label. Open-shell references
are out of scope.

q = D_ij² / D_ii is the primary number: a density statement — the
charge transferred between the two orbitals — and the number JANPA
prints in its own charge-transfer table; it keeps its meaning under
DFT. E(2) = n_i F_ij² / (F_jj − F_ii) kcal/mol is the companion
estimate: well-defined at HF; under DFT the operator is the auxiliary
Kohn–Sham matrix, so E(2) is indicative. JANPA itself has no E(2)
feature; the table is this tool’s computation.

Pairs with |F_ij| / (F_jj − F_ii) ≥ 0.25 are marked `*`. The two
orbitals are strongly mixed and the second-order estimate is not
meaningful. The bridged ethyl cation surfaces as one such pair,
not as a hyperconjugation energy.

PNAO has no pair analysis.

The NBO 3.1 record in [VALIDATION.md](VALIDATION.md) is HF/cc-pVDZ
only, and its NBO side comes from Gaussian's non-authentic NBO 3.1
module ([caveat](VALIDATION.md)). It does not validate DFT, the Molden
viewer layer, or a modern NBO-grade reference.
