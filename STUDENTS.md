# If your homework said “NBO”

You are in the right place. This page is the survival note.

Your instructor said “look at the NBOs.” You have an ORCA file, a script
named `orca_to_janpa.py`, and a sudden need to know what any of that means. That
is a reasonable emergency. The rest of this repo can wait.

## Two pictures of the same molecule

Gen chem drew bonds as lines between atoms. That Lewis picture is still
the one your homework wants: a pair of electrons in a bond, a lone pair on an atom.

The ORCA calculation does not hand you that picture. It solves for a
"wavefunction" (whatever that is lol) and
then reports "canonical" (whatever that means) molecular orbitals that are
spread over the entire molecule. They are not wrong. They just aren't “the C–H bond.”

Localization is the step that cuts the wavefunction into objects that
look like bonds and lone pairs. That is the whole point of the
assignment.

## Your instructor said NBO. This tool’s matching set is CLPO.

Natural Bond Orbitals (NBO) are Weinhold’s Lewis-like set and are a
specific way to perform the localization. This pipeline
does not run that program. It exports "chemist’s localized
property-optimized orbitals" (CLPO) from JANPA.

CLPO does the same *job* as NBO:

- **BD** — a bonding pair (two atoms)
- **LP** — a lone pair (one atom)
- **NB** — an antibond (empty; the name is “antibonding,” not
  “nonbonding,” which is unhelpful and not your fault)
- **RY** — Rydberg (also empty; leftover on one atom)

Occupancy near 2 means “this orbital holds a pair.” Occupancy near 0
means “this one is empty.” That is the occupancy lesson.

Same job is not the same program. Do not copy a CLPO number into a
report as “NBO” unless your instructor said that is acceptable.

```bash
python orca_to_janpa.py path/to/molecule --clpo --avogadro
```

Open `<base>_CLPO.molden`. That is the viewer file.

## If they asked for charges

Pictures of bonds are CLPO. Atomic charges and Wiberg bond indices
(how much bond exists between two atoms) come from natural atomic
orbitals (NAO).

```bash
python orca_to_janpa.py path/to/molecule --nao
```

The numbers live in `<base>.JANPA`. You do not need the orbital pictures
for a charge table.

## What not to touch

You do not need six flags. You need `--clpo`. Add `--nao` if the
homework asked for charges.

`--pnao`, `--lho`, `--aho`, `--lpo`, and the E(2) pair table are real.
They are also not this assignment. If the prompt did not name them,
skip them. [REVIEW.md](REVIEW.md) is the next page, not this one.

If something still looks insane, it may be. Ask your instructor before
you invent a story for it.
