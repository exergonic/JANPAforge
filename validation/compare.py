"""Compare G09/NBO 3.1 against this pipeline's JANPA/CLPO analysis.

NOTE: "NBO 3.1" is the module bundled with Gaussian -- not authentic NBO;
see the caveat at the top of ../VALIDATION.md (Weinhold, J. Comput. Chem.
47 (2026) e70374).

Recomputes the tables in ../VALIDATION.md from the inputs and outputs in
the per-molecule folders next to this script:

    python compare.py [molecule ...]      # default: all molecules

Expected layout per molecule folder (see VALIDATION.md):
    g09_<mol>_hf_ccpvdz.out      Gaussian output (incl. NBO sections)
    orca_<mol>_hf_ccpvdz.out     ORCA output
    janpa_<mol>_hf_ccpvdz.JANPA  janpa stdout from the pipeline run
    clpo_<mol>_hf_ccpvdz_E2.txt  pair-interaction table from --e2
"""
import re
import sys
from pathlib import Path

BASE = Path(__file__).resolve().parent
MOLS = ("water", "formaldehyde", "isobutene")

# Each NBO section is a fixed-format block: the scans below start a few
# lines past the section header and stop at the section's own end marker;
# the large slice bounds are only safety caps against format drift.
#
# Charge row: element, index, then five columns -- charge, core, valence,
# rydberg, total.  Only the charge is used here.
RE_G09_CHARGE = re.compile(r"^\s*([A-Za-z]{1,2})\s+(\d+)\s+(-?\d+\.\d+)")
# Wiberg matrix row: index. element  values...
RE_G09_WIBERG = re.compile(r"^\s*(\d+)\.\s+([A-Za-z]{1,2})\s+(.*)$")
# E(2) row: idx. donor / idx. acceptor   E(2)  dE  F(i,j)
RE_G09_E2 = re.compile(r"^\s*\d+\.\s+(\S.*?)\s+/\s+\d+\.\s+(.*?)\s{2,}"
                       r"(\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$")
# Occupied-orbital row of the NBO summary: idx. CR|BD|LP (n) label  occ  E
RE_G09_OCC = re.compile(r"^\s*(\d+)\.\s+(CR|BD|LP)\s+\(\s*\d+\)\s+"
                        r"(.*?)\s{2,}(\d+\.\d+)")
# Pipeline atom labels are element+index with no separator: C1, H11.
RE_PIPE_LABEL = re.compile(r"^[A-Za-z]+\d+$")
# Pipeline pair row: # donor -> acceptor   F_ij  dE  q  E2; only E2 is used.
RE_PIPE_E2 = re.compile(r"^\s*\d+\s+(\S+)\s+->\s+(\S+)\s+-?[\d.]+\s+"
                        r"[\d.]+\s+[\d.]+\s+([\d.]+)\s*$")


def read_g09(path):
    lines = path.read_text(errors="replace").splitlines()
    out = {"energy": None, "charges": [], "wiberg": [], "wiberg_totals": [],
           "occupancies": [], "e2": []}
    for i, line in enumerate(lines):
        if "SCF Done" in line:
            out["energy"] = float(line.split("=")[1].split()[0])
        if "Summary of Natural Population Analysis" in line:
            for l2 in lines[i + 2:i + 60]:
                m = RE_G09_CHARGE.match(l2)
                if m:
                    out["charges"].append((m.group(1), int(m.group(2)),
                                           float(m.group(3))))
                elif out["charges"] and l2.strip().startswith("="):
                    break
        if "Wiberg bond index matrix" in line:
            # The matrix prints in column blocks (9 wide in this G09 build),
            # so an atom's row label appears once per block.
            cols = []
            for l2 in lines[i + 1:]:
                if "Wiberg bond index, Totals" in l2:
                    break
                header = re.match(r"^Atom\s+([\d\s]+)$", l2.strip())
                if header:
                    cols = [int(x) for x in header.group(1).split()]
                    continue
                m = RE_G09_WIBERG.match(l2)
                if m and cols:
                    row = int(m.group(1)) - 1
                    while len(out["wiberg"]) <= row:
                        out["wiberg"].append({})
                    for col, value in zip(cols, m.group(3).split()):
                        out["wiberg"][row][col - 1] = float(value)
        if "Wiberg bond index, Totals" in line:
            for l2 in lines[i + 4:i + 4 + 60]:
                m = RE_G09_WIBERG.match(l2)
                if m:
                    out["wiberg_totals"].append(float(m.group(3).split()[0]))
                elif out["wiberg_totals"]:
                    break
        if "Natural Bond Orbitals (Summary)" in line:
            for l2 in lines[i + 5:i + 5 + 300]:
                m = RE_G09_OCC.match(l2)
                if m:
                    out["occupancies"].append((m.group(2), int(m.group(1)),
                                               m.group(3).strip(),
                                               float(m.group(4))))
                elif "RY*" in l2 and out["occupancies"]:
                    break
        if "Second Order Perturbation" in line:
            for l2 in lines[i + 6:i + 6 + 600]:
                m = RE_G09_E2.match(l2)
                if m:
                    out["e2"].append((m.group(1), m.group(2),
                                      float(m.group(3))))
                elif "Natural Bond Orbitals" in l2:
                    break
    return out


def read_pipeline(folder, mol):
    out = {"energy": None, "charges": [], "wiberg": [], "wiberg_totals": [],
           "occupancies": [], "e2": [], "e2_totals": ""}
    orca_lines = (folder / f"orca_{mol}_hf_ccpvdz.out").read_text(
        errors="replace").splitlines()
    for line in orca_lines:
        if "FINAL SINGLE POINT ENERGY" in line:
            out["energy"] = float(line.split()[-1])
    janpa_lines = (folder / f"janpa_{mol}_hf_ccpvdz.JANPA").read_text(
        errors="replace").splitlines()
    for i, line in enumerate(janpa_lines):
        if "Final electron populations and NPA charges" in line:
            # Row: center, nuclear charge, electron population, NMB
            # population, NPA charge -- the charge is the last field.
            for l2 in janpa_lines[i + 3:i + 3 + 80]:
                parts = l2.split()
                if len(parts) == 5 and RE_PIPE_LABEL.match(parts[0]):
                    out["charges"].append((parts[0], float(parts[4])))
                elif out["charges"]:
                    break
        if "Wiberg-Mayer bond indices" in line:
            # Row i lists columns i..n (upper triangle); the diagonal
            # itself prints in parentheses:  ( 1.5723)  0.7862  0.7862
            for l2 in janpa_lines[i + 2:i + 2 + 80]:
                parts = l2.replace("(", " ").replace(")", " ").split()
                if parts and parts[0].isdigit():
                    row = int(parts[0])
                    while len(out["wiberg"]) < row:
                        out["wiberg"].append([])
                    out["wiberg"][row - 1].extend(float(v) for v in parts[1:])
                elif out["wiberg"]:
                    break
        if "Summary of CLPO results" in line:
            for l2 in janpa_lines[i + 2:i + 2 + 200]:
                cols = l2.split("\t")
                if len(cols) >= 3 and cols[0].strip().isdigit():
                    try:
                        occ = float(cols[2].strip())
                    except ValueError:
                        continue
                    if occ > 0.5:       # occupied CLPOs are ~2, virtuals ~0
                        out["occupancies"].append((int(cols[0]),
                                                   cols[1].strip(), occ))
            break
    for line in (folder / f"clpo_{mol}_hf_ccpvdz_E2.txt").read_text(
            errors="replace").splitlines():
        m = RE_PIPE_E2.match(line)
        if m:
            out["e2"].append((m.group(1), m.group(2), float(m.group(3))))
        elif line.startswith("totals"):
            out["e2_totals"] = line
    return out


def main():
    which = sys.argv[1:] or list(MOLS)
    for mol in which:
        folder = BASE / mol
        g09 = read_g09(folder / f"g09_{mol}_hf_ccpvdz.out")
        pipeline = read_pipeline(folder, mol)
        print(f"\n{'=' * 72}\n{mol.upper()}  (HF/cc-pVDZ, fixed geometry)\n{'=' * 72}")
        print(f"E(RHF)  G09 {g09['energy']:.10f}   pipeline {pipeline['energy']:.10f}   "
              f"delta {abs(g09['energy'] - pipeline['energy']):.2e} Ha")
        print("\n-- NPA charges (idx: G09 / pipeline / delta) --")
        worst = 0.0
        for (sym, idx, q_g09), (_, q_pipe) in zip(
                sorted(g09["charges"], key=lambda r: r[1]),
                sorted(pipeline["charges"],
                       key=lambda r: int(re.sub(r"\D", "", r[0])))):
            dq = q_pipe - q_g09
            worst = max(worst, abs(dq))
            print(f"   {sym}{idx:<3} {q_g09:+.5f} / {q_pipe:+.5f} / {dq:+.5f}")
        print(f"   max |delta q| = {worst:.4f} e")
        print("\n-- Wiberg bond indices (off-diagonal, G09 > 0.05) --")
        worst = 0.0
        for i in range(len(g09["wiberg"])):
            for j in range(i + 1, len(g09["wiberg"])):
                w_g09 = g09["wiberg"][i].get(j)
                if w_g09 is None:
                    continue
                row = pipeline["wiberg"][i] if i < len(pipeline["wiberg"]) else []
                # JANPA's row i lists columns i..n, so column j is at j - i.
                w_pipe = row[j - i] if 0 <= j - i < len(row) else None
                if w_pipe is None:
                    print(f"   {i + 1}-{j + 1}: {w_g09:.4f} /    n/a")
                else:
                    worst = max(worst, abs(w_pipe - w_g09))
                    if w_g09 > 0.05 or w_pipe > 0.05:
                        print(f"   {i + 1}-{j + 1}: {w_g09:.4f} / "
                              f"{w_pipe:.4f} / {w_pipe - w_g09:+.4f}")
        print(f"   max |delta W| (all pairs) = {worst:.4f}")
        print("   totals: " + " ".join(
            f"{i + 1}:{g09['wiberg_totals'][i]:.4f}/"
            f"{pipeline['wiberg'][i][0]:.4f}"
            for i in range(min(len(g09["wiberg_totals"]),
                               len(pipeline["wiberg"])))))
        print("\n-- occupied orbitals: G09 NBO --")
        for kind, idx, label, occ in g09["occupancies"]:
            print(f"   {idx:>3}. {kind:<3} {label:<28} {occ:.5f}")
        print("-- occupied orbitals: CLPO --")
        for idx, label, occ in pipeline["occupancies"][:14]:
            print(f"   {idx:>3}. {label:<44} {occ:.5f}")
        print(f"\n-- E(2), G09 (first 12 of {len(g09['e2'])} rows, "
              f"not ranked) --")
        for donor, acceptor, e2 in g09["e2"][:12]:
            print(f"   {e2:6.2f}  {donor} -> {acceptor}")
        print(f"-- E(2), pipeline (first 12 table rows, NB-acceptor "
              f"block first; "
              f"{sum(1 for r in pipeline['e2'] if r[2] >= 0.5)} pairs "
              f">= 0.5 kcal/mol) --")
        for donor, acceptor, e2 in pipeline["e2"][:12]:
            print(f"   {e2:6.2f}  {donor} -> {acceptor}")
        print(f"   {pipeline['e2_totals']}")


if __name__ == "__main__":
    main()
