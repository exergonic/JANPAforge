#!/usr/bin/env python3
"""ORCA (.gbw) -> JANPA pipeline for HF/DFT (no .47 file needed).

Standalone script -- stdlib only, no installs, no venv needed:
``python orca_to_janpa.py <base>``.

For closed-shell HF/DFT the JANPA route does NOT use ORCA's NBO
interface at all:

    1. ``orca_2mkl <base> -molden``      -> ``<base>.molden.input``
    2. ``molden2molden -fromorca3bf -orca3signs`` -> ``<base>.PURE``
    3. ``janpa -i <base>.PURE``          -> NPA/CLPO analysis on stdout

The ``.47`` file (with ``-ds47``) is only needed for correlated
densities (MP2 relaxed / QCISD / CCSD orbopt).  On ORCA >= 3.1 the
``NPA``/``NBO`` keywords require the licensed Weinhold NBO6/NBO7
binary (``NBOEXE``); the old ``NBOEXE=pwd`` dummy-exe trick from the
JANPA wiki (ORCA 3.0.x era) no longer produces a ``.47`` file --
ORCA just logs ``Now starting NBO....`` + the dummy output and moves
on.  See ``diagnose_nbo_output()`` and README.md.

Viewer workflow: one pipeline flag per JANPA orbital set (``--clpo``,
``--lho``, ``--aho``, ``--lpo``, ``--nao``, ``--pnao``, or ``--all-sets``)
writes the spherical substrate (``<base>_<SET>_spherical.molden``: JANPA's
export with corrected markers and spin, kept as the analysis input) and
the viewer file (``<base>_<SET>.molden``: cartesian, real Fock energies,
occupied-first order; integer Occup with ``--avogadro`` until the upstream
Avogadro occupancy bug is fixed).  The standalone ``--to-cart`` /
``--fix-markers`` / ``--sort-energy`` / ``--e2`` modes re-process existing
files.  See README.md.
"""

from __future__ import annotations

import argparse
import math
import re
import shutil
import subprocess
import sys
from operator import mul
from pathlib import Path
from typing import NamedTuple

DEFAULT_ORCA_DIR = Path("C:/ORCA_6.1.1")
DEFAULT_JANPA_DIR = Path("C:/Users/mccan/Code/third-party/JANPA")

MOLDEN2MOLDEN_JAR = "molden2molden.jar"
JANPA_JAR = "janpa.jar"


def find_orca_2mkl(orca_dir: Path = DEFAULT_ORCA_DIR) -> Path:
    """Locate orca_2mkl.exe (prefers the given dir, falls back to PATH)."""
    cand = orca_dir / "orca_2mkl.exe"
    if cand.is_file():
        return cand
    found = shutil.which("orca_2mkl")
    if found:
        return Path(found)
    raise FileNotFoundError(
        f"orca_2mkl not found in {orca_dir} nor on PATH. "
        "Pass --orca-dir C:/ORCA_6.1.1"
    )


def find_jars(janpa_dir: Path = DEFAULT_JANPA_DIR) -> tuple[Path, Path]:
    """Return (molden2molden.jar, janpa.jar), raising if missing."""
    m2m = janpa_dir / MOLDEN2MOLDEN_JAR
    jp = janpa_dir / JANPA_JAR
    missing = [str(p) for p in (m2m, jp) if not p.is_file()]
    if missing:
        raise FileNotFoundError(f"JANPA jar(s) missing: {missing}")
    return m2m, jp


def resolve_base(target: str | Path) -> tuple[Path, Path]:
    """Resolve user target to (workdir, basename-without-extension).

    Accepts a basename, a .gbw path, a .mp2nos path, or an .inp path.
    Returns the directory holding the calculation and the stem to use
    as the ORCA BaseName.
    """
    p = Path(str(target))
    if p.suffix.lower() in {".gbw", ".mp2nos", ".inp", ".out"}:
        return p.parent.resolve(), p.stem
    # bare basename, possibly with directory
    if p.parent != Path("."):
        return p.parent.resolve(), p.name
    return Path.cwd().resolve(), p.name


def orbitals_source(workdir: Path, base: str) -> tuple[str, list[str]]:
    """Pick the orbital file for orca_2mkl.

    Returns (label, extra_argv).  Per the JANPA wiki, MP2 natural
    orbitals live in <base>.mp2nos and must be used instead of the
    .gbw for MP2 densities.  orca_2mkl reads those via -anyorbs.
    """
    mp2nos = workdir / f"{base}.mp2nos"
    if mp2nos.is_file():
        return ("mp2nos", ["-anyorbs", str(mp2nos)])
    return ("gbw", [])


def run_orca_2mkl(workdir: Path, base: str, orca_2mkl: Path) -> Path:
    """Run ``orca_2mkl <base> -molden``; return the .molden.input path."""
    label, extra = orbitals_source(workdir, base)
    gbw = workdir / f"{base}.gbw"
    if label == "gbw" and not gbw.is_file():
        raise FileNotFoundError(f"No {gbw} (and no {base}.mp2nos) in {workdir}")
    argv = [str(orca_2mkl), base, "-molden", *extra]
    r = subprocess.run(
        argv, cwd=str(workdir), stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, text=True, timeout=600,
    )
    molden = workdir / f"{base}.molden.input"
    if r.returncode != 0 or not molden.is_file():
        raise RuntimeError(
            f"orca_2mkl failed (rc={r.returncode}) using {label} source.\n"
            f"cmd: {' '.join(argv)}\n{molden} not written.\n{r.stdout[-3000:]}"
        )
    return molden


def run_molden2molden(
    workdir: Path,
    molden_input: Path,
    output: Path,
    m2m_jar: Path,
    dot47: Path | None = None,
) -> Path:
    """Convert ORCA-style Molden to conventional Molden (.PURE).

    Without ``dot47`` this is the HF/DFT path (no density replacement).
    With ``dot47`` (a .47 / .mdcip.47 file) the SCF MOs are replaced by
    natural orbitals from the correlated density (MP2/CC/QCISD route).
    """
    argv = [
        "java", "-jar", str(m2m_jar),
        "-i", str(molden_input),
        "-o", str(output),
        "-fromorca3bf", "-orca3signs",
    ]
    if dot47 is not None:
        if not dot47.is_file():
            raise FileNotFoundError(f".47 file not found: {dot47}")
        argv += ["-ds47", str(dot47)]
    r = subprocess.run(
        argv, cwd=str(workdir), stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, text=True, timeout=600,
    )
    if r.returncode != 0 or not output.is_file():
        raise RuntimeError(
            f"molden2molden failed (rc={r.returncode}).\n"
            f"cmd: {' '.join(argv)}\n{r.stdout[-3000:]}"
        )
    if "Data loaded successfully" not in r.stdout:
        print(
            "WARNING: molden2molden output lacks success marker; "
            "check the log above.",
            file=sys.stderr,
        )
    return output


def run_janpa(
    workdir: Path,
    pure_file: Path,
    janpa_jar: Path,
    out_file: Path | None = None,
    extra_args: list[str] | None = None,
) -> str:
    """Run JANPA on the .PURE file; return stdout (also saved if asked)."""
    argv = ["java", "-jar", str(janpa_jar), "-i", str(pure_file)]
    if extra_args:
        argv += extra_args
    r = subprocess.run(
        argv, cwd=str(workdir), stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, text=True, timeout=1200,
    )
    if r.returncode != 0:
        raise RuntimeError(
            f"janpa failed (rc={r.returncode}).\n"
            f"cmd: {' '.join(argv)}\n{r.stdout[-4000:]}"
        )
    if "Total number of electrons" not in r.stdout:
        raise RuntimeError(
            "janpa ran but its output lacks the electron-count marker "
            "-- the .PURE file may be corrupt.\n" + r.stdout[-4000:]
        )
    if out_file is not None:
        out_file.write_text(r.stdout, encoding="utf-8")
    return r.stdout


def diagnose_nbo_output(out_path: Path) -> str:
    """Explain why an ORCA .out did or did not yield a .47 file.

    Returns a human-readable diagnosis string.
    """
    lines: list[str] = []
    text = out_path.read_text(encoding="utf-8", errors="replace")
    workdir = out_path.parent.resolve()
    stem = out_path.stem
    f47 = sorted(workdir.glob(f"{stem}*.47"))
    lines.append(f"ORCA output : {out_path}")
    lines.append(
        ".47 files   : "
        + (", ".join(p.name for p in f47) if f47 else "(none found)")
    )
    if "Now starting NBO" not in text:
        lines.append(
            "NBO section : absent -- the input had neither !NPA/!NBO nor "
            "a %nbo block, so ORCA never attempted the NBO interface."
        )
        return "\n".join(lines)
    # NBO interface was attempted; classify the outcome.
    idx = text.index("Now starting NBO")
    window = text[idx : idx + 2000]
    if "N A T U R A L" in window or "NATURAL POPULATIONS" in window:
        lines.append(
            "NBO section : REAL NBO output present -- a licensed NBO6/NBO7 "
            "binary ran. A .47/ARCHIVE file should exist alongside the .gbw."
        )
    elif "Stopping NBO" in window:
        # Dummy-exe signature: the only line between the markers is the
        # dummy program's own stdout (e.g. `pwd` prints the directory).
        between = window.split("Now starting NBO", 1)[1].split(
            "Stopping NBO", 1
        )[0].strip()
        lines.append(
            "NBO section : DUMMY-EXE signature -- between 'Now starting "
            "NBO....' and 'Stopping NBO...' ORCA printed only:\n"
            f"    {between!r}\n"
            "That is the stdout of whatever NBOEXE/GENEXE pointed at "
            "(e.g. `pwd`), which exits 0 without doing NBO analysis. "
            "ORCA takes the zero exit code as success and continues, "
            "writing no .47 file and no NBO analysis."
        )
        lines.append(
            "Why        : the NBOEXE=pwd trick dates to the ORCA 3.0.x era, "
            "when ORCA generated the .47 itself (gennbo) and only needed a "
            "zero-exit placeholder. Since ORCA 3.1 the NPA/NBO keywords "
            "drive the licensed Weinhold NBO6/NBO7 binary; a placeholder "
            "cannot produce the .47/ARCHIVE. The ORCA 6.1 manual states "
            "NBOEXE must be the real nbo7 executable (i4 integer build)."
        )
    else:
        lines.append(
            "NBO section : 'Now starting NBO' found but outcome unclear; "
            "search the .out for 'NBO' and error lines manually."
        )
    if not f47:
        lines.append(
            "Next step  : for HF/DFT you do NOT need the .47 -- run this "
            "script's gbw->molden->PURE->janpa pipeline instead. The .47 "
            "(-ds47) is only for correlated densities (MP2 relaxed, "
            "CCSD/QCISD orbopt), which on ORCA 6 require the NBO license "
            "(or use the .mp2nos path for MP2, see README)."
        )
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Spherical <-> cartesian Molden helpers (Molden convention).
#
# JANPA writes spherical MOs ([5D]/[7F] markers, 5-component d and
# 7-component f in Molden order).  Viewers that assume cartesian
# functions (notably Avogadro) misrender these files, exactly as a
# Psi4 ``puream=1`` Molden would.  ``convert_to_cart()`` expands the
# [MO] coefficients to 6-component d / 10-component f (Molden cartesian
# order) and drops the spherical markers -- the analogue of rerunning
# with ``puream=0``.
#
# Convention (verified against molden2molden -cart2pure unit-vector
# responses and the JGints bytecode): both sides use normalized
# functions, pure = normalized solid harmonics, cart = normalized
# cartesian Gaussians, and cart->pure is the overlap projection
# T[r][k] = <phi_r|chi_k>.  Normalized cartesians of one shell are NOT
# mutually orthogonal (e.g. <xx|yy> = 1/3), so textbook orthonormal
# matrices do not apply; the maps below are derived numerically from
# Gaussian moments (exponent-independent).  Every conversion is gated
# on projection error and overlap-norm drift; the script fails loudly
# instead of writing a bad file.
# ---------------------------------------------------------------------------

SPHER_N = {"s": 1, "p": 3, "d": 5, "f": 7, "g": 9}
CART_N = {"s": 1, "p": 3, "d": 6, "f": 10, "g": 15}


# Cartesian monomials per L in Molden order, as (lx, ly, lz) exponents.
CART_MONOMIALS = {
    "d": [(2, 0, 0), (0, 2, 0), (0, 0, 2),
          (1, 1, 0), (1, 0, 1), (0, 1, 1)],
    "f": [(3, 0, 0), (0, 3, 0), (0, 0, 3),
          (1, 2, 0), (2, 1, 0), (2, 0, 1),
          (1, 0, 2), (0, 1, 2), (0, 2, 1), (1, 1, 1)],
}

# Normalized solid harmonics (Molden order) as {monomial: coeff} direction
# combos. Same polynomial set as the Molden manual and MOrbVis's
# evaluator (verified term-by-term); exact normalization and the maps
# below are derived numerically, because normalized cartesian Gaussians
# of one shell are NOT mutually orthogonal (e.g. <xx|yy> = 1/3), so
# textbook orthonormal-basis matrices do not apply directly.
SPHER_COMBOS = {
    "d": [{(2, 0, 0): -1.0, (0, 2, 0): -1.0, (0, 0, 2): 2.0},
          {(1, 0, 1): 1.0},
          {(0, 1, 1): 1.0},
          {(2, 0, 0): 1.0, (0, 2, 0): -1.0},
          {(1, 1, 0): 1.0}],
    "f": [{(0, 0, 3): 2.0, (2, 0, 1): -3.0, (0, 2, 1): -3.0},
          {(1, 0, 2): 4.0, (3, 0, 0): -1.0, (1, 2, 0): -1.0},
          {(0, 1, 2): 4.0, (2, 1, 0): -1.0, (0, 3, 0): -1.0},
          {(2, 0, 1): 1.0, (0, 2, 1): -1.0},
          {(1, 1, 1): 1.0},
          {(3, 0, 0): 1.0, (1, 2, 0): -3.0},
          {(2, 1, 0): 3.0, (0, 3, 0): -1.0}],
}


def _dfact(n: int) -> float:
    r = 1.0
    while n > 1:
        r *= n
        n -= 2
    return r


def _moment(a: int) -> float:
    """<x^a> under a Gaussian weight, up to one common factor.

    Only ratios enter normalized overlaps, so the bare double factorial
    suffices (odd moments vanish)."""
    if a % 2:
        return 0.0
    return _dfact(a - 1) if a > 0 else 1.0


def _mono_overlap(m1: tuple, m2: tuple) -> float:
    return (_moment(m1[0] + m2[0]) * _moment(m1[1] + m2[1])
            * _moment(m1[2] + m2[2]))


def _invert(A: list[list[float]]) -> list[list[float]]:
    n = len(A)
    M = [row[:] + [1.0 if i == j else 0.0 for j in range(n)]
         for i, row in enumerate(A)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(M[r][col]))
        M[col], M[piv] = M[piv], M[col]
        d = M[col][col]
        if abs(d) < 1e-12:
            raise RuntimeError("singular overlap in pure/cart map builder")
        M[col] = [v / d for v in M[col]]
        for r in range(n):
            if r != col and M[r][col] != 0.0:
                f = M[r][col]
                M[r] = [a - f * b for a, b in zip(M[r], M[col])]
    return [row[n:] for row in M]


def _build_L_maps(L: str):
    """(E, T, S): pure->cart expansion, cart->pure projection (overlap),
    cartesian overlap -- all for normalized functions.

    phi_r = sum_k E[r][k] chi_k with E = S^-1 T^T; cartesian MO
    coefficients follow as C_k = sum_r E[r][k] c_r.  Exponent-
    independent, so built once at import.
    """
    monos = CART_MONOMIALS[L]
    n = len(monos)
    combos = SPHER_COMBOS[L]
    m = len(combos)
    nrm = [math.sqrt(_mono_overlap(x, x)) for x in monos]
    S = [[_mono_overlap(a, b) / (na * nb)
          for b, nb in zip(monos, nrm)] for a, na in zip(monos, nrm)]
    T = []
    for combo in combos:
        n2 = sum(ca * cb * _mono_overlap(ma, mb)
                 for ma, ca in combo.items()
                 for mb, cb in combo.items())
        nphi = math.sqrt(n2)
        T.append([sum(ca * _mono_overlap(ma, mk)
                      for ma, ca in combo.items()) / (nphi * nk)
                  for mk, nk in zip(monos, nrm)])
    Sinv = _invert(S)
    E = [[sum(Sinv[k][j] * T[r][j] for j in range(n)) for k in range(n)]
         for r in range(m)]
    return E, T, S


_MAPS = {L: _build_L_maps(L) for L in ("d", "f")}


def spher_d_to_cart(c: list[float]) -> list[float]:
    """5 spherical d -> 6 cartesian (Molden xx,yy,zz,xy,xz,yz)."""
    E, _, _ = _MAPS["d"]
    return [sum(E[r][k] * c[r] for r in range(5)) for k in range(6)]


def spher_f_to_cart(c: list[float]) -> list[float]:
    """7 spherical f -> 10 cartesian (Molden xxx..xyz order)."""
    E, _, _ = _MAPS["f"]
    return [sum(E[r][k] * c[r] for r in range(7)) for k in range(10)]


def _check_block(L: str, c_spher: list[float],
                 C_cart: list[float]) -> tuple[float, float]:
    """Projection error max|T.C - c| and relative norm drift
    |C'SC - c'c| / c'c."""
    _, T, S = _MAPS[L]
    back = [sum(T[r][k] * C_cart[k] for k in range(len(C_cart)))
            for r in range(len(c_spher))]
    dc = max(abs(a - b) for a, b in zip(c_spher, back))
    n0 = sum(a * a for a in c_spher)
    n1 = sum(C_cart[i] * S[i][j] * C_cart[j]
             for i in range(len(C_cart)) for j in range(len(C_cart)))
    dn = abs(n1 - n0) / max(n0, 1e-30)
    return dc, dn


def molden_shell_ls(path: Path) -> list[str]:
    """Angular-momentum letter of each contracted shell in [GTO] order."""
    ls: list[str] = []
    in_gto = False
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if line.startswith("["):
            in_gto = line.upper() == "[GTO]"
            continue
        if not in_gto or not line:
            continue
        m = re.match(r"^([A-Za-z])\s+\d+", line)
        if m and m.group(1).lower() in SPHER_N:
            ls.append(m.group(1).lower())
    if not ls:
        raise ValueError(f"No [GTO] shells parsed from {path}")
    return ls


def _is_marker(line: str) -> bool:
    return line.strip().upper() in ("[5D]", "[7F]", "[9G]")


def convert_to_cart(in_path: Path, out_path: Path,
                    spin: str | None = None,
                    avogadro: bool = False,
                    sort: "SortPlan | None" = None) -> str:
    """Rewrite a spherical Molden file with cartesian d/f MO coefficients.

    The [GTO] section is copied verbatim (same shells/exponents/
    contractions); only [MO] coefficient vectors are expanded and the
    [5D]/[7F]/[9G] markers dropped.  Returns a validation report; raises
    on round-trip failure instead of writing a bad file.

    With avogadro=True, MOs are reordered occupied-first and Occup is
    rewritten as integers (2/0, threshold occ>1.0): Avogadro parses
    Occup as int (1.986 -> 1, undercounting electrons 16 -> 8) and
    fills orbitals positionally, so fractional occupations mislabel
    the occupied/virtual sets.

    With ``sort`` (a SortPlan from ``plan_orbital_order``) the MO blocks
    are written in the plan's order and each ``Ene=`` is replaced by the
    plan's value -- real Fock expectation energies, or the class-order
    sequential numbers when no energy data exists.
    """
    shell_ls = molden_shell_ls(in_path)
    if "g" in shell_ls:
        raise NotImplementedError("g shells not implemented for --to-cart")
    expect = sum(SPHER_N[L] for L in shell_ls)
    before, blocks, after = split_molden_mo_section(in_path)
    worst_dc = 0.0
    worst_dn = 0.0
    mo_blocks: list[tuple[list[str], list[float], float]] = []
    for header, coeff_lines in blocks:
        coeffs = [float(l.split()[1]) for l in coeff_lines]
        if len(coeffs) != expect:
            raise ValueError(
                f"MO #{len(mo_blocks) + 1}: {len(coeffs)} coeffs but "
                f"basis has {expect} spherical functions -- wrong file?")
        cart: list[float] = []
        pos = 0
        for L in shell_ls:
            n = SPHER_N[L]
            c = coeffs[pos:pos + n]
            pos += n
            if L in ("s", "p"):
                cart.extend(c)
                continue
            cc = spher_d_to_cart(c) if L == "d" else spher_f_to_cart(c)
            cart.extend(cc)
            dc, dn = _check_block(L, c, cc)
            worst_dc = max(worst_dc, dc)
            worst_dn = max(worst_dn, dn)
        occup = 0.0
        for h in header:
            hs = h.strip()
            if hs.startswith("Occup="):
                try:
                    occup = float(hs.split("=", 1)[1])
                except ValueError:
                    pass
        mo_blocks.append((list(header), cart, occup))

    if sort is not None:
        if len(sort.perm) != len(mo_blocks):
            raise RuntimeError(
                f"sort plan covers {len(sort.perm)} orbitals but the file "
                f"has {len(mo_blocks)} -- plan/file mismatch"
            )
        ordered = [mo_blocks[i] for i in sort.perm]
    elif avogadro:
        ordered = ([b for b in mo_blocks if b[2] > 1.0]
                   + [b for b in mo_blocks if b[2] <= 1.0])
    else:
        ordered = mo_blocks
    n_mo = len(ordered)
    out: list[str] = list(before)
    out.append("[MO]")
    for pos, (header, coeffs, occup) in enumerate(ordered):
        is_occ = avogadro and occup > 1.0
        for h in header:
            hs = h.strip()
            if sort is not None and hs.startswith("Ene="):
                out.append(f" Ene= {sort.enes[pos]:.14E}")
            elif hs.startswith("Occup=") and avogadro:
                out.append(f" Occup= {2 if is_occ else 0}")
            elif hs.startswith("Spin="):
                out.append(f"Spin= {spin}" if spin else h)
            else:
                out.append(h)
        for i, v in enumerate(coeffs, 1):
            out.append(f"  {i:<6d}{v: .12f}".rstrip())
    out.extend(after)

    report = (f"{n_mo} MOs -> "
              f"{sum(CART_N[L] for L in shell_ls)} cartesian functions; "
              f"projection max|T.C-c|={worst_dc:.2e}, "
              f"norm drift={worst_dn:.2e}")
    if avogadro:
        nocc = sum(1 for b in mo_blocks if b[2] > 1.0)
        kept = [b[2] for b in mo_blocks if b[2] > 1.0]
        drop = [b[2] for b in mo_blocks if b[2] <= 1.0]
        report += (f"; Avogadro layout: {nocc} doubly-occupied + "
                   f"{n_mo - nocc} virtual (threshold occ>1.0")
        if kept and drop:
            report += (f", min kept={min(kept):.4f}, "
                       f"max dropped={max(drop):.4f}")
        report += f"; m_electrons={2 * nocc})"
        true_e = sum(b[2] for b in mo_blocks)
        if abs(2 * nocc - true_e) > 0.5:
            report += (f" [note: sum(Occup)={true_e:.3f} but the integer "
                       f"layout makes Avogadro count {2 * nocc}]")
    if worst_dc > 1e-9 or worst_dn > 1e-9 or n_mo == 0:
        raise RuntimeError(f"VALIDATION FAILED: {report}")
    out_path.write_text("\n".join(out) + "\n", encoding="utf-8", newline="\n")
    return report


def split_molden_mo_section(
    path: Path,
) -> tuple[list[str], list[tuple[list[str], list[str]]], list[str]]:
    """Split a Molden file into (before, MO blocks, after).

    Blocks are (header_lines, coefficient_lines) per MO, in file order;
    lines before/after the [MO] section are returned verbatim.  Spherical
    markers ([5D]/[7F]/[9G]) are dropped everywhere -- callers re-emit
    their own.
    """
    before: list[str] = []
    after: list[str] = []
    blocks: list[tuple[list[str], list[str]]] = []
    header: list[str] = []
    coeffs: list[str] = []
    in_mo = False

    def flush() -> None:
        if header or coeffs:
            blocks.append((list(header), list(coeffs)))
            header.clear()
            coeffs.clear()

    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if _is_marker(line):
            continue
        if not in_mo:
            if line.startswith("[MO]"):
                in_mo = True
            else:
                before.append(raw)
            continue
        if line.startswith("[") or line.startswith("Sym="):
            flush()
            if line.startswith("["):
                in_mo = False
                after.append(raw)
            else:
                header.append(raw)
            continue
        if line == "":
            continue
        if line.startswith(("Ene=", "Occup=", "Spin=")):
            header.append(raw)
            continue
        parts = line.split()
        if len(parts) >= 2:
            try:
                float(parts[1])
                coeffs.append(raw)
                continue
            except ValueError:
                pass
        header.append(raw)
    flush()
    return before, blocks, after


def fix_spherical_markers(in_path: Path, out_path: Path,
                          spin: str | None = None,
                          sort: "SortPlan | None" = None) -> str:
    """Ensure [5D]/[7F]/[9G] markers match the shells actually present.

    JANPA omits [7F] (breaking MOrbVis f shells) and writes a spurious
    [9G] with no g shells (tripping Avogadro).  Markers are rewritten
    from the [GTO] content: [5D] if d or f shells exist, [7F] if f
    shells exist, [9G] only if g shells exist.  With ``sort`` (a
    SortPlan) the MO blocks are reordered and their ``Ene=`` values
    replaced just like in ``convert_to_cart``.
    """
    shell_ls = molden_shell_ls(in_path)
    want = []
    if "d" in shell_ls or "f" in shell_ls:
        want.append("[5D]")
    if "f" in shell_ls:
        want.append("[7F]")
    if "g" in shell_ls:
        want.append("[9G]")
    before, blocks, after = split_molden_mo_section(in_path)
    if sort is not None:
        if len(sort.perm) != len(blocks):
            raise RuntimeError(
                f"sort plan covers {len(sort.perm)} orbitals but the file "
                f"has {len(blocks)} -- plan/file mismatch"
            )
        blocks = [blocks[i] for i in sort.perm]
    out: list[str] = list(before)
    out.extend(want)
    out.append("[MO]")
    for pos, (header, coeffs) in enumerate(blocks):
        for h in header:
            hs = h.strip()
            if hs.startswith("Spin=") and spin:
                out.append(f"Spin= {spin}")
            elif sort is not None and hs.startswith("Ene="):
                out.append(f" Ene= {sort.enes[pos]:.14E}")
            else:
                out.append(h)
        out.extend(coeffs)
    out.extend(after)
    out_path.write_text("\n".join(out) + "\n", encoding="utf-8", newline="\n")
    return (f"markers now {want or '(none)'} "
            f"for shells {sorted(set(shell_ls))}")


# ---------------------------------------------------------------------------
# Viewer orbital ordering for JANPA exports.
#
# JANPA's CLPO/LHO/NAO exports carry no orbital energies: every ``Ene=`` is
# a sequential number and the MO blocks sit in JANPA's internal
# hybrid-pairing order, which makes the Avogadro/MOrbVis orbital list
# unreadable.  The physically meaningful energy of a localized orbital is
# the expectation value of the Fock operator, E_i = <phi_i|F|phi_i> -- the
# diagonal element of the Fock matrix in the localized basis.
#
# Two routes, tried in order by ``plan_orbital_order()``:
#
#   energy route -- needs the Fock and overlap matrices in the AO basis:
#     ``janpa -doFock -Fock_AO_File <base>.fock_ao.txt -S_Matrix_File
#     <base>.S.txt`` dump both (spherical AO basis, [GTO] function order --
#     the same convention as the molden files).  JANPA builds F from the
#     input file's MO coefficients and ``Ene=`` values, so the dumps MUST
#     come from a run on the canonical ``<base>.PURE``; the dumps are gated
#     against that file (F C = S C diag(eps) on the canonical MOs) and
#     against the export (orthonormality under S), and E is cross-checked
#     against sum_k |<phi_i|MO_k>|^2 eps_k -- all before anything is
#     written.  The pipeline set flags (e.g. ``--clpo``) produce every
#     input.
#
#   class route -- needs only the janpa stdout saved as ``<base>.JANPA``:
#     its CLPO summary labels every orbital (BD)/(NB)/(LP)/(RY).  Order =
#     [BD+LP] [NB] [RY], original file order inside each class.  The stdout
#     index -> file position mapping is verified against the occupancies.
#
# Both routes put the occupied block (Occup > 1.0) first, so Avogadro's
# positional occupied/virtual assignment stays correct.
# ---------------------------------------------------------------------------

HARTREE_EV = 27.211386245988

SORT_ORTHO_TOL = 1e-4   # max |C^T S C - I| accepted for either orbital set
SORT_FOCK_TOL = 1e-3    # max |F C - S C eps| accepted (relative to scale)
SORT_DE_TOL = 1e-3      # max |E(Fock) - sum_k w eps_k|, Hartree
LABEL_OCC_TOL = 2e-5    # occupancy agreement when mapping stdout -> file

_EXPORT_SUFFIXES = ("clpo", "lho", "nao", "pnao", "aho", "lpo",
                    "cart", "fixed", "spherical")
_CLASS_RANK = {"BD": 0, "LP": 0, "NB": 1, "RY": 2}

# The JANPA Molden exports this script can carry through the viewer/analysis
# pipeline.
class SetExport(NamedTuple):
    """One JANPA orbital set: tag, janpa option, route-B chain, capabilities.

    ``chain`` is the NAO-space transformation chain that route B of --e2
    uses (an empty tuple means the F_NAO dump itself is the check, i.e.
    the NAO set; None means no route-B check is possible).  ``orthonormal``
    is False for the pre-orthogonalization intermediate PNAO: its orbitals
    are normalized but mutually non-orthogonal, so the orthonormality gate
    is relaxed for it (normalized orbitals still have a well-defined
    E = <phi|F|phi>, so the energy sort stays valid).  ``pair`` is False
    where the pair-interaction analysis (E2, q) is not defined: PNAO only.
    """

    tag: str
    opt: str
    chain: tuple[str, ...] | None
    orthonormal: bool
    pair: bool


_SET_EXPORTS = {
    "clpo": SetExport("CLPO", "-CLPO_Molden_File",
                      ("clpo2lho", "lho2nao"), True, True),
    "lho":  SetExport("LHO", "-LHO_Molden_File", ("lho2nao",), True, True),
    "aho":  SetExport("AHO", "-AHO_Molden_File", ("aho2nao",), True, True),
    "lpo":  SetExport("LPO", "-LPO_Molden_File", ("lpo2aho", "aho2nao"),
                      True, True),
    "nao":  SetExport("NAO", "-NAO_Molden_File", (), True, True),
    "pnao": SetExport("PNAO", "-PNAO_Molden_File", None, False, False),
}


def _set_from_stem(stem: str) -> SetExport | None:
    """The SetExport whose tag appears in a filename stem, if any."""
    low = stem.lower()
    for se in _SET_EXPORTS.values():
        if f"_{se.tag.lower()}" in low:
            return se
    return None

# janpa dump options: key -> (option, filename pattern).  The S / Fock dumps
# are shared by every set; the transformation chains feed the route-B
# cross-check of --e2 (NAO-space chains, each verified against route A).
_DUMP_EXPORTS = {
    "s":        ("-S_Matrix_File", "{b}.S.txt"),
    "fock_ao":  ("-Fock_AO_File", "{b}.fock_ao.txt"),
    "fock_nao": ("-Fock_NAO_File", "{b}.fock_nao.txt"),
    "clpo2lho": ("-CLPO2LHO_File", "{b}.clpo2lho.txt"),
    "lho2nao":  ("-LHO2NAO_File", "{b}.lho2nao.txt"),
    "aho2nao":  ("-AHO2NAO_File", "{b}.aho2nao.txt"),
    "lpo2aho":  ("-LPO2AHO_File", "{b}.lpo2aho.txt"),
}


def requested_sets(args) -> list[str]:
    """The orbital-set keys requested on the command line, in fixed order."""
    if getattr(args, "all_sets", False):
        return list(_SET_EXPORTS)
    return [k for k in _SET_EXPORTS if getattr(args, k)]


class SortPlan(NamedTuple):
    """Output order for one molden orbital set (see plan_orbital_order)."""

    perm: list[int]      # perm[i] = input MO block index -> output slot i
    enes: list[float]    # Ene= value per output slot (Ha, or sequential)
    report: str          # human-readable audit trail


def derive_base(stem: str) -> str:
    """'ethene_CLPO' -> 'ethene'; strips known export/conversion suffixes."""
    parts = stem.split("_")
    while len(parts) > 1 and parts[-1].lower() in _EXPORT_SUFFIXES:
        parts.pop()
    return "_".join(parts)


def molden_spin_labels(path: Path) -> list[str]:
    """The 'Spin=' label of every MO block, in file order."""
    labels: list[str] = []
    for raw in path.read_text(encoding="utf-8",
                             errors="replace").splitlines():
        line = raw.strip()
        if line.startswith("Spin="):
            labels.append(line.split("=", 1)[1].strip())
    return labels


def resolve_spin_label(in_path: Path, pure: Path | None = None) -> str | None:
    """Spin= value to retag a uniform closed-shell set with, or None.

    JANPA writes ``Spin= Beta`` on every orbital even for closed-shell
    sets; the canonical source (``<base>.PURE`` next to the export, or the
    given ``pure`` path) carries the calculation's own labeling.  Retag
    only when both sides are uniform, so open-shell sets are never touched.
    """
    labels = set(molden_spin_labels(in_path))
    if len(labels) != 1:
        return None
    if pure is None:
        pure = in_path.with_name(f"{derive_base(in_path.stem)}.PURE")
    if not pure.is_file():
        return None
    src = set(molden_spin_labels(pure))
    if len(src) != 1:
        return None
    value = src.pop()
    return value if value != labels.pop() else None


def read_matrix_dump(path: Path) -> list[list[float]]:
    """Read a JANPA matrix export (title / 'n m' / column labels / n rows)."""
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    try:
        n = int(lines[1].split()[0])
    except (IndexError, ValueError):
        raise RuntimeError(f"{path}: not a JANPA matrix dump")
    rows: list[list[float]] = []
    for raw in lines[2:]:
        if len(rows) == n:
            break
        toks = raw.rstrip().split("\t")
        if len(toks) < n:
            continue
        try:
            rows.append([float(t) for t in toks[:n]])
        except ValueError:
            continue
    if len(rows) != n:
        raise RuntimeError(f"{path}: parsed {len(rows)} of {n} matrix rows")
    return rows


def molden_mo_vectors(path: Path) -> list[tuple[float, float, list[float]]]:
    """(Ene, Occup, coefficients) per MO block, in file order."""
    _before, blocks, _after = split_molden_mo_section(path)
    out: list[tuple[float, float, list[float]]] = []
    for header, coeff_lines in blocks:
        if not coeff_lines:
            continue
        coeffs = [float(l.split()[1]) for l in coeff_lines]
        ene = 0.0
        occ = 0.0
        for h in header:
            hs = h.strip()
            if hs.startswith("Ene="):
                ene = float(hs.split("=", 1)[1])
            elif hs.startswith("Occup="):
                occ = float(hs.split("=", 1)[1])
        out.append((ene, occ, coeffs))
    if not out:
        raise RuntimeError(f"{path}: no [MO] blocks found")
    return out


def clpo_summary_from_log(log_path: Path) -> dict:
    """{index: (class, display_label, occupancy)} from the CLPO summary.

    ``class`` is 'BD'/'NB'/'LP'/'RY'/'??' and ``display_label`` looks like
    'C1-H3:BD' (same style the CT analysis prints).  Empty dict when the
    log has no CLPO summary block.
    """
    text = log_path.read_text(encoding="utf-8", errors="replace")
    if "*** Summary of CLPO results" not in text:
        return {}
    block = text.split("*** Summary of CLPO results", 1)[1]
    block = block.split("Number of two-center", 1)[0]
    out: dict[int, tuple[str, str, float]] = {}
    for line in block.splitlines():
        cells = line.split("\t")
        if len(cells) < 3:
            continue
        head = cells[0].strip().split()
        if not head or not head[0].isdigit():
            continue
        try:
            occ = float(cells[2].strip())
        except ValueError:
            continue
        desc = cells[1].strip()
        if "antibonding" in desc:
            atoms, cls = desc.split(",")[0].strip(), "NB"
        elif "(BD)" in desc:
            atoms = desc.split("(BD)", 1)[1].split(",")[0].strip()
            cls = "BD"
        elif "(LP)" in desc:
            atoms = desc.split("(LP)", 1)[1].strip()
            cls = "LP"
        elif "(RY)" in desc:
            atoms = desc.replace("(RY)", "").strip()
            cls = "RY"
        else:
            atoms, cls = desc, "??"
        label = f"{atoms}:{cls}" if cls != "??" else atoms
        out[int(head[0])] = (cls, label, occ)
    return out


def clpo_classes_from_log(log_path: Path):
    """({index: 'BD'|'LP'|'NB'|'RY'}, {index: occupancy}) from janpa stdout.

    Returns (None, None) when the log has no CLPO summary block.
    """
    summary = clpo_summary_from_log(log_path)
    if not summary:
        return None, None
    labels = {i: v[0] for i, v in summary.items()}
    occs = {i: v[2] for i, v in summary.items()}
    return labels, occs


def clpo_ct_pairs_from_log(log_path: Path) -> list[tuple[int, float, int]]:
    """JANPA's own charge-transfer table, for cross-checking ``--e2``.

    Returns [(donor_index, charge, acceptor_index), ...], empty when the
    CT section is absent.
    """
    text = log_path.read_text(encoding="utf-8", errors="replace")
    if "Approximate charge transfer" not in text:
        return []
    sec = text.split("Approximate charge transfer", 1)[1]
    sec = sec.split("IntErfragment", 1)[0]
    pat = re.compile(r"^\s*(\d+)\s+.+?\s+[\d.-]+\s+-->\s+([\d.]+)\s+-->"
                     r"\s+[\d.]+\s+.+?\s+(\d+)\s*$")
    pairs = []
    for line in sec.splitlines():
        mt = pat.match(line)
        if mt:
            pairs.append((int(mt.group(1)), float(mt.group(2)),
                          int(mt.group(3))))
    return pairs


def _matmul(A: list[list[float]], B: list[list[float]]) -> list[list[float]]:
    return [[sum(map(mul, row, col)) for col in zip(*B)] for row in A]


def _max_dev_from_identity(M: list[list[float]]) -> float:
    worst = 0.0
    for i, row in enumerate(M):
        for j, v in enumerate(row):
            worst = max(worst, abs(v - (1.0 if i == j else 0.0)))
    return worst


def _load_pair_inputs(target: Path, pure: Path, s_path: Path, f_path: Path):
    """Parse and validate the shared inputs of the energy analyses.

    Returns (S, F, eps, Cm, Ct, occ_mo, occ_tgt) with Cm/Ct the canonical
    and target MO block columns in AO order.
    """
    S = read_matrix_dump(s_path)
    F = read_matrix_dump(f_path)
    n = len(S)
    if len(F) != n or any(len(row) != n for row in S + F):
        raise RuntimeError(f"{s_path} / {f_path}: not square {n}x{n} dumps")
    tgt = molden_mo_vectors(target)
    src = molden_mo_vectors(pure)
    if len(tgt) != n or len(src) != n:
        raise RuntimeError(
            f"dimension mismatch: {target.name} has {len(tgt)} MOs, "
            f"{pure.name} has {len(src)}, dumps are {n}x{n} -- dumps and "
            "files belong to different calculations?"
        )
    for name, blocks in ((target.name, tgt), (pure.name, src)):
        for k, (_, _, c) in enumerate(blocks):
            if len(c) != n:
                raise RuntimeError(
                    f"{name}: MO {k + 1} has {len(c)} coefficients, "
                    f"expected {n}"
                )
    eps = [s[0] for s in src]
    Cm = [[src[k][2][i] for k in range(n)] for i in range(n)]  # [ao][mo]
    Ct = [[tgt[k][2][i] for k in range(n)] for i in range(n)]
    occ_mo = [s[1] for s in src]
    occ_tgt = [t[1] for t in tgt]
    return S, F, eps, Cm, Ct, occ_mo, occ_tgt


def _fock_gates(S, F, eps, Cm, Ct):
    """Orthonormality + canonical-residual numbers, shared by the analyses.

    Gates both orbital sets under S and the canonical residual
    F C = S C diag(eps) -- the dumps must describe the .PURE.  Returns
    (S*Ct, max|Cm'SCm - I|, max|Ct'SCt - I|, max|F C - S C eps|,
    max|F C|); each caller applies its own tolerances so its error text
    can name the analysis.
    """
    n = len(S)
    SCm = _matmul(S, Cm)
    SCt = _matmul(S, Ct)
    ortho_mo = _max_dev_from_identity(_matmul(list(zip(*Cm)), SCm))
    ortho_tgt = _max_dev_from_identity(_matmul(list(zip(*Ct)), SCt))
    FCm = _matmul(F, Cm)
    resid = 0.0
    scale = 0.0
    for a in range(n):
        for k in range(n):
            resid = max(resid, abs(FCm[a][k] - SCm[a][k] * eps[k]))
            scale = max(scale, abs(FCm[a][k]))
    return SCt, ortho_mo, ortho_tgt, resid, scale


def fock_orbital_energies(target: Path, pure: Path, s_path: Path,
                          f_path: Path,
                          require_ortho: bool = True) -> tuple[list[float], str]:
    """E_i = <phi_i|F|phi_i> for every MO block of ``target``, fully gated.

    S and F come from JANPA's ``-doFock`` dumps (spherical AO basis, [GTO]
    function order).  Gates: orthonormality of both orbital sets under S
    (the target's is skipped when ``require_ortho`` is False -- PNAO is
    normalized but not mutually orthogonal by construction), the canonical
    residual F C = S C diag(eps) on the .PURE MOs, and an independent
    cross-check of E via the canonical-MO expansion weights.  Raises
    instead of returning untrustworthy numbers.
    """
    S, F, eps, Cm, Ct, _, _ = _load_pair_inputs(target, pure, s_path, f_path)
    n = len(S)

    SCt, ortho_mo, ortho_tgt, resid, scale = _fock_gates(S, F, eps, Cm, Ct)
    M = _matmul(list(zip(*Ct)), _matmul(F, Ct))  # F in the target basis
    energies = [M[i][i] for i in range(n)]

    # Independent cross-check: E_i = sum_k |<phi_i|MO_k>|^2 eps_k.
    U = _matmul(list(zip(*Cm)), SCt)           # [canonical k][target i]
    worst_de = 0.0
    for i in range(n):
        e2 = sum(U[k][i] * U[k][i] * eps[k] for k in range(n))
        worst_de = max(worst_de, abs(e2 - energies[i]))

    se = _set_from_stem(target.stem)
    tag = se.tag if se else "target"
    if require_ortho:
        oth = f"{tag} orthonormality {ortho_tgt:.2e}"
    else:
        oth = (f"{tag} non-orthogonality {ortho_tgt:.2e} (expected; "
               "normalized but mutually non-orthogonal, not gated)")
    report = (f"gates: {oth} | "
              f"MO orthonormality {ortho_mo:.2e} | "
              f"F C - S C eps {resid:.2e} (scale {scale:.2f}) | "
              f"E cross-check {worst_de:.2e} Ha")
    if ((require_ortho and ortho_tgt > SORT_ORTHO_TOL)
            or ortho_mo > SORT_ORTHO_TOL
            or resid > SORT_FOCK_TOL * max(scale, 1.0)
            or worst_de > SORT_DE_TOL):
        raise RuntimeError("sort-by-energy validation FAILED: " + report)
    return energies, report


def plan_orbital_order(molden_path: Path, pure: Path | None = None,
                       s_matrix: Path | None = None,
                       fock_ao: Path | None = None) -> SortPlan:
    """Decide the output order for a JANPA orbital export.

    Energy route when the ``-doFock`` dumps are available, else the CLPO
    class order from the janpa stdout; raises with instructions when
    neither is possible.  Output order: occupied (Occup > 1.0) first, then
    ascending E -- or, in the class route, [BD+LP] [NB] [RY] in file order.
    """
    base = derive_base(molden_path.stem)
    d = molden_path.parent
    se = _set_from_stem(molden_path.stem)
    require_ortho = se is None or se.orthonormal
    tgt = molden_mo_vectors(molden_path)
    nblk = len(tgt)
    occs = [t[1] for t in tgt]
    log_path = d / f"{base}.JANPA"
    labels = None
    lab_occ = None
    worst_lab = None
    label_why = f"no {log_path.name}"
    if log_path.is_file():
        labels, lab_occ = clpo_classes_from_log(log_path)
        if not labels:
            label_why = f"{log_path.name} has no CLPO summary block"
        elif len(labels) != nblk:
            hint = ""
            if 0 < len(labels) < nblk:
                hint = (" -- partial label list: rerun janpa with "
                        "-RyOccPrintThreshold -1 (or use the --clpo flag)")
            label_why = (f"{log_path.name}: {len(labels)} labels vs "
                         f"{nblk} MO blocks{hint}")
        else:
            worst_lab = max(abs(occs[i] - lab_occ[i + 1])
                            for i in range(nblk))
            if worst_lab <= LABEL_OCC_TOL:
                label_why = ""
            else:
                label_why = (f"{log_path.name}: occupancy mismatch vs "
                             f"{molden_path.name} (max {worst_lab:.1e}) -- "
                             "the log does not describe this file")
    label_ok = label_why == ""
    if not label_ok:
        labels = None  # never attach labels that belong to another set

    pure = pure or d / f"{base}.PURE"
    s_matrix = s_matrix or d / f"{base}.S.txt"
    fock_ao = fock_ao or d / f"{base}.fock_ao.txt"
    missing = [p.name for p in (pure, s_matrix, fock_ao) if not p.is_file()]

    if not missing:
        energies, gates = fock_orbital_energies(molden_path, pure, s_matrix,
                                                fock_ao, require_ortho)
        order = sorted(range(nblk),
                       key=lambda i: (occs[i] <= 1.0, energies[i]))
        enes = [energies[i] for i in order]
        mode = "energy"
        lines = [f"sort  : energy (<phi|F|phi>; {s_matrix.name}, "
                 f"{fock_ao.name}, {pure.name})",
                 f"        {gates}"]
        if not require_ortho:
            lines.append(
                f"        note: pre-orthogonalization set -- orbitals are "
                f"normalized but mutually non-orthogonal, and Occup sums "
                f"to {sum(occs):.3f} e, not the electron count; energies "
                f"are single-orbital expectation values, pair analysis is "
                f"not defined for it")
    elif label_ok:
        order = sorted(range(nblk),
                       key=lambda i: _CLASS_RANK.get(labels[i + 1], 3))
        enes = [float(i) for i in range(nblk)]  # sequential, as JANPA writes
        mode = "class"
        lines = [f"sort  : class (energy data missing: {', '.join(missing)})",
                 f"        labels from {log_path.name}; stdout->file mapping "
                 f"verified, max |dOccup| {worst_lab:.1e}"]
    else:
        raise RuntimeError(
            f"cannot order orbitals: no energy inputs ({', '.join(missing)})"
            f" and no usable CLPO labels ({label_why}).\n"
            f"  produce them with:  python orca_to_janpa.py {base} --clpo "
            f"(or the set flag matching this export: --lho/--aho/--lpo/"
            f"--nao/--pnao)"
        )

    n_occ = sum(1 for o in occs if o > 1.0)
    if labels is None:
        lines.append(f"        labels: unused ({label_why})")
    else:
        counts = {}
        for i in range(nblk):
            lab = labels[i + 1]
            counts[lab] = counts.get(lab, 0) + 1
        lines.append("        classes: " + ", ".join(
            f"{c} {counts[c]}" for c in ("BD", "LP", "NB", "RY")
            if c in counts))
        odd = [i + 1 for i in range(nblk)
               if (occs[i] > 1.0) != (_CLASS_RANK.get(labels[i + 1], 3) == 0)]
        if odd:
            lines.append(
                f"        note: Occup>1.0 disagrees with BD/LP class for "
                f"{len(odd)} orbital(s): {odd[:6]}")
    rule = "ascending E" if mode == "energy" else "[BD+LP] [NB] [RY]"
    lines.append(f"        order: occupied (Occup>1.0: {n_occ}) first, "
                 f"then {rule}")
    if mode == "energy":
        e_occ = enes[:n_occ]
        if e_occ:
            lines.append(f"        occupied block: {min(e_occ):+.6f} .. "
                         f"{max(e_occ):+.6f} Ha "
                         f"({min(e_occ) * HARTREE_EV:+.2f} .. "
                         f"{max(e_occ) * HARTREE_EV:+.2f} eV)")
            for pos in range(min(n_occ, 12)):
                lab = f" {labels[order[pos] + 1]}" if labels else ""
                lines.append(
                    f"          {pos + 1:>2}  E={enes[pos]:+.6f} Ha "
                    f"({enes[pos] * HARTREE_EV:+8.3f} eV)  "
                    f"occ={occs[order[pos]]:.5f}{lab}")
            if n_occ > 12:
                lines.append(f"          ... ({n_occ - 12} more)")
        lines.append(f"        E range: {min(enes):+.6f} .. "
                     f"{max(enes):+.6f} Ha")
    return SortPlan(order, enes, "\n".join(lines))


# ---------------------------------------------------------------------------
# Pairwise donor -> acceptor interaction analysis (E(2)-like + charge
# transfer).
#
# JANPA does NOT print NBO-style second-order interaction energies: its wiki
# ("E2_pert") explains that the ingredients are the Fock matrix in a
# localized basis, and warns that E(2) only has a well-defined meaning for
# true Hartree-Fock wavefunctions -- under DFT the Fock matrix belongs to
# the auxiliary Kohn-Sham system, so the numbers are indicative at best.
# What JANPA does print is the charge transferred between localized orbitals
# (its experimental "charge transfer analysis", fixed 0.01 e print
# threshold).
#
# Both quantities are computable from the dumps this script already gates:
#
#   E2(i -> j) = n_i F_ij^2 / (F_jj - F_ii)     [kcal/mol, donor i occupied]
#   q(i -> j)  = D_ij^2 / D_ii                 [electrons]
#
# with F the Fock matrix in the localized basis (F_ab = <phi_a|F|phi_b>) and
# D the 1-RDM in the same basis (D_ab = <phi_a|D|phi_b>; built from the
# .PURE MOs, whose Occup= values carry the electrons).  q is exactly the
# number JANPA's own CT table prints -- verified against every printed pair
# on the test molecules -- so this mode computes it for ALL pairs, not just
# those above JANPA's fixed print threshold.
# ---------------------------------------------------------------------------

E2_KCAL = 627.5094740631     # Hartree -> kcal/mol
E2_ROUTE_TOL = 1e-3          # route A vs B Fock agreement (relative)
CT_REPRO_TOL = 2e-5          # |q_ours - q_janpa| vs the printed values
E2_STRONG_RATIO = 0.25       # |F_ij|/dE above which second order fails


def pair_interaction_analysis(
        target: Path,
        pure: Path | None = None,
        s_matrix: Path | None = None,
        fock_ao: Path | None = None,
        set_key: str | None = None,
) -> tuple[list[str], list[str], str]:
    """E(2)-like + charge-transfer pair table for a JANPA orbital export.

    ``set_key`` selects the transformation chain for the route-B cross-check
    (CLPO only); without it the CLPO tag in the filename decides.

    Returns (header_lines, table_rows, totals_line).
    """
    base = derive_base(target.stem)
    d = target.parent
    se = _SET_EXPORTS.get(set_key) if set_key else _set_from_stem(target.stem)
    if se is not None and not se.pair:
        raise RuntimeError(
            f"{se.tag} is a pre-orthogonalization intermediate set: its "
            "orbitals are normalized but not mutually orthogonal, and its "
            "occupancies do not sum to the electron count, so the "
            "pair-interaction analysis (E2, q) is not defined for it (the "
            "viewer file and its energy sort are unaffected)."
        )
    pure = pure or d / f"{base}.PURE"
    s_matrix = s_matrix or d / f"{base}.S.txt"
    fock_ao = fock_ao or d / f"{base}.fock_ao.txt"
    missing = [p.name for p in (pure, s_matrix, fock_ao) if not p.is_file()]
    if missing:
        raise RuntimeError(
            f"pair analysis needs {', '.join(missing)} next to "
            f"{target.name}\n"
            f"  produce them with:  python orca_to_janpa.py {base} --clpo"
        )
    S, F, eps, Cm, Ct, occ_mo, occ_tgt = _load_pair_inputs(
        target, pure, s_matrix, fock_ao)
    n = len(S)
    SCt, ortho_mo, ortho_tgt, resid, scale = _fock_gates(S, F, eps, Cm, Ct)
    F_loc = _matmul(list(zip(*Ct)), _matmul(F, Ct))
    U = _matmul(list(zip(*Cm)), SCt)      # [canonical k][target i]

    # 1-RDM in the target basis: D = sum_k occ_k |u_k><u_k|
    Dm = [[0.0] * n for _ in range(n)]
    for k in range(n):
        ok = occ_mo[k]
        if ok == 0.0:
            continue
        uk = U[k]
        for i in range(n):
            uki = uk[i]
            if uki == 0.0:
                continue
            row = Dm[i]
            f = ok * uki
            for j in range(n):
                row[j] += f * uk[j]

    checks = (f"gates  : {target.name} orthonormality {ortho_tgt:.2e} | "
              f"{pure.name} orthonormality {ortho_mo:.2e} | "
              f"F C - S C eps {resid:.2e} (scale {scale:.2f})")
    if (ortho_tgt > SORT_ORTHO_TOL or ortho_mo > SORT_ORTHO_TOL
            or resid > SORT_FOCK_TOL * max(scale, 1.0)):
        raise RuntimeError("pair analysis validation FAILED: " + checks)

    # optional independent route B from JANPA's own NAO data (wiki E2_pert
    # recipe).  The path is set-specific: a NAO-space chain product
    # (CLPO/LHO/AHO/LPO) is checked for orthonormality on the way; for NAO
    # the F_NAO dump itself is the independent side; PNAO has no NAO-space
    # path and skips the check with a note.
    nao = d / f"{base}.fock_nao.txt"
    if se is None or se.chain is None:
        checks += (" | route B skipped (no NAO-space transformation chain "
                   "for this set)")
    elif not nao.is_file():
        checks += " | route B skipped (F_NAO dump missing)"
    else:
        chain_paths = [d / _DUMP_EXPORTS[c][1].format(b=base)
                       for c in se.chain]
        if not all(p.is_file() for p in chain_paths):
            checks += (" | route B skipped (chain dumps missing: "
                       + ", ".join(p.name for p in chain_paths
                                   if not p.is_file()) + ")")
        else:
            F_NAO = read_matrix_dump(nao)
            dev_b = None
            if se.chain:
                T = read_matrix_dump(chain_paths[0])
                for p in chain_paths[1:]:
                    T = _matmul(T, read_matrix_dump(p))
                F_b = _matmul(T, _matmul(F_NAO, list(zip(*T))))
                dev_b = _max_dev_from_identity(_matmul(T, list(zip(*T))))
                what = "chain " + "+".join(se.chain)
            else:
                F_b = F_NAO          # NAO: the dump is the other path
                what = "direct F_NAO dump"
            diff_b = max(abs(F_loc[i][j] - F_b[i][j])
                         for i in range(n) for j in range(n))
            scale_b = max(abs(F_loc[i][j]) for i in range(n) for j in range(n))
            checks += (f" | route A vs B {diff_b:.2e} ({what}, "
                       f"max|F| {scale_b:.1f}"
                       + (f", basis {dev_b:.1e}" if dev_b is not None else "")
                       + ")")
            if ((dev_b is not None and dev_b > SORT_ORTHO_TOL)
                    or diff_b > E2_ROUTE_TOL * max(scale_b, 1.0)):
                raise RuntimeError("pair analysis validation FAILED: " + checks)

    # labels and JANPA's own CT numbers -- only when the log describes
    # this exact export (occupancy fingerprint)
    labels = None
    ct_note = ""
    log_path = d / f"{base}.JANPA"
    if not log_path.is_file():
        ct_note = (f"no {log_path.name} -- no orbital labels and no CT "
                   "cross-check")
    else:
        summary = clpo_summary_from_log(log_path)
        if not summary or len(summary) != len(occ_tgt):
            ct_note = (f"{log_path.name} does not describe this export "
                       f"({len(summary)} labels vs {len(occ_tgt)} MOs)")
        else:
            worst_lab = max(abs(occ_tgt[i] - summary[i + 1][2])
                            for i in range(len(occ_tgt)))
            if worst_lab > LABEL_OCC_TOL:
                ct_note = (f"{log_path.name} occupancy mismatch "
                           f"({worst_lab:.1e}) -- labels/CT check skipped")
            else:
                labels = {i: v[1] for i, v in summary.items()}
                pairs = clpo_ct_pairs_from_log(log_path)
                if pairs:
                    worst_ct = max(abs(Dm[i - 1][j - 1] ** 2
                                       / Dm[i - 1][i - 1] - q)
                                   for i, q, j in pairs)
                    ct_note = (f"JANPA CT reproduced {len(pairs)}/"
                               f"{len(pairs)} printed pairs "
                               f"(max |dq| {worst_ct:.1e})")
                    if worst_ct > CT_REPRO_TOL:
                        raise RuntimeError(
                            "pair analysis validation FAILED: " + checks
                            + " | " + ct_note)
                else:
                    ct_note = (f"{log_path.name}: no CT pairs above "
                               "JANPA's 0.01 e print threshold")

    rows = []
    n_skip = 0
    for i in range(n):
        if occ_tgt[i] <= 1.0:
            continue
        for j in range(n):
            if occ_tgt[j] > 1.0:
                continue
            de = F_loc[j][j] - F_loc[i][i]
            dii = Dm[i][i]
            q = (Dm[i][j] ** 2 / dii) if dii > 1e-12 else 0.0
            e2 = None
            strong = False
            if de > 1e-9:
                e2 = occ_tgt[i] * F_loc[i][j] ** 2 / de * E2_KCAL
                strong = abs(F_loc[i][j]) / de >= E2_STRONG_RATIO
            else:
                n_skip += 1
            rows.append((e2, q, i, j, F_loc[i][j], de, strong))
    rows.sort(key=lambda r: (r[0] is None, -(r[0] or 0.0)))

    fallback = se.tag if se else "MO"

    def name(idx: int) -> str:
        return labels[idx + 1] if labels else f"{fallback} {idx + 1}"

    table_rows = []
    for k, (e2, q, i, j, fij, de, strong) in enumerate(rows, 1):
        e2s = f"{e2:11.2f}" if e2 is not None else "          -"
        table_rows.append(
            f"{k:>4}  {name(i):>18} -> {name(j):<18} "
            f"{fij:+.6f} {de:9.5f} {e2s} {q:9.5f}"
            + (" *" if strong else ""))

    sum_e2 = sum(r[0] for r in rows if r[0] is not None)
    sum_q = sum(r[1] for r in rows)
    n_strong = sum(1 for r in rows if r[6])
    n_don = sum(1 for o in occ_tgt if o > 1.0)
    n_acc = len(occ_tgt) - n_don
    totals = (f"totals : sum E2 = {sum_e2:.1f} kcal/mol | "
              f"sum q = {sum_q:.5f} e | {n_don} donors x {n_acc} acceptors "
              f"= {len(rows)} pairs"
              + (f" | {n_strong} strongly mixed (*)" if n_strong else "")
              + (f" | {n_skip} skipped: dE <= 0" if n_skip else ""))

    header = [
        f"pair-interaction analysis  [{target.name}]",
        f"inputs : {s_matrix.name}, {fock_ao.name}, {pure.name}",
        f"check  : {ct_note}",
        checks,
        "note   : E2 = n_i F_ij^2/(F_jj-F_ii) is a perturbation estimate "
        "with a well-defined",
        "         meaning in Hartree-Fock; under DFT the 'Fock' operator "
        "belongs to the",
        "         auxiliary Kohn-Sham system, so read E2 with care there "
        "and prefer the",
        "         charge q = D_ij^2/D_ii (the quantity JANPA's CT analysis "
        "prints; see its",
        "         wiki page 'E2_pert' and Nikolaienko et al., J. Comput. "
        "Chem. 39 (2018) 1090).",
        f"         rows marked * have |F_ij|/(F_jj-F_ii) >= "
        f"{E2_STRONG_RATIO}: the two orbitals are strongly mixed, the "
        "second-order",
        "         estimate is not meaningful for them (often a sign the "
        "Lewis-like",
        "         reference itself is inadequate -- e.g. 3c-2e bonding).",
        f"{'#':>4}  {'donor':>18} -> {'acceptor':<18} {'F_ij':>9} "
        f"{'dE/Ha':>9} {'E2(kcal/mol)':>11} {'q(e)':>9}",
    ]
    return header, table_rows, totals


def _emit_e2(target: Path, pure: Path | None, s_matrix: Path | None,
             fock_ao: Path | None, out: Path | None = None,
             set_key: str | None = None) -> Path:
    """Print the pair-interaction report and write the table file."""
    header, table, totals = pair_interaction_analysis(target, pure, s_matrix,
                                                      fock_ao,
                                                      set_key=set_key)
    for line in header:
        print(line)
    for line in table[:25]:
        print(line)
    if len(table) > 25:
        print(f"      ... ({len(table) - 25} more rows in the file)")
    print(totals)
    out = out or target.with_name(target.stem + "_E2.txt")
    out.write_text("\n".join(header + table + [totals]) + "\n",
                   encoding="utf-8", newline="\n")
    print(f"Wrote {out}")
    return out


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="orca_to_janpa",
        description=(
            "ORCA .gbw -> Molden -> JANPA (.PURE) pipeline for HF/DFT. "
            "No .47 file needed."
        ),
    )
    ap.add_argument(
        "target",
        nargs="?",
        help="Basename or path to .gbw/.mp2nos/.inp (e.g. ethene_NBO)",
    )
    ap.add_argument("--orca-dir", type=Path, default=DEFAULT_ORCA_DIR)
    ap.add_argument("--janpa-dir", type=Path, default=DEFAULT_JANPA_DIR)
    ap.add_argument(
        "--dot47",
        type=Path,
        default=None,
        help="Optional .47/.mdcip.47 for the correlated-density (-ds47) route",
    )
    ap.add_argument(
        "--janpa-out",
        type=Path,
        default=None,
        help="Save JANPA stdout here (default: <base>.JANPA next to input)",
    )
    ap.add_argument(
        "--janpa-args",
        nargs=argparse.REMAINDER,
        default=[],
        help="Extra args passed through to janpa.jar after --",
    )
    _set_blurb = {
        "clpo": "the NBO-analog Lewis-like set (BD/NB/LP/RY); the only set "
                "with CLPO labels in <base>.JANPA and a verified route-B "
                "cross-check for --e2",
        "lho": "localized hybrids (the atom-centered hybrids CLPOs are "
               "built from)",
        "aho": "atomic hybrids (the LPO-construction intermediate)",
        "lpo": "property-optimized localized orbitals",
        "nao": "natural atomic orbitals (the NPA/Wiberg set)",
        "pnao": "pre-orthogonalization NAOs",
    }
    for _key, _se in _SET_EXPORTS.items():
        ap.add_argument(
            f"--{_key}",
            action="store_true",
            help=f"Pipeline mode: export the {_se.tag} set ({_set_blurb[_key]})."
            f"  Writes <base>_{_se.tag}_spherical.molden (the substrate: "
            f"JANPA's export with corrected markers/spin, kept as the "
            f"analysis input) and <base>_{_se.tag}.molden (the viewer file: "
            f"cartesian, real energies, occupied-first; integer Occup with "
            f"--avogadro), plus the shared dumps (<base>.S.txt, "
            f"<base>.fock_ao.txt, <base>.fock_nao.txt) and its "
            f"transformation chain",
        )
    ap.add_argument(
        "--all-sets",
        action="store_true",
        help="Pipeline mode for every JANPA orbital set at once: --clpo "
        "--lho --aho --lpo --nao --pnao",
    )
    ap.add_argument(
        "--diagnose",
        type=Path,
        default=None,
        metavar="ORCA.OUT",
        help="Only diagnose why an ORCA .out has no .47; do not run pipeline",
    )
    ap.add_argument(
        "--to-cart",
        type=Path,
        default=None,
        metavar="SPHER.MOLDEN",
        help="Spherical->cartesian Molden conversion (d 5->6, f 7->10, "
        "markers dropped) for cartesian-only viewers e.g. Avogadro",
    )
    ap.add_argument(
        "--cart-out",
        type=Path,
        default=None,
        help="Output for --to-cart (default: <stem>_cart.molden)",
    )
    ap.add_argument(
        "--fix-markers",
        type=Path,
        default=None,
        metavar="MOLDEN",
        help="Rewrite [5D]/[7F]/[9G] markers from the actual shells and "
        "default-fix the spin label (this is also the substrate sanitizer "
        "used by the set flags)",
    )
    ap.add_argument(
        "--markers-out",
        type=Path,
        default=None,
        help="Output for --fix-markers (default: <stem>_fixed.molden)",
    )
    ap.add_argument(
        "--spin",
        choices=("Alpha", "Beta"),
        default=None,
        help="Retag Spin= lines (overrides the automatic correction: a "
        "uniform closed-shell set gets the canonical source's own label)",
    )
    ap.add_argument(
        "--sort-energy",
        action="store_true",
        help="Standalone --to-cart/--fix-markers modifier (the set flags "
        "apply this by default): write the MOs occupied-first then in "
        "ascending Fock energy E=<phi|F|phi> (reads the <base>.S.txt / "
        "<base>.fock_ao.txt dumps and <base>.PURE from a set-flag run); "
        "falls back to CLPO class order (BD+LP, NB, RY) from <base>.JANPA "
        "when those inputs are missing",
    )
    ap.add_argument(
        "--pure",
        type=Path,
        default=None,
        help="--sort-energy/--e2 input override: canonical MOs + energies "
        "(default: <base>.PURE next to the molden file)",
    )
    ap.add_argument(
        "--s-matrix",
        type=Path,
        default=None,
        help="--sort-energy/--e2 input override: JANPA overlap dump "
        "(default: <base>.S.txt)",
    )
    ap.add_argument(
        "--fock-ao",
        type=Path,
        default=None,
        help="--sort-energy/--e2 input override: JANPA Fock dump in the AO "
        "basis (default: <base>.fock_ao.txt)",
    )
    ap.add_argument(
        "--e2",
        nargs="?",
        const="",
        default=None,
        metavar="MOLDEN",
        help="Pairwise donor->acceptor interaction table for a JANPA "
        "export: E2 = n_i F_ij^2/(F_jj-F_ii) in kcal/mol (Fock matrix in "
        "the localized basis) plus the charge transfer q = D_ij^2/D_ii "
        "that JANPA's own CT analysis prints.  Pass the file (usually "
        "<base>_CLPO_spherical.molden), or use bare --e2 with a set flag "
        "to analyze the export from this run.  Cross-checks JANPA's "
        "printed CT values when the log describes the export",
    )
    ap.add_argument(
        "--e2-out",
        type=Path,
        default=None,
        help="Output for --e2 (default: <stem>_E2.txt)",
    )
    ap.add_argument(
        "--avogadro",
        action="store_true",
        help="--to-cart/set-flag modifier: write integer Occup= 2/0 "
        "(threshold occ>1.0) instead of the true fractional occupancies. "
        "Avogadro parses Occup as int, so fractional occupations miscount "
        "electrons and mislabel virtuals; this is the workaround until the "
        "upstream bug is fixed.  The other fix it used to carry (Alpha "
        "spin) is now the default.",
    )
    return ap


def main(argv: list[str] | None = None) -> None:
    ap = build_parser()
    args = ap.parse_args(argv)

    sets = requested_sets(args)
    e2_inline = args.e2 == ""
    if e2_inline and not sets:
        ap.error("--e2 without a file requires a set flag (--clpo/--lho/"
                 "--aho/--lpo/--nao/--pnao/--all-sets); it analyzes the "
                 "export from this run")
    if args.e2 not in (None, "") and sets:
        ap.error("--e2 takes no file with a set flag (use bare --e2)")
    if args.avogadro and args.to_cart is None and not sets:
        ap.error("--avogadro applies to --to-cart and the set flags")
    if (args.sort_energy and args.to_cart is None
            and args.fix_markers is None and not sets):
        ap.error("--sort-energy applies to --to-cart/--fix-markers and the "
                 "set flags")
    if ((args.pure or args.s_matrix or args.fock_ao)
            and not (args.sort_energy or args.e2 is not None)):
        ap.error("--pure/--s-matrix/--fock-ao are --sort-energy/--e2 inputs")

    if args.diagnose is not None:
        print(diagnose_nbo_output(args.diagnose))
        return
    if args.e2:
        _emit_e2(Path(args.e2), args.pure, args.s_matrix, args.fock_ao,
                 args.e2_out)
        return
    if args.to_cart is not None:
        out = args.cart_out or args.to_cart.with_name(
            args.to_cart.stem + "_cart.molden")
        spin = (args.spin or resolve_spin_label(args.to_cart)
                or ("Alpha" if args.avogadro else None))
        sort = None
        if args.sort_energy:
            sort = plan_orbital_order(args.to_cart, args.pure,
                                      args.s_matrix, args.fock_ao)
            print(sort.report)
        print(convert_to_cart(args.to_cart, out, spin, args.avogadro, sort))
        print(f"Wrote {out}")
        return
    if args.fix_markers is not None:
        out = args.markers_out or args.fix_markers.with_name(
            args.fix_markers.stem + "_fixed.molden")
        spin = args.spin or resolve_spin_label(args.fix_markers)
        sort = None
        if args.sort_energy:
            sort = plan_orbital_order(args.fix_markers, args.pure,
                                      args.s_matrix, args.fock_ao)
            print(sort.report)
        print(fix_spherical_markers(args.fix_markers, out, spin, sort))
        print(f"Wrote {out}")
        return
    if not args.target:
        ap.error("target is required (or use --diagnose ORCA.OUT)")

    workdir, base = resolve_base(args.target)
    orca_2mkl = find_orca_2mkl(args.orca_dir)
    m2m_jar, janpa_jar = find_jars(args.janpa_dir)

    print(f"[1/3] orca_2mkl  : {base}  (cwd={workdir})")
    molden = run_orca_2mkl(workdir, base, orca_2mkl)
    print(f"      wrote {molden.name}")

    pure = workdir / f"{base}.PURE"
    tag = f" + -ds47 {args.dot47.name}" if args.dot47 else " (HF/DFT, no .47)"
    print(f"[2/3] molden2molden{tag} -> {pure.name}")
    run_molden2molden(workdir, molden, pure, m2m_jar, args.dot47)
    print("      wrote", pure.name)

    out_file = args.janpa_out or (workdir / f"{base}.JANPA")
    print(f"[3/3] janpa -> {out_file.name}")
    passthrough = (
        args.janpa_args[1:]
        if args.janpa_args[:1] == ["--"]
        else args.janpa_args
    )
    janpa_extra = passthrough
    if sets:
        # One Molden export per requested set plus the dumps the analysis
        # modes need: -doFock builds the Fock matrix from the .PURE orbital
        # energies; the Fock_NAO dump and the per-set transformation chains
        # feed the route-B cross-check of --e2 (NAO-space chains for
        # CLPO/LHO/AHO/LPO, the F_NAO dump directly for NAO; JANPA wiki
        # "E2_pert" recipe).  janpa writes every substrate; the viewer
        # files are converted from them.
        chain_keys: list[str] = []
        for k in sets:
            for c in (_SET_EXPORTS[k].chain or ()):
                if c not in chain_keys:
                    chain_keys.append(c)
        set_args: list[str] = []
        for k in sets:
            se = _SET_EXPORTS[k]
            set_args += [se.opt, f"{base}_{se.tag}_spherical.molden"]
        set_args += ["-doFock"]
        for dk in ("fock_ao", "fock_nao", "s", *chain_keys):
            opt, pat = _DUMP_EXPORTS[dk]
            set_args += [opt, pat.format(b=base)]
        set_args += ["-MatrixFloatNumberFormat", "%.9f",
                     "-RyOccPrintThreshold", "-1"]
        janpa_extra = set_args + passthrough  # user args last: they win
    stdout = run_janpa(workdir, pure, janpa_jar, out_file, janpa_extra)
    # Print the electron-count + NPA summary lines as a quick receipt.
    for line in stdout.splitlines():
        if "Total number of electrons" in line or "Sum of electrons" in line:
            print("      " + line.strip())
    print(f"Done. JANPA output saved to {out_file}")
    if sets:
        for k in sets:
            se = _SET_EXPORTS[k]
            substrate = workdir / f"{base}_{se.tag}_spherical.molden"
            viewer = workdir / f"{base}_{se.tag}.molden"
            # Substrate: minimal sanitation only -- markers corrected to
            # match the actual shells and the spin label corrected from the
            # canonical source; JANPA's order, Ene placeholders and
            # fractional Occup are kept.  Analysis input (same spherical
            # basis as the dumps).
            spin = args.spin or resolve_spin_label(substrate,
                                                   workdir / f"{base}.PURE")
            print("      " + fix_spherical_markers(substrate, substrate, spin))
            # Viewer file: cartesian, real Fock energies, occupied-first;
            # --avogadro adds the integer-Occup workaround for Avogadro's
            # electron-count bug (kept behind the flag until upstream fixes
            # it).  PNAO sorts like every other set: its orbitals are
            # normalized, so E = <phi|F|phi> is well-defined -- the sort
            # report states the non-orthogonality and the relaxed gate.
            sort = plan_orbital_order(substrate, args.pure, args.s_matrix,
                                      args.fock_ao)
            print(sort.report)
            print("      " + convert_to_cart(
                substrate, viewer,
                spin or ("Alpha" if args.avogadro else None),
                args.avogadro, sort))
            print(f"Viewer file  : {viewer.name}  (cartesian, real energies, "
                  "occupied-first"
                  + (", integer Occup" if args.avogadro else "") + ")")
            print(f"Substrate    : {substrate.name}  (spherical -- analysis "
                  "input)")
            if e2_inline:
                if se.pair:
                    _emit_e2(substrate, args.pure, args.s_matrix,
                             args.fock_ao,
                             workdir / f"{base}_{se.tag}_E2.txt", set_key=k)
                else:
                    print(f"      {se.tag}_E2: skipped (pair analysis is "
                          "not defined for this set)")
        print(f"Data         : {base}.S.txt, {base}.fock_ao.txt, "
              f"{base}.fock_nao.txt"
              + "".join(f", {base}.{c}.txt" for c in chain_keys))
        _first = next((_SET_EXPORTS[k] for k in sets
                       if _SET_EXPORTS[k].pair), None)
        if _first is not None:
            print(f"Next         : python {Path(sys.argv[0]).name} --e2 "
                  f"{base}_{_first.tag}_spherical.molden")


if __name__ == "__main__":
    main()
