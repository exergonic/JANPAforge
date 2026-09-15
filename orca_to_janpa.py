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

Viewer workflow: ``--to-cart`` / ``--fix-markers`` rewrite a JANPA export
for Avogadro / MOrbVis, and ``--sort-energy`` orders the MOs by their Fock
expectation value -- or, when that data is absent, by CLPO class (bonding,
antibonding, Rydberg).  The pipeline ``--clpo`` flag exports the CLPOs and
the energy data that ``--sort-energy`` needs.  See README.md.
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

# Molden component orders (match MOrbVis's evaluator + Molden manual).
CART_D_ORDER = ("xx", "yy", "zz", "xy", "xz", "yz")
CART_F_ORDER = ("xxx", "yyy", "zzz", "xyy", "xxy", "xxz",
                "xzz", "yzz", "yyz", "xyz")

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
    """Projection error max|T.C - c| and overlap-norm drift |C'SC - c'c|."""
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
    fills orbitals positionally by energy order, so fractional
    occupations with all-zero energies mislabel virtuals as occupied.

    With ``sort`` (a SortPlan from ``plan_orbital_order``) the MO blocks
    are written in the plan's order and each ``Ene=`` is replaced by the
    plan's value -- real Fock expectation energies, or the class-order
    sequential numbers when no energy data exists.
    """
    shell_ls = molden_shell_ls(in_path)
    if "g" in shell_ls:
        raise NotImplementedError("g shells not implemented for --to-cart")
    lines = in_path.read_text(encoding="utf-8", errors="replace").splitlines()
    out: list[str] = []
    worst_dc = 0.0
    worst_dn = 0.0
    mo_blocks: list[tuple[list[str], list[float], float]] = []
    cur_header: list[str] = []
    cur_coeffs: list[float] = []
    in_mo = False

    def flush_mo() -> None:
        nonlocal worst_dc, worst_dn
        if not cur_header and not cur_coeffs:
            return
        expect = sum(SPHER_N[L] for L in shell_ls)
        if len(cur_coeffs) != expect:
            raise ValueError(
                f"MO #{len(mo_blocks) + 1}: {len(cur_coeffs)} coeffs but "
                f"basis has {expect} spherical functions -- wrong file?")
        cart: list[float] = []
        pos = 0
        for L in shell_ls:
            n = SPHER_N[L]
            c = cur_coeffs[pos:pos + n]
            pos += n
            if L in ("s", "p"):
                cart.extend(c)
            elif L == "d":
                cc = spher_d_to_cart(c)
                cart.extend(cc)
                dc, dn = _check_block("d", c, cc)
                worst_dc = max(worst_dc, dc)
                worst_dn = max(worst_dn, dn)
            elif L == "f":
                cc = spher_f_to_cart(c)
                cart.extend(cc)
                dc, dn = _check_block("f", c, cc)
                worst_dc = max(worst_dc, dc)
                worst_dn = max(worst_dn, dn)
        occup = 0.0
        for h in cur_header:
            if h.strip().startswith("Occup="):
                try:
                    occup = float(h.split("=", 1)[1])
                except ValueError:
                    pass
        mo_blocks.append((list(cur_header), cart, occup))
        cur_header.clear()
        cur_coeffs.clear()

    for raw in lines:
        line = raw.strip()
        if line.startswith("[MO]"):
            in_mo = True
            out.append(raw)
            continue
        if in_mo and line.startswith("["):
            flush_mo()
            in_mo = False
            out.append(raw)
            continue
        if _is_marker(line):
            continue  # drop spherical markers in cartesian output
        if in_mo:
            if line.startswith("Sym="):
                flush_mo()
                cur_header.append(raw)
            elif line.startswith(("Ene=", "Occup=")):
                cur_header.append(raw)
            elif line.startswith("Spin="):
                cur_header.append(
                    f"Spin= {spin}" if spin else raw)
            elif line == "":
                continue
            else:
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        cur_coeffs.append(float(parts[1]))
                        continue
                    except ValueError:
                        pass
                cur_header.append(raw)
        else:
            out.append(raw)
    flush_mo()

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
    for pos, (header, coeffs, occup) in enumerate(ordered):
        is_occ = avogadro and occup > 1.0
        for h in header:
            hs = h.strip()
            if sort is not None and hs.startswith("Ene="):
                out.append(f" Ene= {sort.enes[pos]:.14E}")
            elif hs.startswith("Occup=") and avogadro:
                out.append(f" Occup= {2 if is_occ else 0}")
            else:
                out.append(h)
        for i, v in enumerate(coeffs, 1):
            out.append(f"  {i:<6d}{v: .12f}".rstrip())

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
#     against sum_k |<CLPO_i|MO_k>|^2 eps_k -- all before anything is
#     written.  The pipeline ``--clpo`` flag produces every input.
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
                    "cart", "fixed")
_CLASS_RANK = {"BD": 0, "LP": 0, "NB": 1, "RY": 2}


class SortPlan(NamedTuple):
    """Output order for one molden orbital set (see plan_orbital_order)."""

    mode: str            # "energy" or "class"
    perm: list[int]      # perm[i] = input MO block index -> output slot i
    enes: list[float]    # Ene= value per output slot (Ha, or sequential)
    report: str          # human-readable audit trail


def derive_base(stem: str) -> str:
    """'ethene_CLPO' -> 'ethene'; strips known export/conversion suffixes."""
    parts = stem.split("_")
    while len(parts) > 1 and parts[-1].lower() in _EXPORT_SUFFIXES:
        parts.pop()
    return "_".join(parts)


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
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    blocks: list[tuple[float, float, list[float]]] = []
    ene = 0.0
    occ = 0.0
    coeffs: list[float] = []
    in_mo = False

    def flush() -> None:
        nonlocal ene, occ
        if coeffs:
            blocks.append((ene, occ, list(coeffs)))
            coeffs.clear()
        ene = occ = 0.0

    for raw in lines:
        line = raw.strip()
        if not in_mo:
            if line.startswith("[MO]"):
                in_mo = True
            continue
        if line.startswith("[") or line.startswith("Sym="):
            flush()
            if line.startswith("["):
                in_mo = False
            continue
        if line.startswith("Ene="):
            ene = float(line.split("=", 1)[1])
        elif line.startswith("Occup="):
            occ = float(line.split("=", 1)[1])
        elif line.startswith("Spin=") or line == "":
            continue
        else:
            parts = line.split()
            if len(parts) >= 2:
                try:
                    coeffs.append(float(parts[1]))
                except ValueError:
                    pass
    flush()
    if not blocks:
        raise RuntimeError(f"{path}: no [MO] blocks found")
    return blocks


def clpo_classes_from_log(log_path: Path):
    """({index: 'BD'|'LP'|'NB'|'RY'}, {index: occupancy}) from janpa stdout.

    Returns (None, None) when the log has no CLPO summary block.
    """
    text = log_path.read_text(encoding="utf-8", errors="replace")
    if "*** Summary of CLPO results" not in text:
        return None, None
    block = text.split("*** Summary of CLPO results", 1)[1]
    block = block.split("Number of two-center", 1)[0]
    labels: dict[int, str] = {}
    occs: dict[int, float] = {}
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
        idx = int(head[0])
        desc = cells[1]
        if "antibonding" in desc:
            lab = "NB"
        elif "(BD)" in desc:
            lab = "BD"
        elif "(LP)" in desc:
            lab = "LP"
        elif "(RY)" in desc:
            lab = "RY"
        else:
            lab = "??"
        labels[idx] = lab
        occs[idx] = occ
    return labels, occs


def _matmul(A: list[list[float]], B: list[list[float]]) -> list[list[float]]:
    return [[sum(map(mul, row, col)) for col in zip(*B)] for row in A]


def _max_dev_from_identity(M: list[list[float]]) -> float:
    worst = 0.0
    for i, row in enumerate(M):
        for j, v in enumerate(row):
            worst = max(worst, abs(v - (1.0 if i == j else 0.0)))
    return worst


def fock_orbital_energies(target: Path, pure: Path, s_path: Path,
                          f_path: Path) -> tuple[list[float], str]:
    """E_i = <phi_i|F|phi_i> for every MO block of ``target``, fully gated.

    S and F come from JANPA's ``-doFock`` dumps (spherical AO basis, [GTO]
    function order).  Gates: orthonormality of both orbital sets under S,
    the canonical residual F C = S C diag(eps) on the .PURE MOs, and an
    independent cross-check of E via the canonical-MO expansion weights.
    Raises instead of returning untrustworthy numbers.
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

    SCm = _matmul(S, Cm)
    SCt = _matmul(S, Ct)
    ortho_mo = _max_dev_from_identity(_matmul(list(zip(*Cm)), SCm))
    ortho_clpo = _max_dev_from_identity(_matmul(list(zip(*Ct)), SCt))
    FCm = _matmul(F, Cm)
    resid = 0.0
    scale = 0.0
    for a in range(n):
        for k in range(n):
            resid = max(resid, abs(FCm[a][k] - SCm[a][k] * eps[k]))
            scale = max(scale, abs(FCm[a][k]))
    F_ct = _matmul(F, Ct)
    M = _matmul(list(zip(*Ct)), F_ct)          # F in the target basis
    energies = [M[i][i] for i in range(n)]

    # Independent cross-check: E_i = sum_k |<phi_i|MO_k>|^2 eps_k.
    U = _matmul(list(zip(*Cm)), SCt)           # [canonical k][target i]
    worst_de = 0.0
    for i in range(n):
        e2 = sum(U[k][i] * U[k][i] * eps[k] for k in range(n))
        worst_de = max(worst_de, abs(e2 - energies[i]))

    report = (f"gates: CLPO orthonormality {ortho_clpo:.2e} | "
              f"MO orthonormality {ortho_mo:.2e} | "
              f"F C - S C eps {resid:.2e} (scale {scale:.2f}) | "
              f"E cross-check {worst_de:.2e} Ha")
    if (ortho_clpo > SORT_ORTHO_TOL or ortho_mo > SORT_ORTHO_TOL
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
                                                fock_ao)
        order = sorted(range(nblk),
                       key=lambda i: (occs[i] <= 1.0, energies[i]))
        enes = [energies[i] for i in order]
        mode = "energy"
        lines = [f"sort  : energy (<phi|F|phi>; {s_matrix.name}, "
                 f"{fock_ao.name}, {pure.name})",
                 f"        {gates}"]
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
            f"  produce them with:  python orca_to_janpa.py {base} --clpo"
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
    return SortPlan(mode, order, enes, "\n".join(lines))


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
    ap.add_argument(
        "--clpo",
        action="store_true",
        help="Pipeline mode: export CLPOs plus the viewer-order energy "
        "artifacts -- <base>_CLPO.molden, <base>.S.txt, <base>.fock_ao.txt "
        "and the full CLPO label list in <base>.JANPA",
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
        help="Rewrite [5D]/[7F]/[9G] markers to match shells present "
        "(JANPA omits [7F], writes spurious [9G])",
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
        help="Retag Spin= lines (JANPA labels closed-shell MOs Beta; "
        "Avogadro prefers Alpha)",
    )
    ap.add_argument(
        "--sort-energy",
        action="store_true",
        help="--to-cart/--fix-markers modifier: write the MOs occupied-"
        "first then in ascending Fock energy E=<phi|F|phi> (reads the "
        "<base>.S.txt / <base>.fock_ao.txt dumps and <base>.PURE, all "
        "written by a '--clpo' janpa run); falls back to CLPO class order "
        "(BD+LP, NB, RY) from <base>.JANPA when those inputs are missing",
    )
    ap.add_argument(
        "--pure",
        type=Path,
        default=None,
        help="--sort-energy input override: canonical MOs + energies "
        "(default: <base>.PURE next to the molden file)",
    )
    ap.add_argument(
        "--s-matrix",
        type=Path,
        default=None,
        help="--sort-energy input override: JANPA overlap dump "
        "(default: <base>.S.txt)",
    )
    ap.add_argument(
        "--fock-ao",
        type=Path,
        default=None,
        help="--sort-energy input override: JANPA Fock dump in the AO "
        "basis (default: <base>.fock_ao.txt)",
    )
    ap.add_argument(
        "--avogadro",
        action="store_true",
        help="Avogadro-view layout for --to-cart: occupied-first MO order "
        "plus integer Occup= 2/0 (threshold occ>1.0). Avogadro parses "
        "Occup as int and fills orbitals positionally, so fractional "
        "occupations miscount electrons and mislabel virtuals. Implies "
        "--spin Alpha unless --spin is given.",
    )
    return ap


def main(argv: list[str] | None = None) -> None:
    ap = build_parser()
    args = ap.parse_args(argv)

    if args.avogadro and args.to_cart is None:
        ap.error("--avogadro is a --to-cart modifier")
    if (args.sort_energy and args.to_cart is None
            and args.fix_markers is None):
        ap.error("--sort-energy is a --to-cart/--fix-markers modifier")
    if ((args.pure or args.s_matrix or args.fock_ao)
            and not args.sort_energy):
        ap.error("--pure/--s-matrix/--fock-ao are --sort-energy inputs")

    if args.diagnose is not None:
        print(diagnose_nbo_output(args.diagnose))
        return
    if args.to_cart is not None:
        out = args.cart_out or args.to_cart.with_name(
            args.to_cart.stem + "_cart.molden")
        spin = args.spin or ("Alpha" if args.avogadro else None)
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
        sort = None
        if args.sort_energy:
            sort = plan_orbital_order(args.fix_markers, args.pure,
                                      args.s_matrix, args.fock_ao)
            print(sort.report)
        print(fix_spherical_markers(args.fix_markers, out, args.spin, sort))
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
    extra = (
        args.janpa_args[1:]
        if args.janpa_args[:1] == ["--"]
        else args.janpa_args
    )
    if args.clpo:
        # CLPO export plus the dumps that make --sort-energy possible
        # (-doFock builds the Fock matrix from the .PURE orbital energies).
        extra = [
            "-CLPO_Molden_File", f"{base}_CLPO.molden",
            "-doFock",
            "-Fock_AO_File", f"{base}.fock_ao.txt",
            "-S_Matrix_File", f"{base}.S.txt",
            "-MatrixFloatNumberFormat", "%.9f",
            "-RyOccPrintThreshold", "-1",
        ] + extra
    stdout = run_janpa(workdir, pure, janpa_jar, out_file, extra)
    # Print the electron-count + NPA summary lines as a quick receipt.
    for line in stdout.splitlines():
        if "Total number of electrons" in line or "Sum of electrons" in line:
            print("      " + line.strip())
    print(f"Done. JANPA output saved to {out_file}")
    if args.clpo:
        print(f"CLPO export  : {base}_CLPO.molden  (+ {base}.S.txt, "
              f"{base}.fock_ao.txt)")
        print(f"Next         : python {Path(sys.argv[0]).name} --to-cart "
              f"{base}_CLPO.molden --avogadro --sort-energy")


if __name__ == "__main__":
    main()
