# HF/cc-pVDZ cross-check record

Per-molecule evidence for the G09/NBO cross-checks published in
[../VALIDATION.md](../VALIDATION.md).

- `g09_<mol>_hf_ccpvdz.gjf` / `.out` — the Gaussian run (NPA charges,
  Wiberg bond indices, NBO E(2) table)
- `orca_<mol>_hf_ccpvdz.inp` / `.out` — the ORCA run
- `janpa_<mol>_hf_ccpvdz.JANPA` — janpa stdout from the pipeline run
- `clpo_<mol>_hf_ccpvdz_E2.txt` — the pair-interaction table (`--e2`)
- the JANPA work files behind that table, under the workdir names its
  header references: `<base>.PURE`, `<base>.S.txt`,
  `<base>.fock_ao.txt`, `<base>.fock_nao.txt`, `<base>.clpo2lho.txt`,
  `<base>.lho2nao.txt`, `<base>.JANPA` and
  `<base>_CLPO_spherical.molden`

`<base>` is `water_ccpvdz_hf`, `formaldehyde_hf`, or `isobutene_hf`
(its name is quoted in the table's own header lines).

`compare.py` recomputes the tables in `../VALIDATION.md` from these
files (`python compare.py [molecule ...]`; `compare_all.txt` is the
saved full run).

To reproduce a pair-interaction table byte-for-byte, from inside its
folder:

```bash
python ../../orca_to_janpa.py --e2 <base>_CLPO_spherical.molden \
    --e2-out clpo_<mol>_hf_ccpvdz_E2.txt
```

The substrate, the S/Fock/NAO dumps, the CLPO transformation chain and
the janpa log are picked up by their workdir names; nothing else is
required.
