# HEO Mixing Thermodynamics Plotter

An interactive Python tool that computes the two central mixing descriptors of
multicomponent (high-entropy) oxide ceramics —

- **ΔH_mix** — regular-solution mixing enthalpy built from binary pair
  enthalpies, evaluated **per crystallographic sublattice** and recombined with
  site-stoichiometry weights,
- **ΔS_mix^N** — site-resolved configurational entropy normalized per mole of
  atoms in the formula unit —

for three structure families — **Perovskite ABO₃**, **Spinel AB₂O₄**, and
**Rock Salt AO** — and plots all entered compositions on a single
**color-coded ΔH–ΔS classification map** to test whether the structure
families separate in this two-descriptor space.

A semi-empirical binary enthalpy database (~70 elements, kJ/mol) is **embedded
in the script** — no external data files are needed.

---

## 1. Executive Summary

| Item | Description |
|---|---|
| Language / dependencies | Python 3.8+, NumPy, Pandas, Matplotlib |
| Interface | Interactive console (structure menu → site compositions → plot) |
| Structures | Perovskite ABO₃, Spinel AB₂O₄, Rock Salt AO |
| Session mode | Multiple compositions per structure, multiple structures per run, one combined figure |
| Enthalpy database | Embedded symmetric ΔH_ij matrix (kJ/mol), ~70 elements |
| Outputs | Per-composition ΔH_mix (kJ/mol) + ΔS_mix^N (J/mol·K) printout, color-coded ΔH–ΔS scatter map with legend and labels |

## 2. Structures and site weights

| Structure | Formula | Sublattice counts (C) | Site weights w = C/ΣC |
|---|---|---|---|
| Perovskite | ABO₃ | A=1, B=1, O=3 | w_A = w_B = 1/5 |
| Spinel | AB₂O₄ | A=1, B=2, O=4 | w_A = 1/7, w_B = 2/7 |
| Rock Salt | AO | A=1, O=1 | w_A = 1/2 |

Oxygen is fixed to pure O (single species → zero configurational contribution).
The same site weights scale both descriptors — this is what makes the ΔH–ΔS
map discriminative between structure families (§7).

## 3. Descriptors

**Configurational entropy (normalized):**

```
ΔS_mix^N = −R · Σ_sites (C_site / ΣC) · Σ_i x_i · ln(x_i)
```

**Mixing enthalpy (regular-solution form):**

```
ΔH_mix = Σ_sites (C_site / ΣC) · Σ_{i<j} 4·h_ij·x_i·x_j
```

where h_ij is the binary mixing enthalpy of the pair (embedded database,
kJ/mol) and x_i are the site-level cation fractions. A site contributes to
either descriptor only when it hosts ≥ 2 distinct species.

## 4. Binary enthalpy database

- Embedded, symmetric ΔH_ij table in kJ/mol (~70 metallic elements, H … Pu),
  in the spirit of the Takeuchi–Inoue semi-empirical compilation.
- Regular-solution interpolation: a 50/50 pair on a site contributes exactly
  h_ij; a pure site contributes zero.
- Unknown/unparseable pairs default to **0 kJ/mol** (silent ideal-like
  contribution) — see Limitations §9.
- Placeholder entries: the Se, Te and Pa rows/columns carry dummy values
  (0.1) — do not use these elements.

## 5. Interactive workflow

```
=== Multiple Structure Mode ===
1. Perovskite ABO3
2. Spinel AB2O4
3. Rock Salt AO
0. Finish and plot
```

Per structure, the program asks for the A-site composition (and B-site for
Perovskite/Spinel) in the form:

```
Na 0.2 Bi 0.2 Ba 0.2 Sr 0.2 Ca 0.2
```

Fractions are auto-normalized to sum = 1. Each composition is printed with its
ΔH_mix and ΔS_mix^N, accumulated across the session, and finally plotted on
the classification map (fixed color per structure family, one legend entry per
family, labels on the points).

## 6. Sanity checks (equimolar anchors)

These values follow analytically from the descriptors and match reference
values in the HEO literature:

| Test case | Expected ΔS_mix^N |
|---|---|
| Rock Salt, 5 equimolar A cations | (R/2)·ln5 = **6.69 J/mol·K** |
| Perovskite, 5 equimolar A cations + single B | (R/5)·ln5 = **2.68 J/mol·K** |
| Spinel, 5+5 equimolar A/B cations | (3R/7)·ln5 = **5.73 J/mol·K** |

Structural checks: a 50/50 pair contributes exactly h_ij kJ/mol to its site
term; single-species sites contribute zero to both descriptors.

## 7. Reading the map

The site weights (C_site/ΣC) rescale both descriptors per structure: a
perovskite keeps ~20% weight on each cation sublattice, a spinel compresses the
two cation-site contributions to 1/7 and 2/7, and a rock salt dilutes its
single cation site to 1/2. The same local chemistry therefore lands in
different regions of the ΔH–ΔS plane depending on the host structure — compare
*relative positions and cluster separation* of the families, not absolute phase
stability.

## 8. Installation & usage

```bash
pip install -r requirements.txt
python heo_mixing_plotter.py
```

## 9. Limitations & scope

1. **Placeholder data** — Se/Te/Pa rows/columns (0.1 kJ/mol) and the actinide
   columns are placeholders; exclude these elements from real analyses.
2. **Silent zeros** — missing or unparseable binary pairs default to
   0 kJ/mol; for publication-grade numbers, spot-check one known composition
   against the source table.
3. **Regular-solution level only** — no elastic/size term, no temperature
   dependence, no DFT/CALPHAD; descriptors are indicators, not phase verdicts.
4. **No charge-neutrality / valence check** — the tool accepts any element on
   any site.
5. **Entropy convention** — normalized per formula unit (C_site/ΣC); the
   literature also uses per-mole-of-cations conventions; state yours when
   comparing numbers.
6. **Fixed stoichiometry** — no oxygen non-stoichiometry, no spinel inversion.

## 10. Program structure

| Block | Purpose |
|---|---|
| `csv_data` (embedded) | symmetric ΔH_ij matrix, kJ/mol, ~70 elements |
| robust loader | `pd.read_csv(..., engine='python', on_bad_lines='warn')` → numeric coercion → missing → 0 |
| structures registry | ABO₃ / AB₂O₄ / AO with per-site counts {A, B, O} |
| `get_site_composition()` | interactive parser (`"Na 0.2 Bi 0.2 …"` → normalized dict) |
| ΔH_mix block | per-site regular solution Σ 4·h_ij·x_i·x_j × site weight |
| ΔS_mix^N block | −R · Σ_sites w·Σ x·ln x |
| results & plot | accumulate (ΔH, ΔS, structure) → color-coded scatter + legend + labels |

## 11. References

1. A. Takeuchi, A. Inoue, *Classification of bulk metallic glasses…*, Materials
   Transactions 46, 2817–2829 (2005) — source of the binary mixing enthalpies.
2. C. M. Rost et al., *Entropy-stabilized oxides*, Nature Communications 6,
   8485 (2015).
3. A. F. Manchón-Gordón, J. S. Blázquez, et al., *Descriptors for Predicting
   Single- and Multi-Phase Formation in High-Entropy Oxides: A Unified
   Framework Approach*, Materials 2025, 18, 3862.
   https://doi.org/10.3390/ma18163862 — the ΔS_mix^N / size-misfit descriptor
   framework.

## 12. License

MIT — see `LICENSE`.