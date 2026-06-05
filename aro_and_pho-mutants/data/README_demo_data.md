##Demo Data — `aro_expts_demo.csv` 

This file is a **synthetic** stand-in for the real experimental data. It is
structurally identical to the original and are intended for testing and sharing
the analysis script (`aro_pho-analysis.R`).

---

## File structure

| Column | Description |
|---|---|
| `Sample` | Strain and treatment label, formatted to match the raw instrument export (including intentional inconsistencies such as trailing spaces, `\n` separators, and spacing variants — preserved so the cleaning steps in the R script are properly exercised) |
| `Replicate` | Replicate number (1–3) |
| `0`, `0.5`, `1`, … `30` | OD readings at each time point (hours); 61 columns covering t = 0 to 30 h in 0.5 h steps |

---

## How the values were generated

Each row is a simulated growth curve drawn from a logistic sigmoid:

$$y(t) = \frac{P}{1 + e^{-r'(t - L')}} - y_0$$

| Parameter | Meaning |
|---|---|
| *P* (plateau) | Approximate saturation OD |
| *r* (rate) | Nominal steepness of exponential growth phase |
| *L* (lag) | Nominal inflection point — mid-lag in hours |
| *r'*, *L'* | Per-replicate draws: *L'* ~ N(*L*, 1.2²), *r'* ~ N(*r*, (0.06*r*)²) |
| *y₀* | Value at *t* = 0 (subtracted so curves start near zero) |

Replicate-to-replicate variation is introduced by independently jittering the
lag (~±1.2 h) and growth rate (~±6%) for each replicate. This produces
visible SD ribbons in the averaged plots while keeping individual curves
smooth. No additional within-curve noise is added.

The nominal parameters were chosen by eye to reproduce the broad
characteristics of the real data (lag duration, growth rate, plateau OD) for
each strain/treatment combination. Biologically meaningful differences are
preserved directionally — for example, aro4 +8 mM Phe has a substantially
delayed lag (~34 h) and reduced plateau (~0.9) relative to WT, consistent with
the expected growth defect. The values are **not** calibrated to any specific
real measurement and should not be used for biological inference.


