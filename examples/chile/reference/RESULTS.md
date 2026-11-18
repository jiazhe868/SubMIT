# Subevent inversion results: 2024-07-19-mww74-chile-argentina-border-region

- Selected subevents: **4** | best misfit **0.8124**
- Total moment: **Mw 7.41** (catalog Mw 7.4)

## L-curve (misfit vs subevents)
```
n  best_misfit  reduction_%  delta_aic  chains  chains_within_1%  island_bests
1       1.1950            -       94.4      32                32  1.1951,1.1954,1.1957,1.1957
2       1.0269         14.1       53.7      32                31  1.0266,1.0270,1.0283,1.0300
3       0.8733         15.0        8.9      32                18  0.8814,1.0100,1.0109,1.0139
4       0.8124          7.0        0.0      32                27  0.8106,0.8114,0.8116,0.8117
5       0.7840          3.5        5.8      32                20  0.7806,0.7858,0.8102,0.8160
misfit components of each best model: penalized = data_misfit x clvd_factor x stress_factor (evaluate_misfit_components.sh)
n    penalized  data_misfit  clvd_factor  stress_factor  worst_axis_deg  p_misfit  sh_misfit  regional_misfit
N=1     1.1950       1.0310       1.1571         1.0018             0.1     1.338      0.663            1.450
N=2     1.0269       0.8883       1.1511         1.0042             3.8     1.144      0.660            1.023
N=3     0.8733       0.7533       1.1506         1.0076             5.0     1.049      0.565            0.729
N=4     0.8124       0.7023       1.1450         1.0102             7.2     0.950      0.545            0.689
N=5     0.7840       0.6800       1.1366         1.0144            17.0     0.877      0.551            0.695
selection: AIC = ln(misfit) / sigma^2 + 2 k N, sigma = 0.05, k = 10 parameters per subevent, so a subevent costs 5% of misfit (ln(misfit) has to fall by 0.05, the misfit by 4.9%, to tie); delta_aic = AIC above the lowest; every N with delta_aic below 2 is tied, and the smallest of them is selected
lowest AIC: N=4
SUGGESTED number of subevents (smallest N with delta_aic below 2): 4
```

## Best model

| # | centroid (s) | X east (km) | Y north (km) | duration (s) | depth (km) | Mw | from E1 (km) | speed since E1 began (km/s) | nodal planes (strike/dip/rake) | angle from the total mechanism (deg) |
|---|---|---|---|---|---|---|---|---|---|---|
| E1 | 1.79 | 0.00 | 0.00 | 3.55 | 122.9 | 6.49 | 0 | - | 336/67/-97 ; 173/24/-74 | 7 |
| E2 | 4.57 | -4.40 | -1.54 | 8.83 | 135.5 | 7.08 | 13 | 2.95 | 346/68/-88 ; 162/22/-94 | 8 |
| E3 | 7.64 | 4.51 | 0.65 | 3.46 | 148.4 | 6.91 | 26 | 3.39 | 336/69/-101 ; 184/24/-64 | 14 |
| E4 | 10.37 | 8.99 | 1.47 | 9.12 | 168.6 | 7.19 | 47 | 4.49 | 340/72/-95 ; 175/19/-76 | 10 |

E1 sits at the hypocentre when there are two or more subevents. The speed is the distance from E1 divided by the time since E1 began (its centroid time minus half its duration); the search caps it (`vmax_cau` in `search_par.file`). The angle is measured between each moment tensor and the double couple of the catalogue tensor; the mechanism prior acts on it (`MechanismSigma` in `Par.file`).

## Figures

- `subevents.pdf/png` — map (lon/lat, 95% error bars), depth section, STFs
- `fits_P.pdf`, `fits_Pvel.pdf`, `fits_SH.pdf`, `fits_rayl.pdf` — waveform fits
- `misfit_evolution.pdf` — misfit vs sample number, all chains
- `histoplot_py.pdf` — posterior histograms (gated ensemble)
- `lcurve.pdf` — model selection: best misfit per number of subevents (star: the selected number, the smallest of those tied with the lowest AIC; dashed line: a misfit 5% lower per added subevent, the cost of a subevent; hollow points: numbers that cannot be selected, their last subevent lies on the bound of the centroid range; small open circles: island bests) and, to its right, the factors of the misfit

## Convergence of the chains
```
chains: 32; best misfit 0.8124
chains whose best is within 1% of the global best: 27
chains whose best is within 2% of the global best: 31
chains whose best is within 5% of the global best: 32
split R-hat of the gated chains (1 = converged, > 1.1 = chains disagree):
  E1: centroid time 1.35, duration 1.35
  E2: centroid time 1.32, x east 1.17, y north 1.05, duration 1.20, depth 1.10
  E3: centroid time 1.24, x east 1.13, y north 1.15, duration 1.33, depth 1.26
  E4: centroid time 1.24, x east 1.13, y north 1.14, duration 1.29, depth 1.23
largest split R-hat: 1.35
```

