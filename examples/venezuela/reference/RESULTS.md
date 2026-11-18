# Subevent inversion results: 2026-06-24-mww72-near-coast-of-venezuela

- Selected subevents: **5** | best misfit **2.2433**
- Total moment: **Mw 7.59** (catalog Mw 7.2)

## L-curve (misfit vs subevents)
```
n  best_misfit  reduction_%  delta_aic  chains  chains_within_1%  island_bests
1       4.0896            -      160.2      32                27  4.0937,4.2518,4.7395,4.8220
2       3.6323         11.2      132.8      32                32  3.6351,3.6356,3.6380,3.8942
3       2.9781         18.0       73.3      32                20  2.9573,2.9714,2.9788,2.9815
4       2.5879         13.1       37.2      32                25  2.5719,2.5797,2.5859,2.5867
5       2.2433         13.3        0.0      32                10  2.2474,2.2509,2.2519,2.2675
6       2.2200          1.0       15.8      32                 3  2.2049,2.2157,2.2159,2.2187
misfit components of each best model: penalized = data_misfit x clvd_factor x stress_factor (evaluate_misfit_components.sh)
n    penalized  data_misfit  clvd_factor  stress_factor  worst_axis_deg  p_misfit  sh_misfit  regional_misfit
N=1     4.0896       3.9579       1.0331         1.0002             0.1     3.625      6.569            3.126
N=2     3.6323       3.3994       1.0492         1.0184            19.4     5.013      2.973            2.297
N=3     2.9781       2.7518       1.0522         1.0285            24.7     3.697      3.188            1.673
N=4     2.5879       2.4397       1.0381         1.0218            21.7     3.044      2.908            1.588
N=5     2.2433       2.0986       1.0416         1.0262            26.3     2.513      2.768            0.842
N=6     2.2200       2.0576       1.0472         1.0303            27.1     2.518      2.650            0.919
selection: AIC = ln(misfit) / sigma^2 + 2 k N, sigma = 0.05, k = 10 parameters per subevent, so a subevent costs 5% of misfit (ln(misfit) has to fall by 0.05, the misfit by 4.9%, to tie); delta_aic = AIC above the lowest; every N with delta_aic below 2 is tied, and the smallest of them is selected
lowest AIC: N=5
SUGGESTED number of subevents (smallest N with delta_aic below 2): 5
```

## Best model

| # | centroid (s) | X east (km) | Y north (km) | duration (s) | depth (km) | Mw | from E1 (km) | speed since E1 began (km/s) | nodal planes (strike/dip/rake) | angle from the total mechanism (deg) |
|---|---|---|---|---|---|---|---|---|---|---|
| E1 | 3.94 | 0.00 | 0.00 | 6.59 | 20.3 | 6.45 | 0 | - | 86/82/-176 ; 356/86/-8 | 6 |
| E2 | 32.17 | 119.19 | 17.88 | 8.77 | 22.3 | 7.21 | 121 | 3.82 | 355/84/-12 ; 86/78/-174 | 6 |
| E3 | 40.79 | 128.10 | 15.58 | 3.72 | 20.4 | 7.30 | 129 | 3.21 | 87/83/177 ; 177/87/7 | 4 |
| E4 | 51.13 | 137.22 | 25.99 | 3.58 | 16.7 | 7.19 | 140 | 2.77 | 86/85/-171 ; 355/81/-5 | 12 |
| E5 | 108.83 | 153.04 | 45.83 | 7.40 | 17.6 | 6.96 | 160 | 1.48 | 177/67/22 ; 78/70/155 | 28 |

E1 sits at the hypocentre when there are two or more subevents. The speed is the distance from E1 divided by the time since E1 began (its centroid time minus half its duration); the search caps it (`vmax_cau` in `search_par.file`). The angle is measured between each moment tensor and the double couple of the catalogue tensor; the mechanism prior acts on it (`MechanismSigma` in `Par.file`).

## Figures

- `subevents.pdf/png` — map (lon/lat, 95% error bars), depth section, STFs
- `fits_P.pdf`, `fits_Pvel.pdf`, `fits_SH.pdf`, `fits_rayl.pdf` — waveform fits
- `misfit_evolution.pdf` — misfit vs sample number, all chains
- `histoplot_py.pdf` — posterior histograms (gated ensemble)
- `lcurve.pdf` — model selection: best misfit per number of subevents (star: the selected number, the smallest of those tied with the lowest AIC; dashed line: a misfit 5% lower per added subevent, the cost of a subevent; hollow points: numbers that cannot be selected, their last subevent lies on the bound of the centroid range; small open circles: island bests) and, to its right, the factors of the misfit

## Convergence of the chains
```
chains: 32; best misfit 2.2433
chains whose best is within 1% of the global best: 10
chains whose best is within 2% of the global best: 23
chains whose best is within 5% of the global best: 32
split R-hat of the gated chains (1 = converged, > 1.1 = chains disagree):
  E1: centroid time 1.23, duration 1.23
  E2: centroid time 1.47, x east 1.36, y north 1.24, duration 1.42, depth 1.30
  E3: centroid time 1.22, x east 1.72, y north 1.24, duration 1.06, depth 1.32
  E4: centroid time 1.27, x east 1.66, y north 1.25, duration 1.07, depth 1.24
  E5: centroid time 1.28, x east 1.29, y north 1.28, duration 1.07, depth 1.18
largest split R-hat: 1.72
```

