# Subevent inversion results: 2025-07-29-mww88-off-east-coast-of-kamchatka

- Selected subevents: **7** | best misfit **0.8823**
- Total moment: **Mw 8.83** (catalog Mw 8.8)

## L-curve (misfit vs subevents)
```
n  best_misfit  reduction_%  delta_aic  chains  chains_within_1%  island_bests
1       2.9170            -      358.3      32                32  2.9157,3.1092,3.1111,3.1124
2       2.3524         19.4      292.3      32                32  2.3551,2.3559,2.3567,2.3591
3       1.9097         18.8      228.9      32                30  1.9108,1.9134,1.9498,2.0402
4       1.4560         23.8      140.4      32                30  1.4556,1.4576,1.5416,1.5439
5       1.1801         18.9       76.3      32                 7  1.1809,1.1849,1.1859,1.1909
6       0.9742         17.4       19.6      32                 1  0.9835,0.9835,0.9862,0.9864
7       0.8823          9.4        0.0      32                 6  0.8813,0.8855,0.8891,0.9394
misfit components of each best model: penalized = data_misfit x clvd_factor x stress_factor (evaluate_misfit_components.sh)
n    penalized  data_misfit  clvd_factor  stress_factor  worst_axis_deg  p_misfit  sh_misfit  regional_misfit
N=1     2.9170       2.9122       1.0004         1.0012             0.1     2.900      2.925            6.205
N=2     2.3524       2.3458       1.0008         1.0020             3.6     2.341      2.349           13.957
N=3     1.9097       1.7425       1.0934         1.0023             2.9     1.752      1.732            9.496
N=4     1.4560       1.3489       1.0665         1.0121            12.6     1.440      1.245            9.549
N=5     1.1801       1.0863       1.0588         1.0261            13.7     1.049      1.131           11.395
N=6     0.9742       0.8899       1.0616         1.0312            15.8     0.894      0.885            7.457
N=7     0.8823       0.7999       1.0682         1.0325            16.2     0.778      0.828            7.446
selection: AIC = ln(misfit) / sigma^2 + 2 k N, sigma = 0.05, k = 10 parameters per subevent, so a subevent costs 5% of misfit (ln(misfit) has to fall by 0.05, the misfit by 4.9%, to tie); delta_aic = AIC above the lowest; every N with delta_aic below 2 is tied, and the smallest of them is selected
NOTE: the AIC is lowest at N=7, the largest number of subevents tried; more subevents might fit better, so read 7 as a lower bound
lowest AIC: N=7
SUGGESTED number of subevents (smallest N with delta_aic below 2): 7
```

## Best model

| # | centroid (s) | X east (km) | Y north (km) | duration (s) | depth (km) | Mw | from E1 (km) | speed since E1 began (km/s) | nodal planes (strike/dip/rake) | angle from the total mechanism (deg) |
|---|---|---|---|---|---|---|---|---|---|---|
| E1 | 32.61 | 0.00 | 0.00 | 30.32 | 20.7 | 8.17 | 0 | - | 58/75/102 ; 200/19/53 | 7 |
| E2 | 61.35 | -66.80 | -62.42 | 24.26 | 39.4 | 7.87 | 93 | 2.13 | 205/32/67 ; 51/61/103 | 28 |
| E3 | 97.88 | -138.71 | -159.97 | 34.51 | 35.3 | 8.20 | 212 | 2.64 | 52/71/100 ; 205/22/64 | 13 |
| E4 | 114.19 | -137.62 | -225.34 | 33.25 | 17.5 | 8.45 | 264 | 2.73 | 59/77/101 ; 197/17/49 | 3 |
| E5 | 134.45 | -196.30 | -258.44 | 34.96 | 27.9 | 8.35 | 325 | 2.77 | 58/75/98 ; 208/17/62 | 6 |
| E6 | 153.22 | -225.16 | -281.08 | 23.20 | 39.9 | 8.06 | 361 | 2.66 | 56/71/103 ; 200/23/56 | 10 |
| E7 | 164.03 | -229.33 | -306.12 | 34.61 | 13.7 | 8.42 | 383 | 2.61 | 61/81/104 ; 184/17/34 | 10 |

E1 sits at the hypocentre when there are two or more subevents. The speed is the distance from E1 divided by the time since E1 began (its centroid time minus half its duration); the search caps it (`vmax_cau` in `search_par.file`). The angle is measured between each moment tensor and the double couple of the catalogue tensor; the mechanism prior acts on it (`MechanismSigma` in `Par.file`).

## Figures

- `subevents.pdf/png` — map (lon/lat, 95% error bars), depth section, STFs
- `fits_P.pdf`, `fits_Pvel.pdf`, `fits_SH.pdf`, `fits_rayl.pdf` — waveform fits
- `misfit_evolution.pdf` — misfit vs sample number, all chains
- `histoplot_py.pdf` — posterior histograms (gated ensemble)
- `lcurve.pdf` — model selection: best misfit per number of subevents (star: the selected number, the smallest of those tied with the lowest AIC; dashed line: a misfit 5% lower per added subevent, the cost of a subevent; hollow points: numbers that cannot be selected, their last subevent lies on the bound of the centroid range; small open circles: island bests) and, to its right, the factors of the misfit

## Convergence of the chains
```
chains: 32; best misfit 0.8823
chains whose best is within 1% of the global best: 6
chains whose best is within 2% of the global best: 19
chains whose best is within 5% of the global best: 32
split R-hat of the gated chains (1 = converged, > 1.1 = chains disagree):
  E1: centroid time 1.06, duration 1.21
  E2: centroid time 1.14, x east 1.10, y north 1.08, duration 1.20, depth 1.33
  E3: centroid time 1.13, x east 1.46, y north 1.13, duration 1.10, depth 1.57
  E4: centroid time 1.18, x east 1.54, y north 1.15, duration 1.08, depth 2.18
  E5: centroid time 1.15, x east 1.10, y north 1.16, duration 1.16, depth 1.25
  E6: centroid time 1.18, x east 1.28, y north 1.22, duration 1.30, depth 1.31
  E7: centroid time 1.20, x east 1.19, y north 1.13, duration 1.25, depth 1.16
largest split R-hat: 2.18
```

