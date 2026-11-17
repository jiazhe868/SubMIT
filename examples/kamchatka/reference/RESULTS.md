# Subevent inversion results: 2025-07-29-mww88-off-east-coast-of-kamchatka

- Selected subevents: **7** | best misfit **0.8193**
- Total moment: **Mw 8.84** (catalog Mw 8.8)

## L-curve (misfit vs subevents)
```
n  best_misfit  reduction_%  delta_aic  chains  chains_within_1%  island_bests
1       2.9090            -      386.8      32                32  2.9085,3.0840,3.0861,3.0863
2       2.2403         23.0      302.4      32                32  2.2404,2.2404,2.2420,2.2442
3       1.7973         19.8      234.2      32                23  1.7991,1.8002,1.8022,1.8036
4       1.4391         19.9      165.3      32                15  1.4408,1.4571,1.4575,1.4642
5       1.1170         22.4       84.0      32                14  1.1202,1.2252,1.2298,1.2300
6       0.9065         18.8       20.5      32                10  0.9066,0.9086,0.9199,0.9321
7       0.8193          9.6        0.0      32                 6  0.8139,0.8216,0.8500,0.8707
misfit components of each best model: penalized = data_misfit x clvd_factor x stress_factor (evaluate_misfit_components.sh)
n    penalized  data_misfit  clvd_factor  stress_factor  worst_axis_deg  p_misfit  sh_misfit  regional_misfit
N=1     2.9090       2.9082       1.0003         1.0000             0.0     2.872      2.948            5.786
N=2     2.2403       2.0863       1.0738         1.0000            15.7     2.175      1.994           10.733
N=3     1.7973       1.6440       1.0933         1.0000            10.4     1.699      1.585            9.832
N=4     1.4391       1.3491       1.0667         1.0000            18.8     1.309      1.393           10.512
N=5     1.1170       1.0449       1.0690         1.0000            16.5     1.034      1.057            9.430
N=6     0.9065       0.8533       1.0623         1.0000            17.6     0.849      0.858            8.087
N=7     0.8193       0.7771       1.0542         1.0000            22.9     0.777      0.775            8.199
selection: AIC = ln(misfit) / sigma^2 + 2 k N, sigma = 0.05, k = 10 parameters per subevent, so a subevent costs 5% of misfit (ln(misfit) has to fall by 0.05, the misfit by 4.9%, to tie); delta_aic = AIC above the lowest; every N with delta_aic below 2 is tied, and the smallest of them is selected
NOTE: the AIC is lowest at N=7, the largest number of subevents tried; more subevents might fit better, so read 7 as a lower bound
lowest AIC: N=7
SUGGESTED number of subevents (smallest N with delta_aic below 2): 7
```

## Best model

| # | centroid (s) | X east (km) | Y north (km) | duration (s) | depth (km) | Mw |
|---|---|---|---|---|---|---|
| E1 | 32.26 | 0.00 | 0.00 | 28.16 | 21.1 | 8.12 |
| E2 | 62.53 | -48.88 | -36.93 | 26.81 | 36.3 | 7.98 |
| E3 | 96.18 | -106.13 | -179.04 | 34.41 | 35.4 | 8.18 |
| E4 | 112.24 | -101.69 | -198.71 | 31.81 | 10.3 | 8.44 |
| E5 | 115.65 | -203.77 | -246.60 | 31.66 | 9.7 | 8.40 |
| E6 | 133.39 | -209.33 | -248.14 | 32.39 | 20.1 | 8.37 |
| E7 | 156.20 | -235.11 | -283.37 | 34.58 | 29.7 | 8.27 |

## Figures

- `subevents.pdf/png` — map (lon/lat, 95% error bars), depth section, STFs
- `fits_P.pdf`, `fits_Pvel.pdf`, `fits_SH.pdf`, `fits_rayl.pdf` — waveform fits
- `misfit_evolution.pdf` — misfit vs sample number, all chains
- `histoplot_py.pdf` — posterior histograms (gated ensemble)
- `lcurve.pdf` — model selection: best misfit per number of subevents (star: the selected number, the smallest of those tied with the lowest AIC; dashed line: a misfit 5% lower per added subevent, the cost of a subevent; hollow points: numbers that cannot be selected, their last subevent lies on the bound of the centroid range; small open circles: island bests) and, to its right, the factors of the misfit

## Convergence of the chains
```
chains: 32; best misfit 0.8193
chains whose best is within 1% of the global best: 6
chains whose best is within 2% of the global best: 21
chains whose best is within 5% of the global best: 31
split R-hat of the gated chains (1 = converged, > 1.1 = chains disagree):
  E1: centroid time 1.10, duration 1.20, depth 1.37
  E2: centroid time 1.09, x east 1.07, y north 1.05, duration 1.13, depth 1.31
  E3: centroid time 1.15, x east 1.58, y north 1.15, duration 1.09, depth 1.49
  E4: centroid time 1.22, x east 1.64, y north 1.20, duration 1.36, depth 1.82
  E5: centroid time 4.47, x east 1.13, y north 1.20, duration 1.27, depth 1.47
  E6: centroid time 7.16, x east 1.29, y north 1.28, duration 1.63, depth 2.10
  E7: centroid time 1.86, x east 1.21, y north 1.17, duration 1.58, depth 1.54
largest split R-hat: 7.16
```

