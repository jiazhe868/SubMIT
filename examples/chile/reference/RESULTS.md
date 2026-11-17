# Subevent inversion results: 2024-07-19-mww74-chile-argentina-border-region

- Selected subevents: **4** | best misfit **0.7318**
- Total moment: **Mw 7.41** (catalog Mw 7.4)

## L-curve (misfit vs subevents)
```
n  best_misfit  reduction_%  delta_aic  chains  chains_within_1%  island_bests
1       1.1946            -      136.0      32                32  1.1946,1.1948,1.1951,1.1956
2       0.9920         17.0       81.7      32                32  0.9928,0.9930,0.9937,0.9939
3       0.8396         15.4       35.0      32                16  0.8400,0.8415,0.8474,0.9663
4       0.7318         12.8        0.0      32                16  0.7369,0.8415,0.8445,0.8498
5       0.7371         -0.7       22.9      32                 9  0.7382,0.7415,0.7417,0.7436
misfit components of each best model: penalized = data_misfit x clvd_factor x stress_factor (evaluate_misfit_components.sh)
n    penalized  data_misfit  clvd_factor  stress_factor  worst_axis_deg  p_misfit  sh_misfit  regional_misfit
N=1     1.1946       1.0319       1.1577         1.0000             0.1     1.339      0.669            1.445
N=2     0.9920       0.8656       1.1460         1.0000             6.2     1.196      0.575            1.136
N=3     0.8396       0.8354       1.0051         1.0000            29.9     1.131      0.591            1.046
N=4     0.7318       0.7094       1.0316         1.0000            29.2     0.929      0.553            0.767
N=5     0.7371       0.7149       1.0310         1.0000            29.3     0.956      0.548            0.773
selection: AIC = ln(misfit) / sigma^2 + 2 k N, sigma = 0.05, k = 10 parameters per subevent, so a subevent costs 5% of misfit (ln(misfit) has to fall by 0.05, the misfit by 4.9%, to tie); delta_aic = AIC above the lowest; every N with delta_aic below 2 is tied, and the smallest of them is selected
lowest AIC: N=4
SUGGESTED number of subevents (smallest N with delta_aic below 2): 4
```

## Best model

| # | centroid (s) | X east (km) | Y north (km) | duration (s) | depth (km) | Mw |
|---|---|---|---|---|---|---|
| E1 | 4.24 | 0.00 | 0.00 | 8.37 | 142.1 | 7.07 |
| E2 | 7.83 | -2.45 | -9.46 | 4.28 | 152.8 | 6.93 |
| E3 | 10.62 | 1.97 | 2.23 | 7.34 | 170.4 | 7.02 |
| E4 | 13.41 | 19.77 | 12.06 | 10.29 | 162.9 | 7.07 |

## Figures

- `subevents.pdf/png` — map (lon/lat, 95% error bars), depth section, STFs
- `fits_P.pdf`, `fits_Pvel.pdf`, `fits_SH.pdf`, `fits_rayl.pdf` — waveform fits
- `misfit_evolution.pdf` — misfit vs sample number, all chains
- `histoplot_py.pdf` — posterior histograms (gated ensemble)
- `lcurve.pdf` — model selection: best misfit per number of subevents (star: the selected number, the smallest of those tied with the lowest AIC; dashed line: a misfit 5% lower per added subevent, the cost of a subevent; hollow points: numbers that cannot be selected, their last subevent lies on the bound of the centroid range; small open circles: island bests) and, to its right, the factors of the misfit

## Convergence of the chains
```
chains: 32; best misfit 0.7318
chains whose best is within 1% of the global best: 16
chains whose best is within 2% of the global best: 26
chains whose best is within 5% of the global best: 32
split R-hat of the gated chains (1 = converged, > 1.1 = chains disagree):
  E1: centroid time 1.51, duration 1.40, depth 1.56
  E2: centroid time 1.66, x east 1.57, y north 1.53, duration 1.27, depth 1.61
  E3: centroid time 1.75, x east 2.41, y north 1.96, duration 1.55, depth 1.42
  E4: centroid time 1.79, x east 1.28, y north 1.82, duration 1.66, depth 1.60
largest split R-hat: 2.41
```

