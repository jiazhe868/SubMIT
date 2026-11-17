# Subevent inversion results: 2026-06-24-mww72-near-coast-of-venezuela

- Selected subevents: **5** | best misfit **2.0562**
- Total moment: **Mw 7.59** (catalog Mw 7.2)

## L-curve (misfit vs subevents)
```
n  best_misfit  reduction_%  delta_aic  chains  chains_within_1%  island_bests
1       4.1238            -      198.4      32                22  4.1482,4.3095,4.7555,4.8563
2       3.2735         20.6      126.0      32                32  3.2764,3.2805,3.3925,3.4261
3       2.5324         22.6       43.3      32                 4  2.5556,2.6558,2.7738,2.7805
4       2.3417          7.5       32.0      32                16  2.3364,2.3453,2.3478,2.3528
5       2.0562         12.2        0.0      32                 5  2.0642,2.1526,2.1845,2.2381
6       1.9567          4.8        0.2      32                 6  1.9585,2.0209,2.0215,2.0420
misfit components of each best model: penalized = data_misfit x clvd_factor x stress_factor (evaluate_misfit_components.sh)
n    penalized  data_misfit  clvd_factor  stress_factor  worst_axis_deg  p_misfit  sh_misfit  regional_misfit
N=1     4.1238       3.9933       1.0327         1.0000             0.0     3.601      6.557            3.166
N=2     3.2735       3.0977       1.0568         1.0000            24.6     4.147      2.845            1.995
N=3     2.5324       2.3787       1.0646         1.0000            29.9     2.794      2.954            1.397
N=4     2.3417       2.1791       1.0746         1.0000            29.9     2.313      2.937            1.373
N=5     2.0562       1.9728       1.0423         1.0000            29.9     2.293      2.794            0.927
N=6     1.9567       1.8905       1.0350         1.0000            30.0     2.034      2.890            1.019
selection: AIC = ln(misfit) / sigma^2 + 2 k N, sigma = 0.05, k = 10 parameters per subevent, so a subevent costs 5% of misfit (ln(misfit) has to fall by 0.05, the misfit by 4.9%, to tie); delta_aic = AIC above the lowest; every N with delta_aic below 2 is tied, and the smallest of them is selected
NOTE: a tie - N=5 (delta_aic 0.0), N=6 (delta_aic 0.2) are equally supported (delta_aic below 2); N=5 is selected as the smallest: report both
lowest AIC: N=5
SUGGESTED number of subevents (smallest N with delta_aic below 2): 5
```

## Best model

| # | centroid (s) | X east (km) | Y north (km) | duration (s) | depth (km) | Mw |
|---|---|---|---|---|---|---|
| E1 | 15.77 | 0.00 | 0.00 | 4.53 | 21.0 | 6.50 |
| E2 | 38.72 | 119.62 | 43.52 | 24.24 | 25.9 | 7.43 |
| E3 | 51.71 | 130.02 | 53.85 | 9.44 | 10.7 | 7.19 |
| E4 | 59.59 | 135.66 | 52.01 | 5.17 | 9.8 | 6.91 |
| E5 | 108.16 | 168.60 | 52.83 | 6.03 | 15.1 | 7.00 |

## Figures

- `subevents.pdf/png` — map (lon/lat, 95% error bars), depth section, STFs
- `fits_P.pdf`, `fits_Pvel.pdf`, `fits_SH.pdf`, `fits_rayl.pdf` — waveform fits
- `misfit_evolution.pdf` — misfit vs sample number, all chains
- `histoplot_py.pdf` — posterior histograms (gated ensemble)
- `lcurve.pdf` — model selection: best misfit per number of subevents (star: the selected number, the smallest of those tied with the lowest AIC; dashed line: a misfit 5% lower per added subevent, the cost of a subevent; hollow points: numbers that cannot be selected, their last subevent lies on the bound of the centroid range; small open circles: island bests) and, to its right, the factors of the misfit

## Convergence of the chains
```
chains: 32; best misfit 2.0562
chains whose best is within 1% of the global best: 5
chains whose best is within 2% of the global best: 17
chains whose best is within 5% of the global best: 32
split R-hat of the gated chains (1 = converged, > 1.1 = chains disagree):
  E1: centroid time 1.86, duration 1.57, depth 2.17
  E2: centroid time 1.69, x east 1.71, y north 1.56, duration 1.55, depth 1.66
  E3: centroid time 1.95, x east 1.81, y north 2.75, duration 1.17, depth 1.66
  E4: centroid time 1.84, x east 1.57, y north 2.58, duration 1.18, depth 1.21
  E5: centroid time 1.69, x east 1.35, y north 1.81, duration 1.11, depth 1.27
largest split R-hat: 2.75
```

