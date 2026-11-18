# Subevent inversion results: 2024-12-05-mww70-off-coast-of-northern

- Selected subevents: **4** | best misfit **3.7900**
- Total moment: **Mw 7.04** (catalog Mw 7.0)

## L-curve (misfit vs subevents)
```
n  best_misfit  reduction_%  delta_aic  chains  chains_within_1%  island_bests
1       4.7984            -       34.4      32                32  4.7954,4.8006,4.8011,4.8044
2       4.2606         11.2        6.8      32                32  4.2541,4.2544,4.2600,4.2614
3       4.3285         -1.6       33.1      32                23  4.3219,4.3235,4.3255,4.3285
4       3.7900         12.4        0.0      32                10  3.8367,3.8766,4.3779,4.3811
5       3.8806         -2.4       29.4      32                 9  3.8862,3.8959,3.9132,3.9170
misfit components of each best model: penalized = data_misfit x clvd_factor x stress_factor (evaluate_misfit_components.sh)
n    penalized  data_misfit  clvd_factor  stress_factor  worst_axis_deg  p_misfit  sh_misfit  regional_misfit
N=1     4.7984       4.5538       1.0535         1.0002             0.0    16.715      2.583            3.492
N=2     4.2606       3.4572       1.1864         1.0388            29.4     5.107      2.812            2.735
N=3     4.3285       3.5007       1.1888         1.0401            29.8     5.208      2.830            2.732
N=4     3.7900       3.0792       1.1629         1.0584            29.9     3.605      2.820            2.784
N=5     3.8806       3.1593       1.1591         1.0597            30.0     3.978      2.789            2.773
selection: AIC = ln(misfit) / sigma^2 + 2 k N, sigma = 0.05, k = 10 parameters per subevent, so a subevent costs 5% of misfit (ln(misfit) has to fall by 0.05, the misfit by 4.9%, to tie); delta_aic = AIC above the lowest; every N with delta_aic below 2 is tied, and the smallest of them is selected
lowest AIC: N=4
SUGGESTED number of subevents (smallest N with delta_aic below 2): 4
```

## Best model

| # | centroid (s) | X east (km) | Y north (km) | duration (s) | depth (km) | Mw | from E1 (km) | speed since E1 began (km/s) | nodal planes (strike/dip/rake) | angle from the total mechanism (deg) |
|---|---|---|---|---|---|---|---|---|---|---|
| E1 | 5.16 | 0.00 | 0.00 | 5.27 | 4.0 | 6.20 | 0 | - | 9/82/11 ; 277/79/172 | 15 |
| E2 | 7.71 | 13.71 | 11.23 | 3.02 | 12.3 | 6.89 | 20 | 3.78 | 192/76/23 ; 96/67/165 | 29 |
| E3 | 11.92 | 29.31 | 12.61 | 3.07 | 13.5 | 6.83 | 33 | 3.54 | 185/82/-40 ; 282/50/-169 | 33 |
| E4 | 17.42 | -4.56 | -4.71 | 3.50 | 4.0 | 5.64 | 7 | 0.44 | 281/75/-179 ; 191/89/-15 | 9 |

E1 sits at the hypocentre when there are two or more subevents. The speed is the distance from E1 divided by the time since E1 began (its centroid time minus half its duration); the search caps it (`vmax_cau` in `search_par.file`). The angle is measured between each moment tensor and the double couple of the catalogue tensor; the mechanism prior acts on it (`MechanismSigma` in `Par.file`).

## Figures

- `subevents.pdf/png` — map (lon/lat, 95% error bars), depth section, STFs
- `fits_P.pdf`, `fits_Pvel.pdf`, `fits_SH.pdf`, `fits_rayl.pdf` — waveform fits
- `misfit_evolution.pdf` — misfit vs sample number, all chains
- `histoplot_py.pdf` — posterior histograms (gated ensemble)
- `lcurve.pdf` — model selection: best misfit per number of subevents (star: the selected number, the smallest of those tied with the lowest AIC; dashed line: a misfit 5% lower per added subevent, the cost of a subevent; hollow points: numbers that cannot be selected, their last subevent lies on the bound of the centroid range; small open circles: island bests) and, to its right, the factors of the misfit

## Convergence of the chains
```
chains: 32; best misfit 3.7900
chains whose best is within 1% of the global best: 10
chains whose best is within 2% of the global best: 23
chains whose best is within 5% of the global best: 32
split R-hat of the gated chains (1 = converged, > 1.1 = chains disagree):
  E1: centroid time 1.65, duration 1.73
  E2: centroid time 2.33, x east 2.57, y north 2.09, duration 2.02, depth 2.92
  E3: centroid time 2.63, x east 2.53, y north 2.03, duration 1.62, depth 3.81
  E4: centroid time 2.05, x east 4.38, y north 1.45, duration 1.52, depth 3.36
largest split R-hat: 4.38
```

