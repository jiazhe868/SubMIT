# SubMIT — Bayesian multiple-subevent inversion of large earthquakes

SubMIT describes a large earthquake as a **small number of subevents**: point sources that each
have their own time, place, depth, duration and moment tensor (the "type" of faulting). It finds the
subevents that best explain hundreds of seismograms recorded around the world and nearby,
**estimates their uncertainties with Markov-chain Monte Carlo (MCMC)**, and **chooses how many
subevents the data actually support**.

This version is a fully automated pipeline: starting from downloaded waveforms it screens the
data, computes Green's functions, builds the priors, runs the inversions for 1, 2, 3, … subevents,
selects the number of subevents and produces publication-style figures — with no hand tuning.
Four validated example earthquakes are included, and each can be reproduced with one command.

The method was introduced and applied in Jia et al. (2020a, 2020b, 2022a, 2025a); the full
references are in [Section 12](#12-credits-citation-and-licenses).

![How SubMIT sees an earthquake](docs/figures/concept.png)

*How the method works, using the 2025 Kamchatka M 8.8 result. **(1)** The earthquake is described
as seven subevents. **(2)** Each releases its seismic moment as a smooth pulse; their sum is the
earthquake's moment-rate function. **(3)** Because the subevents are at different places, a station
the rupture runs toward sees the pulses arrive bunched together, while a station on the other side
sees them spread out. Matching these differences across ~100–200 stations is what locates the
subevents in space and time.*

---

## Contents

1. [Quick start: reproduce an example in three commands](#1-quick-start)
2. [The four example earthquakes](#2-the-four-example-earthquakes)
3. [Installation](#3-installation)
4. [How the method works (with equations)](#4-how-the-method-works)
5. [The automated pipeline](#5-the-automated-pipeline)
6. [Results of the examples](#6-results-of-the-examples)
7. [Running your own earthquake](#7-running-your-own-earthquake)
8. [Reading the outputs](#8-reading-the-outputs)
9. [Using an AI coding assistant](#9-using-an-ai-coding-assistant)
10. [Troubleshooting](#10-troubleshooting)
11. [Repository layout](#11-repository-layout)
12. [Credits, citation and licenses](#12-credits-citation-and-licenses)

---

## 1. Quick start

On a Linux machine with SAC, an MPI library, Fortran/C compilers, Java and conda
([details](#3-installation)):

```bash
git clone https://github.com/jiazhe868/SubMIT.git && cd SubMIT
conda env create -f environment.yml && conda activate submit

./submit check                         # are all dependencies installed?
./submit build                         # compile the bundled Green's-function codes
./submit example california            # reproduce the 2024 California M7.0 result (~10 min)
./submit compare california            # check it against the published numbers
```

The figures and a summary (`RESULTS.md`) appear in `work/california/IRIS/figs_and_results/`.
`./submit compare` prints a table like this one (from a test run on a fresh install):

```
L-curve (best penalized misfit for each number of subevents)
  n  published  reproduced     diff
  1     2.7122      2.7122    +0.0%
  2     1.9677      1.9677    +0.0%
  3     1.7905      1.7905    +0.0%
  ...
selected number of subevents: published 3, reproduced 3

selected 3-subevent model (published -> reproduced)
          time (s)               x,y (km)     depth (km)      dur (s)          Mw  MT sim
E1     4.5->4.5        0,0   ->    0,0       3.1->3.1      3.3->3.3   6.55->6.55 +1.00
E2     7.7->7.7        2,-0  ->    2,-0      3.0->3.0      8.7->8.7   6.91->6.91 +1.00
E3    14.9->14.9      49,-18 ->   49,-18     3.7->3.7     11.7->11.7  6.62->6.62 +1.00
total Mw: published 7.04, reproduced 7.04
```

### Three ways to reproduce an example

| mode | what it does | time on 32 cores |
|------|--------------|------------------|
| `--mode figures` (default) | rebuilds the Green's functions and regenerates every figure from the published best models | 5–60 min |
| `--mode frozen` | re-runs all the MCMC inversions with the exact published configuration files, then selects the number of subevents | 1–14 h |
| `--mode full` | runs the whole automated pipeline from the processed waveforms: screening, weighting, priors, inversions, selection, figures | 1.5–16 h |

```bash
./submit example chile --mode frozen      # re-run the inversions of the Chile example
./submit example chile --mode full        # re-derive everything from the waveforms
./submit example chile --mode full --quick   # 10-minute smoke test (numbers not meaningful)
./submit status chile                     # progress of a running example
./submit list                             # all bundled examples
```

Long runs keep going after you log out if you start them with
`nohup ./submit example kamchatka --mode frozen > kam.log 2>&1 &`.

Measured inversion times for the whole subevent grid (32 MPI ranks, AMD EPYC 7532):
California ≈ 3 h, Chile ≈ 5 h, Venezuela ≈ 5 h, Kamchatka ≈ 14 h.
Green's functions add 5–60 min (more for deep events, which need a wide depth range).

`frozen` and `full` runs are not bit-for-bit identical to the published ones: the 32 chains
exchange states through MPI at times that depend on the machine's speed. Expect the same selected
number of subevents and best misfits within a few percent. `full` mode also reflects the current
code, so it can differ more: it builds Green's functions with the current velocity-model rules
(the published Chile run used an earlier regional model), and, for example, the Chile
station AU.MAW, dropped by hand in the published run, is now handled automatically.

---

## 2. The four example earthquakes

![Moment-rate functions of the four examples](docs/figures/results_overview.png)

| example | earthquake | depth | selected subevents | total Mw (catalog) | best misfit |
|---------|-----------|-------|:------------------:|:------------------:|:-----------:|
| `chile` | 2024-07-19 Chile–Argentina border | 123 km (intraslab) | **4** | 7.41 (7.4) | 0.869 |
| `california` | 2024-12-05 offshore Cape Mendocino | shallow strike-slip | **3** | 7.04 (7.0) | 1.791 |
| `venezuela` | 2026-06-24 near the coast of Venezuela | shallow | **6** | 7.58 (7.2) | 3.072 |
| `kamchatka` | 2025-07-29 Kamchatka megathrust | shallow subduction | **7** | 8.83 (8.8) | 1.107 |

They were chosen to span the problems the method has to handle: a deep event whose depth phases
arrive late (Chile), a strike-slip event with a dense regional network (California), a
complex event that appears to combine smooth, slower slip with short, sharp ruptures (Venezuela) and a 450-km-long
megathrust rupture (Kamchatka). Misfit values are only comparable between runs of the *same*
event ([why](#misfit)). Detailed results are in [Section 6](#6-results-of-the-examples).

---

## 3. Installation

SubMIT runs on Linux (tested on Ubuntu 22.04). You need:

| requirement | why | how to get it |
|-------------|-----|---------------|
| **SAC** (with `libsac.a`, `libsacio.a`) | reading/writing seismograms, data processing | request it from [IRIS/EarthScope](https://ds.iris.edu/ds/nodes/dmc/forms/sac/); put `sac` on your `PATH` and set `SACAUX` as in its install notes |
| **MPI** + **C and Fortran compilers** | the inversion runs one Markov chain per MPI rank | Open MPI or MPICH with `gcc`/`gfortran` (`sudo apt install libopenmpi-dev gfortran`), or Intel oneAPI (`mpiicx`/`mpiifx`, faster) |
| **Java** | TauP travel times (bundled) | `sudo apt install default-jre` |
| **Perl, gawk** | Green's-function and data scripts | usually preinstalled |
| **Python 3 + packages** | preparation, screening, figures | `conda env create -f environment.yml` |
| ~2–4 GB disk per example | Green's functions | |

Then:

```bash
conda activate submit
./submit check     # lists anything missing
./submit build     # compiles fk, mtel3 and a test build of the inversion code
```

The inversion and forward codes (`finv`, `ffwd`) are compiled automatically inside every run
directory, so they always match the configuration of that run. The build uses Intel oneAPI if it
finds it and GNU MPI wrappers otherwise; add `--gnu` to `./submit example` to force the GNU
toolchain (all checks in this README were run with GNU). With Intel, point `INTEL_SETVARS` at your
`setvars.sh` if it is not in `/opt/intel/oneapi/`.

Useful environment variables: `SUBMIT_PYTHON` (which Python to use), `SUBMIT_NCORE` (MPI ranks;
default = physical cores, at most 32), `SACHOME` (if `sac` is not on `PATH`).

---

## 4. How the method works

### 4.1 The subevent model

This parameterization follows Jia et al. (2020a, 2020b, 2022a). Each subevent $`i = 1 \dots N`$ has

| symbol | meaning | type |
|--------|---------|------|
| $`t_i`$ | centroid time after the origin time | nonlinear (searched) |
| $`x_i, y_i`$ | horizontal position east / north of the epicenter (with two or more subevents, subevent 1 is fixed at the epicenter; a single subevent can move) | nonlinear |
| $`z_i`$ | centroid depth | nonlinear |
| $`T_i`$ | duration | nonlinear |
| $`\mathbf m_i`$ | five independent components of the deviatoric moment tensor | linear (solved exactly) |

So $`N`$ subevents have $`5N`$ nonlinear and $`5N`$ linear unknowns. Rupture directivity *within* a
subevent can be switched on, but is off by default (`vr = 0.01 km/s` in `Input.model` makes each
subevent a point source); directivity of the whole earthquake emerges from the arrangement of the
subevents.

Each subevent releases its moment with a Gaussian moment-rate function whose full width at 10% of
the peak is its duration. As seen at station $`j`$ it is shifted by a time $`\Delta t_{ij}`$ that
depends on where the subevent is (Section 4.2):

```math
s_{ij}(t) = \frac{1}{\sqrt{2\pi}\,\sigma_i}\exp\!\left[-\frac{(t-t_i-\Delta t_{ij})^2}{2\sigma_i^2}\right],
\qquad \sigma_i = \frac{T_i}{2\sqrt{2\ln 10}} .
```

### 4.2 Forward problem: predicting a seismogram

The displacement predicted at station $`j`$ is the sum over subevents of moment-tensor components
times Green's functions, convolved with the subevent's moment-rate function:

```math
u_j(t) \;=\; \sum_{i=1}^{N}\sum_{k=1}^{5} m_{ik}\;\big[\,G_{jk}(z_i,\,t) * s_{ij}(t)\,\big].
```

For **teleseismic P and SH waves** (stations at about 40°–90°), Green's functions $`G_{jk}`$ are computed at the
station's distance on a 2-km depth grid (linearly interpolated to $`z_i`$) with the ray-theory code
`mtel3`, using a CRUST1.0 source-side crust. The subevent's position enters through a time shift
relative to the first subevent,

```math
\Delta t_{ij} \;=\; -\,\frac{r_i}{c_j}\cos(\phi_j-\psi_i)\;-\;(z_i-z_1)\,\eta_j,
\qquad \eta_j=\sqrt{v^{-2}-c_j^{-2}},
```

where $`(r_i,\psi_i)`$ are the distance and azimuth of the subevent from the epicenter, $`\phi_j`$
the station azimuth, $`c_j`$ the apparent horizontal velocity of the ray (from TauP) and $`\eta_j`$
the vertical slowness at the source (P or S velocity $`v`$). This is the term that produces the
bunching/stretching in the figure at the top.

For **regional three-component records** (up to ~950 km), Green's functions are computed with the
frequency–wavenumber code `fk` in a 1-D CRUST1.0 model on a distance–depth grid, and each subevent
uses the Green's function for its *own* distance and azimuth to the station (bilinear
interpolation), so the timing is exact rather than approximated.

All records are band-pass filtered identically to the data (teleseismic 0.005–0.2 Hz, regional
0.02–0.15 Hz by default; narrower for some extended ruptures).

### 4.3 Linear step: the moment tensors

For every trial set of nonlinear parameters $`\theta`$, the moment tensors are found by damped
least squares, with an extra constraint that their sum stays close to the catalog moment tensor
of the whole earthquake (GCMT or USGS):

```math
\hat{\mathbf m}(\theta) = \arg\min_{\mathbf m}\;
\big\|\mathbf W\big(\mathbf G(\theta)\,\mathbf m-\mathbf d\big)\big\|^2
\;+\;\alpha^2\|\mathbf m\|^2
\;+\;\lambda^2\Big\|\sum_{i=1}^N \mathbf m_i-\mathbf m_{\mathrm{CMT}}\Big\|^2 ,
\qquad \lambda = \kappa_{\mathrm{CMT}}\,\frac{\|\mathbf W\mathbf d\|}{M_0^{\mathrm{CMT}}} .
```

$`\mathbf W`$ holds the weights of the three wave types (P, SH, regional; [Section 5](#5-the-automated-pipeline)),
$`\alpha`$ is a small Tikhonov damping (`Tikhonov_alpha`) and $`\kappa_{\mathrm{CMT}}`$
(`CMTscaling`, default 3) sets how firmly the total is anchored; it is scaled by the data norm so
the same value works for an M 6.9 and an M 8.8. This is solved exactly (Cholesky) in a fraction of a second,
which is why the moment tensors do not need to be sampled.

### 4.4 Misfit and physical penalties <a name="misfit"></a>

Each predicted trace is aligned with the data by cross-correlation (shift
$`|\tau_j|\le\tau_{\max}`$: 2 s for P, 5 s for SH, distance-scaled for regional records) and the
misfit combines amplitude and waveform-shape agreement:

```math
\chi(\theta) \;=\; \frac{\displaystyle\sum_j \big\|d_j(t+\tau_j)-u_j(t)\big\|^2\; e^{\,|1-C_j|}}
                        {\displaystyle\sum_j \|d_j\|\,\|u_j\|},
```

where $`C_j`$ is the correlation coefficient. Values near 1 mean "residual energy comparable to the
signal", and they depend on the station set, windows and weights — so compare misfits only between
runs of the same event.

Two penalties keep the solution physical:

```math
\Phi(\theta) \;=\; \chi(\theta)\;
\exp\!\Big(\frac{\bar\epsilon}{\kappa}\Big)\;
\exp\!\Big(\max\!\big(0,\tfrac{\Theta_P}{30^\circ}-1\big)+\max\!\big(0,\tfrac{\Theta_T}{30^\circ}-1\big)\Big).
```

* **Double-couple penalty.** $`\bar\epsilon=\sum_i \sqrt{\epsilon_i}\,M_{0,i}/\sum_i M_{0,i}`$ is the
  moment-weighted non-double-couple (CLVD) content; $`\kappa`$ is `DCconstrain` (default 1).
* **Stress-consistency penalty.** $`\Theta_P`$ ($`\Theta_T`$) is the largest angle, over *all*
  subevents, between a subevent's pressure (tension) axis and that of the summed moment tensor.
  Mechanisms may rotate freely within 30° of the overall mechanism; beyond that the misfit grows
  exponentially. Every subevent counts, however small, so a tiny subevent cannot "soak up"
  unexplained energy with an implausible, opposite mechanism.

### 4.5 Bayesian sampling

The posterior probability of the nonlinear parameters is

```math
p(\theta\mid\mathbf d)\;\propto\;\exp\!\left[-\frac{\Phi(\theta)}{2\sigma^2\,\Phi_{\min}}\right]\pi(\theta),
\qquad \sigma = 0.05 ,
```

i.e. models whose misfit is a few percent above the best one found ($`\Phi_{\min}`$) remain
plausible. The prior $`\pi(\theta)`$ is

```math
\pi(\theta) \;\propto\; \prod_{i=1}^N \rho(x_i,y_i)\;\times\;
\mathbb 1\!\left[\;t_{i+1}-t_i \ge 2.5\,\mathrm s,\;\; t_i \ge T_i/2,\;\;
|\mathbf x_i-\mathbf x_1| \le v_{\max}\,(t_i-t_1),\;\; \theta\in\text{bounds}\right],
```

where $`\rho`$ is a smoothed map of the aftershocks (subevents are more likely where aftershocks
occurred), subevents are numbered in time order, and the causality condition (default
$`v_{\max}=1.5\,V_S`$) forbids a subevent that would require the rupture to travel faster than
plausible.

The posterior is sampled with **32 parallel Markov chains** (one per MPI rank) using
Metropolis–Hastings with single-parameter moves, adaptive joint moves along the learned posterior
covariance (Haario et al., 2001, tuned to 23% acceptance) and one-dimensional moves along its
principal directions. After burn-in, chains stuck far above the best solution are restarted from
good states ("hybrid" consolidation), and the posterior is built from the chains that converged.
Chain lengths grow with $`N`$ (1000·N burn-in + 1000·N samples by default, capped at 6000).

### 4.6 How many subevents?

More subevents always fit the data at least as well, so the inversion is repeated for
$`N = 1, 2, 3, \dots`$ and the answer is the **smallest $`N`$ whose best misfit is within 5% of the
best misfit over all $`N`$**. If that is the largest $`N`$ tried, the grid is extended automatically.

![Model selection for the four examples](docs/figures/model_selection.png)

*Misfit (relative to each event's best) against the number of subevents. Stars mark the selection.
Thanks to the penalties, the curves flatten or rise once extra subevents stop being supported
by the data.*

---

## 5. The automated pipeline

```mermaid
flowchart TD
    A["Waveforms downloaded from IRIS/EarthScope<br/>(teleseismic + regional)"] --> B["Step 1: instrument response removal,<br/>rotation, picks, 1 sample/s"]
    B --> C["Step 2a: screening, stage 1<br/>(amplitudes, noise, dead channels, drift)"]
    C --> D["Step 2b: Green's functions<br/>(mtel3 teleseismic, fk regional; CRUST1.0)"]
    D --> E["Step 2c: priors and configuration<br/>(aftershock map, catalog moment tensor, bounds, windows)"]
    E --> F["Step 3a: 1-subevent inversion"]
    F --> G["calibrate from the 1-subevent model:<br/>time bounds, windows, wave-type weights,<br/>screening stage 2 (flipped / mis-scaled stations)"]
    G --> G2["re-run N = 1 with the calibrated settings<br/>(so every N is scored on the same data)"]
    G2 --> H["Step 3b: inversions for N = 2, 3, ...<br/>(32 MCMC chains each)"]
    H --> I{"best N inside the<br/>5% band and below<br/>the largest N tried?"}
    I -- no --> J["add one more subevent"] --> H
    I -- yes --> K["Step 4: forward model, uncertainties,<br/>figures and RESULTS.md"]
```

What each automatic decision does, and why it matters:

| stage | what happens | why |
|-------|--------------|-----|
| **Integration** | teleseismic velocity records become displacement | the inversion fits displacement |
| **Screening, stage 1** | rejects records whose amplitude is >8× off the network median (robust median/MAD), whose signal-to-noise ratio is too low (5.5 teleseismic, 3 regional), whose channel is dead, or whose regional record is dominated by long-period drift | one broken sensor can bias a whole solution |
| **Green's functions** | depth grid centered on the catalog depth, wider for deep events (half-width from the Wells & Coppersmith, 1994, rupture length) | subevents of deep ruptures can be 50 km apart in depth |
| **Aftershock prior** | USGS aftershocks above a magnitude floor; background seismicity removed with a Poisson test against the year before; small clusters dropped; gaps between large clusters bridged with weak smooth corridors | keeps subevents where the fault actually slipped |
| **Windows** | long enough for the whole rupture plus, for deep events, the depth phases (pP, sP) | a 120-km-deep event's sS arrives ~60 s after S |
| **1-subevent calibration** | the 1-subevent solution sets the time range searched for N ≥ 2 and the window lengths | stops late subevents from fitting noise |
| **Wave-type weights** | P, SH and regional weights are iterated until their contributions to the misfit are about 1.2 : 1 : 1 | otherwise the nearest, largest records drown the rest |
| **Screening, stage 2** | stations that anti-correlate with the 1-subevent prediction (flipped polarity) or are ×4 off in amplitude (wrong gain) are removed | these errors only show up against a model |
| **Consistent scoring** | the 1-subevent inversion is repeated with the calibrated windows, weights and stations | the number of subevents is chosen by comparing misfits, which is only fair if every N is scored on the same data |

![Data processing](docs/figures/data_processing.png)

*Examples of what the processing does to real records from these events: drift removal, alignment
on the arrival, rejection of impossible amplitudes and of a sensor with reversed polarity, and
weighting of the three wave types.*

The pipeline is a set of shell steps that you can also run by hand from an `IRIS/` folder (this is
what `./submit` does):

```bash
bash ../programs/step2_prepare_inversions.sh        # screening, GFs, priors, configuration
bash ../programs/step3_do_inversions.sh             # 1-subevent calibration + N = 2..5
bash ../programs/step3_do_inversions_tempnsub.sh 6  # one more N, without the calibration hooks
python ../programs/prepare_SubMIT/plot_lcurve.py    # L-curve + selected N
bash ../programs/step4_fwd_tempnsub.sh 4            # forward model + figures for N = 4
```

---

## 6. Results of the examples

For each example, the main figure shows the subevents in map view (beachballs = moment tensors,
error bars = 95% intervals from the MCMC ensemble, ★ = epicenter), the moment-rate functions, the
total moment tensor and a depth section. Below it, a collapsible block holds the diagnostics:

* **Posterior histograms** — one row per subevent (E1, E2, …), one column per parameter (centroid
  time, duration, east and north location, depth). Narrow, single-peaked histograms mean a
  well-resolved parameter; broad or multi-peaked ones show the trade-offs the data allow. With two
  or more subevents, E1's horizontal position is fixed at the epicenter.
* **Waveform fits** — recorded (black) and predicted (red) seismograms for every station, labeled
  with the station name, epicentral distance and azimuth. The full plates (all stations, all
  pages) and the convergence of every Markov chain are in `examples/<name>/reference/` and are
  regenerated by `./submit example <name>`.

### 2024 Chile–Argentina border, Mw 7.4 (123 km deep)

![Chile](docs/figures/chile_subevents.png)

Four subevents within ~17 s. Horizontally they stay within 15 km of the epicenter, but their
depths step from 126 km to 170 km: the rupture grew mainly **downward** through the subducting
slab. This event (the "2024 Mw 7.4 Calama earthquake") was studied in detail with SubMIT by
Jia et al. (2025a), who found five subevents over ~20 s at 125–174 km depth. Here the 5-subevent
model fits only 1.4% worse than the 4-subevent one, and the 5% rule prefers the simpler model.

<details>
<summary><b>Uncertainties and waveform fits (Chile)</b></summary>

**Posterior distributions** of each subevent's parameters from the MCMC ensemble:

![Chile posterior histograms](docs/figures/chile_hist.png)

**Teleseismic P waves** (displacement; first of 2 pages):

![Chile P-wave fits](docs/figures/chile_fits_P.png)

**Teleseismic SH waves** (displacement; first of 2 pages):

![Chile SH-wave fits](docs/figures/chile_fits_SH.png)

**Regional three-component records** (first of 4 pages):

![Chile regional fits](docs/figures/chile_fits_rayl.png)

</details>

### 2025 Kamchatka megathrust, Mw 8.8

![Kamchatka](docs/figures/kamchatka_subevents.png)

Seven thrust subevents (Mw 8.1–8.5) over ~230 s. After a first subevent at the epicenter and a
second ~80 km to the north, the rupture ran **~450 km to the south-southwest** along the
trench. All seven mechanisms agree with the overall thrust mechanism (tensor similarity ≥ 0.94).

<details>
<summary><b>Uncertainties and waveform fits (Kamchatka)</b></summary>

**Posterior distributions** of each subevent's parameters from the MCMC ensemble:

![Kamchatka posterior histograms](docs/figures/kamchatka_hist.png)

**Teleseismic P waves** (displacement; first of 2 pages):

![Kamchatka P-wave fits](docs/figures/kamchatka_fits_P.png)

**Teleseismic SH waves** (displacement; first of 2 pages):

![Kamchatka SH-wave fits](docs/figures/kamchatka_fits_SH.png)

**Regional three-component records**; only two regional stations, with a small weight, so this event relies mainly on teleseismic data:

![Kamchatka regional fits](docs/figures/kamchatka_fits_rayl.png)

</details>

### 2026 near the coast of Venezuela, Mw 7.2

![Venezuela](docs/figures/venezuela_subevents.png)

Six subevents over ~110 s propagating **~170 km east-northeast** along the coast, combining
short, strong bursts with a longer, smoother subevent. The total moment (Mw 7.58) exceeds the
catalog Mw 7.2 — the one result in this set that we flag for further study.

<details>
<summary><b>Uncertainties and waveform fits (Venezuela)</b></summary>

**Posterior distributions** of each subevent's parameters from the MCMC ensemble:

![Venezuela posterior histograms](docs/figures/venezuela_hist.png)

**Teleseismic P waves** (displacement; first of 2 pages):

![Venezuela P-wave fits](docs/figures/venezuela_fits_P.png)

**Teleseismic SH waves** (displacement; first of 2 pages):

![Venezuela SH-wave fits](docs/figures/venezuela_fits_SH.png)

**Regional three-component records**:

![Venezuela regional fits](docs/figures/venezuela_fits_rayl.png)

</details>

### 2024 offshore Cape Mendocino, California, Mw 7.0

![California](docs/figures/california_subevents.png)

Three strike-slip subevents within ~21 s; the third, ~15 s after the origin, lies ~50 km to the
east-southeast. This matches the number of subevents of an earlier hand-tuned inversion.

<details>
<summary><b>Uncertainties and waveform fits (California)</b></summary>

**Posterior distributions** of each subevent's parameters from the MCMC ensemble:

![California posterior histograms](docs/figures/california_hist.png)

**Teleseismic P waves** (displacement; first of 2 pages):

![California P-wave fits](docs/figures/california_fits_P.png)

**Teleseismic SH waves** (displacement; first of 2 pages):

![California SH-wave fits](docs/figures/california_fits_SH.png)

**Regional three-component records** (first of 7 pages):

![California regional fits](docs/figures/california_fits_rayl.png)

</details>

---

## 7. Running your own earthquake

1. **Download waveforms** with [Wilber 3](https://ds.iris.edu/wilber3/) as SAC files with
   response (SACPZ) files: teleseismic broadband stations at 30°–90° (`BH?` channels), and, in a
   separate request, regional stations within ~10° if there are any. Keep Wilber's event-folder
   name (`YYYY-MM-DD-mwXX-region`, e.g. `2024-07-19-mww74-chile-argentina-border-region`): the
   magnitude is read from it.
2. **Put the downloads in a new folder** (short path, see [Troubleshooting](#10-troubleshooting)):

   ```text
   /data/myevent/
   ├── IRIS/       the teleseismic .tar file(s) from Wilber
   └── IRISloc/    the regional .tar file(s) from Wilber (optional)
   ```

3. **Run it:**

   ```bash
   ./submit prepare /data/myevent     # step 1: responses, rotation, picks, 1 sample/s, SNR selection
   ./submit run /data/myevent         # steps 2-4: everything else, including choosing N
   ```

   The aftershock catalog and the catalog moment tensor are downloaded from the USGS during
   step 2, so the machine needs internet access. Results appear in
   `/data/myevent/IRIS/figs_and_results/`. Add `--quick` to `run` first if you want a 10-minute
   check that everything works.

Settings you may want to change afterwards are in `Par.file` (windows, frequency bands,
weights, regularization) and `search_par.file` (search bounds, chain lengths, and
`neq_min`/`neq_max`, which set the number of subevents) inside each `IRIS/inv_<event>_<N>sub/`
folder; `programs/submit.conf` holds the sampling interval. After editing, re-run one
subevent count with `bash ../programs/step3_do_inversions_tempnsub.sh N` from `IRIS/`, then
`step4_fwd_tempnsub.sh N`.

---

## 8. Reading the outputs

`IRIS/figs_and_results/` collects everything for the selected model:

| file | content |
|------|---------|
| `RESULTS.md` | summary: selected N, misfits, subevent table with uncertainties, total moment |
| `subevents.png/.pdf` | the figure shown in Section 6 |
| `lcurve.pdf/.txt` | misfit against the number of subevents and the selection |
| `fits_P.pdf`, `fits_SH.pdf`, `fits_Pvel.pdf`, `fits_rayl.pdf` | data (black) vs prediction (red) at every station, with correlation coefficients |
| `histoplot_py.pdf` | posterior distributions of every subevent parameter |
| `misfit_evolution.pdf` | misfit along each Markov chain (convergence check) |

Per-N details are in `IRIS/inv_<event>_<N>sub/` (`best_model_hybrid.dat`, chains, logs) and
`IRIS/fwd_<event>_<N>sub/` (`fm.dat` = moment tensors in 10²⁷ dyne·cm as Mxx Mxy Mxz Myy Myz Mzz
with x north, y east, z down; `Input.model` = time, x east, y north, duration, vr, θ, depth).

---

## 9. Using an AI coding assistant

Everything above is driven by the `./submit` command, so an AI coding assistant that can run
shell commands can install, run and check SubMIT for you by following this README. Open the
repository in the assistant and ask, for example:

* *"Check whether my machine has everything SubMIT needs, and help me install what is missing."*
* *"Reproduce the California example and tell me whether it matches the published result."*
* *"Re-run the Chile inversions in frozen mode in the background and report when they finish."*
* *"Explain the L-curve and the selected model of the Kamchatka example."*
* *"Set up a new run for the earthquake whose Wilber downloads are in ~/data/myevent."*

---

## 10. Troubleshooting

| symptom | fix |
|---------|-----|
| SAC reports files it cannot read, or integrated records are unchanged | SAC cannot handle file names longer than ~120 characters: use a short work path (`--workdir /data/w`) |
| `COMPILE FAILED` | run `./submit check`; with GNU compilers try `./submit example … --gnu` |
| `mpirun` fails right away after building with GNU | an Intel `mpirun` was picked up; use `--gnu` (sets `SUBMIT_NO_INTEL=1`) |
| `ERROR: missing Green's function file` | a station has no Green's function; rerun the Green's-function stage or remove the station from `stations*.info` and update the `num_sta_*` counts in `Par.file` |
| Python errors about `numpy`/`obspy` | `conda activate submit`, or set `SUBMIT_PYTHON` to that environment's `python3` |
| need to stop a run | `sh programs/stop_step3.sh work/<name>/IRIS` (stops the step-3 orchestrator *and* its MPI jobs; killing only the shell leaves them running) |

---

## 11. Repository layout

```text
submit                  command-line driver (check, build, example, status, compare)
environment.yml         Python environment
examples/<name>/        four validated earthquakes
  inputs.tar.gz           processed waveforms (output of step 1)
  catalog/                archived aftershock catalog and moment tensor (exact reproduction)
  gf_inputs/              velocity models + grid of the published Green's functions
  frozen/<N>sub/          exact configuration + best model of every published inversion
  reference/              published figures, L-curve, moment tensors, best model
programs/               the package (paths inside are relative: keep this layout)
  step1_* ... step4_*     pipeline steps
  code_SubMIT/            MCMC inversion (finv.f90 + C forward/linear solver)
  fwd_SubMIT/             forward modeling and figures
  prepare_SubMIT/         screening, priors, weights, bounds, selection
  gf_SubMIT/              Green's-function drivers
  fk3.2/ mtel3/ crust1.0/ TauP-2.0/   bundled third-party codes and models
  hist/                   posterior histograms and convergence plots
tools/                  comparison and README-figure scripts
docs/figures/           figures used in this README
```

---

## 12. Credits, citation and licenses

SubMIT was developed by Zhe Jia. If you use it, please cite the method papers (Jia et al., 2020a,
2020b, 2022a) and Jia et al. (2025a), which used this code for the 2024 Mw 7.4 Calama earthquake
(the `chile` example).

### Method and applications of the multiple-subevent inversion

* Jia, Z., Shen, Z., Zhan, Z., Li, C., Peng, Z., & Gurnis, M. (2020a). The 2018 Fiji Mw 8.2 and
  7.9 deep earthquakes: One doublet in two slabs. *Earth and Planetary Science Letters*, 531,
  115997. https://doi.org/10.1016/j.epsl.2019.115997
* Jia, Z., Wang, X., & Zhan, Z. (2020b). Multifault models of the 2019 Ridgecrest sequence
  highlight complementary slip and fault junction instability. *Geophysical Research Letters*,
  47(17), e2020GL089802. https://doi.org/10.1029/2020GL089802
* Jia, Z., Zhan, Z., & Kanamori, H. (2022a). The 2021 South Sandwich Island Mw 8.2 earthquake: A
  slow event sandwiched between regular ruptures. *Geophysical Research Letters*, 49(3),
  e2021GL097104. https://doi.org/10.1029/2021GL097104
* Jia, Z., Mao, W., Flores, M. C., Barra, S., Ruiz, S., Potin, B., Becker, T. W., Moreno, M.,
  Baez, J. C., Ceroni, D., & Cabrera, L. (2025a). Deep intra-slab rupture and mechanism transition
  of the 2024 Mw 7.4 Calama earthquake. *Nature Communications*, 16, 8140.
  https://doi.org/10.1038/s41467-025-63480-5
* Kutschera, F., Jia, Z., Oryan, B., Wong, J. W. C., Fan, W., & Gabriel, A.-A. (2024). The
  multi-segment complexity of the 2024 Mw 7.5 Noto Peninsula earthquake governs tsunami
  generation. *Geophysical Research Letters*, 51(21), e2024GL109790.
  https://doi.org/10.1029/2024GL109790

### Related studies

* Jia, Z., Zhan, Z., & Helmberger, D. (2022b). Bayesian differential moment tensor inversion:
  Theory and application to the North Korea nuclear tests. *Geophysical Journal International*,
  229(3), 2034–2046. https://doi.org/10.1093/gji/ggac053
* Jia, Z., Jin, Z., Marchandon, M., Ulrich, T., Gabriel, A.-A., Fan, W., Shearer, P., Zou, X.,
  Rekoske, J., Bulut, F., Garagon, A., & Fialko, Y. (2023). The complex dynamics of the 2023
  Kahramanmaraş, Turkey, Mw 7.8–7.7 earthquake doublet. *Science*, 381(6661), 985–990.
  https://doi.org/10.1126/science.adi0685
* Jia, Z., Fan, W., Mao, W., Shearer, P. M., & May, D. A. (2025b). Dual mechanism transition
  controls rupture development of large deep earthquakes. *AGU Advances*, 6(3), e2025AV001701.
  https://doi.org/10.1029/2025AV001701 (subevent models of 40 deep earthquakes with a related,
  grid-search subevent method)
* Li, J., & Jia, Z. (2026). The 2025 Mw 8.8 Kamchatka megathrust: A rapid recurrence with complex
  heterogeneous rupture. *Geophysical Research Letters*, 53(9), e2026GL121923.
  https://doi.org/10.1029/2026GL121923 (finite-fault model of the `kamchatka` example event)

### Other references used in this README

* Haario, H., Saksman, E., & Tamminen, J. (2001). An adaptive Metropolis algorithm. *Bernoulli*,
  7(2), 223–242. https://doi.org/10.2307/3318737
* Wells, D. L., & Coppersmith, K. J. (1994). New empirical relationships among magnitude, rupture
  length, rupture width, rupture area, and surface displacement. *Bulletin of the Seismological
  Society of America*, 84(4), 974–1002.

### Bundled third-party software and models (with their own terms)

* **fk** — Lupei Zhu (Saint Louis University), frequency–wavenumber synthetics: Zhu, L., &
  Rivera, L. A. (2002). A note on the dynamic and static displacements from a point source in
  multilayered media. *Geophysical Journal International*, 148(3), 619–627.
  https://doi.org/10.1046/j.1365-246X.2002.01610.x. Copyright notice in `programs/fk3.2/README`
  (redistribution permitted with the notice; the author asks to be informed of redistribution).
* **TauP 2.0** — Crotwell, H. P., Owens, T. J., & Ritsema, J. (1999). The TauP Toolkit: Flexible
  seismic travel-time and ray-path utilities. *Seismological Research Letters*, 70(2), 154–160.
  https://doi.org/10.1785/gssrl.70.2.154. GPL-3.0 (`programs/TauP-2.0/gpl-3.0.txt`, source
  included); the bundled copy omits the documentation and the optional Jython console.
* **CRUST1.0** — Laske, G., Masters, G., Ma, Z., & Pasyanos, M. (2013). Update on CRUST1.0 — A
  1-degree global model of Earth's crust. *Geophysical Research Abstracts*, 15, EGU2013-2658.
* **r8lib** — John Burkardt (LGPL), linear algebra.
* SAC is required but **not** included; obtain it from IRIS/EarthScope.

Earthquake data are from IRIS/EarthScope and the networks that operate the stations; the aftershock catalogs
and catalog moment tensors are from the USGS ComCat and the Global CMT project.
