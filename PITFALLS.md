# Failure atlas

How this pipeline failed silently, how each failure was found, and what now prevents it from
coming back.

Every entry is a defect of SubMIT v2 (GitHub commit 4ed2f61), of a fix proposed for one (A4), or of
the v3 rewrite itself (marked where it occurs), caught in November 2026: by a cohort bug report on
the 2026 Mw 7.4 Colombia earthquake, by clean-environment runs (fresh clone, fresh conda
environment) on three earthquakes the code had never seen (2026-10-09 Panama Mw 7.7, 2026-03-24
Tonga Mw 7.5 at 234 km depth, 2025-12-08 Aomori Mw 7.6), and by an independent code review (F). Nothing below is a generic best-practice list. Where nobody is known
to have believed a wrong explanation, "The false story" says so and names what made the defect easy
to miss.

**The organising observation: only two of the fourteen incidents stopped the pipeline (E1, E2). In
ten, the pipeline ran to the end and drew plausible figures, and one of those (A1) was explained
away twice instead of being found: its unrotated SH records were removed as "polarity flips" and
"gain errors", and the poor SH fit they caused was read as a smooth background rupture.**

Companion document: [README.md](README.md) explains the method and the pipeline.

---

## Index

| id | incident | category | cost if uncaught |
|---|---|---|---|
| A1 | horizontals of unequal length left unrotated | data integrity | 26 of 56 SH records of the published Venezuela model were north or east components |
| A2 | broadband-only regional stations dropped, co-located sensors counted twice | data integrity | Panama kept 4 of 19 regional stations; IU.MAJO counted 3 times |
| A3 | regional signal-to-noise measured against the P coda | data integrity | every regional station of a deep event rejected |
| A4 | clip rules carried over from one network | data integrity | 10 of 22 good Venezuela stations, and IU.PET, flagged as clipped |
| A5 | a record at another sampling interval in the bundled inputs | data integrity | a 50-samples/s record read as if sampled every second |
| B1 | teleseismic timing referenced to the depth of E1 | forward model | all teleseismic arrivals shifted by (z1 - evdp) x vertical slowness |
| B2 | sP shifted twice | forward model | sP of a subevent 40 km from the reference about 3 s early |
| B3 | 0/0 correlation and NaN misfits under -Ofast | forward model | chains frozen at NaN, memory read outside the trace |
| C1 | prior corridors re-attached to zeroed speckle | prior | Kamchatka support 40% instead of 16%; two published subevents on corridors |
| D1 | one lucky chain decides the number of subevents | sampling & selection | selected N not reproducible; ensembles that mix modes |
| D2 | the 1-subevent consolidation stored a misfit of another model | sampling & selection | calibration runs whose recorded best did not belong to its model |
| E1 | fresh clone has no greenFuncDir_disp; the clean-up then deletes every Green's function | systems | run aborts after the longest stage |
| E2 | fixed array limits | systems | step 3 aborted for 108 and 232 teleseismic stations |
| E3 | silent crashes hidden by unchecked exit codes | systems | Qt backend abort in caldens.py, orphaned MPI jobs, clobbered arguments |
| F1-F8 | defects found by an independent review of the v3 code, each checked before it was fixed | review | see section F |

Each incident uses the same template: **Symptom → The false story → Root cause → Fix →
Permanent regression test → Invalidated artifacts → Generalisable lesson.**

---

# A · Data integrity

## A1 — Horizontals of unequal length left unrotated

**Symptom.** In the processed teleseismic data of the examples, some `.t` (transverse) records carry
`kcmpnm = BHN` and `cmpaz = 0`, and their `.r` partner `kcmpnm = BHE`. CN.DRLN (Chile): BHE 1462 and
BHN 1475 samples. The step-1 log has one `ERROR 2010: Number of points in pair of files are not
equal` per such pair, and the bundled inputs hold such stations for every example: 31 for Venezuela
(most of the GE network), 5 for Chile, 4 each for California and Kamchatka.

**The false story.** The README read the Venezuela result as a rupture with a smooth, long-lasting
background component, because SH fitted poorly (mean correlation 0.26 against 0.59 for P). The
stage-2 screen listed CN.DRLN and GE.SNAA as polarity flips and CN.GAC as a gain error.

**Root cause.** SAC's `rot to gcp` refuses two traces of different length, and the next command of
the same SAC session, `w rtr/....r rtr/....t`, writes the unrotated east and north traces under the
names of the rotated ones. Wilber delivers E and N with a different number of samples whenever one
component has a short gap. The selected published models used these records as SH data: Venezuela
26 of 56 stations (published correlation 0.17 against 0.34 for the correctly rotated ones),
California 3 of 67 (0.06 against 0.74), Kamchatka 2 of 45 (0.23 against 0.63), Chile 2 of 51
(0.75 against 0.91).

**Fix.** `programs/MyExtractSeed.sh` cuts a horizontal pair whose lengths differ to the common
window before rotating. Pairs of equal length are untouched. Every other teleseismic difference in
the four examples follows from the rotation: a few stations cross the signal-to-noise threshold or
change the duplicate-sensor choice with their rotated transverse component (Venezuela GE.KBS in and
GE.MATE out; Kamchatka G.CAN in, G.ECH.00 instead of .10).

**Permanent regression test.** Step 1 on the original Wilber downloads of the four examples,
compared with `examples/*/inputs.tar.gz` by `tools/compare_step1_outputs.py`, which also lists the
radial and transverse records that are still unrotated: 62, 10, 8 and 8 records in the bundled
Venezuela, Chile, California and Kamchatka inputs, none in the v3 output.

**Invalidated artifacts.** The SH data of the bundled example inputs for the stations above, and
every published result that used them, above all the Venezuela interpretation. The examples have to
be reprocessed and re-inverted.

**Generalisable lesson.** **A tool that fails inside a batch session does not stop the session: check
the output for the property the step was supposed to create (here, cmpaz of the transverse
component), not only that a file appeared.**

## A2 — Broadband-only regional stations dropped, co-located sensors counted twice

**Symptom.** Panama (HN?, HH? and BH? downloaded from 19 stations): 4 regional stations survived
step 1, all with strong-motion records. Colombia (cohort report): the nearest kept station was at
496 km while HEL (169 km) and RUS (369 km) had broadband records. Aomori: 11 records from 7 stations
entered the inversion (IU.MAJO.00, .10 and .60; II.ERM and IU.YSS twice).

**The false story.** "Strong-motion records are the clean ones; broadband sensors clip near a large
earthquake", which is why v2 used broadband channels only when no strong-motion channel existed.

**Root cause.** The regional step processed one channel family per event (HN?, else BH?; HH? never)
and never chose between sensors of one site.

**Fix.** `programs/MyExtractSeed_loc.sh` processes every HN?, HH? and BH? family;
`programs/select/choose_regional_instruments.py` keeps one record per station: the broadband record
unless it is clipped or disagrees with a co-located strong-motion record at 0.02-0.1 Hz (cc >= 0.9,
amplitude ratio 0.8-1.25; thresholds from the cohort report). Panama now keeps 17 of 19 stations and
Aomori one record for each of its 7; CU.BCIP's broadband record disagreed with its accelerometer
(cc 0.06, amplitude x0.12) and was replaced by it.

**Permanent regression test.** `tests/test_regressions.py::test_clip_rule_*`, and the step-1 example
comparison.

**Invalidated artifacts.** Regional data sets of downloads with several channel families.

**Generalisable lesson.** **A rule that selects data by type ("use strong motion") silently deletes
every station that lacks that type; select by station, then by quality.**

## A3 — Regional signal-to-noise measured against the P coda

**Symptom.** Tonga (234 km deep): 1 of 5 regional stations kept (v2), 2 of 5 after the fix of A2.
The log shows vertical components at 10-12 dB and horizontals at -1.7 to 1.9 dB.

**The false story.** None documented. Low horizontal signal-to-noise at island stations looks like
ordinary site noise, so the rejections drew no attention.

**Root cause.** `programs/select/selectSNRloc.pl` compared the 50 s after the S pick with the 50 s
before it. For a deep source at 5-7 degrees S-P is about 63 s, so the "noise" window was the P coda.

**Fix.** Every component is measured from the P pick (header t2 = t1 of the vertical, written in
`programs/MyExtractSeed_loc.sh`): signal after P against pre-event noise. Tonga keeps all 5 regional
stations.

**Permanent regression test.** Step-1 example comparison: California and Chile each gain one station,
Venezuela is unchanged.

**Invalidated artifacts.** Regional data sets of intermediate-depth and deep events.

**Generalisable lesson.** A noise window must lie before every arrival of the event, not before the
phase being scored.

## A4 — Clip rules carried over from one network

**Symptom.** While fixing A2: 10 of the 22 strong-motion stations of the Venezuela example were
flagged as clipped, and IU.PET (Kamchatka) was dropped as clipped.

**The false story.** The proposed rules assumed that any repeated sample near the peak, or any peak
near 2^23 counts, means a saturated sensor.

**Root cause.** Two rules of the cohort's Colombian test did not transfer. "Three equal samples above
half the peak" fires on smooth, oversampled 200-sps accelerograms (CM.ARGC HNE: six equal samples at
0.6% of full scale). "Peak above 95% of 2^23 counts" assumes 24-bit digitizers; IU.PET peaks at
4.8-6.6 x 2^23 without any flat top.

**Fix.** Clipped = a plateau of at least 3 samples AT the peak value with the peak above 10% of 2^23.
Sensor saturation is left to the agreement test of A2.

**Permanent regression test.** `tests/test_regressions.py::test_clip_rule_ignores_repeated_samples_of_a_smooth_strong_motion_peak`,
`::test_clip_rule_flags_a_flat_topped_broadband_record`.

**Invalidated artifacts.** None released.

**Generalisable lesson.** A quality threshold measured on one network is a hypothesis for the next
one; test it on data that is known to be good.

## A5 — A record at another sampling interval in the bundled inputs

**Symptom.** The first v3 run of `./submit example california --mode full` stopped in step 3:
`sampling interval differs from Par.file deltat in data/G.DZM.00.z`.

**The false story.** None; v2 never complained.

**Root cause.** v2's step 1 left a few records at their original sampling rate, and the bundled
California inputs hold G.DZM.00 at 0.02 s. No published configuration fits it, but step 2 of a
full-mode run lists every record, and v2's forward code then read the 50-samples/s trace as if it
were sampled every second. The v3 forward code checks the sampling interval of every record and
stops instead.

**Fix.** `programs/prepare_SubMIT/screen_station_amplitudes.py` sets aside, with a logged reason,
every record whose sampling interval is not DELTAT of `programs/submit.conf`, before the amplitude
screens; the v3 step 1 resamples such records to 1 sample/s in the first place.

**Permanent regression test.** `tests/test_regressions.py::test_records_at_another_sampling_interval_are_set_aside`,
and `./submit example california --mode full --quick` (9 min 55 s on 32 ranks, exit 0).

**Invalidated artifacts.** v2 full-mode runs on the bundled California inputs, if G.DZM.00 passed
their screens.

**Generalisable lesson.** A reader that assumes a property of its input (here, the sampling
interval) must check it; the earliest stage that can check it should set the record aside.

---

# B · Forward model

## B1 — Teleseismic timing referenced to the depth of E1

**Symptom.** Colombia (cohort report): E1 at 79.8 km for a 110-km hypocentre, and every subevent
time 3.2 s later than the hand-picked onset and the back-projection.

**The false story.** None documented. The shift is invisible inside one run: every subevent moves
together, and the cohort found it only because it disagreed with a hand-picked onset and a
back-projection.

**Root cause.** `sub_forward.c` shifted every teleseismic Green's function by
`-(z_i - z_1) * eta`, with z_1 the depth of E1. The data are aligned on the arrival from the catalogue
hypocentre (SAC t1 at evdp), and mtel3 puts the direct P of every grid depth at the same time, so the
reference must be evdp. With z_1 free, P moved by (z_1 - evdp) eta_p and SH by (z_1 - evdp) eta_s,
two different amounts, against regional Green's functions whose timing is absolute.

**Fix.** `programs/code_SubMIT/forward_model.c` uses the evdp of the data headers as reference depth.
On the published Chile 4-subevent model the data misfit falls from 0.8537 to 0.8452 with this fix
alone.

**Permanent regression test.** The equivalence run of the rewritten forward code (README section 9).

**Invalidated artifacts.** Every published misfit and best model was computed with the old reference;
`./submit example ... --mode figures` no longer reproduces them to the digit.

**Generalisable lesson.** **A synthetic and the data must be referenced to the same pick; when a
parameter (here the depth of E1) also moves the reference, the fit absorbs the error into time.**

## B2 — sP shifted twice

**Symptom.** In the Chile Green's functions the sP peak moves 0.318 s per km of source depth, close to
the expected eta_p + eta_s = 0.324 s/km of the P ray parameter; the forward code moved the S leg by a
further (z_i - z_ref)(eta_s - eta_p).

**The false story.** None was told; the depth phases were simply mistimed.

**Root cause.** The P synthetic convolved the P leg (P, pP) and the S leg (sP) of each Green's
function with two pulses, the S-leg pulse shifted with the vertical S slowness of the direct-S ray
parameter. mtel3 already delays sP correctly for each grid depth, so the whole Green's function
needs one shift. For a subevent 40 km from the reference depth, sP was about 3 s early.

**Fix.** One pulse for the summed P and S legs (`forward_model.c`), which also halves the
convolutions. Published Chile 4-subevent data misfit: 0.8436 with this fix alone, 0.8244 with B1 and
B2 together (-3.4%).

**Permanent regression test.** As B1.

**Invalidated artifacts.** As B1; most for deep events.

**Generalisable lesson.** When a table (here Green's functions per depth) already contains an effect,
the code that uses it must not model the effect again.

## B3 — 0/0 correlation and NaN misfits under -Ofast

**Symptom.** Cohort report: 20 of 96 one-subevent chains ended at NaN.

**The false story.** None documented; the cohort report identified the NaN guard directly. A
frozen chain looks like an ordinary non-converged one in the misfit plots.

**Root cause.** `-Ofast` assumes finite numbers, so NaN comparisons may be compiled away and a NaN
proposal accepted. Separately, the cross-correlation of an all-zero synthetic (a zero component
weight, a pulse outside the window) divided 0 by 0 and left the lag uninitialized, which then
indexed memory outside the trace.

**Fix.** `forward_model.c` returns lag 0 and correlation 0 for an all-zero trace and tests the misfit
on its IEEE bits; `finv.f90` tests it again and rejects a proposal with zero prior density
explicitly; the Makefile adds `-fno-finite-math-only`.

**Permanent regression test.** None automatic; the NaN guard is two lines in each file.

**Invalidated artifacts.** Runs with NaN chains.

**Generalisable lesson.** Fast-math flags remove exactly the checks that catch degenerate input.

---

# C · Prior

## C1 — Prior corridors re-attached to zeroed speckle

**Symptom.** Kamchatka catalogue: caldens.py zeroed 6 minor density components and then drew 141
bridges; the prior covered 40% of the grid instead of 16%. Colombia (cohort report): 30% grew to 52%,
and one subevent of the 5-subevent model sat on a corridor without seismicity.

**The false story.** None was told about the corridors, which no figure showed. The published
Kamchatka path, including the second subevent ~80 km north, was presented as resolved by the data.

**Root cause.** The minimum spanning tree spanned every aftershock, and the bridge test skipped only
edges inside one surviving component, so every edge to a zeroed point (label 0) became a corridor.

**Fix.** `programs/prepare_SubMIT/caldens.py` builds the tree over aftershocks on surviving components
only. The cohort's fix (skip every edge touching label 0) would also cut the wanted bridge of a
doublet whose tree path passes through sparse in-between aftershocks.

**Permanent regression test.** `tests/test_regressions.py::test_prior_does_not_bridge_zeroed_minor_components`,
`::test_prior_keeps_the_bridge_between_doublet_clusters`.

**Invalidated artifacts.** The Kamchatka published model: E2 (-35, 77 km), the "second subevent ~80 km
to the north", and E6 lie where the corrected prior is zero (0.24 and 0.21 before).

**Generalisable lesson.** **When a rule removes points and a later rule connects points, the
connection must run over the points that survived.**

---

# D · Sampling and selection

## D1 — One lucky chain decides the number of subevents

**Symptom.** Panama, v2 code, 16 chains: N = 2 selected with a best misfit of 1.5303 reached by ONE
chain; the median chain's best was 16.6% worse, and a bootstrap over the chains selects N = 2 in 67%
of resamples and N = 3 in 33%. The published runs: typically 2-6 of 32 chains within 2% of the best
and split R-hat up to 11 (California E3 x). Colombia (user runs): 5, 4, 3, 2 and 4 subevents under
small configuration changes.

**The false story.** None documented. Because the Colombia runs changed a setting between runs,
each change of the selected N could be read as a response to that setting rather than to the luck of
the search.

**Root cause.** Independent chains of 2000 x N iterations each settle in a local mode of a very
peaked posterior (a model 1% worse has e^-2 of the density), so the best misfit of each N is the best
of a few random modes, with a spread of ~5% between seeds, the size of the selection band. The v2
consolidation copied a random state within 10% of the best at the end of the burn-in, with equal
weight, which kept separate modes in one ensemble.

**Fix.** `programs/code_SubMIT/finv.f90`: half of the chains start from the best (N-1) model plus
one subevent drawn from the prior; 5% of the proposals relocate one subevent to a position drawn from
the prior; the population is resampled with weights exp(-dphi / (2 sigma^2 phi_ref)) at 25, 50 and
75% of a burn-in of 3000 x N (islands of 8 chains, the last round global). On the California
3-subevent configuration, 10-16 of 16 chains then ended within 2% of the best in seven of eight
runs. The best misfit still differed between runs (1.840 and 1.844 for two seeds of the final
sampler, 1.761, 1.840 and 1.847 for an earlier version), and after the last round the chains are
copies of one island, so their agreement says nothing about that spread. finv therefore writes each island's best before the last round
(`island_best_misfits.txt`), and `programs/prepare_SubMIT/plot_lcurve.py` bootstraps those
independent replicates (the chains, for v2 runs), reports the share of resamples that select each
N and warns below 80%. `programs/hist/posterior_summary.py` writes the chains within 1/2/5% of the
best and the split R-hat of every parameter (`convergence.txt`, RESULTS.md). The first version of
the warm start had its own defect, which the Panama test log showed ("warm" on 8 of the 16 chains of
the N = 1 run): without `warm_start.dat` the warm model has zero subevents, which equals N - 1 for
N = 1, so half of the 1-subevent chains started at the epicentre instead of from the prior.

**Permanent regression test.** `tests/test_regressions.py::test_lcurve_flags_a_selection_that_rests_on_one_chain`,
`::test_lcurve_uses_islands_when_resampled_chains_agree`.

**Invalidated artifacts.** Selected N of runs without a convergence report.

**Generalisable lesson.** **A model-selection rule applied to "the best of the chains" is only as
reliable as the search: report how often the rule picks the same answer across independent
searches, and make sure the searches you resample are still independent.**

## D2 — The 1-subevent consolidation stored a misfit of another model

**Symptom.** Calama 1-subevent run: split R-hat of E1's y 4.7, for a five-parameter problem.

**The false story.** None documented. A high R-hat for one parameter reads as a poorly resolved
parameter, not as a bookkeeping error.

**Root cause.** The adaptive moves and the burn-in consolidation covered the active parameters of the
adaptive map, which excluded x and y of E1 for every N, although a single subevent may move. The
consolidation copied a peer's depth, duration, time and misfit, and kept the chain's own x and y: the
recorded misfit did not belong to the recorded model. The single-parameter moves also chose an
inactive slot for half of the x/y moves of one subevent, and the first and last of three or more
subevents half as often as the others.

**Fix.** `finv.f90`: x and y of E1 are adaptive parameters when it is the only subevent; the subevent
of a move is uniform; the population state copied between chains is the complete state.

**Permanent regression test.** None automatic.

**Invalidated artifacts.** 1-subevent calibration runs (they set the windows of N >= 2).

**Generalisable lesson.** A copy of "the state" must copy all of it, including what the bookkeeping
structure was not built for.

---

# E · Systems

## E1 — Fresh clone has no greenFuncDir_disp; the clean-up then deletes every Green's function

**Symptom.** Fresh clone: `cp -r greenFuncDir/vmodel_* greenFuncDir_disp/` fails, and step 2 stops
with "no teleseismic Green's functions" after the longest stage.

**The false story.** None documented. The error names the missing Green's functions, not the
missing folder whose absence led to their deletion.

**Root cause.** git does not track empty folders. After the failed copy, `dotel.sh` went on and
deleted the computed teleseismic and regional Green's functions (`rm loc/*.grn.*`,
`rm -rf greenFuncDir/vmodel_*`).

**Fix.** `programs/gf_SubMIT/dotel.sh` creates both folders and stops on a failed copy before any
deletion.

**Permanent regression test.** The clean-clone run of `./submit build` and step 2 (README section 3).

**Invalidated artifacts.** None.

**Generalisable lesson.** A clean-up must run only after the step it cleans up after has been
verified.

## E2 — Fixed array limits

**Symptom.** Tonga (108 teleseismic stations) and Aomori (232) aborted in step 3: `num_sta_P/SH
(216/98) exceeds maxstanum=200`.

**The false story.** None; the run simply stopped after two hours of preparation.

**Root cause.** The inversion arrays were sized at compile time for 200 P traces (displacement plus
velocity, i.e. 100 stations), while the README advised downloading 100-150 stations and the screens
kept more than 100 good ones.

**Fix.** `forward_model.c` sizes every array from the actual counts, and
`screen_station_amplitudes.py` thins a network denser than 100 stations by removing the most
redundant station in azimuth and distance first.

**Permanent regression test.** `tests/test_regressions.py::test_dense_teleseismic_networks_are_thinned_to_the_limit`.

**Invalidated artifacts.** None.

**Generalisable lesson.** A limit that the documentation invites users to exceed is a bug.

## E3 — Silent crashes hidden by unchecked exit codes

**Symptom.** In a fresh conda environment with `DISPLAY` set, caldens.py aborted ("Could not find the
Qt platform plugin xcb") after writing its outputs, and v2's step 2, which ran its stages through
`awk ... | sh` without exit codes, went on. Killing step 3 orphaned the awk that started the mpirun
jobs. Sourcing Intel's setvars.sh rewrote the positional parameters of the calling script. On the
test machine, the OpenMPI launcher hung even for `hostname`.

**The false story.** None; each looked like a successful run or a transient hang.

**Root cause.** Exit codes and process trees that no stage checked.

**Fix.** Every plotting script selects the Agg backend; steps 2-4 are plain loops that check each
stage; arguments are captured before sourcing compiler environments
(`programs/compiler_environment.sh`); `./submit check` starts one process through mpirun with a hard
timeout.

A fourth hazard showed up in the test runs themselves: bash reads a script while it runs it. When I
copied a corrected step-3 driver into the test installation seconds after one run had started it,
that run finished all five inversions and then resumed reading at its old byte offset in the new
file ("syntax error: unexpected end of file"), so it stopped before the L-curve. Step 3 is now one
`{ ...; exit; }` block, which bash parses completely before running it, so editing the script
during a multi-hour run, which the README invites, no longer corrupts the run.

**Permanent regression test.** `./submit check`, and the clean-environment runs. The first
clean run of `./submit example california` with the v3 code caught a defect of the rewrite itself:
steps 3 and 4 run under `set -u`, and `compiler_environment.sh` read the unset `INTEL_SETVARS`
(and sources Intel scripts that read unset variables), which stopped step 4 before any figure.

**Invalidated artifacts.** None.

**Generalisable lesson.** **A pipeline that does not check exit codes reports the last stage that
printed something, not the last stage that worked.**

---

# F · Found by review

During the clean-environment runs, an independent review of the new and rewritten code (a reviewer
that was given the code, not my conclusions) reported eleven defects. All eleven held up, and they
are grouped here into eight entries. Each was checked before it was fixed: reproduced where a test
could show it (F3-F7), and by following the code path otherwise (F1, F2, F8 and the mpirun case of
F7, which depends on the MPI launcher). They are listed more briefly than A-E because only F4 changed
a test run (Aomori); the provenance of each fix is in its code comment.

| id | defect | introduced | reproduction or measurement | fix |
|---|---|---|---|---|
| F1 | proposals redrawn until every hard rule held, so the chains sampled pi(x) Z(x), Z = share of feasible proposals from x: states near the bounds and the time-gap and causality rules were under-weighted | v2 | follows from the proposal density q 1_F / Z(x); 23-27% of the proposals break a rule on the California 3-subevent configuration | `finv.f90`: an infeasible random-walk proposal is a rejection (no forward run); relocations keep their redraw, which is exact for an independence proposal. Benchmark: best 1.840 and 1.844 (seeds 0, 500) against 1.761-1.847 before, 16% faster |
| F2 | the density at the epicentre entered the prior although E1 is pinned there for N >= 2; where the aftershock map is zero at the epicentre every proposal was rejected and all chains froze | v2 | the examples have density 0.41-1.0 at the epicentre, so none was hit | `finv.f90`: the pinned E1 is left out of the prior product (its factor cancelled in every ratio anyway) |
| F3 | prior values written at the lower-left corner of their 4-km bins, while finv reads them as cell centres: the whole prior shifted 2 km south-west; the search box ended at the lower edge of the last nonzero bin; an all-zero map was written as NaN | v2 | Kamchatka catalogue: identical values, coordinates +2.00 km in x and y after the fix | `caldens.py`: centres, upper edge, guard before the division |
| F4 | wave-type weights balanced on the step-2 windows, not on the shortened windows of the N >= 2 runs | v2 | Aomori: P:SH:regional 2.67 : 1 : 0.054 as run, 2.44 : 1 : 0.045 on the 55-s windows actually used | `calibrate_from_1sub.sh`: the forward template takes the N = 2 Par.file |
| F5 | regional 1/2 horizontals of unequal length written unrotated, the regional form of A1 | v2 | II.MSVF.00.BH2 (Tonga) shortened by 30 samples: SAC ERROR 2010 before the fix, rotated (cmpaz 0 and 90) after it, the other 13 records byte-identical | `MyExtractSeed_loc.sh`: cut to the common window, then check that every rotated north record has cmpaz 0 |
| F6 | Green's-function key int(dist) against the saclst print that names the files in dotel.sh | v3 rewrite | 3 of 309 test distances mismatched (float32 6202.9995: key 6202, file 6203.grn.*), all matched after the fix | `prepare_stationinfo.py`, `test_green_function_key_follows_the_saclst_print_used_by_dotel` |
| F7 | lost exit codes: the calibration returned 1 for a filler regional set (<= 2 stations) whenever the screen removed a station, which the new step-3 check turned into an abort; failed inversions were skipped and the L-curve built from the N that worked; two step-2 stages unchecked; mpirun forwarding the event list to rank 0 | v3 rewrite and v2 | filler case run on a copy of the Aomori folders with two regional stations: exit 0 after the fix | explicit statuses, step 3 stops on a failed N, `mpirun ... < /dev/null` |
| F8 | with 8 or fewer ranks finv has one island, and plot_lcurve.py then bootstrapped the chains, which descend from that one population and always agree | v3 rewrite | | one island counts as one replicate, so the reproducibility is reported as not measured; `test_lcurve_does_not_trust_the_chains_of_a_single_island` |

**Generalisable lesson.** **A review by someone who has not seen your reasoning finds the defects your
reasoning hides: give the reviewer the code, ask for concrete failure scenarios, and check each one
before fixing it.**
