"""Regression tests for the defects fixed in SubMIT v3 (PITFALLS.md); each test names what it guards.

They need only the conda environment (numpy, scipy, obspy, matplotlib) - no SAC, MPI or Green's
functions - and run in a few seconds on the catalogues bundled with the examples:

    python tests/test_regressions.py        (or: python -m pytest tests/test_regressions.py)
"""
import glob
import os
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np

repository_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
programs = os.path.join(repository_root, "programs")


def run_in(folder, command):
    """run a pipeline script inside folder; return its stdout (fail loudly with its output)"""
    result = subprocess.run(command, cwd=folder, capture_output=True, text=True)
    assert result.returncode == 0, f"{' '.join(command)} failed:\n{result.stdout}\n{result.stderr}"
    return result.stdout


# ---------------------------------------------------------------------------
# aftershock prior (PITFALLS C1)
# ---------------------------------------------------------------------------
def density_support_after_caldens(example_name):
    """support fraction of seisdens.dat and the caldens log for one bundled catalogue"""
    catalogue = os.path.join(repository_root, "examples", example_name, "catalog")
    longitude, latitude = open(os.path.join(catalogue, "mainshock.dat")).read().split()[:2]
    with tempfile.TemporaryDirectory() as folder:
        for name in ("lola.dat", "lola_bg.dat"):
            shutil.copy(os.path.join(catalogue, name), folder)
        log = run_in(folder, [sys.executable, os.path.join(programs, "prepare_SubMIT", "caldens.py"), longitude, latitude])
        density = np.loadtxt(os.path.join(folder, "seisdens.dat"))
    return float(np.mean(density[:, 2] > 0)), log


def test_prior_does_not_bridge_zeroed_minor_components():
    # v2 zeroed 6 minor components of the Kamchatka catalogue and re-attached them with 141 corridors
    support, log = density_support_after_caldens("kamchatka")
    assert "inter-cluster bridge" not in log, log
    assert support < 0.25, f"prior support {support:.2f} of the grid (v2 with the spurious corridors: 0.40)"


def test_prior_keeps_the_bridge_between_doublet_clusters():
    # the 2026 Venezuela doublet needs one corridor between its two rupture zones
    support, log = density_support_after_caldens("venezuela")
    assert "1 inter-cluster bridge" in log, log


# ---------------------------------------------------------------------------
# number of subevents (PITFALLS D1, D3)
# ---------------------------------------------------------------------------
def write_inversion_folder(root, subevent_count, chain_best_misfits):
    folder = os.path.join(root, f"inv_test_{subevent_count}sub")
    os.makedirs(folder)
    for rank, misfit in enumerate(chain_best_misfits):
        with open(os.path.join(folder, f"{rank}best.dat"), "w") as best_file:
            best_file.write(f"{subevent_count:2d} <neq " + "".join(f"{misfit:10.4f}" + "    0.0000" * 7 for _ in range(subevent_count)) + "\n")


def write_island_bests(root, subevent_count, island_best_misfits):
    folder = os.path.join(root, f"inv_test_{subevent_count}sub")
    with open(os.path.join(folder, "island_best_misfits.txt"), "w") as island_file:
        for island, misfit in enumerate(island_best_misfits, start=1):
            island_file.write(f"{island:4d}{8:4d}{misfit:12.5f}\n")


def lcurve_report_of(best_misfit_per_count, island_bests_per_count=None):
    """output of plot_lcurve.py for an L-curve given as {subevent count: best misfit}, every chain at the best"""
    with tempfile.TemporaryDirectory() as folder:
        for subevent_count, best_misfit in best_misfit_per_count.items():
            write_inversion_folder(folder, subevent_count, [best_misfit] * 4)
        for subevent_count, island_bests in (island_bests_per_count or {}).items():
            write_island_bests(folder, subevent_count, island_bests)
        return run_in(folder, [sys.executable, os.path.join(programs, "prepare_SubMIT", "plot_lcurve.py")])


def selected_count_in(output):
    """the selected number of subevents: the last line of lcurve.txt, which ./submit reads"""
    match = re.search(r"^SUGGESTED number of subevents \(smallest N with delta_aic below 2\): (\d+)$", output, re.MULTILINE)
    assert match, f"no selection line in:\n{output}"
    return int(match.group(1))


def lcurve_row_of(output, subevent_count):
    """fields of the row of one N in the table of plot_lcurve.py: n, best misfit, reduction (%), AIC above the
    lowest, chains, chains within 1%, island bests"""
    for line in output.splitlines():
        fields = line.split()
        if len(fields) == 7 and fields[0] == str(subevent_count):
            return fields
    raise AssertionError(f"no row for N={subevent_count} in:\n{output}")


def test_every_added_subevent_has_to_lower_the_misfit_by_five_percent():
    # published v2 Venezuela run: the rule "smallest N within 5% of the lowest misfit" rejected 5
    # subevents because 7 fit 6.5% better, 3.3% per added subevent, and selected 6. The same rule gave
    # 5 and 6 in the two v3 runs (16 and 32 ranks); the AIC gives 5 in all three (in the 32-rank run 6
    # subevents are tied with 5, 0.2 above it, and the smaller is selected).
    venezuela_runs = {
        "published v2": {1: 5.5472, 2: 4.6327, 3: 3.8535, 4: 3.5992, 5: 3.2028, 6: 3.0719, 7: 2.9938},
        "v3, 16 ranks": {1: 4.1304, 2: 3.2668, 3: 2.6448, 4: 2.3678, 5: 2.1477, 6: 2.0976},
        "v3, 32 ranks": {1: 4.1238, 2: 3.2735, 3: 2.5324, 4: 2.3417, 5: 2.0562, 6: 1.9567, 7: 1.9224},
    }
    for run_name, best_misfit_per_count in venezuela_runs.items():
        output = lcurve_report_of(best_misfit_per_count)
        assert selected_count_in(output) == 5, (run_name, output)


def test_three_searches_of_one_event_select_the_same_number_and_a_tie_is_named():
    # California v3: the 5% band selected 4 with 16 ranks and 5 with 32 ranks, where 4 subevents missed
    # the band by 0.05%. Two 32-rank searches that differ only in their chain seeds then put the lowest
    # AIC at 4 and at 5: the fifth subevent lowered the misfit by 4.8% and by 5.0%, and it costs 4.9%.
    # Both name the tie of 4 and 5 and both select 4, the smaller of the tied numbers.
    output_16_ranks = lcurve_report_of({1: 2.3985, 2: 1.8079, 3: 1.6501, 4: 1.5350, 5: 1.5069})
    output_32_ranks = lcurve_report_of({1: 2.3857, 2: 1.8898, 3: 1.6189, 4: 1.4992, 5: 1.4271, 6: 1.4594})
    output_32_ranks_other_seeds = lcurve_report_of({1: 2.3906, 2: 1.7922, 3: 1.6305, 4: 1.4782, 5: 1.4044, 6: 1.3572})
    for output in (output_16_ranks, output_32_ranks, output_32_ranks_other_seeds):
        assert selected_count_in(output) == 4, output
    assert "a tie" not in output_16_ranks, output_16_ranks
    assert ("NOTE: a tie - N=4 (delta_aic 0.0), N=5 (delta_aic 0.3) are equally supported (delta_aic below 2); "
            "N=4 is selected as the smallest: report both") in output_32_ranks, output_32_ranks
    assert ("NOTE: a tie - N=4 (delta_aic 0.5), N=5 (delta_aic 0.0) are equally supported (delta_aic below 2); "
            "N=4 is selected as the smallest: report both") in output_32_ranks_other_seeds, output_32_ranks_other_seeds
    assert lcurve_row_of(output_32_ranks, 5)[1:4] == ["1.4271", "4.8", "0.3"], output_32_ranks
    assert lcurve_row_of(output_32_ranks_other_seeds, 5)[1:4] == ["1.4044", "5.0", "0.0"], output_32_ranks_other_seeds


def test_grid_extension_follows_the_lowest_aic_and_not_the_selection():
    # the second California search on its first grid (1 to 5 subevents): 4 is selected, and the lowest
    # AIC is at 5, the largest N tried. ./submit has to add 6 subevents, which might be clearly better
    # than both; it reads the two numbers with the two functions tested here.
    with tempfile.TemporaryDirectory() as folder:
        for subevent_count, best_misfit in {1: 2.3906, 2: 1.7922, 3: 1.6305, 4: 1.4782, 5: 1.4044}.items():
            write_inversion_folder(folder, subevent_count, [best_misfit] * 4)
        output = run_in(folder, [sys.executable, os.path.join(programs, "prepare_SubMIT", "plot_lcurve.py")])
        reader_functions = [line for line in open(os.path.join(repository_root, "submit"))
                            if line.startswith(("selected_subevent_count()", "lowest_aic_subevent_count()"))]
        assert len(reader_functions) == 2, reader_functions
        numbers_read_by_submit = run_in(folder, ["bash", "-c", "".join(reader_functions)
                                                 + "echo $(selected_subevent_count) $(lowest_aic_subevent_count)"])
    assert numbers_read_by_submit.split() == ["4", "5"], numbers_read_by_submit
    assert "the AIC is lowest at N=5, the largest number of subevents tried" in output, output
    assert "read 4 as a lower bound" in output, output
    submit_source = open(os.path.join(repository_root, "submit")).read()
    assert 'while [ "$(lowest_aic_subevent_count)" = "$n" ]' in submit_source, "the grid extension no longer follows the lowest AIC"


def write_last_centroid_and_range(root, subevent_count, best_misfit, last_centroid_s, centroid_range_s):
    """give one fake inversion folder what the bound test of plot_lcurve.py reads: the best model (its
    subevents 5 s apart, the last one at last_centroid_s) and the centroid range of search_par.file"""
    folder = os.path.join(root, f"inv_test_{subevent_count}sub")
    centroid_times = [last_centroid_s - 5.0 * (subevent_count - 1 - subevent) for subevent in range(subevent_count)]
    with open(os.path.join(folder, "best_model_hybrid.dat"), "w") as model_file:
        model_file.write(f"{subevent_count:2d} <neq " + "".join(
            f"{best_misfit:10.4f}{time:10.4f}" + "    0.0000" * 6 for time in centroid_times) + "\n")
    for best_file in glob.glob(os.path.join(folder, "[0-9]*best.dat")):
        shutil.copy(os.path.join(folder, "best_model_hybrid.dat"), best_file)
    with open(os.path.join(folder, "search_par.file"), "w") as search_file:
        search_file.write(f"neq_min {subevent_count}\ncen_min {centroid_range_s[0]}\ncen_max {centroid_range_s[1]}\n")


def test_a_best_model_that_ends_on_the_bound_of_its_centroid_range_is_not_selected():
    # 2026 Vanuatu Mw 7.3, 121 km deep, v3 with the record screen, two sets of chain seeds: three
    # subevents fit 11% better than one, and the third sat at 28.0 s of a centroid range that ends at
    # 28 s (53 km east in a box that ends at 53 km), where it fitted the three regional stations (with
    # a wider range it settles where its P wave arrives with the pP phases of the first two subevents).
    # The lowest AIC selected 3 in both searches; a subevent in the corner of its search range is not a
    # subevent to report.
    vanuatu_searches = {
        "seeds 0":   {1: (1.7243, 9.3), 2: (1.6653, 10.5), 3: (1.5341, 28.0), 4: (1.4955, 28.0), 5: (1.5449, 27.9)},
        "seeds 500": {1: (1.7138, 9.3), 2: (1.7534, 10.4), 3: (1.5197, 27.97), 4: (1.4842, 28.0), 5: (1.7611, 28.0)},
    }
    for search_name, best_misfit_and_last_centroid in vanuatu_searches.items():
        with tempfile.TemporaryDirectory() as folder:
            for subevent_count, (best_misfit, last_centroid) in best_misfit_and_last_centroid.items():
                write_inversion_folder(folder, subevent_count, [best_misfit] * 4)
                write_last_centroid_and_range(folder, subevent_count, best_misfit, last_centroid, (1.5, 28))
            output = run_in(folder, [sys.executable, os.path.join(programs, "prepare_SubMIT", "plot_lcurve.py")])
        assert selected_count_in(output) == 1, (search_name, output)
        assert "NOTE: N=3,4,5 cannot be selected" in output, (search_name, output)
        assert lcurve_row_of(output, 3)[3] == "bound" and lcurve_row_of(output, 2)[3] != "bound", (search_name, output)
        assert "lowest AIC: N=1" in output, (search_name, output)


def test_a_last_subevent_inside_its_centroid_range_keeps_its_place_in_the_selection():
    # California, 32 ranks: the last subevent of the 4-subevent model lies at 15.37 s of a range from
    # 1.5 to 16 s, 0.63 s from its end, and that of the other seeds at 15.67 s: outside the last 1% of
    # the range (0.145 s), where the pressed models of all v3 runs lie (0.75% at most; the nearest free
    # one, the 3-subevent model listed here, is at 2.07%). The 6-subevent model ends at 15.91 s, inside it.
    with tempfile.TemporaryDirectory() as folder:
        for subevent_count, (best_misfit, last_centroid) in {3: (1.6305, 15.70), 4: (1.4782, 15.67), 5: (1.4044, 15.51),
                                                             6: (1.3572, 15.91)}.items():
            write_inversion_folder(folder, subevent_count, [best_misfit] * 4)
            write_last_centroid_and_range(folder, subevent_count, best_misfit, last_centroid, (1.5, 16))
        output = run_in(folder, [sys.executable, os.path.join(programs, "prepare_SubMIT", "plot_lcurve.py")])
    assert selected_count_in(output) == 4, output
    assert "NOTE: N=6 cannot be selected" in output, output
    assert lcurve_row_of(output, 4)[3] == "0.5" and lcurve_row_of(output, 5)[3] == "0.0", output


def test_step_3_and_the_selection_call_the_same_share_of_the_centroid_range_its_bound():
    # "on the bound" is one definition used in two files: step 3 widens the range of the 1-subevent run
    # with it, the selection excludes an N with it
    step3_source = open(os.path.join(programs, "step3_do_inversions.sh")).read()
    selection_source = open(os.path.join(programs, "prepare_SubMIT", "plot_lcurve.py")).read()
    step3_share = float(re.search(r"u - c < ([0-9.]+) \* \(u - l\)", step3_source).group(1))
    selection_share = float(re.search(r"^bound_share_of_range\s*=\s*([0-9.]+)", selection_source, re.MULTILINE).group(1))
    assert step3_share == selection_share == 0.01, (step3_share, selection_share)


def test_when_every_best_model_ends_on_the_bound_nothing_is_excluded_and_the_selection_is_not_trusted():
    # 2026 Colombia before the centroid range was widened: every N had a subevent at 30.8-31.0 s of a
    # range that ends at 31 s
    with tempfile.TemporaryDirectory() as folder:
        for subevent_count, (best_misfit, last_centroid) in {1: (3.2649, 31.0), 2: (3.7805, 31.0), 3: (3.5659, 30.9)}.items():
            write_inversion_folder(folder, subevent_count, [best_misfit] * 4)
            write_last_centroid_and_range(folder, subevent_count, best_misfit, last_centroid, (1.5, 31))
        output = run_in(folder, [sys.executable, os.path.join(programs, "prepare_SubMIT", "plot_lcurve.py")])
    assert selected_count_in(output) == 1, output
    assert "No N could be excluded: do not trust the selection below" in output, output
    assert "bound" not in lcurve_row_of(output, 1), output


def test_a_misfit_that_rises_with_more_subevents_selects_one_subevent():
    # 2026 Tonga, 234 km deep, v3 with 16 ranks: no second subevent lowered the misfit
    output = lcurve_report_of({1: 0.9084, 2: 0.9199, 3: 0.9394, 4: 0.8988, 5: 0.8912})
    assert selected_count_in(output) == 1, output
    assert lcurve_row_of(output, 2)[1:3] == ["0.9199", "-1.3"], output
    assert "a tie" not in output, output


def test_one_failed_search_in_the_middle_does_not_end_the_selection():
    # published v2 Kamchatka run: the 5-subevent search ended above the 4-subevent one (1.5965 against
    # 1.5580) and 6 subevents then fit 21% better. The lowest AIC of the whole grid is at 7; a rule
    # that stops at the first N that fails to improve would have stopped at 4.
    output = lcurve_report_of({1: 3.0806, 2: 2.3534, 3: 1.9398, 4: 1.5580, 5: 1.5965, 6: 1.2614, 7: 1.1070, 8: 1.1205})
    assert selected_count_in(output) == 7, output


def test_island_bests_are_listed_and_do_not_vote():
    # the islands of a v3 run stop at 75% of the burn-in: bootstrapped as if they were finished
    # searches, they voted for 6 subevents in 71% of the resamples of the 32-rank Venezuela run, whose
    # finished chains select 5
    output = lcurve_report_of({5: 2.0562, 6: 1.9567, 7: 1.9224},
                              {5: [2.0642, 2.1526, 2.1845, 2.2381], 6: [1.9585, 2.0209, 2.0215, 2.0420],
                               7: [1.9172, 1.9246, 1.9246, 1.9424]})
    assert selected_count_in(output) == 5, output
    assert lcurve_row_of(output, 5) == ["5", "2.0562", "-", "0.0", "4", "4", "2.0642,2.1526,2.1845,2.2381"], output
    assert "WARNING" not in output and "bootstrap" not in output, output


def test_selection_and_sampler_use_the_same_posterior_width():
    # the 5% of the rule is 2 k sigma^2 with the sigma of the posterior that finv samples: the two
    # constants live in different files and must not drift apart
    sampler_source = open(os.path.join(programs, "code_SubMIT", "finv.f90")).read()
    selection_source = open(os.path.join(programs, "prepare_SubMIT", "plot_lcurve.py")).read()
    forward_source = open(os.path.join(programs, "code_SubMIT", "forward_model.c")).read()
    sampler_sigma = float(re.search(r"sigma_relative\s*=\s*([0-9.]+)", sampler_source).group(1))
    selection_sigma = float(re.search(r"^posterior_sigma\s*=\s*([0-9.]+)", selection_source, re.MULTILINE).group(1))
    # the mechanism prior of the forward code weighs an angle against a change of misfit with the same width
    forward_sigma = float(re.search(r"const double posterior_sigma\s*=\s*([0-9.]+)", forward_source).group(1))
    parameters_per_subevent = int(re.search(r"^parameters_per_subevent\s*=\s*(\d+)", selection_source, re.MULTILINE).group(1))
    assert sampler_sigma == selection_sigma == forward_sigma, (sampler_sigma, selection_sigma, forward_sigma)
    assert abs(2 * parameters_per_subevent * selection_sigma**2 - 0.05) < 1e-12


def test_lcurve_reports_a_selection_at_the_largest_n_tried_as_a_lower_bound():
    # ./submit stops extending the grid at 7 subevents (2026-11-16); in the v3 Kamchatka rerun the
    # misfit still fell 10% from 7 to 8, so a selection at the edge must not read as an optimum
    output = lcurve_report_of({5: 1.12, 6: 0.90, 7: 0.82})
    assert selected_count_in(output) == 7, output
    assert "read 7 as a lower bound" in output, output


def test_submit_refuses_more_subevents_than_the_limit():
    with tempfile.TemporaryDirectory() as folder:
        result = subprocess.run(["bash", os.path.join(repository_root, "submit"), "run", folder, "--nmax", "8"],
                                capture_output=True, text=True)
    assert result.returncode != 0, result
    assert "outside 1..7" in result.stderr, result.stderr


def test_submit_refuses_to_start_before_the_bundled_codes_are_built():
    # PITFALLS E4: an example started in a fresh clone without ./submit build computed the regional
    # Green's functions with another fk found on the PATH and failed in step 4. The copy of the driver
    # below sees an examples folder and an empty programs folder, as a clone does before the build.
    with tempfile.TemporaryDirectory() as folder:
        shutil.copy(os.path.join(repository_root, "submit"), folder)
        os.makedirs(os.path.join(folder, "examples", "california"))
        os.makedirs(os.path.join(folder, "programs"))
        os.makedirs(os.path.join(folder, "own_event", "IRIS"))
        for arguments in (["example", "california"], ["run", os.path.join(folder, "own_event")]):
            result = subprocess.run(["bash", os.path.join(folder, "submit")] + arguments, capture_output=True, text=True)
            assert result.returncode != 0, result
            assert "is not built - run ./submit build first" in result.stderr, result.stderr
            assert not os.path.exists(os.path.join(folder, "work")), "the example was set up before the check"


# ---------------------------------------------------------------------------
# depth grid and Green's functions (PITFALLS B5)
# ---------------------------------------------------------------------------
def first_rock_depth_of(model_rows):
    """the first depth (km) that programs/gf_SubMIT/loc/doloc.sh allows for a CRUST1.0 column, given as
    rows of (Vp, Vs, density, thickness): the two lines of the script are run as they stand"""
    script = open(os.path.join(programs, "gf_SubMIT", "loc", "doloc.sh")).read()
    rule_lines = [line for line in script.splitlines() if line.startswith(("soft_column_km=", "first_rock_depth="))]
    assert len(rule_lines) == 2, rule_lines
    with tempfile.TemporaryDirectory() as folder:
        with open(os.path.join(folder, "vmodel.txt"), "w") as model:
            for row in model_rows:
                model.write(" ".join(str(value) for value in row) + "\n")
        output = run_in(folder, ["bash", "-c", "\n".join(rule_lines) + '\necho "$soft_column_km $first_rock_depth"'])
    soft_column, first_depth = output.split()
    return float(soft_column), int(first_depth)


def test_depth_grid_starts_below_water_and_soft_sediments():
    # 2024 California: 2.75 km of water over 0.63 km of sediment with Vs 0.34 km/s; the 3-km node of
    # the old grid was a source in that sediment, and every shallow subevent of the v3 models sat on it
    california = [(1.50, 0.01, 1.02, 2.75), (1.75, 0.34, 1.82, 0.63), (5.00, 2.70, 2.55, 0.64), (6.50, 3.70, 2.85, 1.41),
                  (7.10, 4.05, 3.05, 4.38), (7.90, 4.40, 3.26, 0)]
    assert first_rock_depth_of(california) == (3.38, 4)
    # 2025 Kamchatka: 2.45 km of water and 5 km of sediments below 2 km/s
    kamchatka = [(1.50, 0.01, 1.02, 2.45), (2.00, 0.55, 1.93, 2.00), (3.50, 1.79, 2.31, 3.00), (5.00, 2.70, 2.55, 0.86),
                 (6.50, 3.70, 2.85, 1.88), (7.77, 4.33, 3.20, 0)]
    assert first_rock_depth_of(kamchatka) == (7.45, 8)
    # 2024 Chile: 60 m of sediment; the 3-km minimum of the grid stays the limit
    chile = [(2.50, 1.08, 2.11, 0.06), (6.00, 3.50, 2.72, 14.30), (8.18, 4.54, 3.37, 0)]
    assert first_rock_depth_of(chile) == (0.06, 1)
    # a soft layer below rock is not part of the column that is skipped
    buried_soft_layer = [(6.00, 3.50, 2.72, 5.0), (3.00, 1.50, 2.30, 2.0), (8.00, 4.50, 3.30, 0)]
    assert first_rock_depth_of(buried_soft_layer) == (0.0, 1)


def test_green_function_scripts_that_count_their_jobs_are_started_with_bash():
    # "jobs -rp" caps the parallel fk and mtel3 jobs; dash, which is sh on Debian and Ubuntu, answers
    # "Illegal option -r", so under "sh doloc.sh" the cap never held
    for script in ("gf_SubMIT/loc/doloc.sh", "gf_SubMIT/dotel.sh"):
        assert "jobs -rp" in open(os.path.join(programs, script)).read(), script
    for caller in (os.path.join(programs, "step2_prepare_inversions.sh"), os.path.join(repository_root, "submit")):
        source = open(caller).read()
        assert "bash doloc.sh" in source and "bash dotel.sh" in source, caller
        assert not re.search(r"(?<!ba)sh do(loc|tel)\.sh", source), caller


def test_teleseismic_green_functions_use_the_depth_grid_and_every_record_of_the_run():
    # dotel.sh derived the depth grid with its own copy of doloc.sh's rule, and listed the station
    # distances from vel_*.z and vel_*.t: there are no velocity transverse records, so a rebuild after
    # the model-based screen left six SH-only California stations without Green's functions
    script = open(os.path.join(programs, "gf_SubMIT", "dotel.sh")).read()
    assert "< loc/distdep.dat" in script and "halfw" not in script, "dotel.sh has to read the grid that doloc.sh wrote"
    distance_line = next(line for line in script.splitlines() if line.startswith("gcdistlst="))
    assert "data/*.z" in distance_line and "data/*.t" in distance_line and "vel_" not in distance_line, distance_line


# ---------------------------------------------------------------------------
# the first subevent at the hypocentre, and the speed cap counted from the time it begins (PITFALLS D6)
# ---------------------------------------------------------------------------
def test_causality_cap_is_the_mean_s_velocity_of_the_searched_depths():
    # the cap was the S velocity of iasp91 at the catalogue depth (3.4 km/s for the California example,
    # whose subevents lie in oceanic mantle with 4.4 km/s), times 1.5 for shallow strike-slip events
    california_model = """#vmodel
10 1.50 0.01 1.02 2.75
1.75 0.34 1.82 0.63
5.00 2.70 2.55 0.64
6.50 3.70 2.85 1.41
7.10 4.05 3.05 4.38
7.90 4.40 3.26 25
8.2 4.5 3.4 250
8.8 4.77 3.51 100
9.7 5.28 3.88 250
10.9 6.1 4.4 0
3 5.8 3.36 2.72 20
6.5 3.75 2.92 15
8.04 4.47 3.32 0
"""
    with tempfile.TemporaryDirectory() as folder:
        os.makedirs(os.path.join(folder, "gf_test"))
        os.makedirs(os.path.join(folder, "inv_test_2sub"))
        open(os.path.join(folder, "gf_test", "vmodel"), "w").write(california_model)
        script = ("import sys; sys.path.insert(0, sys.argv[1]); import get_search_par as g; "
                  "print(g.mean_s_velocity_of_source_model('test', 4, 20), g.mean_s_velocity_of_source_model('test', 300, 700), "
                  "g.mean_s_velocity_of_source_model('test', 20, 4))")
        output = run_in(os.path.join(folder, "inv_test_2sub"), [sys.executable, "-c", script, os.path.join(programs, "prepare_SubMIT")])
    rock_grid_mean, deep_mean, empty_range = output.split()
    # 0.02 km at 2.70, 1.41 km at 3.70, 4.38 km at 4.05 and 10.19 km at 4.40 km/s
    assert abs(float(rock_grid_mean) - 4.24) < 0.005, output
    # 300-384.81 km in the 4.77 layer, 384.81-634.81 in the 5.28 layer, the rest in the half-space
    expected_deep = (84.81 * 4.77 + 250 * 5.28 + 65.19 * 6.1) / 400
    assert abs(float(deep_mean) - expected_deep) < 0.005, output
    assert empty_range == "None", output
    source = open(os.path.join(programs, "prepare_SubMIT", "get_search_par.py")).read()
    assert "vfac = -1.5" not in source, "no automatic supershear headroom (PITFALLS D6)"


def test_results_table_lists_the_speed_since_e1_began_and_the_angle_from_the_total_mechanism():
    # the figure shows centroids; the causality rule acts on the distance from E1 over the time since
    # E1 began, and the mechanism prior on the angle from the double couple of the catalogue tensor
    with tempfile.TemporaryDirectory() as folder:
        inversion, forward = os.path.join(folder, "inv_test_2sub"), os.path.join(folder, "fwd_test_2sub")
        os.makedirs(inversion)
        os.makedirs(forward)
        # E1 lasts 6 s around 4 s, so it begins at 1 s; E2 lies 50 km away 10 s after that
        open(os.path.join(forward, "Input.model"), "w").write("4.0 0.0 0.0 6.0 0.01 0.0 10.0\n11.0 30.0 40.0 5.0 0.01 0.0 10.0\n")
        # fm.dat: Mxx Mxy Mxz Myy Myz Mzz. E1 is a vertical strike-slip fault striking north or east; E2 is
        # the same mechanism with its strike turned by 15 degrees
        turned = np.radians(2 * 15.0)
        open(os.path.join(forward, "fm.dat"), "w").write(
            "1e27 0.0 1.0 0.0 0.0 0.0 0.0\n" + f"1e27 {-np.sin(turned):.6f} {np.cos(turned):.6f} 0.0 {np.sin(turned):.6f} 0.0 0.0\n")
        open(os.path.join(forward, "misfit.dat"), "w").write("1.0000\n")
        open(os.path.join(inversion, "totalmt.dat"), "w").write("0.0\n0.0\n1.0\n0.0\n0.0\n")
        run_in(folder, [sys.executable, os.path.join(programs, "fwd_SubMIT", "make_results_doc.py"), "test", "2"])
        table = [line for line in open(os.path.join(folder, "figs_and_results", "RESULTS.md")) if line.startswith("| E")]
    first, second = ([cell.strip() for cell in line.strip().strip("|").split("|")] for line in table)
    assert first[7:9] == ["0", "-"] and first[10] == "0", first
    assert second[7:9] == ["50", "5.00"], second
    assert second[10] == "30", second
    strikes = sorted(round(float(plane.split("/")[0])) % 180 for plane in first[9].split(";"))
    assert strikes == [0, 90], first[9]
    turned_strikes = sorted(round(float(plane.split("/")[0])) % 90 for plane in second[9].split(";"))
    assert turned_strikes in ([15, 15], [75, 75]), second[9]


def write_chain_file(folder, chain, rows_of_subevents):
    """one chain file in finv's format: per sample, 8 numbers per subevent in 10-character fields
    (misfit, centroid time, x, y, duration, rupture velocity, rupture direction, depth)"""
    with open(os.path.join(folder, f"{chain}chain.dat"), "w") as chain_file:
        for subevents in rows_of_subevents:
            chain_file.write(f"{len(subevents):2d} <neq " + "".join(f"{value:10.4f}" for subevent in subevents for value in subevent) + "\n")


def test_convergence_report_lists_only_the_sampled_parameters_of_the_first_subevent():
    # E1 of two or more subevents sits at the hypocentre: finv samples its centroid time and its
    # duration only, and a split R-hat of its constant position and depth is 0/0
    generator = np.random.default_rng(7)
    with tempfile.TemporaryDirectory() as folder:
        for chain in range(4):
            samples = []
            for _ in range(60):
                first = (1.0, 4 + generator.normal(0, 0.3), 0.0, 0.0, 6 + generator.normal(0, 0.3), 0.01, 0.0, 122.9)
                second = (1.0, 9 + generator.normal(0, 0.3), 5 + generator.normal(0, 1), -3 + generator.normal(0, 1),
                          8 + generator.normal(0, 0.3), 0.01, 0.0, 150 + generator.normal(0, 2))
                samples.append((first, second))
            write_chain_file(folder, chain, samples)
        run_in(folder, [sys.executable, os.path.join(programs, "hist", "posterior_summary.py")])
        report = open(os.path.join(folder, "convergence.txt")).read()
    first_line = next(line for line in report.splitlines() if line.strip().startswith("E1:"))
    second_line = next(line for line in report.splitlines() if line.strip().startswith("E2:"))
    assert re.fullmatch(r"\s*E1: centroid time [0-9.]+, duration [0-9.]+", first_line), first_line
    assert all(name in second_line for name in ("centroid time", "x east", "y north", "duration", "depth")), second_line
    assert "nan" not in report, report


def test_bundled_models_keep_the_first_subevent_at_the_hypocentre_and_obey_the_speed_cap():
    # every best model of every bundled example with two or more subevents: E1 at the epicentre and at the
    # catalogue depth (moved into the depth grid), and no centroid ahead of a front that leaves E1 with
    # vmax_cau when E1 begins. The first v3 examples had E1 14-23 km below the Chile hypocentre and
    # Venezuela's E1 -> E2 at 5.5 km/s (PITFALLS D6); their Par.file had no mechanism prior (B4).
    checked_model_count = 0
    for example in sorted(glob.glob(os.path.join(repository_root, "examples", "*", "frozen"))):
        hypocentre_depth = float(open(os.path.join(os.path.dirname(example), "event_list.dat")).read().split()[-1])
        for folder in sorted(glob.glob(os.path.join(example, "*sub"))):
            subevent_count = int(os.path.basename(folder)[:-3])
            if subevent_count < 2:
                continue
            search = dict(line.split()[:2] for line in open(os.path.join(folder, "search_par.file")) if len(line.split()) >= 2)
            numbers = open(os.path.join(folder, "best_model_hybrid.dat")).read().split("<neq")[1].split()
            # per subevent: misfit, centroid time, x east, y north, duration, rupture velocity, rupture direction, depth
            model = np.array(numbers[:8 * subevent_count], dtype=float).reshape(subevent_count, 8)
            expected_depth = min(max(hypocentre_depth, float(search["z_min"])), float(search["z_max"]))
            assert abs(model[0, 2]) < 1e-3 and abs(model[0, 3]) < 1e-3, f"{folder}: E1 is not at the epicentre"
            assert abs(model[0, 7] - expected_depth) < 0.01, f"{folder}: E1 at {model[0, 7]} km, hypocentre at {expected_depth} km"
            speed_cap = float(search["vmax_cau"])
            assert 2.0 < speed_cap < 6.0, f"{folder}: vmax_cau {speed_cap} is not a speed in km/s"
            first_begins = model[0, 1] - model[0, 4] / 2
            for subevent in range(1, subevent_count):
                distance = np.sqrt((model[subevent, 2] - model[0, 2])**2 + (model[subevent, 3] - model[0, 3])**2
                                   + (model[subevent, 7] - model[0, 7])**2)
                speed = distance / (model[subevent, 1] - first_begins)
                assert speed <= speed_cap * 1.001, f"{folder}: E{subevent + 1} at {speed:.2f} km/s, the cap is {speed_cap} km/s"
            assert "MechanismSigma=" in open(os.path.join(folder, "Par.file")).read(), f"{folder}: no mechanism prior in Par.file"
            checked_model_count += 1
    assert checked_model_count >= 12, checked_model_count


def test_bundled_chile_search_from_a_model_without_a_free_slot_ends_below_it():
    # PITFALLS D7: the four subevents of the Chile model follow each other closer than twice the time gap,
    # so a fifth one drawn from the prior fits only after the last. That search began at 0.8716, above the
    # 0.8599 of four subevents, and ended at 0.8714, and the example selected 3 subevents with one set of
    # chain seeds and 5 with another. With the insertion of finv.f90 (draw_warm_start_model) the
    # 5-subevent search starts between the subevents and ends below the 4-subevent one.
    sampler_source = open(os.path.join(programs, "code_SubMIT", "finv.f90")).read()
    time_gap = float(re.search(r"t_gap_nominal\s*=\s*([0-9.]+)", sampler_source).group(1))
    assert "warm_model_is_crowded" in sampler_source, "finv.f90 has no draw for a warm model without a free slot"

    def best_model(subevent_count):
        folder = os.path.join(repository_root, "examples", "chile", "frozen", f"{subevent_count}sub")
        numbers = open(os.path.join(folder, "best_model_hybrid.dat")).read().split("<neq")[1].split()
        # per subevent: misfit, centroid time, x east, y north, duration, rupture velocity, rupture direction, depth
        return np.array(numbers[:8 * subevent_count], dtype=float).reshape(subevent_count, 8)

    four, five = best_model(4), best_model(5)
    gaps_of_four = np.diff(four[:, 1])
    assert (gaps_of_four < 2 * time_gap).all(), f"the 4-subevent model has a free slot (gaps {gaps_of_four}): pick another crowded model"
    assert five[0, 0] < four[0, 0], f"5 subevents ended at {five[0, 0]}, not below the {four[0, 0]} of 4 subevents"


def test_bundled_depth_grids_start_in_rock():
    # the depth grid of every bundled example starts at least 0.5 km below the water and the soft
    # sediments of its source-side model (California's 3-km node was a source in 0.63 km of sediment,
    # PITFALLS B5)
    checked_example_count = 0
    for inputs in sorted(glob.glob(os.path.join(repository_root, "examples", "*", "gf_inputs"))):
        numbers = open(os.path.join(inputs, "vmodel_tel")).read().split("\n", 1)[1].split()
        layer_count = int(numbers[0])
        soft_column = 0.0
        for index in range(layer_count):
            _vp, vs, _density, thickness = (float(value) for value in numbers[1 + 4 * index:5 + 4 * index])
            if vs >= 2.0:
                break
            soft_column += thickness
        first_depth = float(open(os.path.join(inputs, "distdep.dat")).read().split()[3])
        assert first_depth >= soft_column + 0.5, f"{inputs}: grid from {first_depth} km, soft column {soft_column} km"
        checked_example_count += 1
    assert checked_example_count == 4, checked_example_count


# ---------------------------------------------------------------------------
# misfit components next to the L-curve (README section 4.6)
# ---------------------------------------------------------------------------
misfit_components_header = ("chain sampled_misfit penalized_misfit data_misfit clvd_factor stress_factor worst_p_deg "
                            "worst_t_deg p_misfit sh_misfit regional_misfit\n")


def write_misfit_components(root, subevent_count, rows):
    """misfit_components.txt of one fake inversion folder: rows of (chain, sampled, data, clvd, stress, worst_p);
    the re-evaluated penalized misfit is written equal to the sampled one, as a consistent re-evaluation gives"""
    with open(os.path.join(root, f"inv_test_{subevent_count}sub", "misfit_components.txt"), "w") as table:
        table.write(misfit_components_header)
        for chain, sampled, data, clvd, stress, worst_p in rows:
            table.write(f"{chain} {sampled} {sampled} {data} {clvd} {stress} {worst_p} 10.0 1.1 1.2 2.5\n")


def test_lcurve_lists_the_misfit_factors_of_each_best_model():
    # California v3 (32 ranks): 4 -> 5 subevents lowered the data misfit 6.6% and raised the CLVD
    # factor 1.9%; the L-curve alone showed only the product
    with tempfile.TemporaryDirectory() as folder:
        write_inversion_folder(folder, 4, [1.4992, 1.52])
        write_inversion_folder(folder, 5, [1.4271, 1.46])
        write_misfit_components(folder, 4, [(0, 1.4992, 1.4348, 1.0449, 1.0, 29.8), (1, 1.52, 1.45, 1.0483, 1.0, 28.0)])
        write_misfit_components(folder, 5, [(0, 1.4271, 1.3398, 1.0651, 1.0, 29.6), (1, 1.46, 1.37, 1.0657, 1.0, 29.0)])
        output = run_in(folder, [sys.executable, os.path.join(programs, "prepare_SubMIT", "plot_lcurve.py")])
        lcurve_text = open(os.path.join(folder, "lcurve.txt")).read()
    component_rows = {line.split()[0]: line.split()[1:6] for line in lcurve_text.splitlines() if line.startswith("N=")}
    assert component_rows["N=4"] == ["1.4992", "1.4348", "1.0449", "1.0000", "29.8"], lcurve_text
    assert component_rows["N=5"] == ["1.4271", "1.3398", "1.0651", "1.0000", "29.6"], lcurve_text
    assert "no misfit components" not in output, output
    assert lcurve_text.rstrip().endswith("SUGGESTED number of subevents (smallest N with delta_aic below 2): 4"), lcurve_text


def test_lcurve_ignores_misfit_components_left_over_from_an_earlier_run():
    # a table whose best sampled misfit is not the folder's best describes other models
    with tempfile.TemporaryDirectory() as folder:
        write_inversion_folder(folder, 1, [2.0, 2.1])
        write_inversion_folder(folder, 2, [1.6, 1.7])
        write_misfit_components(folder, 1, [(0, 2.0, 1.9, 1.0526, 1.0, 0.1)])
        write_misfit_components(folder, 2, [(0, 1.65, 1.5, 1.1, 1.0, 20.0)])
        output = run_in(folder, [sys.executable, os.path.join(programs, "prepare_SubMIT", "plot_lcurve.py")])
    component_rows = {line.split()[0]: line.split()[1] for line in output.splitlines() if line.startswith("N=")}
    assert component_rows == {"N=1": "2.0000"}, output
    assert "NOTE: no misfit components for N=2" in output, output


def test_misfit_components_come_from_one_ffwd_run_and_must_match_the_sampled_misfits():
    # evaluate_misfit_components.sh with a stand-in for ffwd that scores a model by the centroid time
    # of its first subevent: the table must keep the chain order and the sampled misfits. One chain
    # that re-evaluates off is counted and kept (32-rank Venezuela 1sub: 1 of 32 chains 0.145% off,
    # a 4-decimal sample next to a jump of the misfit, once suppressed the table of that N); misfits
    # that are off for every chain (another code version) suppress the table
    fake_ffwd = ("#!/bin/bash\n"
                 "gawk 'NF > 0 && !in_model { time = $1; in_model = 1 }\n"
                 "     NF == 0 && in_model { print \"bandresid: 2e-6 1e-6 3e-7 bandnorm: 1e-6 1e-6 1e-7\";\n"
                 "                       printf \"misfit_components: %f 1.050000 1.000000 %f 29.5 12.0\\n\", time / 1.05, time;\n"
                 "                       print \"penalized misfit\", time; in_model = 0 }' \"$1\"\n")
    script = os.path.join(programs, "prepare_SubMIT", "evaluate_misfit_components.sh")
    with tempfile.TemporaryDirectory() as folder:
        fake_ffwd_path = os.path.join(folder, "fake_ffwd")
        open(fake_ffwd_path, "w").write(fake_ffwd)
        os.chmod(fake_ffwd_path, 0o755)
        inversion_dir = os.path.join(folder, "inv_test_2sub")
        os.makedirs(inversion_dir)
        for chain, (sampled, first_time) in enumerate(((1.5022, 1.5022), (1.4992, 1.4992), (1.6, 1.6))):
            subevent_blocks = (f"{sampled:10.4f}{first_time:10.4f}" + "    0.0000" * 6 +
                               f"{sampled:10.4f}" + "    5.0000" + "    0.0000" * 6)
            open(os.path.join(inversion_dir, f"{chain}best.dat"), "w").write(" 2 <neq " + subevent_blocks + "\n")
        environment = dict(os.environ, FFWD_PROGRAM=fake_ffwd_path)
        result = subprocess.run(["bash", script, "inv_test_2sub"], cwd=folder, env=environment, capture_output=True, text=True)
        table_lines = open(os.path.join(inversion_dir, "misfit_components.txt")).read().splitlines()
        # chain 2 alone re-evaluated 6% away from its sampled misfit
        open(os.path.join(inversion_dir, "2best.dat"), "w").write(
            " 2 <neq " + f"{1.6:10.4f}{1.7:10.4f}" + "    0.0000" * 6 + f"{1.6:10.4f}" + "    5.0000" + "    0.0000" * 6 + "\n")
        one_chain_off = subprocess.run(["bash", script, "inv_test_2sub"], cwd=folder, env=environment, capture_output=True, text=True)
        one_chain_off_table = open(os.path.join(inversion_dir, "misfit_components.txt")).read().splitlines()
        # every chain re-evaluated 5% away
        for chain, sampled in enumerate((1.5022, 1.4992, 1.6)):
            open(os.path.join(inversion_dir, f"{chain}best.dat"), "w").write(
                " 2 <neq " + f"{sampled:10.4f}{1.05 * sampled:10.4f}" + "    0.0000" * 6 + f"{sampled:10.4f}" + "    5.0000" + "    0.0000" * 6 + "\n")
        mismatch = subprocess.run(["bash", script, "inv_test_2sub"], cwd=folder, env=environment, capture_output=True, text=True)
        mismatch_left_a_table = os.path.exists(os.path.join(inversion_dir, "misfit_components.txt"))
    assert "misfit components: inv_test_2sub, 3 chains; best data 1.4278 x CLVD 1.0500 x stress 1.0000" in result.stdout, result
    assert table_lines[0] == misfit_components_header.strip(), table_lines
    assert [line.split()[:3] for line in table_lines[1:]] == [["0", "1.5022", "1.502200"], ["1", "1.4992", "1.499200"],
                                                              ["2", "1.6000", "1.600000"]], table_lines
    assert table_lines[1].split()[8:] == ["2.0000", "1.0000", "3.0000"], table_lines
    assert "1 chain(s) re-evaluated more than 0.1% from the sampled misfit" in one_chain_off.stdout, one_chain_off
    assert one_chain_off_table[3].split()[:3] == ["2", "1.6000", "1.700000"], one_chain_off_table
    assert "differ from the sampled ones by 5.00% for the median chain" in mismatch.stdout, mismatch
    assert not mismatch_left_a_table


# ---------------------------------------------------------------------------
# regional instruments (PITFALLS A2, A4)
# ---------------------------------------------------------------------------
def load_instrument_chooser():
    sys.path.insert(0, os.path.join(programs, "select"))
    import choose_regional_instruments
    return choose_regional_instruments


def test_clip_rule_ignores_repeated_samples_of_a_smooth_strong_motion_peak():
    # 2026 Venezuela CM.ARGC HNE: 6 equal samples near a peak at 0.6% of full scale, not clipped
    chooser = load_instrument_chooser()
    time = np.arange(0, 60, 0.005)
    smooth = np.round(53696 * np.sin(2 * np.pi * 0.05 * time))
    assert chooser.longest_flat_run_at_peak(smooth) >= 3
    peak = np.max(np.abs(smooth))
    assert not (chooser.longest_flat_run_at_peak(smooth) >= chooser.flat_run_minimum
                and peak >= chooser.plateau_minimum_fraction * chooser.full_scale_counts)


def test_clip_rule_flags_a_flat_topped_broadband_record():
    chooser = load_instrument_chooser()
    clipped = np.clip(2e7 * np.sin(np.linspace(0, 20, 4000)), -8388607, 8388607)
    assert chooser.longest_flat_run_at_peak(clipped) >= chooser.flat_run_minimum


# ---------------------------------------------------------------------------
# Green's-function file names (station lists must name the files dotel.sh wrote)
# ---------------------------------------------------------------------------
def load_function_from_script(script_path, function_name):
    """one function of a script that runs on import (prepare_stationinfo.py writes the station lists)"""
    import ast
    tree = ast.parse(open(script_path).read())
    function_node = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == function_name)
    namespace = {}
    exec(compile(ast.Module(body=[function_node], type_ignores=[]), script_path, "exec"), namespace)
    return namespace[function_name]


def test_green_function_key_follows_the_saclst_print_used_by_dotel():
    # dotel.sh names the files after saclst's 7-digit print of dist; int(dist) gave 6202 for a
    # 6203.grn.* file (float32 dist 6202.9995), and step 3 stopped on a missing Green's function
    key = load_function_from_script(os.path.join(programs, "code_SubMIT", "prepare_stationinfo.py"),
                                    "green_function_distance_key")
    assert key(6202.99951171875) == 6203
    assert key(10999.9990234375) == 11000
    assert key(12345.6789) == 12345
    assert key(8713.0) == 8713


# ---------------------------------------------------------------------------
# records at another sampling interval (bundled California inputs: G.DZM.00 at 0.02 s)
# ---------------------------------------------------------------------------
def test_records_at_another_sampling_interval_are_set_aside():
    from obspy import Trace
    from obspy.io.sac import SACTrace
    sys.path.insert(0, os.path.join(programs, "prepare_SubMIT"))
    import screen_station_amplitudes as screen
    with tempfile.TemporaryDirectory() as folder:
        for name, interval in (("XX.ONE.00.z", 1.0), ("G.DZM.00.z", 0.02)):
            trace = SACTrace.from_obspy_trace(Trace(np.zeros(100, dtype=np.float32)))
            trace.delta = interval
            trace.write(os.path.join(folder, name))
        report = []
        flagged = screen.sampling_interval_screen(sorted(glob.glob(os.path.join(folder, "*.z"))), report, 1.0)
    assert [os.path.basename(f) for f in flagged] == ["G.DZM.00.z"], flagged


# ---------------------------------------------------------------------------
# dense teleseismic networks (PITFALLS E2)
# ---------------------------------------------------------------------------
def test_dense_teleseismic_networks_are_thinned_to_the_limit():
    # 2025 Aomori: 232 stations aborted v2 in step 3 (num_sta_P/SH exceeds maxstanum)
    from obspy import Trace
    from obspy.io.sac import SACTrace
    sys.path.insert(0, os.path.join(programs, "prepare_SubMIT"))
    import screen_station_amplitudes as screen
    random_generator = np.random.default_rng(3)
    with tempfile.TemporaryDirectory() as folder:
        os.makedirs(os.path.join(folder, "data", "excluded"))
        for index in range(150):
            for component in "zrt":
                trace = SACTrace.from_obspy_trace(Trace(random_generator.normal(size=600).astype(np.float32)))
                trace.delta, trace.b, trace.t1 = 1.0, -300.0, 0.0
                trace.az, trace.gcarc = float(random_generator.uniform(0, 30)), float(random_generator.uniform(40, 90))
                trace.write(os.path.join(folder, "data", f"XX.S{index:03d}.00.{component}"))
        previous_folder = os.getcwd()
        os.chdir(folder)
        try:
            report = []
            screen.thin_dense_teleseismic_networks(report)
        finally:
            os.chdir(previous_folder)
        remaining = glob.glob(os.path.join(folder, "data", "*.z"))
    assert len(remaining) == 100, len(remaining)


# ---------------------------------------------------------------------------
# bundled examples rebuilt from a run (tools/make_example_from_run.py)
# ---------------------------------------------------------------------------
def test_frozen_example_keeps_configuration_but_never_code_chains_or_data():
    # ./submit example copies frozen/<N>sub over a fresh programs/code_SubMIT: a stale finv.f90 or
    # caldens.py in frozen/ would silently replace the current code; island bests and misfit
    # components are kept so that figures mode measures reproducibility and draws the misfit factors
    # at every N
    import importlib.util
    specification = importlib.util.spec_from_file_location(
        "make_example_from_run", os.path.join(repository_root, "tools", "make_example_from_run.py"))
    packager = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(packager)
    with tempfile.TemporaryDirectory() as folder:
        inversion_dir = os.path.join(folder, "inv_test_2sub")
        os.makedirs(os.path.join(inversion_dir, "data"))
        for name in ("Par.file", "search_par.file", "stations.info", "best_model_hybrid.dat", "island_best_misfits.txt",
                     "misfit_components.txt", "misfit_components.log", "make_ffwd.log", "ffwd", "ffwd.o",
                     "finv.f90", "forward_model.c", "caldens.py", "0best.dat", "0chain.dat", "mpi.out", "warm_start.dat"):
            open(os.path.join(inversion_dir, name), "w").write("1 2 3\n")
        frozen_dir = os.path.join(folder, "frozen_2sub")
        packager.write_frozen_folder(inversion_dir, frozen_dir, is_selected=False)
        frozen_names = set(os.listdir(frozen_dir))
    assert frozen_names == {"Par.file", "search_par.file", "stations.info", "best_model_hybrid.dat",
                            "island_best_misfits.txt", "misfit_components.txt"}, frozen_names



# ---------------------------------------------------------------------------
# what the first 1-subevent run decides: centroid range (PITFALLS D4) and records (PITFALLS A6)
# ---------------------------------------------------------------------------
def run_centroid_range_widening(first_centroid, centroid_after_each_rerun):
    """widen_centroid_range_while_on_bound of step3_do_inversions.sh on two fake inversion folders (range
    1.5-31 s, windows of the 2026 Colombia run), with a stand-in for the inversion that puts the 1-subevent
    centroid where centroid_after_each_rerun says: a number, or "bound" for just below the current upper bound.
    Returns the output, the folder contents {relative path: text} and the number of reruns."""
    step3_text = open(os.path.join(programs, "step3_do_inversions.sh")).read()
    function_text = re.search(r"^widen_centroid_range_while_on_bound\(\) \{.*?^\}", step3_text, re.DOTALL | re.MULTILINE).group(0)
    stand_in = ('run_inversion() { echo rerun >> reruns.txt\n'
                '    if [ "$new_centroid" = bound ]; then value=$(gawk \'$1=="cen_max"{print $2 - 0.01}\' "$1/search_par.file"); else value=$new_centroid; fi\n'
                '    gawk -v c="$value" \'{$4 = c; print}\' "$1/best_model_hybrid.dat" > "$1/best.tmp" && mv "$1/best.tmp" "$1/best_model_hybrid.dat"; }\n')
    with tempfile.TemporaryDirectory() as folder:
        for subevent_count in (1, 2):
            inversion_dir = os.path.join(folder, f"inv_test_{subevent_count}sub")
            os.makedirs(inversion_dir)
            open(os.path.join(inversion_dir, "search_par.file"), "w").write("burn_in 3000\ncen_min 1.5\ncen_max 31\ndura_max 30\n")
            open(os.path.join(inversion_dir, "Par.file"), "w").write(
                "nd_timeP= 123 #ending time of inverse time window (P)\nnd_timeSH= 123 #ending time (SH)\n"
                "bg_timeRayl= -20 #beginning (regional)\nnd_timeRayl= 261 #ending (regional)\n"
                "stf_btime= -18 #beginning of the moment-rate window\nstf_etime= 55 #ending of the moment-rate window\n"
                "deltat= 1.0 #sampling interval\n")
        open(os.path.join(folder, "inv_test_1sub", "best_model_hybrid.dat"), "w").write(
            f" 1 <neq 3.4508 {first_centroid} -39.1 -33.0 3.83 0.01 0.0 61.3\n")
        open(os.path.join(folder, "harness.sh"), "w").write(
            f"new_centroid={centroid_after_each_rerun}\n{stand_in}{function_text}\nwiden_centroid_range_while_on_bound inv_test_1sub test\n")
        output = run_in(folder, ["bash", "harness.sh"])
        contents = {os.path.relpath(os.path.join(root, name), folder): open(os.path.join(root, name)).read()
                    for root, _, names in os.walk(folder) for name in names}
    return output, contents, contents.get("reruns.txt", "").count("rerun")


def test_centroid_range_is_doubled_when_the_first_one_subevent_run_ends_on_its_bound():
    # 2026 Colombia Mw 7.4 (v3, two seeds): every chain of the 1-subevent run ended at 30.9-31.0 s of a
    # range that ends at 31 s, the main P pulse being about 22 s late, and the misfit then rose from 1 to
    # 2 subevents. The range and the windows of every N grow by the length of the range.
    output, contents, rerun_count = run_centroid_range_widening(30.9548, 34.0)
    assert "lengthened by 30 s" in output, output
    assert rerun_count == 1
    for subevent_count in (1, 2):
        assert "cen_max 61\n" in contents[f"inv_test_{subevent_count}sub/search_par.file"], contents
        par_text = contents[f"inv_test_{subevent_count}sub/Par.file"]
        for expected_line in ("nd_timeP= 153 #ending time of inverse time window (P)", "nd_timeSH= 153 #", "nd_timeRayl= 291 #",
                              "stf_etime= 85 #", "stf_btime= -18 #", "bg_timeRayl= -20 #"):
            assert expected_line in par_text, par_text


def test_centroid_range_is_left_alone_when_the_centroid_lies_inside_it():
    # the other ten v3 runs: the 1-subevent centroid lies in the first 8-92% of its range (Venezuela,
    # a doublet whose single subevent takes the second event: 114 s of 124)
    output, contents, rerun_count = run_centroid_range_widening(28.5, 28.5)
    assert output.strip() == "" and rerun_count == 0, output
    assert "cen_max 31\n" in contents["inv_test_2sub/search_par.file"] and "nd_timeP= 123 #" in contents["inv_test_2sub/Par.file"]


def test_centroid_range_widening_gives_up_after_two_doublings():
    # a centroid that follows the bound (a wrong origin time, a record that drifts): the loop must end
    output, contents, rerun_count = run_centroid_range_widening(30.99, "bound")
    assert rerun_count == 2 and "still on the upper bound" in output, output
    assert "cen_max 121\n" in contents["inv_test_1sub/search_par.file"], contents
    assert "nd_timeP= 213 #" in contents["inv_test_1sub/Par.file"], contents      # 123 + 30 + 60
    assert "nd_timeRayl= 351 #" in contents["inv_test_1sub/Par.file"], contents   # 261 + 30 + 60, below the 479-s cap


def write_forward_run_for_the_station_screen(folder, amplitude_of_station_3):
    """fwd_test/ as ffwd leaves it for screen_stations_model.py: ten P displacement records (and the velocity
    twin of station 3) with observed traces of unit amplitude, except station 3; inv_test_2sub/data holds the
    same record files"""
    from obspy import Trace
    forward_dir = os.path.join(folder, "fwd_test")
    os.makedirs(os.path.join(forward_dir, "waveforms"))
    random_generator = np.random.default_rng(3)
    record_names = [f"data/XX.S{station:02d}..z" for station in range(10)] + ["data/vel_XX.S03..z"]
    for data_dir in (os.path.join(forward_dir, "data"), os.path.join(folder, "inv_test_2sub", "data")):
        os.makedirs(data_dir)
        for record_name in record_names:
            open(os.path.join(data_dir, os.path.basename(record_name)), "w").write("record\n")
    with open(os.path.join(forward_dir, "stations.info"), "w") as station_list:
        for record_name in record_names[:10]:
            station_list.write(f"{record_name} 60.0 120.0 16.0 9.0 6672\n")
    open(os.path.join(forward_dir, "stationsSH.info"), "w").write("")
    with open(os.path.join(forward_dir, "sta_residual.dat"), "w") as residual_file:
        residual_file.write("# P waves\n")
        for index in range(10):
            samples = random_generator.normal(0, 1, 140)
            if index == 3:
                samples = samples * amplitude_of_station_3
            for prefix in ("P_obs_", "P_syn_"):
                Trace(data=samples.astype(np.float32)).write(os.path.join(forward_dir, "waveforms", f"{prefix}{index:04d}.sac"), format="SAC")
            residual_file.write(f"{index + 1} 1.0 0.8 0.0\n")
        residual_file.write("# SH waves\n")


def test_a_record_far_above_its_group_in_the_inversion_band_is_excluded_without_a_model():
    # 2026 Vanuatu Mw 7.3: the drifting vertical record of IM.IL31 had 27 times the median rms of the P
    # displacement records in the band and window of the inversion and carried 68% of the P misfit. The
    # screen that would have removed it was skipped as a whole, because the 1-subevent result was
    # unrepresentative; the test on the records needs no model and runs in that case too.
    screen_script = os.path.join(programs, "prepare_SubMIT", "screen_stations_model.py")
    with tempfile.TemporaryDirectory() as folder:
        write_forward_run_for_the_station_screen(folder, amplitude_of_station_3=27.0)
        output = run_in(folder, [sys.executable, screen_script, "test"])
        report = open(os.path.join(folder, "excluded_stations_model.txt")).read()
        kept = sorted(os.listdir(os.path.join(folder, "inv_test_2sub", "data")))
        excluded = sorted(os.listdir(os.path.join(folder, "inv_test_2sub", "data", "excluded")))
    assert "DATA  P waves data/XX.S03..z rms x" in report and "MODEL" not in report, report
    assert excluded == ["XX.S03..z", "vel_XX.S03..z"], excluded        # the velocity twin leaves with it
    assert "XX.S03..z" not in kept and "XX.S04..z" in kept, kept
    assert "2 trace(s) excluded" in output, output
    with tempfile.TemporaryDirectory() as folder:
        write_forward_run_for_the_station_screen(folder, amplitude_of_station_3=5.8)       # the largest of nine events
        result = subprocess.run([sys.executable, screen_script, "test"], cwd=folder, capture_output=True, text=True,
                                env=dict(os.environ, SCREEN_DATA_ONLY="1"))
    assert "nothing flagged" in result.stdout, result


if __name__ == "__main__":
    failures = 0
    for name, function in sorted(globals().items()):
        if name.startswith("test_") and callable(function):
            try:
                function()
                print(f"ok      {name}")
            except Exception as error:
                failures += 1
                print(f"FAILED  {name}: {error}")
    sys.exit(1 if failures else 0)
