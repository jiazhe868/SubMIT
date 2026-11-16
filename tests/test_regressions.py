"""Regression tests for the defects fixed in SubMIT v3 (PITFALLS.md); each test names what it guards.

They need only the conda environment (numpy, scipy, obspy, matplotlib) - no SAC, MPI or Green's
functions - and run in a few seconds on the catalogues bundled with the examples:

    python tests/test_regressions.py        (or: python -m pytest tests/test_regressions.py)
"""
import glob
import os
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
# number of subevents (PITFALLS D1)
# ---------------------------------------------------------------------------
def write_inversion_folder(root, subevent_count, chain_best_misfits):
    folder = os.path.join(root, f"inv_test_{subevent_count}sub")
    os.makedirs(folder)
    for rank, misfit in enumerate(chain_best_misfits):
        with open(os.path.join(folder, f"{rank}best.dat"), "w") as best_file:
            best_file.write(f"{subevent_count:2d} <neq " + "".join(f"{misfit:10.4f}" + "    0.0000" * 7 for _ in range(subevent_count)) + "\n")


def test_lcurve_flags_a_selection_that_rests_on_one_chain():
    # 2026 Panama HEAD run: N=2 was selected because ONE of 16 chains reached 1.5303 (median 17% worse)
    with tempfile.TemporaryDirectory() as folder:
        write_inversion_folder(folder, 1, [1.79] * 16)
        write_inversion_folder(folder, 2, [1.5303] + [1.78] * 15)
        write_inversion_folder(folder, 3, [1.5264] * 5 + [1.55] * 11)
        write_inversion_folder(folder, 4, [1.5076] * 2 + [1.56] * 14)
        write_inversion_folder(folder, 5, [1.4662] + [1.60] * 15)
        output = run_in(folder, [sys.executable, os.path.join(programs, "prepare_SubMIT", "plot_lcurve.py")])
    assert "SUGGESTED number of subevents (elbow): 2" in output, output
    assert "WARNING: unstable selection" in output, output


def write_island_bests(root, subevent_count, island_best_misfits):
    folder = os.path.join(root, f"inv_test_{subevent_count}sub")
    with open(os.path.join(folder, "island_best_misfits.txt"), "w") as island_file:
        for island, misfit in enumerate(island_best_misfits, start=1):
            island_file.write(f"{island:4d}{8:4d}{misfit:12.5f}\n")


def test_lcurve_uses_islands_when_resampled_chains_agree():
    # v3 chains descend from one island after the global resampling round, so they agree even when
    # independent searches do not (California 3sub seeds 0/500: 1.761 vs 1.847); the islands decide
    with tempfile.TemporaryDirectory() as folder:
        write_inversion_folder(folder, 1, [2.00] * 16)
        write_inversion_folder(folder, 2, [1.80] * 16)
        write_inversion_folder(folder, 3, [1.70] * 16)
        write_island_bests(folder, 1, [2.00, 2.01])
        write_island_bests(folder, 2, [1.80, 1.81])
        write_island_bests(folder, 3, [1.70, 1.90])
        output = run_in(folder, [sys.executable, os.path.join(programs, "prepare_SubMIT", "plot_lcurve.py")])
    assert "SUGGESTED number of subevents (elbow): 3" in output, output
    assert "WARNING: unstable selection" in output, output


def test_lcurve_does_not_trust_the_chains_of_a_single_island():
    # with 8 or fewer ranks finv has one island, whose resampled chains always agree: no replicates
    with tempfile.TemporaryDirectory() as folder:
        for count, misfit in ((1, 2.0), (2, 1.8), (3, 1.7)):
            write_inversion_folder(folder, count, [misfit] * 8)
            write_island_bests(folder, count, [misfit])
        output = run_in(folder, [sys.executable, os.path.join(programs, "prepare_SubMIT", "plot_lcurve.py")])
    assert "reproducibility not measured" in output, output
    assert "WARNING" not in output, output


def test_lcurve_reports_a_selection_at_the_largest_n_tried_as_a_lower_bound():
    # ./submit stops extending the grid at 7 subevents (2026-11-16); in the v3 Kamchatka rerun the
    # misfit still fell 10% from 7 to 8, so a selection at the edge must not read as an elbow
    with tempfile.TemporaryDirectory() as folder:
        write_inversion_folder(folder, 5, [1.12] * 16)
        write_inversion_folder(folder, 6, [0.90] * 16)
        write_inversion_folder(folder, 7, [0.82] * 16)
        output = run_in(folder, [sys.executable, os.path.join(programs, "prepare_SubMIT", "plot_lcurve.py")])
    assert "SUGGESTED number of subevents (elbow): 7" in output, output
    assert "read 7 as a lower bound" in output, output


def test_submit_refuses_more_subevents_than_the_limit():
    with tempfile.TemporaryDirectory() as folder:
        result = subprocess.run(["bash", os.path.join(repository_root, "submit"), "run", folder, "--nmax", "8"],
                                capture_output=True, text=True)
    assert result.returncode != 0, result
    assert "outside 1..7" in result.stderr, result.stderr


def test_lcurve_stays_quiet_when_every_chain_agrees():
    with tempfile.TemporaryDirectory() as folder:
        write_inversion_folder(folder, 1, [2.0] * 16)
        write_inversion_folder(folder, 2, [1.60] * 16)
        write_inversion_folder(folder, 3, [1.58] * 16)
        output = run_in(folder, [sys.executable, os.path.join(programs, "prepare_SubMIT", "plot_lcurve.py")])
    assert "SUGGESTED number of subevents (elbow): 2" in output, output
    assert "WARNING" not in output, output


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
    # caldens.py in frozen/ would silently replace the current code; island bests are kept so that
    # figures mode measures reproducibility at every N
    import importlib.util
    specification = importlib.util.spec_from_file_location(
        "make_example_from_run", os.path.join(repository_root, "tools", "make_example_from_run.py"))
    packager = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(packager)
    with tempfile.TemporaryDirectory() as folder:
        inversion_dir = os.path.join(folder, "inv_test_2sub")
        os.makedirs(os.path.join(inversion_dir, "data"))
        for name in ("Par.file", "search_par.file", "stations.info", "best_model_hybrid.dat", "island_best_misfits.txt",
                     "finv.f90", "forward_model.c", "caldens.py", "0best.dat", "0chain.dat", "mpi.out", "warm_start.dat"):
            open(os.path.join(inversion_dir, name), "w").write("1 2 3\n")
        frozen_dir = os.path.join(folder, "frozen_2sub")
        packager.write_frozen_folder(inversion_dir, frozen_dir, is_selected=False)
        frozen_names = set(os.listdir(frozen_dir))
    assert frozen_names == {"Par.file", "search_par.file", "stations.info", "best_model_hybrid.dat",
                            "island_best_misfits.txt"}, frozen_names


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
