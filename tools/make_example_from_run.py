#!/usr/bin/env python3
"""Rebuild one bundled example (examples/NAME/) from a finished ./submit prepare + ./submit run folder.

Why. ./submit example reproduces an example from what examples/NAME/ holds: the step-1 records,
the configuration of every inversion folder and the published figures. The v2 examples were put
together by hand, and their inputs still carried horizontals that SAC had failed to rotate
(PITFALLS A1), so they could not simply be refreshed. This script takes everything from one run, so
an example always matches a run that the current code made from the raw downloads.

What goes where (the layout ./submit example expects):
  inputs.tar.gz   IRIS/<event>/data and IRISloc/<event>_loc/data: the step-1 records
  event_list.dat  IRIS/event_list.dat of the run
  selected_n      the selection in IRIS/lcurve.txt
  frozen/<N>sub/  the final configuration of every inversion folder (after the 1-subevent
                  calibration), its best model, its island bests and its misfit components, so
                  figures mode can list the island bests and draw the factors of the misfit at every
                  N; samples.tar.gz (all chains) for the selected N only.
                  Code and scripts are left out on purpose: figures and frozen modes copy frozen/
                  over a fresh copy of programs/code_SubMIT, and stale copies would override it.
  gf_inputs/      velocity models (vmodel_tel, vmodel_loc) and distance/depth grid (distdep.dat)
                  of the Green's functions
  reference/      figures of the selected model (first pages of the PDFs as PNG, at the sizes of the
                  v2 examples), its model files (Input.model, fm.dat, sta_residual.dat), lcurve.txt
                  and RESULTS.md
  catalog/        not touched: the archived aftershock catalog and moment tensor the run replayed
frozen/, reference/ and gf_inputs/ are replaced as a whole, so no file of an older run survives.

History: written 2026-11-16 to rebuild the four examples from v3 runs at 32 MPI ranks; the examples
bundled since 2026-11-17 were made with it.

usage (from the repository root):
    python tools/make_example_from_run.py RUN_FOLDER examples/NAME
then check the result with  ./submit example NAME --workdir <short path>  and  ./submit compare NAME <same path>
"""
import fnmatch
import glob
import os
import re
import shutil
import subprocess
import sys
import tarfile

# files of an inversion folder that make up its configuration, kept in frozen/<N>sub/
configuration_files = [
    "Par.file", "search_par.file", "rank.dat", "stations.info", "stationsSH.info", "stationsloc.info",
    "vel_stations.info", "weights.dat", "distdep.dat", "edges.dat", "seisdens.dat", "seisgrids.dat",
    "lola.dat", "lola_bg.dat", "mainshock.dat", "totalmt.dat", "totalmt_source.txt", "extended.flag",
    "extended_source.txt", "excluded_stations.txt", "Input.model", "best_model_hybrid.dat", "misfit.dat",
    "island_best_misfits.txt", "misfit_components.txt",
]
# (Input.model and misfit.dat are written by step 4 for the selected N only. misfit_components.txt is a
# result, not configuration: it is kept so that figures mode, which ships the chains of the selected N
# only, still draws the factors of every N next to the L-curve; frozen mode recomputes it.)
required_configuration_files = ["Par.file", "search_par.file", "stations.info", "best_model_hybrid.dat"]
# files that are deliberately not frozen: data, chains (samples.tar.gz instead), outputs, logs, figures,
# and code or scripts that ./submit example copies fresh from programs/
left_out_patterns = [
    "data", "dataloc", "__pycache__", "[0-9]*best.dat", "[0-9]*chain.dat", "allsamples.dat", "pooled_ensemble.dat", "mpi.out",
    "*.o", "finv", "ffwd", "*.log", "*.pdf", "*.png", "*.c", "*.h", "*.f", "*.f90", "Makefile", "*.py",
    "convergence.txt", "warm_start.dat", "catalog.txt", "calibration_pass_model.txt",
]
# reference figures: PDF -> first page as PNG, at the resolution that gives the v2 example sizes
# (fits: letter page at 120 dpi = 1020 x 1320; histograms 100 dpi; L-curve 6 x 4.5 in at 110 dpi = 660 x 495,
# or about 990 x 485 since it carries the misfit components, 9 x 4.4 in)
figure_pages = {"fits_P": 120, "fits_SH": 120, "fits_rayl": 120, "histoplot_py": 100, "lcurve": 110}


def subevent_count_of(folder_name):
    return int(re.search(r"_(\d+)sub$", folder_name).group(1))


def matches_any(name, patterns):
    return any(fnmatch.fnmatch(name, pattern) for pattern in patterns)


def selected_subevent_count(lcurve_path):
    for line in open(lcurve_path):
        if line.startswith("SUGGESTED"):
            return int(line.split()[-1])
    sys.exit(f"no SUGGESTED line in {lcurve_path}: the run has not finished its selection")


def write_frozen_folder(inversion_dir, frozen_dir, is_selected):
    """copy the configuration of one inversion folder; for the selected N also pack its chains"""
    os.makedirs(frozen_dir)
    for name in configuration_files:
        source = os.path.join(inversion_dir, name)
        if os.path.exists(source):
            shutil.copy2(source, frozen_dir)
    unknown = [name for name in sorted(os.listdir(inversion_dir))
               if name not in configuration_files and not matches_any(name, left_out_patterns)]
    if unknown:
        print(f"  NOTE {os.path.basename(inversion_dir)}: not frozen (neither configuration nor known output): {' '.join(unknown)}")
    # a configuration that names a path on this machine would break everywhere else
    for name in os.listdir(frozen_dir):
        text = open(os.path.join(frozen_dir, name), errors="replace").read()
        assert not re.search(r"(/disk/|/home/|/tmp/|/Users/)", text), f"{frozen_dir}/{name} contains an absolute path"
    if is_selected:
        chain_files = sorted(glob.glob(os.path.join(inversion_dir, "[0-9]*best.dat")) +
                             glob.glob(os.path.join(inversion_dir, "[0-9]*chain.dat")))
        chain_files += [os.path.join(inversion_dir, "allsamples.dat")]
        with tarfile.open(os.path.join(frozen_dir, "samples.tar.gz"), "w:gz") as archive:
            for path in chain_files:
                archive.add(path, arcname=os.path.basename(path))


def write_reference_folder(iris, event_name, selected, reference_dir):
    results_dir = os.path.join(iris, "figs_and_results")
    forward_dir = os.path.join(iris, f"fwd_{event_name}_{selected}sub")   # step 4's forward folder
    os.makedirs(reference_dir)
    for name in ("subevents.png", "stations.png", "lcurve.txt", "RESULTS.md"):
        shutil.copy2(os.path.join(results_dir, name), reference_dir)
    for name in ("Input.model", "fm.dat", "sta_residual.dat", "misfit_evolution.png"):
        shutil.copy2(os.path.join(forward_dir, name), reference_dir)
    for stem, dots_per_inch in figure_pages.items():
        subprocess.run(["pdftoppm", "-png", "-r", str(dots_per_inch), "-f", "1", "-l", "1", "-singlefile",
                        os.path.join(results_dir, f"{stem}.pdf"), os.path.join(reference_dir, stem)], check=True)


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    run_folder, example_dir = (os.path.abspath(path) for path in sys.argv[1:3])
    iris = os.path.join(run_folder, "IRIS")
    event_name = open(os.path.join(iris, "event_list.dat")).read().split()[0]
    example_event = open(os.path.join(example_dir, "event_list.dat")).read().split()[0]
    assert event_name == example_event, f"the run is {event_name}, the example is {example_event}"
    selected = selected_subevent_count(os.path.join(iris, "lcurve.txt"))
    inversion_dirs = sorted(glob.glob(os.path.join(iris, f"inv_{event_name}_*sub")), key=subevent_count_of)
    counts = [subevent_count_of(folder) for folder in inversion_dirs]
    assert selected in counts, f"selected N={selected} has no inversion folder"
    print(f"{event_name}: subevent counts {counts}, selected {selected}")

    # check every input before anything in the example is deleted
    gf_dir = os.path.join(iris, f"gf_{event_name}")
    forward_dir = os.path.join(iris, f"fwd_{event_name}_{selected}sub")
    needed = [os.path.join(folder, name) for folder in inversion_dirs for name in required_configuration_files]
    needed += [os.path.join(iris, f"inv_{event_name}_{selected}sub", "allsamples.dat")]
    needed += [os.path.join(gf_dir, "loc", "distdep.dat"), os.path.join(gf_dir, "loc", "vmodel"), os.path.join(gf_dir, "vmodel")]
    needed += [os.path.join(iris, "figs_and_results", name) for name in ("subevents.png", "stations.png", "lcurve.txt", "RESULTS.md")]
    needed += [os.path.join(iris, "figs_and_results", f"{stem}.pdf") for stem in figure_pages]
    needed += [os.path.join(forward_dir, name) for name in ("Input.model", "fm.dat", "sta_residual.dat", "misfit_evolution.png")]
    missing = [path for path in needed if not os.path.exists(path)]
    assert not missing, "the run is incomplete, nothing was changed; missing:\n  " + "\n  ".join(missing)

    for name in ("frozen", "reference", "gf_inputs"):
        shutil.rmtree(os.path.join(example_dir, name), ignore_errors=True)
    for inversion_dir, count in zip(inversion_dirs, counts):
        write_frozen_folder(inversion_dir, os.path.join(example_dir, "frozen", f"{count}sub"), count == selected)

    os.makedirs(os.path.join(example_dir, "gf_inputs"))
    shutil.copy2(os.path.join(gf_dir, "loc", "distdep.dat"), os.path.join(example_dir, "gf_inputs", "distdep.dat"))
    shutil.copy2(os.path.join(gf_dir, "loc", "vmodel"), os.path.join(example_dir, "gf_inputs", "vmodel_loc"))
    shutil.copy2(os.path.join(gf_dir, "vmodel"), os.path.join(example_dir, "gf_inputs", "vmodel_tel"))
    grid_of_gf = open(os.path.join(gf_dir, "loc", "distdep.dat")).read().split()
    grid_of_run = open(os.path.join(inversion_dirs[0], "distdep.dat")).read().split()
    assert grid_of_gf == grid_of_run, f"Green's-function grid {grid_of_gf} differs from the inversion grid {grid_of_run}"

    write_reference_folder(iris, event_name, selected, os.path.join(example_dir, "reference"))
    shutil.copy2(os.path.join(iris, "event_list.dat"), os.path.join(example_dir, "event_list.dat"))
    open(os.path.join(example_dir, "selected_n"), "w").write(f"{selected}\n")

    # step-1 records, sorted so that the archive lists them in a fixed order
    with tarfile.open(os.path.join(example_dir, "inputs.tar.gz"), "w:gz") as archive:
        for data_dir in (os.path.join("IRIS", event_name, "data"), os.path.join("IRISloc", f"{event_name}_loc", "data")):
            archive.add(os.path.join(run_folder, data_dir), arcname=data_dir, recursive=False)
            for name in sorted(os.listdir(os.path.join(run_folder, data_dir))):
                archive.add(os.path.join(run_folder, data_dir, name), arcname=os.path.join(data_dir, name))
    print(f"wrote {example_dir}: inputs.tar.gz, frozen/ ({len(counts)} folders), gf_inputs/, reference/, selected_n")


if __name__ == "__main__":
    main()
