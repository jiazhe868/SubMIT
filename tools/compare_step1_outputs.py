#!/usr/bin/env python3
"""Compare the step-1 output of an event folder with the processed inputs bundled with an example.

Why it exists. Step 1 (./submit prepare) turns the original Wilber downloads into the records that
every later step reads, and two of its defects were silent: horizontal pairs of unequal length were
written unrotated under the names of the radial and transverse components (PITFALLS A1; 26 of the 56
SH records of the published Venezuela model), and regional records were scored against the P coda
(A3). This script is the regression check for step 1: it lists the records that appear or disappear,
the records whose samples or length changed, and, in both data sets, the radial and transverse
records that are still unrotated (their component name ends in E, N, 1 or 2 instead of R or T).

usage:
    ./submit prepare /data/chile_rerun          # on the original Wilber downloads of an example
    python tools/compare_step1_outputs.py /data/chile_rerun examples/chile/inputs.tar.gz
"""
import glob
import os
import sys
import tarfile
import tempfile
import warnings

import numpy as np
from obspy import read

# obspy rounds the 1-s sample spacing that SAC stores as 0.99999994 and warns about it on every file
warnings.filterwarnings("ignore", message="Sample spacing read from SAC file")

relative_sample_tolerance = 1e-5       # float32 SAC samples written by two runs of the same code agree to ~1e-7


def data_folder(root, kind):
    """the data/ folder of the teleseismic (IRIS/<event>/data) or regional (IRISloc/<event>_loc/data) records"""
    pattern = os.path.join(root, "IRIS", "[12]*", "data") if kind == "teleseismic" \
        else os.path.join(root, "IRISloc", "[12]*_loc", "data")
    folders = [folder for folder in glob.glob(pattern) if os.path.isdir(folder)]
    return folders[0] if folders else None


def unrotated_horizontal_records(folder):
    """radial/transverse files whose SAC component name shows they were never rotated (PITFALLS A1)"""
    unrotated = []
    for path in sorted(glob.glob(os.path.join(folder, "*.[rt]"))):
        component_name = read(path, headonly=True)[0].stats.sac.get("kcmpnm", "").strip()
        if component_name and component_name[-1] in "EN12":
            unrotated.append(os.path.basename(path))
    return unrotated


def compare_folders(new_folder, reference_folder):
    """print what changed between two data folders, file by file"""
    new_files = {name for name in os.listdir(new_folder) if not name.startswith(".")}
    reference_files = {name for name in os.listdir(reference_folder) if not name.startswith(".")}
    length_changed, samples_changed, largest_change = [], [], (0.0, None)
    for name in sorted(new_files & reference_files):
        new_trace = read(os.path.join(new_folder, name))[0]
        reference_trace = read(os.path.join(reference_folder, name))[0]
        if new_trace.stats.npts != reference_trace.stats.npts or abs(new_trace.stats.delta - reference_trace.stats.delta) > 1e-6:
            length_changed.append(name)
            continue
        reference_peak = max(np.abs(reference_trace.data).max(), 1e-30)
        relative_change = float(np.abs(new_trace.data - reference_trace.data).max() / reference_peak)
        if relative_change > relative_sample_tolerance:
            samples_changed.append(name)
        if relative_change > largest_change[0]:
            largest_change = (relative_change, name)
    print(f"  records: {len(new_files)} new, {len(reference_files)} reference, {len(new_files & reference_files)} in both")
    print(f"  only in new ({len(new_files - reference_files)}): {' '.join(sorted(new_files - reference_files))}")
    print(f"  only in reference ({len(reference_files - new_files)}): {' '.join(sorted(reference_files - new_files))}")
    print(f"  different length ({len(length_changed)}): {' '.join(length_changed)}")
    print(f"  different samples ({len(samples_changed)}, largest {largest_change[0]:.2e} in {largest_change[1]}): "
          f"{' '.join(samples_changed)}")


def main(event_folder, reference_tarball):
    with tempfile.TemporaryDirectory() as reference_root:
        with tarfile.open(reference_tarball) as archive:
            archive.extractall(reference_root, filter="data")
        for kind in ("teleseismic", "regional"):
            new_folder, reference_folder = data_folder(event_folder, kind), data_folder(reference_root, kind)
            print(f"{kind} records")
            if new_folder is None or reference_folder is None:
                print(f"  missing data folder (new: {new_folder}, reference: {reference_folder})")
                continue
            compare_folders(new_folder, reference_folder)
            if kind == "teleseismic":
                for label, folder in (("new", new_folder), ("reference", reference_folder)):
                    unrotated = unrotated_horizontal_records(folder)
                    print(f"  unrotated radial/transverse records in {label} ({len(unrotated)}): {' '.join(unrotated)}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
