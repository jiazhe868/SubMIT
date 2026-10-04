#!/usr/bin/env python3
"""Compare a reproduced example run with the published reference result.

usage: python tools/compare_with_reference.py examples/<name> work/<name>/IRIS
Prints the L-curve side by side, the selected number of subevents, and for the
selected model each subevent's timing, location, duration, magnitude and
moment-tensor similarity to the published one.
"""
import os
import re
import sys
import numpy as np


def lcurve(path):
    """{n: best misfit} and the suggested n from an lcurve.txt."""
    vals, sel = {}, None
    if not os.path.exists(path):
        return vals, sel
    for line in open(path):
        m = re.match(r"\s*(\d+)\s+([0-9.]+)\s", line)
        if m:
            vals[int(m.group(1))] = float(m.group(2))
        m = re.search(r"SUGGESTED.*:\s*(\d+)", line)
        if m:
            sel = int(m.group(1))
    return vals, sel


def tensors(path):
    """rows of (Mxx Mxy Mxz Myy Myz Mzz) x 1e27 dyne-cm from fm.dat"""
    return np.loadtxt(path, ndmin=2)[:, 1:7]


def full(m):
    xx, xy, xz, yy, yz, zz = m
    return np.array([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]])


def mw(m):
    M = full(m)
    m0 = np.sqrt((M * M).sum() / 2) * 1e27
    return 2 / 3 * np.log10(m0) - 10.7


def similarity(a, b):
    """normalized tensor inner product: 1 identical mechanism, -1 opposite"""
    A, B = full(a), full(b)
    return float((A * B).sum() / np.linalg.norm(A) / np.linalg.norm(B))


def main(example, iris):
    ref_l, ref_sel = lcurve(os.path.join(example, "reference", "lcurve.txt"))
    new_l, new_sel = lcurve(os.path.join(iris, "lcurve.txt"))
    print("\nL-curve (best penalized misfit for each number of subevents)")
    print(f"{'n':>3} {'published':>10} {'reproduced':>11} {'diff':>8}")
    for n in sorted(set(ref_l) | set(new_l)):
        a, b = ref_l.get(n), new_l.get(n)
        d = f"{100 * (b - a) / a:+.1f}%" if a and b else ""
        print(f"{n:>3} {a if a else '-':>10} {b if b else '-':>11} {d:>8}")
    print(f"\nselected number of subevents: published {ref_sel}, reproduced {new_sel}")

    ev = open(os.path.join(example, "event_list.dat")).read().split()[0]
    n = new_sel or ref_sel
    fwd = os.path.join(iris, f"fwd_{ev}_{n}sub")
    if n != ref_sel or not os.path.exists(os.path.join(fwd, "fm.dat")):
        print("(subevent-level comparison needs the same selected n and a finished step 4)")
        return
    ra, rb = tensors(os.path.join(example, "reference", "fm.dat")), tensors(os.path.join(fwd, "fm.dat"))
    ma = np.loadtxt(os.path.join(example, "reference", "Input.model"), ndmin=2)
    mb = np.loadtxt(os.path.join(fwd, "Input.model"), ndmin=2)
    print(f"\nselected {n}-subevent model (published -> reproduced)")
    print(f"{'':>3} {'time (s)':>14} {'x,y (km)':>22} {'depth (km)':>14} {'dur (s)':>12} {'Mw':>11} {'MT sim':>7}")
    for i in range(len(ra)):
        print(f"E{i + 1:<2} {ma[i, 0]:6.1f}->{mb[i, 0]:<6.1f} "
              f"{ma[i, 1]:5.0f},{ma[i, 2]:<4.0f}->{mb[i, 1]:5.0f},{mb[i, 2]:<4.0f} "
              f"{ma[i, 6]:6.1f}->{mb[i, 6]:<6.1f} {ma[i, 3]:5.1f}->{mb[i, 3]:<5.1f} "
              f"{mw(ra[i]):4.2f}->{mw(rb[i]):<4.2f} {similarity(ra[i], rb[i]):+.2f}")
    print(f"total Mw: published {mw(ra.sum(0)):.2f}, reproduced {mw(rb.sum(0)):.2f}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
