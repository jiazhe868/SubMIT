"""Pool post-burn-in MCMC samples across chains WITH a convergence gate.

Chains whose best misfit plateaus above the global best never reached the
stationary distribution; their samples are not posterior draws and would inject
spurious mass into near-zero-probability basins (the acceptance temperature
2*0.05^2*like_min means a chain stuck 10% above the global best has relative
density ~exp(-20)). Gate: include a chain iff chain_best <= (1+tol)*global_best.
Included chains are exchangeable draws from the same posterior -> pooled
UNWEIGHTED (weighting by misfit would double-count the likelihood).

usage: python pool_chains.py [tol]   (run in an inv_* dir; default tol 0.05)
writes pooled_ensemble.dat (gated concatenation of *chain.dat) + a gate report.
"""
import glob
import sys

tol = float(sys.argv[1]) if len(sys.argv) > 1 else 0.05
best = {}
for f in glob.glob("*best.dat"):
    c = f.replace("best.dat", "")
    with open(f) as fp:
        line = fp.readline().split()
        if len(line) > 2:
            best[c] = float(line[2])
gbest = min(best.values())
keep = {c for c, b in best.items() if b <= (1 + tol) * gbest}
nlines = 0
with open("pooled_ensemble.dat", "w") as out:
    for c in sorted(keep):
        with open(c + "chain.dat") as fp:
            for line in fp:
                out.write(line)
                nlines += 1
print(f"pool_chains: global best {gbest:.4f}; gate <= {(1+tol)*gbest:.4f}; "
      f"kept {len(keep)}/{len(best)} chains, {nlines} samples -> pooled_ensemble.dat")
drop = sorted((b, c) for c, b in best.items() if c not in keep)
if drop:
    print("excluded (chain: best):", ", ".join(f"{c}:{b:.3f}" for b, c in drop[:10]),
          "..." if len(drop) > 10 else "")
