#!/usr/bin/env python3
import glob
import re
import os
import numpy as np
import matplotlib.pyplot as plt

def process_and_plot(top_fraction=0.25, save_pdf='histoplot_py.pdf', save_samples='allsamples.dat'):
    # 1. Discover all *best.dat files and find matching *chain.dat files
    best_files = glob.glob('*best.dat')
    
    chain_scores = []
    
    for bfile in best_files:
        # Match chain prefix
        match = re.match(r'^(.*)best\.dat$', bfile)
        if not match:
            continue
        prefix = match.group(1)
        cfile = f"{prefix}chain.dat"
        
        if not os.path.exists(cfile):
            continue
            
        # Read the last non-empty line of *best.dat to get the min misfit score (1st numerical parameter)
        with open(bfile, 'r') as f:
            lines = [line.strip() for line in f if line.strip()]
        
        if not lines:
            continue
            
        last_line = lines[-1]
        # Remove '<neq' if present
        clean_line = last_line.replace('<neq', '')
        tokens = clean_line.split()
        
        # In best.dat, tokens are: [num_subevents, val1, val2, val3, ...]
        # 2nd token (index 1) is the first value (misfit score)
        if len(tokens) >= 2:
            score = float(tokens[1])
            chain_scores.append((score, prefix, cfile))

    num_chains = len(chain_scores)
    if num_chains == 0:
        print("No matching *best.dat and *chain.dat files found.")
        return

    # Sort chains by score ascending (lowest misfit score is best)
    chain_scores.sort(key=lambda x: x[0])
    
    num_selected = max(1, int(np.round(num_chains * top_fraction)))
    selected_chains = chain_scores[:num_selected]
    
    print(f"Total chains found: {num_chains}")
    print(f"Selecting best {top_fraction*100:.0f}% ({num_selected} chains):")
    for score, prefix, cfile in selected_chains:
        print(f"  Chain: {prefix} | Best score: {score}")

    # 2. Read selected *chain.dat files and collect parameter samples
    # Based on prepare.sh:
    # columns extracted (1-indexed):
    # 2, 3, 4, 5, 6, 9, 11, 12, 13, 14, 17, 19, 20, 21, 22, 25, 27, 28, 29, 30, 33, 2
    # 0-indexed indices (when <neq is removed):
    # Columns are built dynamically from the subevent count (the original list
    # was hardcoded for exactly 4 subevents and silently dropped every row for
    # any other nsub). Token layout after removing '<neq':
    # 0=neq, then per subevent k (0-based): misfit,cen,x,y,dura,vr,theta,depth
    # at 8k+1..8k+8. Extracted: [misfit, (cen,x,y,dura,depth) x nsub].
    # If the convergence-gated pool exists (pool_chains.py), use it so ALL
    # posterior products share one ensemble definition; else top-fraction chains
    import os as _os
    if _os.path.exists('pooled_ensemble.dat'):
        print('Using convergence-gated pooled_ensemble.dat')
        selected_chains = [(0.0, 'pooled', 'pooled_ensemble.dat')]
    all_samples = []
    nsub = None
    for _, prefix, cfile in selected_chains:
        with open(cfile, 'r') as f:
            for line in f:
                line_str = line.strip()
                if not line_str:
                    continue
                tokens = line_str.replace('<neq', '').split()
                if nsub is None:
                    nsub = int(float(tokens[0]))
                    cols_to_extract = [1]
                    for k in range(nsub):
                        cols_to_extract += [8*k+2, 8*k+3, 8*k+4, 8*k+5, 8*k+8]
                if len(tokens) > max(cols_to_extract):
                    row = [float(tokens[idx]) for idx in cols_to_extract]
                    all_samples.append(row)

    samples = np.array(all_samples)
    # enforce the causal labeling convention per SAMPLE: subevent blocks are
    # sorted by centroid time, so E2 is by definition the 2nd pulse in time
    # (ensembles from runs predating the sampler's hard ordering rule mix
    # labels - e.g. E2 samples later than E4 - which smears every marginal)
    if samples.size and samples.shape[1] >= 1 + 5 * nsub:
        for r in range(samples.shape[0]):
            blocks = samples[r, 1:1 + 5 * nsub].reshape(nsub, 5)
            samples[r, 1:1 + 5 * nsub] = blocks[np.argsort(blocks[:, 0])].ravel()
    print(f"Total samples collected: {samples.shape[0]}")

    if save_samples:
        np.savetxt(save_samples, samples, fmt='%.6f')
        print(f"Saved extracted samples to {save_samples}")

    # 3. Plotting histograms for 4 subevents
    # Each subevent has 5 columns in samples array corresponding to:
    # 1. Centroid time (s): col index (k-1)*5 + 1
    # 2. Duration (s): col index (k-1)*5 + 4
    # 3. WE location (km): col index (k-1)*5 + 2
    # 4. NS location (km): col index (k-1)*5 + 3
    # 5. Depth (km): col index (k-1)*5 + 5

    # xlim=None -> data-driven limits (hardcoded event-specific limits made the
    # script single-event); sample column layout: 0=misfit, then 5 per subevent
    param_configs = [
        {'title': 'Centroid time (s)', 'offset': 1, 'bin_width': 0.5, 'xlim': None},
        {'title': 'Duration (s)',      'offset': 4, 'bin_width': 0.3, 'xlim': None},
        {'title': 'WE location (km)',  'offset': 2, 'bin_width': 3.0, 'xlim': None},
        {'title': 'NS location (km)',  'offset': 3, 'bin_width': 3.0, 'xlim': None},
        {'title': 'Depth (km)',        'offset': 5, 'bin_width': 1.0, 'xlim': None},
    ]

    plt.rcParams.update({"font.size": 12})
    fig, axes = plt.subplots(max(nsub, 2), 5, figsize=(15, 2.8 * max(nsub, 2)))
    plt.subplots_adjust(wspace=0.3, hspace=0.4)

    for k in range(nsub):
        for col_idx, config in enumerate(param_configs):
            ax = axes[k, col_idx]
            data_col = (k * 5) + config['offset']
            data = samples[:, data_col]

            xlim = config['xlim']
            if xlim is None:
                # shared per-column range across ALL subevents (comparability)
                colall = np.concatenate([samples[:, (kk * 5) + config['offset']]
                                         for kk in range(nsub)])
                pad = 0.06 * (colall.max() - colall.min() + 1e-9)
                xlim = (colall.min() - pad, colall.max() + pad)
            # ~25 bins over the panel's OWN posterior span (binning over the
            # shared column range left narrow marginals - e.g. E1 centroid -
            # with a single bar); the shared xlim below keeps comparability
            dspan = data.max() - data.min()
            # pinned/degenerate parameters (e.g. E1 x=y=0) would produce a
            # 1e6-density single-bin spike: annotate instead of a histogram
            if dspan < 1e-3:
                ax.axvline(data[0], color='steelblue', lw=2)
                ax.annotate(f"fixed at {data[0]:.2f}", (0.5, 0.5),
                            xycoords='axes fraction', ha='center', fontsize=10)
                ax.set_xlim(xlim)
                if col_idx == 0:
                    ax.set_ylabel(f"E{k+1} prob.")
                if k == nsub - 1:
                    ax.set_xlabel(config['title'])
                continue
            # ~25 bins, but snapped to a NICE width (1/2/2.5/5 x 10^k) and
            # aligned to multiples of that width: arbitrary widths beat
            # against the sampler's discreteness and render uneven combs
            # width: ~25 bins over the marginal's own span, but never finer
            # than 1/120 of the SHARED axis (a 1-s marginal on a 15-s axis
            # otherwise renders as a solid sub-pixel blob)
            raw = max(dspan / 25.0, (xlim[1] - xlim[0]) / 120.0)
            mag = 10.0 ** np.floor(np.log10(raw))
            bin_width = min((v * mag for v in (1.0, 2.0, 2.5, 5.0, 10.0)),
                            key=lambda v: abs(v - raw) if v >= raw else 1e9)
            lo = np.floor(data.min() / bin_width) * bin_width
            hi = np.ceil(data.max() / bin_width) * bin_width + bin_width
            bins = np.arange(lo, hi + 0.5 * bin_width, bin_width)

            # black edges become solid ink when bars are thin on screen
            ec = 'black' if bin_width >= (xlim[1] - xlim[0]) / 60.0 else 'none'
            counts, bin_edges, patches = ax.hist(
                data, bins=bins, density=True, 
                edgecolor=ec, facecolor='skyblue', alpha=0.85
            )

            # Compute 2.5% to 97.5% cumulative density bounds (matching MATLAB script)
            total_counts = len(data)
            hist_counts, _ = np.histogram(data, bins=bins)
            bin_probs = hist_counts / total_counts if total_counts > 0 else np.zeros_like(hist_counts)
            cum_probs = np.cumsum(bin_probs)
            
            # Find trim1 (where cum prob > 0.025 from left)
            idx1 = np.searchsorted(cum_probs, 0.025)
            if idx1 < len(bin_edges) - 1:
                trim1 = (bin_edges[idx1] + bin_edges[idx1 + 1]) / 2.0
            else:
                trim1 = bin_edges[0]

            # Find trim2 (where cum prob > 0.025 from right)
            cum_probs_rev = np.cumsum(bin_probs[::-1])
            idx2 = np.searchsorted(cum_probs_rev, 0.025)
            if idx2 < len(bin_edges) - 1:
                trim2_idx = len(bin_edges) - 1 - idx2
                trim2 = (bin_edges[trim2_idx - 1] + bin_edges[trim2_idx]) / 2.0
            else:
                trim2 = bin_edges[-1]

            # Draw interval line at y=0
            ax.plot([trim1, trim2], [0, 0], '-ks', linewidth=2, markersize=5, markerfacecolor='gray')

            ax.set_xlim(xlim)
            
            # Labels
            if col_idx == 0:
                ax.set_ylabel(f"E{k+1} prob.", fontsize=12)
            if k == nsub - 1:   # bottom row (was hardcoded k==3 = 4-subevent only)
                ax.set_xlabel(config['title'], fontsize=11)

    plt.tight_layout()
    plt.savefig(save_pdf, dpi=300)
    print(f"Saved histogram plot to {save_pdf}")

if __name__ == '__main__':
    process_and_plot()
