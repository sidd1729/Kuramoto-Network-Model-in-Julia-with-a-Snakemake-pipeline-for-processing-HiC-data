import cooltools
import cooler
import bioframe
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import warnings
import pyarrow.feather as feather
import sys
import os
from tqdm import tqdm
from cytoolz import merge

#Function taken straight from cooltools tutorial
def saddleplot(
    track,
    saddledata,
    n_bins,
    vrange=None,
    qrange=(0.0, 1.0),
    cmap="coolwarm",
    scale="log",
    vmin=0.5,
    vmax=2,
    color=None,
    title=None,
    xlabel=None,
    ylabel=None,
    clabel=None,
    fig=None,
    fig_kws=None,
    heatmap_kws=None,
    margin_kws=None,
    cbar_kws=None,
    subplot_spec=None,
):
    """
    Generate a saddle plot.
    Parameters
    ----------
    track : pd.DataFrame
        See cooltools.digitize() for details.
    saddledata : 2D array-like
        Saddle matrix produced by `make_saddle`. It will include 2 flanking
        rows/columns for outlier signal values, thus the shape should be
        `(n+2, n+2)`.
    cmap : str or matplotlib colormap
        Colormap to use for plotting the saddle heatmap
    scale : str
        Color scaling to use for plotting the saddle heatmap: log or linear
    vmin, vmax : float
        Value limits for coloring the saddle heatmap
    color : matplotlib color value
        Face color for margin bar plots
    fig : matplotlib Figure, optional
        Specified figure to plot on. A new figure is created if none is
        provided.
    fig_kws : dict, optional
        Passed on to `plt.Figure()`
    heatmap_kws : dict, optional
        Passed on to `ax.imshow()`
    margin_kws : dict, optional
        Passed on to `ax.bar()` and `ax.barh()`
    cbar_kws : dict, optional
        Passed on to `plt.colorbar()`
    subplot_spec : GridSpec object
        Specify a subregion of a figure to using a GridSpec.
    Returns
    -------
    Dictionary of axes objects.
    """

#     warnings.warn(
#         "Generating a saddleplot will be deprecated in future versions, "
#         + "please see https://github.com/open2c_examples for examples on how to plot saddles.",
#         DeprecationWarning,
#     )

    from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
    from matplotlib.colors import Normalize, LogNorm
    from matplotlib import ticker
    import matplotlib.pyplot as plt

    class MinOneMaxFormatter(ticker.LogFormatter):
        def set_locs(self, locs=None):
            self._sublabels = set([vmin % 10 * 10, vmax % 10, 1])

        def __call__(self, x, pos=None):
            if x not in [vmin, 1, vmax]:
                return ""
            else:
                return "{x:g}".format(x=x)

    track_value_col = track.columns[3]
    track_values = track[track_value_col].values

    digitized_track, binedges = cooltools.digitize(
        track, n_bins, vrange=vrange, qrange=qrange
    )
    x = digitized_track[digitized_track.columns[3]].values.astype(int).copy()
    x = x[(x > -1) & (x < len(binedges) + 1)]

    # Old version
    # hist = np.bincount(x, minlength=len(binedges) + 1)

    groupmean = track[track.columns[3]].groupby(digitized_track[digitized_track.columns[3]]).mean()

    if qrange is not None:
        lo, hi = qrange
        binedges = np.linspace(lo, hi, n_bins + 1)

    # Barplot of mean values and saddledata are flanked by outlier bins
    n = saddledata.shape[0]
    X, Y = np.meshgrid(binedges, binedges)
    C = saddledata
    if (n - n_bins) == 2:
        C = C[1:-1, 1:-1]
        groupmean = groupmean[1:-1]

    # Layout
    if subplot_spec is not None:
        GridSpec = partial(GridSpecFromSubplotSpec, subplot_spec=subplot_spec)
    grid = {}
    gs = GridSpec(
        nrows=3,
        ncols=3,
        width_ratios=[0.2, 1, 0.1],
        height_ratios=[0.2, 1, 0.1],
        wspace=0.05,
        hspace=0.05,
    )

    # Figure
    if fig is None:
        fig_kws_default = dict(figsize=(5, 5))
        fig_kws = merge(fig_kws_default, fig_kws if fig_kws is not None else {})
        fig = plt.figure(**fig_kws)

    # Heatmap
    if scale == "log":
        norm = LogNorm(vmin=vmin, vmax=vmax)
    elif scale == "linear":
        norm = Normalize(vmin=vmin, vmax=vmax)
    else:
        raise ValueError("Only linear and log color scaling is supported")

    grid["ax_heatmap"] = ax = plt.subplot(gs[4])
    heatmap_kws_default = dict(cmap="coolwarm", rasterized=True)
    heatmap_kws = merge(
        heatmap_kws_default, heatmap_kws if heatmap_kws is not None else {}
    )
    img = ax.pcolormesh(X, Y, C, norm=norm, **heatmap_kws)
    plt.gca().yaxis.set_visible(False)

    # Margins
    margin_kws_default = dict(edgecolor="k", facecolor=color, linewidth=1)
    margin_kws = merge(margin_kws_default, margin_kws if margin_kws is not None else {})
    # left margin hist
    grid["ax_margin_y"] = plt.subplot(gs[3], sharey=grid["ax_heatmap"])

    plt.barh(
        binedges, height=1/len(binedges), width=groupmean, align="edge", **margin_kws
    )

    plt.xlim(plt.xlim()[1], plt.xlim()[0])  # fliplr
    plt.ylim(hi, lo)
    plt.gca().spines["top"].set_visible(False)
    plt.gca().spines["bottom"].set_visible(False)
    plt.gca().spines["left"].set_visible(False)
    plt.gca().xaxis.set_visible(False)
    # top margin hist
    grid["ax_margin_x"] = plt.subplot(gs[1], sharex=grid["ax_heatmap"])

    plt.bar(
        binedges, width=1/len(binedges), height=groupmean, align="edge", **margin_kws
    )

    plt.xlim(lo, hi)
    # plt.ylim(plt.ylim())  # correct
    plt.gca().spines["top"].set_visible(False)
    plt.gca().spines["right"].set_visible(False)
    plt.gca().spines["left"].set_visible(False)
    plt.gca().xaxis.set_visible(False)
    plt.gca().yaxis.set_visible(False)

#     # Colorbar
    grid["ax_cbar"] = plt.subplot(gs[5])
    cbar_kws_default = dict(fraction=0.8, label=clabel or "")
    cbar_kws = merge(cbar_kws_default, cbar_kws if cbar_kws is not None else {})
    if scale == "linear" and vmin is not None and vmax is not None:
        grid["ax_cbar"] = cb = plt.colorbar(img, **cbar_kws)
        # cb.set_ticks(np.arange(vmin, vmax + 0.001, 0.5))
        # # do linspace between vmin and vmax of 5 segments and trunc to 1 decimal:
        decimal = 10
        nsegments = 5
        cd_ticks = np.trunc(np.linspace(vmin, vmax, nsegments) * decimal) / decimal
        cb.set_ticks(cd_ticks)
    else:
        print('cbar')

        cb = plt.colorbar(img, format=MinOneMaxFormatter(), cax=grid["ax_cbar"], **cbar_kws)
        cb.ax.yaxis.set_minor_formatter(MinOneMaxFormatter())

    # extra settings
    grid["ax_heatmap"].set_xlim(lo, hi)
    grid["ax_heatmap"].set_ylim(hi, lo)
    grid['ax_heatmap'].grid(False)
    if title is not None:
        grid["ax_margin_x"].set_title(title)
    if xlabel is not None:
        grid["ax_heatmap"].set_xlabel(xlabel)
    if ylabel is not None:
        grid["ax_margin_y"].set_ylabel(ylabel)

    return grid

#"/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/genome/danRer11.fa"
#System arguments are 1: Cool file path name, 2: Resolution, 3: fasta file path 4:base PC output path , 5: base saddle output path
clr = cooler.Cooler(f"{sys.argv[1]}::/resolutions/{sys.argv[2]}")
bins = clr.bins()[:]
dr_11_genome = bioframe.load_fasta(sys.argv[3])
#print(dr_11_genome.items())
dr_11_genome_renamed = {chrom.replace('chr', '') if chrom.startswith('chr') else chrom: seq for chrom, seq in dr_11_genome.items()}
#print(dr_11_genome.items())
#print(bins)
#print(type(bins['chrom']),type(map(lambda x: f'chr{x}',bins['chrom'])))
#def chrom_converter(x):
    #return f'chr{x}'
#bins['chrom'] = bins['chrom'].map(chrom_converter)
#print(bins)
gc_cov = bioframe.frac_gc(bins[['chrom','start','end']],dr_11_genome_renamed) #genome renamed to be used
#print(gc_cov)
#chrom_name_values = list(map(chrom_converter,clr.chromnames))
#print(type(clr.chromnames))
view_df = pd.DataFrame({'chrom': clr.chromnames,
                        'start': 0,
                        'end': clr.chromsizes.values,
                        'name': clr.chromnames})
#print(view_df)
cis_eigs = cooltools.eigs_cis(clr,gc_cov,view_df=view_df,n_eigs=3)
eigenvector_track = cis_eigs[1][['chrom','start',
                                 'end','E1']]
eigenvector_track_all = cis_eigs[1][['chrom','start','end','E1','E2','E3']]
#print(cis_eigs[0])
#print(cis_eigs[1])

cvd = cooltools.expected_cis(clr=clr,view_df=view_df)
Q_LO = 0.025
Q_HI = 0.975
N_GROUPS = 38
##plt.rcParams['font.size'] = 8
#for i in tqdm(clr.chromnames):
    #print(i,type(i))
    #eigenvector_track_use_all = eigenvector_track_all[eigenvector_track_all['chrom'] == i]
    #eigenvector_track_use_saddle = eigenvector_track[eigenvector_track['chrom'] == i ]
eigenvector_track_use = pd.DataFrame(eigenvector_track_all['E1'])
eigenvector_track_use_2 = pd.DataFrame(eigenvector_track_all['E2'])
eigenvector_track_use_3 = pd.DataFrame(eigenvector_track_all['E3'])
feather.write_feather(eigenvector_track_use,f"{sys.argv[4]}/{os.path.splitext(os.path.basename(sys.argv[1]))[0]}_{sys.argv[2]}_PC1.feather")
feather.write_feather(eigenvector_track_use_2,f"{sys.argv[4]}/{os.path.splitext(os.path.basename(sys.argv[1]))[0]}_{sys.argv[2]}_PC2.feather")
feather.write_feather(eigenvector_track_use_3,f"{sys.argv[4]}/{os.path.splitext(os.path.basename(sys.argv[1]))[0]}_{sys.argv[2]}_PC3.feather")
    #print(f"eigenvector track used is \n {eigenvector_track_use} \n")
    #view_df_use = view_df[view_df['chrom'] == i].copy()
    #print(f"view df use is \n {view_df_use} \n")
    #view_df_use = view_df_use.reset_index(drop=True)#sort_values(['chrom','start','end'])
    #print(f"view df use is \n {view_df_use} \n")
    #cvd_use = cvd[cvd['region1'] == i]
    #print(f"cvd use is \n {cvd_use} \n")
interaction_sum, interaction_count = cooltools.saddle(
    clr,cvd,eigenvector_track,'cis',n_bins = N_GROUPS, qrange = (Q_LO,Q_HI),view_df=view_df)
saddleplot(eigenvector_track,interaction_sum/interaction_count, N_GROUPS,qrange = (Q_LO,Q_HI), cbar_kws ={'label':f"Saddle plot of {os.path.splitext(os.path.basename(sys.argv[1]))[0]} at {sys.argv[2]} resolution"})
    #plt.xlabel(f"Saddle plot for chromosome {i} of {os.path.splitext(os.path.basename(sys.argv[1]))[0]} at {sys.argv[2]} resolution", fontsize = 6)
plt.savefig(f"{sys.argv[5]}/{os.path.splitext(os.path.basename(sys.argv[1]))[0]}_{sys.argv[2]}_saddleplot.png")
plt.show()
    #plt.close()