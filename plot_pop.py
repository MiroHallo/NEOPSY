#!/usr/bin/env python3
# =============================================================================
# PLOT RESULTS OF THE MULTIZONAL TRANSDIMENSIONAL INVERSION (NEOPSY)
#
# Author: Miroslav HALLO, Kyoto University
# E-mail: hallo.miroslav.2a@kyoto-u.ac.jp
# Revision 2026/05: First version
# Tested with: Python 3.12.3, Matplotlib 3.10.9, NumPy 2.4.5, SciPy 1.17.1
# Method:
# Hallo, M., Imperatori, W., Panzera, F., Fäh, D. (2021). Joint multizonal
#      transdimensional Bayesian inversion of surface wave dispersion and
#      ellipticity curves for local near-surface imaging, Geophys. J. Int.,
#      226 (1), 627-659. https://doi.org/10.1093/gji/ggab116
#
# Copyright (C) 2026 Miroslav Hallo
#
# This program is published under the GNU General Public License (GNU GPL).
#
# This program is free software: you can modify it and/or redistribute it
# or any derivative version under the terms of the GNU General Public
# License as published by the Free Software Foundation, either version 3
# of the License, or (at your option) any later version.
#
# This code is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY. We would like to kindly ask you to acknowledge the authors
# and don't remove their names from the code.
#
# You should have received a copy of the GNU General Public License along
# with this program. If not, see <http://www.gnu.org/licenses/>.
#
# =============================================================================

from types import SimpleNamespace
from typing import List, Optional
from pathlib import Path

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.ticker import ScalarFormatter
import numpy as np
from scipy.signal import savgol_filter, find_peaks

from lib import (fscale, dusk, read_pop_headers,
                 read_input_data, read_model_data, read_velocity_model,
                 read_synthetic_model, read_pop_binary_data,
                 read_qwl_headers, read_model_qwl, compute_qwl)


# =============================================================================
# CONSTANTS & PATHS
# =============================================================================

# Get the directory where the current script is located
PROJ_ROOT = Path(__file__).parent.resolve()

# NEOPSY working directory with inversion results
INV_NAME = 'inv'
INV_PATH = PROJ_ROOT / INV_NAME

DATAFILE = INV_PATH / 'in_data.txt'  # Data of DC and ELL curves (file)
MAX_PREFIX = INV_PATH / 'out_modelML'  # ML model (prefix of input files)
MAP_PREFIX = INV_PATH / 'out_modelMAP-SP'  # MAP model (prefix of input files)
POP_PREFIX = INV_PATH / 'out_pop'  # Ensemble stats (prefix of input files)
OUT_FOLDER = INV_PATH / 'results'  # Path to save output figures and results

# Reference rock velocity model for amplification (file)
REF_FILE = PROJ_ROOT / 'data' / 'vs_ref_Swiss.ascii'


# =============================================================================
# PLOTTING PARAMETERS
# =============================================================================

OUT_PREFIX = 'PDF'  # Prefix of output files (Site/Station/Code/Version)
N_FITS = 300  # Number of random data fits from the solution ensemble to plot
PLOT_DEPTH_MAX = 9000  # Maximal depth to be plotted [m]

SYNT_PLOT = False  # Read and plot synthetic model (synthetic test, True/False)
SYNT_FILE = PROJ_ROOT / 'data' / 'data.model'  # Target model (if SYNT_PLOT)

PLOT_IN_DPD = True  # Plot ML and MAP models over PDF (True/False)
PLOT_LEG = True  # Plot Legends (True/False)
FIGURE_REN = 2  # Save figures (0=NO-show, 1=YES-show, 2=YES-close)

USE_DATA_THR = False  # Use joint axes for the data fit plots (True/False)
SLOW_THR = [0.5, 8.0]  # Joint slowness axes threshold [ms/m] (if USE_DATA_THR)
FREQ_THR = [1.0, 30.0]  # Joint freqency axes threshold  [Hz] (if USE_DATA_THR)
MAX_QWL_PROP = 0.5  # Upper threshold of probability for the QWL plots


# =============================================================================
# FUNCTIONS
# =============================================================================

# Plot data fit (slowness or velocity)
def plot_data_fit(indata: SimpleNamespace, maxdata: SimpleNamespace,
                  mapdata: SimpleNamespace, syndata: np.ndarray,
                  out_path: Path, plot_type: str, use_data_thr: bool,
                  slow_thr: List[float], freq_thr: List[float],
                  plot_leg: bool, fig_ren: int) -> None:
    """
    Plots the data fits (observed vs. synthetic curves).

    Args:
        indata (SimpleNamespace): Input (measured) curves.
        maxdata (SimpleNamespace): ML model curves.
        mapdata (SimpleNamespace): MAP model curves.
        syndata (np.ndarray): Ensemble data curves (N_FITS).
        out_path (Path): Path to output directory for results.
        plot_type (str): 'slowness' (ms/m) or plot_type='velocity' (m/s).
        use_data_thr (bool): Use joint axes (True/False).
        slow_thr (List[float]): Joint slowness axes threshold [ms/m].
        freq_thr (List[float]): Joint freqency axes threshold  [Hz].
        plot_leg (bool): Plot Legends (True/False).
        fig_ren (int): Save figures (0=NO-show, 1=YES-show, 2=YES-close).
    Returns:
        None.
    """
    max_mode = len(indata.f_n)
    max_mode_a = int(np.sum(indata.f_n > 0))

    if max_mode_a == 0:
        print(f"[*] No active modes to plot for {plot_type}.")
        return

    # Determine layout grid and figure sizing
    if max_mode_a < 4:
        p_cols = max_mode_a
        p_rows = 1
    elif max_mode_a < 7:
        p_cols = 3
        p_rows = 2
    else:
        p_cols = 3
        p_rows = 3
    px = (p_cols * 7) / 2.54
    py = (p_rows * 7) / 2.54

    # Start figure
    plt.rcParams.update({
        'font.size': 8,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'axes.labelsize': 9,
        'axes.titlesize': 9,
        'legend.fontsize': 7,
    })
    fig, axes = plt.subplots(p_rows, p_cols, figsize=(px, py), squeeze=False)
    fig.patch.set_facecolor('white')
    axes = axes.flatten()

    mode_labels = {
        0: 'Rayleigh (fundamental)', 1: 'Rayleigh (1st higher)',
        2: 'Rayleigh (2nd higher)', 3: 'Rayleigh (3rd higher)',
        4: 'Love (fundamental)', 5: 'Love (1st higher)',
        6: 'Love (2nd higher)', 7: 'Love (3rd higher)',
        8: 'Rayleigh wave ellipticity', 9: 'Rayleigh wave ellipticity angle'
    }

    m_a = 0
    for m in range(max_mode):
        n_pts = indata.f_n[m]
        if n_pts <= 0 or m > 9:
            continue

        ax = axes[m_a]
        m_a += 1
        ftxt = mode_labels[m]
        ci = int(np.sum(indata.f_n[:m]))

        mask = indata.okf[:n_pts, m]
        f_in = indata.data[:n_pts, 0, m]
        f_in_masked = f_in[mask]

        amin, amax, fticks = fscale(np.min(f_in_masked), np.max(f_in_masked))
        frange = [amin, amax]
        lh = [None] * 5

        # -------------------------------------------
        # MODES 1-8: Dispersion curves (Slowness or Velocity)
        if m <= 7:
            # Prepare curves
            if plot_type == 'slowness':
                ens_curves = syndata[:, 0, ci:ci+n_pts].T * 1000.0
                d_obs = indata.data[:n_pts, 1, m] * 1000.0
                d_obs_s = indata.data_s[:n_pts, m] * 1000.0
                y_max = maxdata.data[:n_pts, 1, m] * 1000.0
                y_map = mapdata.data[:n_pts, 1, m] * 1000.0
                y_label = 'Slowness (ms/m)'
                y_err = d_obs_s[mask]
            else:  # velocity
                ens_curves = 1.0 / (syndata[:, 0, ci:ci+n_pts].T)
                d_obs = 1.0 / indata.data[:n_pts, 1, m]
                s_val = indata.data[:n_pts, 1, m]
                s_err = indata.data_s[:n_pts, m]
                err_low = (1.0 / s_val) - (1.0 / (s_val + s_err))
                err_high = -(1.0 / s_val) + (1.0 / (s_val - s_err))
                y_err = np.vstack([err_low[mask], err_high[mask]])
                y_max = 1.0 / maxdata.data[:n_pts, 1, m]
                y_map = 1.0 / mapdata.data[:n_pts, 1, m]
                y_label = 'Velocity (m/s)'

            # Plot curves
            ax.plot(f_in, ens_curves, color=(0.7, 0.7, 0.8), linewidth=0.5,
                    zorder=1)
            lh[4], = ax.plot([], [], color=(0.7, 0.7, 0.8), linewidth=0.5)
            lh[0], = ax.plot(f_in_masked, d_obs[mask],
                             color='k', linestyle=':', linewidth=0.6, zorder=3)
            lh[1] = ax.errorbar(f_in_masked, d_obs[mask], yerr=y_err, fmt='.',
                                color='k', markersize=4, capsize=1,
                                elinewidth=0.6, zorder=4)
            lh[2], = ax.plot(maxdata.data[:n_pts, 0, m], y_max, color='b',
                             zorder=5)
            lh[3], = ax.plot(mapdata.data[:n_pts, 0, m], y_map, color='m',
                             zorder=6)
            ax.set_ylabel(y_label)
            ax.set_xscale('log')

            # Use joint axes
            if use_data_thr:
                if plot_type == 'slowness':
                    ax.set_ylim(slow_thr)
                else:
                    ax.set_ylim([1000.0 / slow_thr[1], 1000.0 / slow_thr[0]])
                ax.set_xlim(freq_thr)
            else:
                ax.set_xlim(frange)
                ax.set_xticks(fticks)

        # -------------------------------------------
        # MODE 9: Ellipticity curve
        elif m == 8:
            ens_curves = 10.0 ** (syndata[:, 0, ci:ci+n_pts].T)
            ax.loglog(f_in, ens_curves, color=(0.7, 0.7, 0.8), linewidth=0.5,
                      zorder=1)
            lh[4], = ax.plot([], [], color=(0.7, 0.7, 0.8), linewidth=0.5)

            d_obs_log = indata.data[:n_pts, 1, m]
            d_obs_s_log = indata.data_s[:n_pts, m]

            y_val = 10.0 ** d_obs_log[mask]
            y_err_lower = y_val - 10.0 ** (d_obs_log[mask] - d_obs_s_log[mask])
            y_err_upper = 10.0 ** (d_obs_log[mask] + d_obs_s_log[mask]) - y_val

            lh[0], = ax.loglog(f_in_masked, y_val, color='k', linestyle=':',
                               linewidth=0.6, zorder=3)
            lh[1] = ax.errorbar(f_in_masked, y_val,
                                yerr=[y_err_lower, y_err_upper], fmt='.',
                                color='k', markersize=4, capsize=1,
                                elinewidth=0.6, zorder=4)

            lh[2], = ax.loglog(maxdata.data[:n_pts, 0, m],
                               10.0 ** maxdata.data[:n_pts, 1, m], color='b',
                               zorder=5)
            lh[3], = ax.loglog(mapdata.data[:n_pts, 0, m],
                               10.0 ** mapdata.data[:n_pts, 1, m], color='m',
                               zorder=6)

            ax.set_ylabel('Ellipticity')
            ax.set_ylim([0.2, 100])
            ax.set_xlim(frange)
            ax.set_xticks(fticks)

        # -------------------------------------------
        # MODE 10: Ellipticity Angle
        elif m == 9:
            def wrapped(x, y_matrix, threshold=90.0):
                diffs = np.abs(np.diff(y_matrix, axis=1))
                jumps = diffs > threshold
                if not np.any(jumps):
                    return x, y_matrix
                y_clean = y_matrix.copy()
                y_clean[:, 1:][jumps] = np.nan
                return x, y_clean

            tmp_ensemble = syndata[:, 0, ci:ci+n_pts]
            _, y_ens_clean = wrapped(f_in, tmp_ensemble)
            ax.plot(f_in, y_ens_clean.T, color=(0.7, 0.7, 0.8), linewidth=0.5,
                    zorder=1)
            lh[4], = ax.plot([], [], color=(0.7, 0.7, 0.8), linewidth=0.5)

            d_obs = indata.data[:n_pts, 1, m]
            _, y_obs_clean = wrapped(f_in_masked, d_obs[mask].reshape(1, -1))
            lh[0], = ax.plot(f_in_masked, y_obs_clean.flatten(), color='k',
                             linestyle=':', linewidth=0.6, zorder=3)

            d_obs_s = indata.data_s[:n_pts, m]
            lh[1] = ax.errorbar(f_in_masked, d_obs[mask], yerr=d_obs_s[mask],
                                fmt='.', color='k', markersize=4, capsize=1,
                                elinewidth=0.6, zorder=4)

            y_max_raw = maxdata.data[:n_pts, 1, m].reshape(1, -1)
            _, y_max_clean = wrapped(maxdata.data[:n_pts, 0, m], y_max_raw)
            lh[2], = ax.plot(maxdata.data[:n_pts, 0, m], y_max_clean.flatten(),
                             color='b', zorder=5)

            y_map_raw = mapdata.data[:n_pts, 1, m].reshape(1, -1)
            _, y_map_clean = wrapped(mapdata.data[:n_pts, 0, m], y_map_raw)
            lh[3], = ax.plot(mapdata.data[:n_pts, 0, m], y_map_clean.flatten(),
                             color='m', zorder=6)

            ax.set_ylabel('Angle (deg)')
            ax.set_ylim([-90, 90])
            ax.set_yticks(np.arange(-90, 91, 30))
            ax.set_xscale('log')
            ax.set_xlim(frange)
            ax.set_xticks(fticks)

        # General subplot tuning
        ax.set_title(ftxt, fontweight='normal')
        ax.set_xlabel('Frequency (Hz)')
        ax.grid(True, which='both', linestyle='--', linewidth=0.3,
                color='gray', alpha=0.3)
        ax.tick_params(axis='both')
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.ticklabel_format(style='plain', axis='x', useOffset=False)

        # Legend
        if m_a == 1 and plot_leg:
            ax.legend([lh[0], lh[1], lh[2], lh[3], lh[4]],
                      ['Data', 'Data errors', 'ML model', 'MAP model',
                       'Predictive dist.'], loc='best')

    for i in range(max_mode_a, len(axes)):
        fig.delaxes(axes[i])

    plt.tight_layout()

    # Save figure
    if fig_ren > 0:
        pdf_out = f"{out_path}_fit_{plot_type}.pdf"
        png_out = f"{out_path}_fit_{plot_type}.png"
        plt.savefig(pdf_out, dpi=600, bbox_inches='tight')
        plt.savefig(png_out, dpi=600, bbox_inches='tight')
        print(f"[+] SUCCESS: Saved datafit to '{pdf_out}' and PNG.")

    if fig_ren > 1:
        plt.close(fig)
    else:
        plt.show()


# -----------------------------------------------------------------------------
# Plot basic histograms (layer interfaces and VR)
def plot_basic_histograms(pop_prefix: Path, out_path: Path,
                          header: SimpleNamespace, maxmodel: SimpleNamespace,
                          mapmodel: SimpleNamespace, fig_ren: int) -> None:
    """
    Reads 1D population statistics, detects major layer interfaces,
    saves results to an ASCII file, and plots basic ensemble histograms.

    Args:
        pop_prefix (Path): Path to ensemble stats (prefix of input files).
        out_path (Path): Path to output directory for results.
        header (SimpleNamespace): POP header metadata.
        maxmodel (SimpleNamespace): ML solution velocity model.
        mapmodel (SimpleNamespace): MAP solution velocity model.
        fig_ren (int): Save figures (0=NO-show, 1=YES-show, 2=YES-close).
    Returns:
        None.
    """
    n_depth = header.n_depth
    n_bins = header.n_bins
    models_count = header.models_count
    d_bins = header.d_bins
    log_d_bins = header.log_d_bins
    all_bins = header.all_bins

    # Read POP data
    dep_path = Path(f"{pop_prefix}_dep1D.txt")
    lay1d_path = Path(f"{pop_prefix}_lay1D.txt")
    vr1d_path = Path(f"{pop_prefix}_vr1D.txt")
    if not (dep_path.exists() and lay1d_path.exists() and vr1d_path.exists()):
        raise FileNotFoundError("One or more 1D POP files are missing.")
    pop_1d = np.loadtxt(dep_path, max_rows=n_depth)
    pop_lay = np.loadtxt(lay1d_path)
    pop_vr = np.loadtxt(vr1d_path)

    # Compute bins
    d_bins_w = d_bins[1] - d_bins[0]
    log_d_bins_w = np.log(log_d_bins[2]) - np.log(log_d_bins[1])
    log_d_bins_c = np.log(log_d_bins[1:n_depth]) + (log_d_bins_w / 2.0)

    vr_bins_w = all_bins[1, 4] - all_bins[0, 4]
    vr_bins_c = all_bins[:n_bins, 4] + (vr_bins_w / 2.0)

    # Find interfaces
    smoo_data = savgol_filter(pop_1d[1:], window_length=9, polyorder=2)
    smoo_am = np.mean(smoo_data)
    locs, properties = find_peaks(smoo_data, distance=9, height=1.5 * smoo_am)
    pks = properties['peak_heights']

    # Save into text file
    ascii_out = f"{out_path}_layers.ascii"
    with open(ascii_out, 'w', encoding='utf-8') as fid:
        fid.write("# The most probable depths of layer interfaces\n")
        fid.write("# Depth[m], Significance(large num. = more significant)\n")
        for k in range(len(locs)):
            depth_val = np.exp(log_d_bins_c[locs[k]])
            significance = pks[k] / smoo_am
            fid.write(f"{depth_val:8.2f} {significance:6.2f}\n")
    print(f"[+] SUCCESS: ASCII saved to '{ascii_out}'.")

    # Start figure
    gs = gridspec.GridSpec(2, 3)
    px = 21 / 2.54
    py = 14 / 2.54
    plt.rcParams.update({
        'font.size': 8,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'axes.labelsize': 9,
        'axes.titlesize': 9,
        'legend.fontsize': 7,
    })
    fig = plt.figure(figsize=(px, py), facecolor='w')

    # Subplot 1: Interfaces (ln(depth))
    ax1 = fig.add_subplot(gs[0:2, 0])
    ax1.barh(log_d_bins_c, pop_1d[1:] / models_count, height=log_d_bins_w,
             color='#666666', edgecolor='none', align='center', zorder=2)
    x_lim_max = ax1.get_xlim()[1]
    for k in range(len(locs)):
        ax1.plot(0.05 * x_lim_max, log_d_bins_c[locs[k]], color='lime',
                 marker='+', markersize=6)
    ax1.set_xlim([0, x_lim_max])
    ax1.set_ylim([np.log(log_d_bins[1]), np.log(log_d_bins[n_depth-1]) +
                  log_d_bins_w])
    ax1.invert_yaxis()
    ax1.set_xlabel('Probability')
    ax1.set_ylabel('ln(depth)')
    ax1.set_title('Interfaces')
    ax1.grid(True, linestyle='--', linewidth=0.3, zorder=1)

    # Subplot 2: Interfaces (Depth m)
    ax2 = fig.add_subplot(gs[0:2, 1])
    for i in range(1, len(pop_1d) - 1):
        height = log_d_bins[i+1] - log_d_bins[i]
        rect = plt.Rectangle((0, log_d_bins[i]), pop_1d[i] / models_count,
                             height, facecolor='#666666', edgecolor='#666666',
                             linewidth=0.1, zorder=2)
        ax2.add_patch(rect)
    ax2.relim()
    ax2.autoscale_view(scalex=True, scaley=False)
    x_lim_max2 = ax2.get_xlim()[1]
    for k in range(len(locs)):
        ax2.plot(0.05 * x_lim_max2, np.exp(log_d_bins_c[locs[k]]),
                 color='lime', marker='+', markersize=6)
    ax2.set_xlim([0, x_lim_max2])
    ax2.set_ylim([0, np.max(d_bins) + d_bins_w])
    ax2.invert_yaxis()
    ax2.set_xlabel('Probability')
    ax2.set_ylabel('Depth (m)')
    ax2.set_title('Interfaces')
    ax2.grid(True, linestyle='--', linewidth=0.3, zorder=1)

    # Subplot 3: Number of layers
    ax3 = fig.add_subplot(gs[0, 2])
    x_lay = np.arange(1, len(pop_lay) + 1)
    ax3.bar(x_lay, pop_lay / models_count, width=1.0, color='#FF6699',
            edgecolor='k', linewidth=0.3, zorder=2)
    y_lim_max3 = ax3.get_ylim()[1]
    info = f" ML: {maxmodel.n_layers} layers\n MAP: {mapmodel.n_layers} layers"
    ax3.text(len(pop_lay), 0.98 * y_lim_max3, info, color='k',
             horizontalalignment='right', verticalalignment='top', fontsize=7)
    ax3.set_xlim([0.5, len(pop_lay) + 0.5])
    ax3.set_xlabel('k layers')
    ax3.set_ylabel('Probability')
    ax3.set_title('Number of layers')
    ax3.grid(True, linestyle='--', linewidth=0.3, zorder=1)

    # Subplot 6: Data variance reduction
    ax6 = fig.add_subplot(gs[1, 2])
    ax6.bar(vr_bins_c, pop_vr, width=vr_bins_w, color='#669999',
            edgecolor='none', zorder=2)
    yl = ax6.get_ylim()
    for x_val in [0, 43.75, 75, 93.75]:
        ax6.plot([x_val, x_val], yl, '-', color='k', linewidth=0.5)
    ax6.text(1, yl[1], 'Fit (within errors) ', color='k', rotation=90,
             horizontalalignment='left', verticalalignment='top', fontsize=6)
    ax6.text(44.75, yl[1], 'Fair fit ', color='k', rotation=90,
             horizontalalignment='left', verticalalignment='top', fontsize=6)
    ax6.text(76, yl[1], 'Good fit ', color='k', rotation=90,
             horizontalalignment='left', verticalalignment='top', fontsize=6)
    ax6.text(94.75, yl[1], 'Perfect fit ', color='k', rotation=90,
             horizontalalignment='left', verticalalignment='top', fontsize=6)
    vr_txt = f" ML: {round(maxmodel.vr)}%\n MAP: {round(mapmodel.vr)}%"
    ax6.text(np.min(all_bins[:n_bins, 4]), 0, vr_txt, color='k',
             horizontalalignment='left', verticalalignment='bottom',
             fontsize=8)
    ax6.set_ylim(yl)
    ax6.set_xlim([np.min(all_bins[:n_bins, 4]), np.max(all_bins[:n_bins, 4]) +
                  vr_bins_w])
    ax6.set_xticks([0, 43.75, 75, 93.75])
    ax6.set_xlabel('(%)')
    ax6.set_ylabel('Number of models')
    ax6.set_title('Data variance reduction')
    ax6.grid(True, linestyle='--', linewidth=0.3, zorder=1)

    plt.tight_layout()

    # Save figure
    if fig_ren > 0:
        pdf_out = f"{out_path}_layers.pdf"
        png_out = f"{out_path}_layers.png"
        plt.savefig(pdf_out, dpi=600, bbox_inches='tight')
        plt.savefig(png_out, dpi=600, bbox_inches='tight')
        print(f"[+] SUCCESS: Saved histograms to '{pdf_out}' and PNG.")

    if fig_ren > 1:
        plt.close(fig)
    else:
        plt.show()


# -----------------------------------------------------------------------------
# Plot posterior marginal PDF
def plot_pdf(pop_prefix: Path, out_path: Path, header: SimpleNamespace,
             maxmodel: SimpleNamespace, mapmodel: SimpleNamespace,
             syntmodel: Optional[SimpleNamespace], plot_in_dpd: bool,
             plot_leg: bool, fig_ren: int, depth_max: float) -> None:
    """
    Plots posterior marginal PDFs and 1D profiles for Vs, Vp, nu, and rho.

    Args:
        pop_prefix (Path): Path to ensemble stats (prefix of input files).
        out_path (Path): Path to output directory for results.
        header (SimpleNamespace): POP header metadata.
        maxmodel (SimpleNamespace): ML solution velocity model.
        mapmodel (SimpleNamespace): MAP solution velocity model.
        syntmodel (SimpleNamespace): Target velocity model (Optional).
        plot_in_dpd: Plot ML and MAP models over PDF (True/False).
        plot_leg: Plot Legends (True/False).
        fig_ren (int): Save figures (0=NO-show, 1=YES-show, 2=YES-close).
        depth_max (float): Maximum depth for the last layer in the plot.
    Returns:
        None.
    """
    n_depth = header.n_depth
    models_count = header.models_count
    d_bins = header.d_bins
    log_d_bins = header.log_d_bins
    all_bins = header.all_bins
    all_bins_p1 = header.all_bins_p1
    d_bins_p1 = header.d_bins_p1

    nmi = np.zeros(4)
    pdf_matrices = []
    ptype_configs = []

    configs = {
        0: {'suffix': 'vs2D.txt', 'label': r'$V_S$',
            'unit': '(m/s)', 'out': 'vs', 'col_idx': 0,
            'max_attr': 'vs', 'map_attr': 'vs', 'syn_attr': 'vs'},
        1: {'suffix': 'vp2D.txt', 'label': r'$V_P$',
            'unit': '(m/s)', 'out': 'vp', 'col_idx': 1,
            'max_attr': 'vp', 'map_attr': 'vp', 'syn_attr': 'vp'},
        2: {'suffix': 'nu2D.txt', 'label': r'$\nu$',
            'unit': '', 'out': 'nu', 'col_idx': 2,
            'max_attr': 'nu', 'map_attr': 'nu', 'syn_attr': 'nu'},
        3: {'suffix': 'rho2D.txt', 'label': r'$\rho$',
            'unit': r'(kg/m$^3$)', 'out': 'rho', 'col_idx': 3,
            'max_attr': 'rho', 'map_attr': 'rho', 'syn_attr': 'rho'}
    }

    # Create DUSK colormap
    color_matrix = dusk()
    cmap_dusk = ListedColormap(color_matrix, name='dusk')

    # Read Depth 1D POP data
    dep_path = Path(f"{pop_prefix}_dep1D.txt")
    if not dep_path.exists():
        raise FileNotFoundError("Depth 1D POP file is missing.")
    pop_1d = np.loadtxt(dep_path, max_rows=n_depth)

    # First run: Read data
    for ptype in range(4):
        cfg = configs[ptype]
        file_2d = Path(f"{pop_prefix}_{cfg['suffix']}")

        if not file_2d.exists():
            print(f"[!] Warning: 2D statistics file missing: {file_2d.name}.")
            pdf_matrices.append(None)
            ptype_configs.append(cfg)
            continue

        pop_2d = np.loadtxt(file_2d, max_rows=n_depth)
        pdf_2d = pop_2d / models_count
        pdf_matrices.append(pdf_2d)
        nmi[ptype] = np.max(pdf_2d)
        ptype_configs.append(cfg)

    # Second run: Compute statistics and plot
    for ptype in range(4):
        pdf_2d = pdf_matrices[ptype]
        if pdf_2d is None:
            continue
        cfg = ptype_configs[ptype]

        # Prepare bins
        bin_width = all_bins[1, cfg['col_idx']] - all_bins[0, cfg['col_idx']]
        tmp_bins = all_bins[:, cfg['col_idx']] + (bin_width / 2.0)
        tmp_bins_p1 = all_bins_p1[:, cfg['col_idx']]
        if np.abs(tmp_bins[-1] - tmp_bins[0]) < 1e-12:
            continue

        # Velocity to slowness
        is_velocity = (ptype == 0 or ptype == 1)
        if is_velocity:
            tmp_bins_calc = 1.0 / tmp_bins
        else:
            tmp_bins_calc = tmp_bins.copy()

        # Prepare MAP/mean profiles
        max_mea_sig = np.zeros((2 * n_depth, 4))
        d_bins2 = np.zeros(2 * n_depth)
        count = 0
        for i in range(n_depth):
            row_pdf = pdf_2d[i, :]
            row_sum = np.sum(row_pdf)
            m_i = np.argmax(row_pdf)
            if row_sum > 0:
                # Arithmetic Mean
                am_val = np.sum(tmp_bins_calc * row_pdf) / row_sum
                # Weighted standard deviation
                var_val = np.sum(row_pdf * (tmp_bins_calc - am_val)**2)/row_sum
                sigma_tmp = np.sqrt(max(0.0, var_val))
            else:
                am_val = tmp_bins_calc[m_i]
                sigma_tmp = 0.0
            max_mea_sig[count, 0] = tmp_bins_calc[m_i]
            max_mea_sig[count, 1] = am_val
            max_mea_sig[count, 2] = am_val - sigma_tmp
            max_mea_sig[count, 3] = am_val + sigma_tmp

            d_bins2[count] = d_bins[i]
            count += 1
            max_mea_sig[count, :] = max_mea_sig[count-1, :]
            if i < n_depth - 1:
                d_bins2[count] = d_bins[i+1]
            else:
                d_bins2[count] = d_bins[-1] + (d_bins[1] - d_bins[0])
            count += 1

        # Slowness to velocity
        if is_velocity:
            max_mea_sig[:, 0:4] = 1.0 / max_mea_sig[:, 0:4]

        # -------------------------------------------
        # Plot the figure with pop statistics
        gs = gridspec.GridSpec(1, 3)
        px = 18 / 2.54
        py = 12 / 2.54
        plt.rcParams.update({
            'font.size': 8,
            'xtick.labelsize': 8,
            'ytick.labelsize': 8,
            'axes.labelsize': 9,
            'axes.titlesize': 9,
            'legend.fontsize': 7,
        })
        fig = plt.figure(figsize=(px, py), facecolor='w')

        x_lim_var = [tmp_bins_p1[0], tmp_bins_p1[-1]]
        y_lim_var = [0, min(depth_max, d_bins_p1[-1])]

        # Subplot 1: Posterior PDF
        ax1 = fig.add_subplot(gs[0, 0:2])
        mesh = ax1.pcolormesh(tmp_bins_p1, d_bins_p1, pdf_2d, cmap=cmap_dusk,
                              shading='flat', zorder=1)
        if plot_in_dpd:
            if syntmodel is not None:
                ax1.plot(getattr(syntmodel, cfg['syn_attr']), syntmodel.depth,
                         color='#4D4D4D', linewidth=0.8, linestyle='-',
                         label='Target', zorder=3)
            ax1.plot(getattr(maxmodel, cfg['max_attr']), maxmodel.depth,
                     color='b', linewidth=0.8, label='ML model', zorder=4)
            ax1.plot(getattr(mapmodel, cfg['map_attr']), mapmodel.depth,
                     color='m', linewidth=0.8, label='MAP model', zorder=5)
        ax1.set_xlim(x_lim_var)
        ax1.set_ylim(y_lim_var)
        ax1.invert_yaxis()

        # Set colormap
        if ptype in [0, 1]:
            mesh.set_clim([0, np.max(nmi[0:2])])
        else:
            delta = np.max(nmi[0:2]) - np.max(nmi[2:4])
            mesh.set_clim([0, np.max(nmi[0:2]) - (delta / 2.0)])

        # Colorbaru
        cbar = fig.colorbar(mesh, ax=ax1, location='left', pad=0.15)
        cbar.set_label('Probability', fontsize=8)
        cbar.ax.tick_params(labelsize=7)

        ax1.set_xlabel(f"{cfg['label']} {cfg['unit']}")
        ax1.set_ylabel('Depth (m)')
        ax1.set_title(f"Posterior marginal PDF of {cfg['label']}")
        ax1.grid(True, linestyle='--', linewidth=0.3, color='lightgray',
                 alpha=0.3)

        # Subplot 2: ML/MAP/AM Profiles
        ax2 = fig.add_subplot(gs[0, 2])
        pop_1d_scaled = ((pop_1d / np.max(pop_1d)) *
                         ((x_lim_var[1] - x_lim_var[0]) / 1.5))
        for i in range(1, len(pop_1d) - 1):
            height = log_d_bins[i+1] - log_d_bins[i]
            rect = plt.Rectangle((x_lim_var[0], log_d_bins[i]),
                                 pop_1d_scaled[i], height, facecolor='#999999',
                                 edgecolor='#999999', linewidth=0.1, zorder=2)
            ax2.add_patch(rect)

        # Plot 1D profiles
        if syntmodel is not None:
            ax2.plot(getattr(syntmodel, cfg['syn_attr']), syntmodel.depth,
                     color='#4D4D4D', linewidth=0.8, label='Target', zorder=3)
        ax2.plot(getattr(maxmodel, cfg['max_attr']), maxmodel.depth, color='b',
                 linewidth=0.8, label='ML model', zorder=4)
        ax2.plot(getattr(mapmodel, cfg['map_attr']), mapmodel.depth, color='m',
                 linewidth=0.8, label='MAP model', zorder=5)
        ax2.plot(max_mea_sig[:, 0], d_bins2, color='r', linewidth=0.8,
                 linestyle='-', label='MAX of PDF', zorder=6)
        ax2.plot(max_mea_sig[:, 1], d_bins2, color='g', linewidth=0.8,
                 linestyle='-', label='AM of PDF', zorder=7)
        ax2.plot(max_mea_sig[:, 2], d_bins2, color='g', linewidth=0.5,
                 linestyle=':', label=r'AM $\pm\sigma$', zorder=8)
        ax2.plot(max_mea_sig[:, 3], d_bins2, color='g', linewidth=0.5,
                 linestyle=':', zorder=8)

        ax2.set_xlim(x_lim_var)
        ax2.set_ylim(y_lim_var)
        ax2.invert_yaxis()

        # Legend
        if plot_leg:
            leg_patch = plt.matplotlib.patches.Patch(facecolor='#999999',
                                                     label='Interfaces')
            handles, labels = ax2.get_legend_handles_labels()
            handles.insert(0, leg_patch)
            labels.insert(0, 'Interfaces')
            ax2.legend(handles, labels, loc='best')

        ax2.set_xlabel(f"{cfg['label']} {cfg['unit']}")
        ax2.set_ylabel('Depth (m)')
        ax2.set_title('Profiles')
        ax2.grid(True, linestyle='--', linewidth=0.3, color='lightgray',
                 alpha=0.3)
        ax2.set_axisbelow(True)

        plt.tight_layout()

        # Save figure
        if fig_ren > 0:
            pdf_out = f"{out_path}_marginal_{cfg['out']}.pdf"
            png_out = f"{out_path}_marginal_{cfg['out']}.png"
            plt.savefig(pdf_out, dpi=600, bbox_inches='tight')
            plt.savefig(png_out, dpi=600, bbox_inches='tight')
            print(f"[+] SUCCESS: Saved marginal PDF to '{pdf_out}' and PNG.")

        if fig_ren > 1:
            plt.close(fig)
        else:
            plt.show()

        # Prepare for MAP uncertainty
        if ptype == 0:
            max_mea_sig_vs = max_mea_sig[::2].copy()
        elif ptype == 1:
            max_mea_sig_vp = max_mea_sig[::2].copy()

        # -------------------------------------------
        # Save ML/MAX models to .asii files
        if ptype == 0:
            # ML model export
            ml_out = f"{out_path}_ML_model.ascii"
            tmp_mod_raw = maxmodel.plot_mod[::2].copy()
            with open(ml_out, 'w', encoding='utf-8') as fid:
                fid.write("# MAXIMUM LIKELIHOOD MODEL. Number of layers:\n")
                fid.write(f"{maxmodel.n_layers}\n")
                fid.write("# Thickness[m]  Vp[m/s]  Vs[m/s] Density[kg/m3]\n")
                for i in range(maxmodel.n_layers - 1):
                    thick = tmp_mod_raw[i + 1, 0] - tmp_mod_raw[i, 0]
                    fid.write(f"{thick:7.1f} "
                              f"{tmp_mod_raw[i, 2]:7.1f} "
                              f"{tmp_mod_raw[i, 1]:7.1f} "
                              f"{tmp_mod_raw[i, 3]:7.1f}\n")
                fid.write("# Last line is the half-space\n")
                last_idx = maxmodel.n_layers - 1
                thick = 0.0
                fid.write(f"{thick:7.1f} "
                          f"{tmp_mod_raw[last_idx, 2]:7.1f} "
                          f"{tmp_mod_raw[last_idx, 1]:7.1f} "
                          f"{tmp_mod_raw[last_idx, 3]:7.1f}\n")
            print(f"[+] SUCCESS: ASCII saved to '{ml_out}'.")

            # MAP model export
            map_out = f"{out_path}_MAP_model.ascii"
            tmp_mod2_raw = mapmodel.plot_mod[::2].copy()
            with open(map_out, 'w', encoding='utf-8') as fid:
                fid.write("# MAXIMUM A POSTERIORI MODEL. Number of layers:\n")
                fid.write(f"{mapmodel.n_layers}\n")
                fid.write("# Thickness[m]  Vp[m/s]  Vs[m/s] Density[kg/m3]\n")
                for i in range(mapmodel.n_layers - 1):
                    thick = tmp_mod2_raw[i + 1, 0] - tmp_mod2_raw[i, 0]
                    fid.write(f"{thick:7.1f} "
                              f"{tmp_mod2_raw[i, 2]:7.1f} "
                              f"{tmp_mod2_raw[i, 1]:7.1f} "
                              f"{tmp_mod2_raw[i, 3]:7.1f}\n")
                fid.write("# Last line is the half-space\n")
                last_idx = mapmodel.n_layers - 1
                thick = 0.0
                fid.write(f"{thick:7.1f} "
                          f"{tmp_mod2_raw[last_idx, 2]:7.1f} "
                          f"{tmp_mod2_raw[last_idx, 1]:7.1f} "
                          f"{tmp_mod2_raw[last_idx, 3]:7.1f}\n")
            print(f"[+] SUCCESS: ASCII saved to '{map_out}'.")

        # -------------------------------------------
        # Save AM profiles to .asii files (Vs and Vp only)
        if ptype == 0 or ptype == 1:
            profile_out = f"{out_path}_AM_profile_{cfg['out']}.ascii"
            with open(profile_out, 'w', encoding='utf-8') as fid:
                fid.write(f"# {cfg['out']} profile from the posterior PDF\n")
                fid.write("# Number of pseudo-layers\n")
                fid.write(f"{n_depth}\n")
                fid.write(f"# Thickness[m],  MAXPDF({cfg['out']})[m/s],  "
                          f"AM({cfg['out']})[m/s],  AM_low({cfg['out']})[m/s] "
                          f" and AM_high({cfg['out']})[m/s]\n")
                for i in range(n_depth):
                    thick = d_bins[i+1] - d_bins[i] if i < n_depth - 1 else 0.0
                    v_tmp = max_mea_sig[2 * i, 0:4]
                    fid.write(f"{thick:7.1f} {v_tmp[0]:7.1f} {v_tmp[1]:7.1f} "
                              f"{v_tmp[3]:7.1f} {v_tmp[2]:7.1f}\n")
            print(f"[+] SUCCESS: ASCII saved to '{profile_out}'.")

    # -------------------------------------------
    # Save MAP model uncertainty to .asii file
    tmp_mod2_raw = mapmodel.plot_mod[::2].copy()
    n_layers_map = mapmodel.n_layers
    thick_map = np.zeros(n_layers_map)
    for i in range(n_layers_map - 1):
        thick_map[i] = tmp_mod2_raw[i + 1, 0] - tmp_mod2_raw[i, 0]
    thick_map[-1] = 0.0

    # Depth Center MAP
    cum_depth_map = np.cumsum(thick_map)
    depth_center_map = cum_depth_map - (thick_map / 2.0)
    depth_center_map[-1] = d_bins_p1[-1]
    d_bins_w = d_bins[1] - d_bins[0]
    depth_center_mms = d_bins + (d_bins_w / 2.0)

    uncertainty_out = f"{out_path}_MAP_model_uncertainty.ascii"
    with open(uncertainty_out, 'w', encoding='utf-8') as fid:
        fid.write("# MAP MODEL - UNCERTAINTY. Number of layers:\n")
        fid.write(f"{n_layers_map}\n")
        fid.write("# Thickness[m],  Vp_low[m/s],  Vp_high[m/s], Vs_low[m/s],  "
                  "Vs_high[m/s]\n")
        for i in range(n_layers_map):
            pos_idx = np.where(depth_center_mms >= depth_center_map[i])[0]
            pos = pos_idx[0] if len(pos_idx) > 0 else len(depth_center_mms) - 1
            correct_vs = max_mea_sig_vs[pos, 1] - tmp_mod2_raw[i, 1]
            thres_vs = max_mea_sig_vs[pos, 2:4] - correct_vs
            thres_vs[thres_vs < 0.0] = 0.0
            correct_vp = max_mea_sig_vp[pos, 1] - tmp_mod2_raw[i, 2]
            thres_vp = max_mea_sig_vp[pos, 2:4] - correct_vp
            thres_vp[thres_vp < 0.0] = 0.0
            if i == n_layers_map - 1:
                fid.write("# Last line is the half-space\n")
            fid.write(
                f"{thick_map[i]:7.1f} {thres_vp[1]:7.1f} {thres_vp[0]:7.1f} "
                f"{thres_vs[1]:7.1f} {thres_vs[0]:7.1f}\n"
            )
    print(f"[+] SUCCESS: ASCII saved to '{uncertainty_out}'.")


# -----------------------------------------------------------------------------
# Plot posterior marginal PDF for QWL
def plot_qwl(pop_prefix: Path, out_path: Path, qwl_head: SimpleNamespace,
             maxqwl: SimpleNamespace, mapqwl: SimpleNamespace,
             syntqwl: Optional[SimpleNamespace], plot_in_dpd: bool,
             plot_leg: bool, fig_ren: int, max_qwl_prop: float) -> None:
    """
    Plots posterior marginal PDFs for QWL parameters and Vs30 distribution.

    Args:
        pop_prefix (Path): Path to ensemble stats (prefix of input files).
        out_path (Path): Path to output directory for results.
        qwl_head (SimpleNamespace): QWL header metadata.
        maxqwl (SimpleNamespace): ML solution QWL data.
        mapqwl (SimpleNamespace): MAP solution  QWL data.
        syntqwl (SimpleNamespace): Target QWL data (Optional).
        plot_in_dpd: Plot ML and MAP data over PDF (True/False).
        plot_leg: Plot Legends (True/False).
        fig_ren (int): Save figures (0=NO-show, 1=YES-show, 2=YES-close).
        max_qwl_prop (float): Upper threshold of probability for the QWL plots.
    Returns:
        None.
    """
    n_freq = qwl_head.n_freq
    models_count = qwl_head.models_count
    f_bins_p1 = qwl_head.f_bins_p1
    qwl_bins_p1 = qwl_head.qwl_bins_p1
    vs30_bins = qwl_head.vs30_bins

    configs = {
        0: {'suffix': 'QWL_dep2D.txt', 'p1_idx': 1, 'qwl_idx': 0,
            'lbl': 'QWL depth (m)', 'attr': 'z_f'},
        1: {'suffix': 'QWL_vs2D.txt', 'p1_idx': 0, 'qwl_idx': 1,
            'lbl': 'QWL velocity (m/s)', 'attr': 'vs_f'},
        2: {'suffix': 'QWL_imp2D.txt', 'p1_idx': 2, 'qwl_idx': 3,
            'lbl': 'QWL impedance', 'attr': 'ic_f'}
    }

    # Create DUSK colormap
    color_matrix = dusk()
    cmap_dusk = ListedColormap(color_matrix, name='dusk')

    # Plot the figure with QWL statistics
    px = 21 / 2.54
    py = 14 / 2.54
    plt.rcParams.update({
        'font.size': 8,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'axes.labelsize': 9,
        'axes.titlesize': 9,
        'legend.fontsize': 7,
    })
    fig = plt.figure(figsize=(px, py), facecolor='w', layout=None)
    gs = gridspec.GridSpec(2, 3, wspace=0.35, hspace=0.35)
    fig.subplots_adjust(left=0.08, right=0.98, top=0.95, bottom=0.08)

    # -------------------------------------------
    # Subplots 1-3: QWL
    for ptype in range(3):
        cfg = configs[ptype]
        file_2d = Path(f"{pop_prefix}_{cfg['suffix']}")
        if not file_2d.exists():
            print(f"[!] Warning: QWL file missing: {file_2d.name}.")
            continue
        pop_qwl2d = np.loadtxt(file_2d, max_rows=n_freq)
        pdf_qwl2d = pop_qwl2d / models_count
        pdf_qwl2d = pdf_qwl2d.T

        ax = fig.add_subplot(gs[0, ptype])
        tmp_bins_p1 = qwl_bins_p1[:, cfg['p1_idx']]
        mesh = ax.pcolormesh(f_bins_p1, tmp_bins_p1, pdf_qwl2d, cmap=cmap_dusk,
                             shading='flat', zorder=1)
        mesh.set_clim([0.0, max_qwl_prop])
        if plot_in_dpd:
            if syntqwl is not None:
                if ptype == 2:
                    syn_y = 1.0 / getattr(syntqwl, cfg['attr'])
                else:
                    syn_y = getattr(syntqwl, cfg['attr'])
                ax.plot(qwl_head.freq, syn_y, color='#4D4D4D', linewidth=0.8,
                        linestyle='-', zorder=3)
            max_y = maxqwl.qwl_data[:, cfg['qwl_idx']]
            map_y = mapqwl.qwl_data[:, cfg['qwl_idx']]
            ax.plot(maxqwl.freq, max_y, color='b', linewidth=0.8, zorder=4)
            ax.plot(mapqwl.freq, map_y, color='m', linewidth=0.8, zorder=5)

        ax.set_xscale('log')
        if ptype == 0:
            ax.set_yscale('log')
            ax.set_yticks(qwl_head.d_tick)
            ax.set_yticklabels(qwl_head.d_tick_label)

        ax.set_xticks(qwl_head.f_tick)
        ax.set_xticklabels(qwl_head.f_tick_label)
        ax.set_xlim([f_bins_p1[0], f_bins_p1[-1]])
        ax.set_ylim([tmp_bins_p1[0], tmp_bins_p1[-1]])
        ax.set_xlabel('Frequency (Hz)')
        ax.set_ylabel(cfg['lbl'])
        ax.set_title('Posterior marginal PDF')
        ax.grid(True, linestyle='--', linewidth=0.3, color='lightgray',
                alpha=0.3)
        ax.set_aspect('auto')

    # -------------------------------------------
    # Read Vs30 statistics
    vs30_file = Path(f"{pop_prefix}_vs30.txt")
    if not vs30_file.exists():
        raise FileNotFoundError(f"Vs30 stats file missing: {vs30_file.name}")
    pop_vs30 = np.loadtxt(vs30_file)

    vs_bins_w_s30 = vs30_bins[1] - vs30_bins[0]
    vs_bins_c_s30 = vs30_bins + (vs_bins_w_s30 / 2.0)

    # Statistical AM / MAP calculation using slowness conversion
    vs_bins_tmp = 1.0 / vs_bins_c_s30
    m_i = np.argmax(pop_vs30)
    row_sum = np.sum(pop_vs30)
    am_val = np.sum(vs_bins_tmp * pop_vs30) / row_sum
    var_val = np.sum(pop_vs30 * (vs_bins_tmp - am_val)**2) / row_sum
    sigma_tmp = np.sqrt(max(0.0, var_val))
    max_mea_sig = np.zeros(6)
    max_mea_sig[0] = 1.0 / vs_bins_tmp[m_i]
    max_mea_sig[1] = 1.0 / am_val
    max_mea_sig[2] = 1.0 / (am_val - sigma_tmp)
    max_mea_sig[3] = 1.0 / (am_val + sigma_tmp)
    max_mea_sig[4] = 1.0 / (am_val - 2.0 * sigma_tmp)
    max_mea_sig[5] = 1.0 / (am_val + 2.0 * sigma_tmp)

    pop_vs30_max = np.max(pop_vs30 / models_count)
    nx_lim = [0.0, 1.2 * pop_vs30_max]

    # -------------------------------------------
    # Subplot 4: Vs30
    ax_hist = fig.add_subplot(gs[1, 0:2])
    hch_lines = []
    hch_lbls = []
    h_bar = ax_hist.bar(vs_bins_c_s30, pop_vs30 / models_count,
                        width=vs_bins_w_s30, facecolor='#666666',
                        edgecolor='none', label='Distribution', zorder=2)
    hch_lines.append(h_bar)
    hch_lbls.append('Distribution')

    if syntqwl is not None:
        l_syn = ax_hist.plot([], [], color='#4D4D4D', linewidth=0.8,
                             label='Target model', zorder=3)[0]
        hch_lines.append(l_syn)
        hch_lbls.append('Target model')
    l_ml = ax_hist.plot([], [], color='b', linewidth=0.8, label='ML model',
                        zorder=4)[0]
    l_map = ax_hist.plot([], [], color='m', linewidth=0.8, label='MAP model',
                         zorder=5)[0]
    hch_lines.extend([l_ml, l_map])
    hch_lbls.extend(['ML model', 'MAP model'])

    if syntqwl is not None:
        p_syn, = ax_hist.plot([syntqwl.vs30], [nx_lim[1]/20.0], 'x',
                              color='#4D4D4D', markersize=5, zorder=6)
        hch_lines.append(p_syn)
        hch_lbls.append(r'$V_{S30}^{\mathrm{Target}}$')

    p_ml, = ax_hist.plot([maxqwl.vs30[1]], [nx_lim[1]/20.0], 'xb',
                         markersize=5, zorder=6)
    p_map, = ax_hist.plot([mapqwl.vs30[1]], [nx_lim[1]/20.0], 'xm',
                          markersize=5, zorder=6)
    p_max, = ax_hist.plot([max_mea_sig[0]], [pop_vs30_max], 'xr', markersize=5,
                          zorder=6)

    # AM mean limits propagation (1 sigma and 2 sigma boundaries)
    l_2sig, = ax_hist.plot([max_mea_sig[4], max_mea_sig[5]],
                           [nx_lim[1]/10.0, nx_lim[1]/10.0], color='lime',
                           linestyle=':', linewidth=1, label='AM (2 sigma)',
                           zorder=7)
    l_1sig = ax_hist.errorbar(max_mea_sig[1], nx_lim[1]/10.0,
                              xerr=[[max_mea_sig[1] - max_mea_sig[3]],
                                    [max_mea_sig[2] - max_mea_sig[1]]],
                              fmt='o', color='lime', ecolor='lime',
                              markersize=3, capsize=3, elinewidth=0.8,
                              label='AM (1 sigma)', zorder=8)

    hch_lines.extend([p_ml, p_map, p_max, l_1sig, l_2sig])
    hch_lbls.extend([r'$V_{S30}^{\mathrm{ML}}$',
                     r'$V_{S30}^{\mathrm{MAP}}$',
                     r'$V_{S30}^{\mathrm{MAX}}$',
                     r'$V_{S30}^{\mathrm{AM}}\ (1\sigma)$',
                     r'$V_{S30}^{\mathrm{AM}}\ (2\sigma)$'])

    ax_hist.set_ylim(nx_lim)
    ax_hist.set_xlim([
        int(max_mea_sig[1] - 5.0 * (max_mea_sig[1] - max_mea_sig[3])),
        int(max_mea_sig[1] + 5.0 * (max_mea_sig[2] - max_mea_sig[1]))
    ])
    ax_hist.set_ylabel('Probability')
    ax_hist.set_xlabel(r'$V_{S30}$ (m/s)')
    ax_hist.set_title(r'Average $V_S$ down to 30 m')
    ax_hist.grid(True, linestyle='--', linewidth=0.3, color='lightgray',
                 alpha=0.3)
    ax_hist.set_axisbelow(True)

    if plot_leg:
        ax_hist.legend(hch_lines, hch_lbls, loc='best', fontsize=7)

    # -------------------------------------------
    # Subplot 6: Text and colorbar
    ax_txt = fig.add_subplot(gs[1, 2])
    ax_txt.axis('off')

    # Construct engineering report matrices
    if syntqwl is not None:
        v_tgt_str = (f"$V_{{S30}}^{{\\mathrm{{Target}}}} = "
                     f"{round(syntqwl.vs30)}$ (m/s)")
        f_tgt_str = (f"$f_{{30}}^{{\\mathrm{{Target}}}} = "
                     f"{syntqwl.vs30 / 120.0:.2f}$ (Hz)")
    else:
        v_tgt_str = ""
        f_tgt_str = ""

    text_block_left = (
        f"{v_tgt_str}\n"
        f"$V_{{S30}}^{{\\mathrm{{ML}}}} = {round(maxqwl.vs30[1])}$ m/s\n"
        f"$V_{{S30}}^{{\\mathrm{{MAP}}}} = {round(mapqwl.vs30[1])}$ m/s\n"
        f"$V_{{S30}}^{{\\mathrm{{MAX}}}} = {round(max_mea_sig[0])}$ m/s\n"
        f"$V_{{S30}}^{{\\mathrm{{AM}}}} = {round(max_mea_sig[1])}$ m/s\n"
        f"({round(max_mea_sig[3])} - {round(max_mea_sig[2])} m/s)"
    )

    text_block_right = (
        f"{f_tgt_str}\n"
        f"$f_{{30}}^{{\\mathrm{{ML}}}} = {maxqwl.vs30[1] / 120.0:.2f}$ Hz\n"
        f"$f_{{30}}^{{\\mathrm{{MAP}}}} = {mapqwl.vs30[1] / 120.0:.2f}$ Hz\n"
        f"$f_{{30}}^{{\\mathrm{{MAX}}}} = {max_mea_sig[0] / 120.0:.2f}$ Hz\n"
        f"$f_{{30}}^{{\\mathrm{{AM}}}} = {max_mea_sig[1] / 120.0:.2f}$ Hz\n"
        f"({max_mea_sig[3] / 120.0:.1f} - {max_mea_sig[2] / 120.0:.1f} Hz)"
    )

    ax_txt.text(0.0, 0.7, text_block_left, va='top', ha='left', fontsize=8)
    ax_txt.text(0.5, 0.7, text_block_right, va='top', ha='left', fontsize=8)

    # Horizontal global colorbar placement matching northoutside constraints
    cbar_ax = fig.add_axes([0.68, 0.42, 0.25, 0.03])
    sm = plt.cm.ScalarMappable(cmap=cmap_dusk,
                               norm=plt.Normalize(vmin=0, vmax=max_qwl_prop))
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation='horizontal')
    cbar.set_label('Probability', fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    # Save figure
    if fig_ren > 0:
        pdf_out = f"{out_path}_QWL.pdf"
        png_out = f"{out_path}_QWL.png"
        plt.savefig(pdf_out, dpi=600, bbox_inches='tight')
        plt.savefig(png_out, dpi=600, bbox_inches='tight')
        print(f"[+] SUCCESS: Saved QWL panel to '{pdf_out}' and PNG.")

    if fig_ren > 1:
        plt.close(fig)
    else:
        plt.show()

    # -------------------------------------------
    # Save Vs30 into text file
    vs30_out = f"{out_path}_Vs30.ascii"
    with open(vs30_out, 'w', encoding='utf-8') as f:
        f.write("# Vs30[m/s],  f30[Hz],  Model Type\n")
        f.write(f"{maxqwl.vs30[1]:8.2f} {maxqwl.vs30[1] / 120.0:8.4f}  "
                "ML model\n")
        f.write(f"{mapqwl.vs30[1]:8.2f} {mapqwl.vs30[1] / 120.0:8.4f}  "
                "MAP model\n")
        f.write(f"{max_mea_sig[0]:8.2f} {max_mea_sig[0] / 120.0:8.4f}  "
                "MAX from PDF\n")
        f.write(f"{max_mea_sig[1]:8.2f} {max_mea_sig[1] / 120.0:8.4f}  "
                "AM from PDF\n")
        f.write(f"{max_mea_sig[3]:8.2f} {max_mea_sig[3] / 120.0:8.4f}  "
                "Low threshold from PDF (-1 sigma)\n"
                )
        f.write(f"{max_mea_sig[2]:8.2f} {max_mea_sig[2] / 120.0:8.4f}  "
                "High threshold from PDF (+1 sigma)\n"
                )
        f.write(f"{max_mea_sig[5]:8.2f} {max_mea_sig[5] / 120.0:8.4f}  "
                "Low threshold from PDF (-2 sigma)\n"
                )
        f.write(f"{max_mea_sig[4]:8.2f} {max_mea_sig[4] / 120.0:8.4f}  "
                "High threshold from PDF (+2 sigma)\n"
                )
    print(f"[+] SUCCESS: ASCII saved to '{vs30_out}'.")


# -----------------------------------------------------------------------------
# Plot posterior marginal PDF for SH amplification
def plot_sh(pop_prefix: Path, out_path: Path, qwl_head: SimpleNamespace,
            maxqwl: SimpleNamespace, mapqwl: SimpleNamespace,
            syntqwl: Optional[SimpleNamespace], plot_leg: bool,
            fig_ren: int) -> None:
    """
    Plots posterior marginal PDFs of SH-wave transfer func. and amplification.

    Args:
        pop_prefix (Path): Path to ensemble stats (prefix of input files).
        out_path (Path): Path to output directory for results.
        qwl_head (SimpleNamespace): QWL header metadata.
        maxqwl (SimpleNamespace): ML solution QWL/SH data.
        mapqwl (SimpleNamespace): MAP solution QWL/SH data.
        syntqwl (SimpleNamespace): Target QWL/SH data (Optional).
        plot_leg (bool): Plot Legends (True/False).
        fig_ren (int): Save figures (0=NO-show, 1=YES-show, 2=YES-close).
    Returns:
        None.
    """
    n_freq = qwl_head.n_freq
    models_count = qwl_head.models_count
    f_bins_p1 = qwl_head.f_bins_p1
    qwl_bins_p1 = qwl_head.qwl_bins_p1
    freq = qwl_head.freq
    amp_vec = qwl_head.amp_vec

    configs = {
        0: {'suffix': 'SH_unr2D.txt', 'idx': 0,
            'title': 'Posterior marginal PDF of\nSH-wave transfer function'},
        1: {'suffix': 'SH_ref2D.txt', 'idx': 1,
            'title': 'Posterior marginal PDF of\namplif. to reference profile'}
    }

    # Create DUSK colormap
    color_matrix = dusk()
    cmap_dusk = ListedColormap(color_matrix, name='dusk')

    # Plot the figure
    px = 20 / 2.54
    py = 14 / 2.54
    plt.rcParams.update({
        'font.size': 8,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'axes.labelsize': 9,
        'axes.titlesize': 9,
        'legend.fontsize': 7,
    })
    fig = plt.figure(figsize=(px, py), facecolor='w', layout=None)
    gs = gridspec.GridSpec(2, 2, wspace=0.20, hspace=0.6)
    fig.subplots_adjust(left=0.08, right=0.96, top=0.92, bottom=0.08)

    log_amp_vec = np.log10(amp_vec)
    d_log = log_amp_vec[1] - log_amp_vec[0]
    tmp_bins_log = log_amp_vec + (d_log / 2.0)
    tmp_bins_p1 = qwl_bins_p1[:, 3]

    for ptype in range(2):
        cfg = configs[ptype]
        file_2d = Path(f"{pop_prefix}_{cfg['suffix']}")
        if not file_2d.exists():
            print(f"[!] Warning: SH file missing: {file_2d}.")
            continue

        pop_qwl2d = np.loadtxt(file_2d, max_rows=n_freq)
        pdf_qwl2d = pop_qwl2d / models_count
        pdf_qwl2d = pdf_qwl2d.T

        # Target model
        if syntqwl is not None:
            if ptype == 0:
                syn_mod_y = np.abs(syntqwl.h3)
            else:
                syn_mod_y = np.abs(syntqwl.h3) / syntqwl.ampcor
        else:
            syn_mod_y = None

        # Extract ML and MAP data
        max_mod_y = maxqwl.sh_data[:, cfg['idx']]
        map_mod_y = mapqwl.sh_data[:, cfg['idx']]

        # Arithmetic Mean & Sigma Profiling (Logarithmic Space)
        max_mea_sig = np.zeros((n_freq, 4))
        for i in range(n_freq):
            r_data = pop_qwl2d[i, :]
            r_sum = np.sum(r_data)
            m_i = np.argmax(r_data)
            if r_sum > 0:
                mean_val = np.sum(tmp_bins_log * r_data) / r_sum
                var_val = np.sum(r_data * (tmp_bins_log - mean_val)**2) / r_sum
                std_val = np.sqrt(max(0.0, var_val))
            else:
                mean_val = tmp_bins_log[m_i]
                std_val = 0.0
            max_mea_sig[i, 0] = tmp_bins_log[m_i]
            max_mea_sig[i, 1] = mean_val
            max_mea_sig[i, 2] = mean_val - std_val
            max_mea_sig[i, 3] = mean_val + std_val

        # Convert back from log space to linear amplification values
        max_mea_sig = 10.0**max_mea_sig

        # -------------------------------------------
        # Subplot 1: PDF plot
        ax_pdf = fig.add_subplot(gs[0, ptype])
        mesh = ax_pdf.pcolormesh(f_bins_p1, tmp_bins_p1, pdf_qwl2d,
                                 cmap=cmap_dusk, shading='flat', zorder=1)
        ax_pdf.set_xscale('log')
        ax_pdf.set_yscale('log')
        ax_pdf.set_xticks(qwl_head.f_tick)
        ax_pdf.set_xticklabels(qwl_head.f_tick_label)
        ax_pdf.set_yticks(qwl_head.a_tick)
        ax_pdf.set_yticklabels(qwl_head.a_tick_label)

        ax_pdf.set_xlim([f_bins_p1[0], f_bins_p1[-1]])
        ax_pdf.set_ylim([tmp_bins_p1[0], tmp_bins_p1[-1]])
        ax_pdf.set_xlabel('Frequency (Hz)')
        ax_pdf.set_ylabel('Amplification')
        ax_pdf.set_title(cfg['title'])
        ax_pdf.grid(True, linestyle='--', linewidth=0.3, color='lightgray',
                    alpha=0.3)
        ax_pdf.set_aspect('auto')

        # Colorbar placement under the active PDF mesh plot
        cbar_ax = fig.add_axes([0.11 + (ptype) * 0.49, 0.49, 0.32, 0.025])
        cbar = fig.colorbar(mesh, cax=cbar_ax, orientation='horizontal')
        cbar.set_label('Probability', fontsize=7)
        cbar.ax.tick_params(labelsize=6)

        # -------------------------------------------
        # Subplot 2: Line profiles
        ax_line = fig.add_subplot(gs[1, ptype])
        hch_lines = []
        hch_lbls = []

        # Target model
        if syntqwl is not None:
            l_syn = ax_line.plot(maxqwl.freq, syn_mod_y, color='#4D4D4D',
                                 linewidth=0.8, linestyle='-',
                                 label='Target model', zorder=3)[0]
            hch_lines.append(l_syn)
            hch_lbls.append('Target model')

        l_ml = ax_line.plot(maxqwl.freq, max_mod_y, 'b-', linewidth=0.8,
                            label='ML model', zorder=4)[0]
        l_map = ax_line.plot(mapqwl.freq, map_mod_y, 'm-', linewidth=0.8,
                             label='MAP model', zorder=5)[0]
        l_am = ax_line.plot(freq[:n_freq], max_mea_sig[:, 1], 'g-',
                            linewidth=0.8, label='AM of PDF', zorder=6)[0]
        l_sig1 = ax_line.plot(freq[:n_freq], max_mea_sig[:, 2], 'g:',
                              linewidth=0.8, label=r'AM $\pm\sigma$',
                              zorder=7)[0]
        ax_line.plot(freq[:n_freq], max_mea_sig[:, 3], 'g:',
                     linewidth=0.8, zorder=7)

        if ptype == 0 or not syntqwl:
            hch_lines.extend([l_ml, l_map, l_am, l_sig1])
            hch_lbls.extend([
                'ML model', 'MAP model', 'AM of PDF', r'AM $\pm\sigma$'
            ])

        ax_line.set_xscale('log')
        ax_line.set_yscale('log')
        ax_line.set_xticks(qwl_head.f_tick)
        ax_line.set_xticklabels(qwl_head.f_tick_label)
        ax_line.set_yticks(qwl_head.a_tick)
        ax_line.set_yticklabels(qwl_head.a_tick_label)

        ax_line.set_xlim([f_bins_p1[0], f_bins_p1[-1]])
        ax_line.set_ylim([tmp_bins_p1[0], tmp_bins_p1[-1]])
        ax_line.set_xlabel('Frequency (Hz)')
        ax_line.set_ylabel('Amplification')
        ax_line.grid(True, linestyle='--', linewidth=0.3, color='lightgray',
                     alpha=0.3)
        ax_line.set_axisbelow(True)
        ax_line.set_aspect('auto')

        if ptype == 0 and plot_leg:
            ax_line.legend(hch_lines, hch_lbls, loc='upper left', fontsize=7)

        # ---------------------------------------------------------------------
        # Save amplification to .asii files
        if ptype == 1:
            # ML model amplification
            ml_out = f"{out_path}_ML_amplification.ascii"
            with open(ml_out, 'w', encoding='utf-8') as f:
                f.write("# MAXIMUM LIKELIHOOD MODEL: Amplification to the "
                        "reference velocity profile\n")
                f.write("# Frequency[Hz],  Amplification[-]\n")
                for i in range(len(maxqwl.freq)):
                    f.write(f"{maxqwl.freq[i]:9.5f} {max_mod_y[i]:9.5f}\n")

            # MAP model amplification
            map_out = f"{out_path}_MAP_amplification.ascii"
            with open(map_out, 'w', encoding='utf-8') as f:
                f.write("# MAXIMUM A POSTERIORI MODEL: Amplification to the "
                        "reference velocity profile\n")
                f.write("# Frequency[Hz],  Amplification[-]\n")
                for i in range(len(mapqwl.freq)):
                    f.write(f"{mapqwl.freq[i]:9.5f} {map_mod_y[i]:9.5f}\n")

            # AM profiles
            am_out = f"{out_path}_AM_profile_amplification.ascii"
            with open(am_out, 'w', encoding='utf-8') as f:
                f.write("# Arithmetic Mean profile: Amplification to the "
                        "reference velocity profile (elastic)\n")
                f.write("# Frequency[Hz],  Amplification[-],  Ampl_low[-],  "
                        "Ampl_high[-]\n")
                for i in range(n_freq):
                    f.write(f"{freq[i]:9.5f} {max_mea_sig[i, 1]:9.5f} "
                            f"{max_mea_sig[i, 2]:9.5f} "
                            f"{max_mea_sig[i, 3]:9.5f}\n")
            print("[+] SUCCESS: 3X ASCII saved to '*amplification.ascii'.")

    # Save figure
    if fig_ren > 0:
        pdf_out = f"{out_path}_SH.pdf"
        png_out = f"{out_path}_SH.png"
        plt.savefig(pdf_out, dpi=600, bbox_inches='tight')
        plt.savefig(png_out, dpi=600, bbox_inches='tight')
        print(f"[+] SUCCESS: Saved marginal PDF for SH '{pdf_out}' and PNG.")

    if fig_ren > 1:
        plt.close(fig)
    else:
        plt.show()


# -----------------------------------------------------------------------------
# Main function
def main():
    """
    Main function - Read and plot results from NEOPSY.

    - Read resulting data (in text files).
    - Plot and save resulting figures.
    """
    print("-" * 50)

    # Create directory for results
    OUT_FOLDER.mkdir(parents=True, exist_ok=True)
    out_path = OUT_FOLDER / OUT_PREFIX

    # -------------------------------------------
    # Load POP headers
    try:
        header = read_pop_headers(POP_PREFIX)
    except Exception as e:
        print(f"[!] Failed to read headers: {e}.")
        return
    print(f"[*] Number of sampling models: {header.models_count}.")

    # -------------------------------------------
    # Read input data (DC and ELL curves)
    try:
        indata = read_input_data(DATAFILE)
    except Exception as e:
        print(f"[!] Failed to read input data: {e}.")
        return

    # -------------------------------------------
    # Read MAX and MAP model data (DC and ELL curves)
    try:
        datafile = Path(f"{MAX_PREFIX}_data.txt")
        maxdata = read_model_data(datafile)
    except Exception as e:
        print(f"[*] MAX model data missing: {e}. Using fallback input data.")
        maxdata = indata

    try:
        datafile = Path(f"{MAP_PREFIX}_data.txt")
        mapdata = read_model_data(datafile)
    except Exception as e:
        print(f"[*] MAP model data missing: {e}. Using fallback input data.")
        mapdata = indata

    # -------------------------------------------
    # Read ML and MAP velocity model
    try:
        model_file = Path(f"{MAX_PREFIX}.txt")
        maxmodel = read_velocity_model(model_file, PLOT_DEPTH_MAX)
    except Exception as e:
        print(f"[!] ERROR: Failed to read ML velocity model: {e}.")
        return

    try:
        model_file = Path(f"{MAP_PREFIX}.txt")
        mapmodel = read_velocity_model(model_file, PLOT_DEPTH_MAX)
    except Exception as e:
        print(f"[!] ERROR: Failed to read MAP velocity model: {e}.")
        return

    # -------------------------------------------
    # Read synthetic target model (if enabled)
    syntmodel = None
    if SYNT_PLOT:
        try:
            syntmodel = read_synthetic_model(SYNT_FILE, PLOT_DEPTH_MAX)
            print("[+] SUCCESS: Loaded target model.")
        except Exception as e:
            print(f"[!] WARNING: Failed to read synthetic model: {e}.")
            syntmodel = None

    # -------------------------------------------
    # Read POP data binary ensemble
    try:
        bin_file = Path(f"{POP_PREFIX}_data.bin")
        syndata = read_pop_binary_data(bin_file, indata.f_n, N_FITS)
    except Exception as e:
        print(f"[!] ERROR: Failed to read ensemble: {e}.")
        return

    # -------------------------------------------
    # Plot data fits
    print("[*] Plot data fits.")
    try:
        plot_data_fit(indata, maxdata, mapdata, syndata, out_path, 'slowness',
                      USE_DATA_THR, SLOW_THR, FREQ_THR, PLOT_LEG, FIGURE_REN)
        plot_data_fit(indata, maxdata, mapdata, syndata, out_path, 'velocity',
                      USE_DATA_THR, SLOW_THR, FREQ_THR, PLOT_LEG, FIGURE_REN)
    except Exception as e:
        print(f"[!] ERROR: Failed to plot results: {e}.")
        return

    # -------------------------------------------
    # Plot basic histograms
    print("[*] Plot basic histograms.")
    try:
        plot_basic_histograms(POP_PREFIX, out_path, header, maxmodel, mapmodel,
                              FIGURE_REN)
    except Exception as e:
        print(f"[!] ERROR: Failed to plot basic histograms: {e}.")
        return

    # -------------------------------------------
    # Plot posterior marginal PDFs, profiles and save ASCII results
    print("[*] Plot posterior marginal PDFs.")
    try:
        plot_pdf(POP_PREFIX, out_path, header, maxmodel, mapmodel, syntmodel,
                 PLOT_IN_DPD, PLOT_LEG, FIGURE_REN, PLOT_DEPTH_MAX)
    except Exception as e:
        print(f"[!] ERROR: Failed to plot posterior PDFs: {e}.")
        return

    # -------------------------------------------
    # Load QWL headers
    try:
        qwl = read_qwl_headers(POP_PREFIX)
    except Exception as e:
        print(f"[!] Failed to read QWL headers: {e}.")
        return

    # -------------------------------------------
    # Load ML and MAP model QWL data
    try:
        qwl_file = Path(f"{MAX_PREFIX}_QWL.txt")
        maxqwl = read_model_qwl(qwl_file)
    except Exception as e:
        print(f"[!] ERROR: Failed to read ML QWL data: {e}.")
        return

    try:
        qwl_file = Path(f"{MAP_PREFIX}_QWL.txt")
        mapqwl = read_model_qwl(qwl_file)
    except Exception as e:
        print(f"[!] ERROR: Failed to read MAP QWL data: {e}.")
        return

    # -------------------------------------------
    # Compute QWL for synthetic target model (if enabled)
    syntqwl = None
    if SYNT_PLOT:
        try:
            syntqwl = compute_qwl(SYNT_FILE, REF_FILE, qwl.freq)
            print("[+] SUCCESS: Computed QWL for target model.")
        except Exception as e:
            print(f"[!] WARNING: Failed to compute QWL for target model: {e}.")
            syntqwl = None

    # -------------------------------------------
    # Plot posterior marginal PDF for QWL representation
    print("[*] Plot posterior marginal PDF for QWL.")
    try:
        plot_qwl(POP_PREFIX, out_path, qwl, maxqwl, mapqwl, syntqwl,
                 PLOT_IN_DPD, PLOT_LEG, FIGURE_REN, MAX_QWL_PROP)
    except Exception as e:
        print(f"[!] ERROR: Failed to plot posterior PDFs: {e}.")
        return

    # -------------------------------------------
    # Plot posterior marginal PDF for SH amplification
    print("[*] Plot posterior marginal PDF for SH amplification.")
    try:
        plot_sh(POP_PREFIX, out_path, qwl, maxqwl, mapqwl, syntqwl,
                PLOT_LEG, FIGURE_REN)
    except Exception as e:
        print(f"[!] ERROR: Failed to plot posterior PDFs: {e}.")
        return

    print("[*] SUCCESS: All done.")


# -----------------------------------------------------------------------------
# Entry point
if __name__ == "__main__":
    main()
