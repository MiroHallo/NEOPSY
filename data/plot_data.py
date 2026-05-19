#!/usr/bin/env python3
# =============================================================================
# PLOT INPUT OBSERVED DATA (dispersion and ellipticity curves)
#
# Author: Miroslav HALLO, Kyoto University
# E-mail: hallo.miroslav.2a@kyoto-u.ac.jp
# Revision 2026/05: First version
# Tested with: Python 3.12.3, Matplotlib 3.10.9, NumPy 2.4.5
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

import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter


# =============================================================================
# INPUT PARAMETERS
# =============================================================================

# Prefix for data files
PREF = ''

# Use modes? (1=YES, 0=NO)
#    R0, R1, R2, R3
#    L0, L1, L2, L3
#    ELL, ELA, None, None
USE = np.array([
    [1, 0, 0, 0],
    [0, 0, 0, 0],
    [1, 0, None, None]
], dtype=float)

# Data filenames
FILENAMES = {
    'R0': f"{PREF}Data_R0.dat",
    'R1': f"{PREF}Data_R1.dat",
    'R2': f"{PREF}Data_R2.dat",
    'R3': f"{PREF}Data_R3.dat",
    'L0': f"{PREF}Data_L0.dat",
    'L1': f"{PREF}Data_L1.dat",
    'L2': f"{PREF}Data_L2.dat",
    'L3': f"{PREF}Data_L3.dat",
    'ELL': f"{PREF}Data_ELL.dat",
    'ELA': f"{PREF}Data_ELA.dat"
}


# =============================================================================
# FUNCTIONS
# =============================================================================

# Automatic frequency axis limits
def fscale(fmin: float, fmax: float) -> tuple[float, float, np.array]:
    """
    Return frequency axis limits and intelligent ticks.

    Args:
        fmin (float): Min frequency [Hz]
        fmax (float): Max frequency [Hz]
    Returns:
        tuple: amin (float): Min anchored frequency [Hz]
               amax (float): Max anchored frequency [Hz]
               fticks: Intelligent ticks [Hz]
    """
    # Define stable anchor points for axis limits
    stable_limits = np.array([0.05, 0.1, 0.2, 0.5, 1, 5, 10, 20, 30, 50, 100])

    # Find closest lower and upper anchor points
    lower_idx = np.where(stable_limits <= fmin)[0]
    upper_idx = np.where(stable_limits >= fmax)[0]

    # If outside the range
    if lower_idx.size > 0:
        amin = stable_limits[lower_idx[-1]]
    else:
        amin = np.min(stable_limits)
    if upper_idx.size > 0:
        amax = stable_limits[upper_idx[0]]
    else:
        amax = np.max(stable_limits)

    # Generate all ticks
    decades = 10.0 ** np.arange(-2, 2)
    major_prio = np.array([1, 2, 5])
    minor_prio = np.array([3, 4, 6, 8])
    detail_prio = np.array([1.5, 2.5, 7, 9])

    ticks_major = np.sort(np.outer(major_prio, decades).flatten())
    ticks_minor = np.sort(np.outer(minor_prio, decades).flatten())
    ticks_detail = np.sort(np.outer(detail_prio, decades).flatten())

    # Select tick set
    ratio = amax / amin
    if ratio > 10:
        all_ticks = ticks_major
    elif ratio > 3:
        all_ticks = np.concatenate([ticks_major, ticks_minor])
    else:
        all_ticks = np.concatenate([ticks_major, ticks_minor, ticks_detail])
    all_ticks = np.sort(all_ticks)

    # Select only ticks that fall within limits
    fticks = all_ticks[(all_ticks >= amin) & (all_ticks <= amax)]

    return float(amin), float(amax), fticks


# -----------------------------------------------------------------------------
# Plot dispersion curves
def plot_dispersion(ftxt: str, filepath: str) -> None:
    """
    Load and plot surface wave dispersion curves

    Args:
        ftxt (str): Code of the data type (R0, R1, R2, R3, L0, L1, L2, L3).
        filepath (str): Path to the file with data.
    Returns:
        None.
    """
    # Load data
    if not os.path.exists(filepath):
        print(f"[!] ERROR: File {filepath} does not exist.")
        return None
    inData = np.loadtxt(filepath, skiprows=1)

    # Convert to velocity
    inData_v = np.zeros((len(inData), 4))
    inData_v[:, 0] = inData[:, 0]
    inData_v[:, 1] = 1.0/inData[:, 1]
    inData_v[:, 2] = (1.0/inData[:, 1]) - (1.0/(inData[:, 1] + inData[:, 2]))
    inData_v[:, 3] = -(1.0/inData[:, 1]) + (1.0/(inData[:, 1] - inData[:, 2]))

    # Get used frequencies
    okf = inData[:, 2] > 0
    if not np.any(okf):
        print(f"[!] WARNING: No valid data to plot in file {filepath}.")
        return None

    # Get frequency limiths
    min_f = np.min(inData[okf, 0])
    max_f = np.max(inData[okf, 0])
    amin, amax, fticks = fscale(min_f, max_f)

    # Start figure
    plt.rcParams.update({
        'font.size': 8,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'axes.labelsize': 9,
        'axes.titlesize': 9,
        'legend.fontsize': 7,
    })
    fig, axs = plt.subplots(2, 1, figsize=(16/2.54, 12/2.54), facecolor='w')

    # Subplot 1: Slowness
    axs[0].semilogx(inData[okf, 0], inData[okf, 1] * 1000,
                    color='k', linestyle=':', label='Data line')
    axs[0].errorbar(inData[okf, 0], inData[okf, 1] * 1000,
                    yerr=inData[okf, 2] * 1000, fmt='.', color='k',
                    markersize=5, capsize=2, elinewidth=1, label='Data errors')
    axs[0].grid(True, which='both', linestyle='-', linewidth=0.5, alpha=0.7)
    axs[0].set_xlim(amin, amax)
    axs[0].set_xticks(fticks)
    axs[0].xaxis.set_major_formatter(ScalarFormatter())
    axs[0].ticklabel_format(style='plain', axis='x', useOffset=False)
    axs[0].set_xlabel('Frequency (Hz)')
    axs[0].set_ylabel('Slowness (ms/m)')
    axs[0].set_title(f"{ftxt} mode", fontweight='normal')
    axs[0].legend(loc='best')

    # Subplot 2: Velocity
    axs[1].semilogx(inData_v[okf, 0], inData_v[okf, 1], color=(0.1, 0.1, 0.9),
                    linestyle=':')
    yerr_asym = np.array([inData_v[okf, 2], inData_v[okf, 3]])
    axs[1].errorbar(inData_v[okf, 0], inData_v[okf, 1], yerr=yerr_asym,
                    fmt='.', color=(0.0, 0.0, 0.9), markersize=5, capsize=2,
                    elinewidth=1)
    axs[1].grid(True, which='both', linestyle='-', linewidth=0.5, alpha=0.7)
    axs[1].set_xlim(amin, amax)
    axs[1].set_xticks(fticks)
    axs[1].xaxis.set_major_formatter(ScalarFormatter())
    axs[1].ticklabel_format(style='plain', axis='x', useOffset=False)
    axs[1].set_xlabel('Frequency (Hz)')
    axs[1].set_ylabel('Velocity (m/s)')

    plt.tight_layout()
    plt.savefig(f"{PREF}Data_{ftxt}.png", dpi=600, bbox_inches='tight')
    plt.close(fig)
    print(f"[+] SUCCESS: Saved to {PREF}Data_{ftxt}.png")


# -----------------------------------------------------------------------------
# Plot ellipticity curve
def plot_ellipticity(ftxt: str, filepath: str) -> None:
    """
    Load and plot Rayleigh wave ellipticity (ELL)

    Args:
        ftxt (str): Figure title.
        filepath (str): Path to the file with data.
    Returns:
        None.
    """
    # Load data
    if not os.path.exists(filepath):
        print(f"[!] ERROR: File {filepath} does not exist.")
        return None
    inData = np.loadtxt(filepath, skiprows=1)

    # Prepare for log-plot
    inData_log = np.zeros((len(inData), 4))
    inData_log[:, 0] = inData[:, 0]
    inData_log[:, 1] = inData[:, 1]
    inData_log[:, 2] = inData[:, 1] - (inData[:, 1] / inData[:, 2])
    inData_log[:, 3] = -inData[:, 1] + (inData[:, 1] * inData[:, 2])

    # Get used frequencies
    okf = inData[:, 2] > 0
    if not np.any(okf):
        print(f"[!] WARNING: No valid data to plot in file {filepath}.")
        return None

    # Get frequency limiths
    min_f = np.min(inData[okf, 0])
    max_f = np.max(inData[okf, 0])
    amin, amax, fticks = fscale(min_f, max_f)

    # Start figure
    plt.rcParams.update({
        'font.size': 8,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'axes.labelsize': 9,
        'axes.titlesize': 9,
        'legend.fontsize': 7,
    })
    fig, ax = plt.subplots(figsize=(16/2.54, 12/2.54), facecolor='w')

    # log-log plot
    ax.loglog(inData_log[okf, 0], inData_log[okf, 1], 'k:', label='Data line')
    yerr = [inData_log[okf, 2], inData_log[okf, 3]]
    ax.errorbar(inData_log[okf, 0], inData_log[okf, 1], yerr=yerr, fmt='.',
                color='k', markersize=5, capsize=2, label='Data errors')
    ax.grid(True, which='both', linestyle='-', linewidth=0.5, alpha=0.7)
    ax.set_xlim(amin, amax)
    ax.set_xticks(fticks)
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.yaxis.set_major_formatter(ScalarFormatter())
    ax.yaxis.set_minor_formatter(FormatStrFormatter("%.1f"))
    ax.ticklabel_format(style='plain', axis='both', useOffset=False)
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('Ellipticity')
    ax.set_title(ftxt, fontweight='normal')
    ax.legend(loc='best')

    plt.tight_layout()
    plt.savefig(f"{PREF}Data_ELL.png", dpi=600, bbox_inches='tight')
    plt.close(fig)
    print(f"[+] SUCCESS: Saved to {PREF}Data_ELL.png")


# -----------------------------------------------------------------------------
# Plot ellipticity angle curve
def plot_ellipticity_angle(ftxt: str, filepath: str) -> None:
    """
    Load and plot Rayleigh wave ellipticity angle (ELA)

    Args:
        ftxt (str): Figure title.
        filepath (str): Path to the file with data.
    Returns:
        None.
    """
    # Load data
    if not os.path.exists(filepath):
        print(f"[!] ERROR: File {filepath} does not exist.")
        return None
    inData = np.loadtxt(filepath, skiprows=1)

    # Convert to Rayleigh wave ellipticity
    inData_e = np.zeros((len(inData), 6))
    inData_e[:, 0] = inData[:, 0]
    inData_e[:, 1] = np.abs(np.tan(np.radians(inData[:, 1])))
    inData_e[:, 2] = np.abs(np.tan(np.radians(np.clip(inData[:, 1]
                            + inData[:, 2], None, 90))))
    inData_e[:, 3] = np.abs(np.tan(np.radians(np.clip(inData[:, 1]
                            - inData[:, 2], -90, None))))
    inData_e[:, 4] = inData_e[:, 1]-np.minimum(inData_e[:, 2], inData_e[:, 3])
    inData_e[:, 5] = -inData_e[:, 1]+np.maximum(inData_e[:, 2], inData_e[:, 3])
    inData_e[:, 5] = np.clip(inData_e[:, 5], None, 100)

    # Get used frequencies
    okf = inData[:, 2] > 0
    if not np.any(okf):
        print(f"[!] WARNING: No valid data to plot in file {filepath}.")
        return None

    # Get frequency limiths
    min_f = np.min(inData[okf, 0])
    max_f = np.max(inData[okf, 0])
    amin, amax, fticks = fscale(min_f, max_f)

    # Start figure
    plt.rcParams.update({
        'font.size': 8,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'axes.labelsize': 9,
        'axes.titlesize': 9,
        'legend.fontsize': 7,
    })
    fig, axs = plt.subplots(2, 1, figsize=(16/2.54, 12/2.54), facecolor='w')

    # Subplot 1: Angle  (semilogx)
    axs[0].semilogx(inData[okf, 0], inData[okf, 1], 'k:', label='Data line')
    axs[0].errorbar(inData[okf, 0], inData[okf, 1], yerr=inData[okf, 2],
                    fmt='.', color='k', markersize=5, capsize=2,
                    label='Data errors')
    axs[0].set_ylim(-90, 90)
    axs[0].set_yticks(np.arange(-90, 91, 30))
    axs[0].grid(True, which='both', linestyle='-', linewidth=0.5, alpha=0.7)
    axs[0].set_xlim(amin, amax)
    axs[0].set_xticks(fticks)
    axs[0].xaxis.set_major_formatter(ScalarFormatter())
    axs[0].ticklabel_format(style='plain', axis='x', useOffset=False)
    axs[0].set_ylabel('Angle (deg)')
    axs[0].set_title(ftxt, fontweight='normal')
    axs[0].legend(loc='best')

    # Subplot 2: Ellipticity (log-log)
    axs[1].loglog(inData_e[okf, 0], inData_e[okf, 1], color=(0.1, 0.1, 0.9),
                  linestyle=':')
    yerr_e = [inData_e[okf, 4], inData_e[okf, 5]]
    axs[1].errorbar(inData_e[okf, 0], inData_e[okf, 1], yerr=yerr_e, fmt='.',
                    color=(0.1, 0.1, 0.9), markersize=5, capsize=2)
    axs[1].grid(True, which='both', linestyle='-', linewidth=0.5, alpha=0.7)
    axs[1].set_xlim(amin, amax)
    axs[1].set_xticks(fticks)
    axs[1].xaxis.set_major_formatter(ScalarFormatter())
    axs[1].ticklabel_format(style='plain', axis='x', useOffset=False)
    axs[1].set_xlabel('Frequency (Hz)')
    axs[1].set_ylabel('Ellipticity')

    plt.tight_layout()
    plt.savefig(f"{PREF}Data_ELA.png", dpi=600, bbox_inches='tight')
    plt.close(fig)
    print(f"[+] SUCCESS: Saved to {PREF}Data_ELA.png")


# -----------------------------------------------------------------------------
# Main function
def main():
    """
    Main function - plot dispersion and ellipticity curves.

    - Read used data (USE key) from FILENAMES.
    - Plot data.
    - Save resulting figures.
    """
    print("-" * 50)
    # Processing Rayleigh and Love dispersion curves
    for w in range(2):  # w=0 (Rayleigh), w=1 (Love)
        for m in range(4):  # m=0..3 (mode 0, 1, 2, 3)
            if USE[w, m] == 1:
                mode_prefix = 'R' if w == 0 else 'L'
                ftxt = f"{mode_prefix}{m}"
                print(f"[*] Plot dispersion curve: {ftxt}")
                plot_dispersion(ftxt, FILENAMES[ftxt])

    # Processing Rayleigh wave ellipticity
    if USE[2, 0] == 1:
        print("[*] Plot Rayleigh wave ellipticity (ELL)")
        plot_ellipticity('Rayleigh wave ellipticity', FILENAMES['ELL'])

    # Processing Rayleigh wave ellipticity angle
    if USE[2, 1] == 1:
        print("[*] Plot Rayleigh wave ellipticity angle (ELA)")
        plot_ellipticity_angle('Rayleigh ellipticity angle', FILENAMES['ELA'])

    print("[*] SUCCESS: All done")


# -----------------------------------------------------------------------------
# Entry point
if __name__ == "__main__":
    main()
