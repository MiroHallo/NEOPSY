#!/usr/bin/env python3
# =============================================================================
# PLOT RESULTS OF THE MULTIZONAL TRANSDIMENSIONAL INVERSION (NEOPSY) - LIBRARY
#
# Author: Miroslav HALLO, Kyoto University
# E-mail: hallo.miroslav.2a@kyoto-u.ac.jp
# Revision 2026/05: First version
# Tested with: Python 3.12.3, NumPy 2.4.5
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
from pathlib import Path

import numpy as np


# =============================================================================
# FUNCTIONS
# =============================================================================

# Automatic frequency axis limits
def fscale(fmin: float, fmax: float) -> tuple[float, float, np.ndarray]:
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
def qwl30(depth: np.ndarray, vs: np.ndarray,
          d30: float = 30.0) -> tuple[float, float]:
    """
    Computes the Vs30 velocity of a 1D layered velocity model.

    Args:
        depth (np.ndarray): Top depth of layers [m].
        vs (np.ndarray): Shear wave velocity [m/s].
        d30 (float): Depth of interest (standard is 30.0) [m].
    Returns:
        tuple[float, float]: f30 (frequency) and Vs30 velocity [m/s].
    """
    n_layers = len(depth)

    # Calculate layer thicknesses
    thicknesses = np.zeros(n_layers)
    thicknesses[:-1] = depth[1:] - depth[:-1]
    thicknesses[-1] = 1e9

    # Find the index of the layer that contains the target depth d30
    idx = np.searchsorted(depth, d30, side='right') - 1
    idx = int(np.clip(idx, 0, n_layers - 1))

    # Compute cumulative vertical travel-times for all layer interfaces
    vs_pro = np.where(vs == 0.0, 1e-9, vs)
    t0 = np.append(0.0, np.cumsum(thicknesses[:-1] / vs_pro[:-1]))

    # Compute total travel time down to exactly d30 meters
    total_time = t0[idx] + (np.abs(d30 - depth[idx]) / vs_pro[idx])

    if total_time <= 0.0:
        return 0.0, 0.0

    # Calculate quarter-wavelength frequency and average Vs30 velocity
    f30 = 1.0 / (4.0 * total_time)
    vs30 = d30 / total_time

    return float(f30), float(vs30)


# -----------------------------------------------------------------------------
def respqwl(depth: np.ndarray, vs: np.ndarray, rho: np.ndarray,
            f: np.ndarray) -> tuple:
    """
    Computes the amplification factor by quarter-wavelength average velocity.

    Args:
        depth (np.ndarray): Top depth of layers [m].
        vs (np.ndarray): Shear wave velocity [m/s].
        rho (np.ndarray): Density [kg/m3].
        f (np.ndarray): Frequency discretization [Hz].
    Returns:
        tuple: z_f - QWL depth (np.ndarray).
               vs_f - QWL velocity (np.ndarray).
               rho_f - QWL density (np.ndarray).
               af - The amplification factor (np.ndarray).
               z2_f - 2nd QWL depth (np.ndarray).
               vs2_f - 2nd QWL velocity (np.ndarray).
               ic_f - QWL Impedance contrast (np.ndarray).
    """
    vs_c = vs[-1]    # Reference (source) shear wave velocity [m/s]
    rho_c = rho[-1]  # Reference (source) density [kg/m3]

    # Depth search discretization [m]
    z = np.logspace(-2.0, 5.0, 10000)

    n_layers = len(depth)
    n_freq = len(f)

    # Calculate layer thicknesses
    thicknesses = np.zeros(n_layers)
    thicknesses[:-1] = depth[1:] - depth[:-1]
    thicknesses[-1] = 1e9

    # Find layer indexes for all discrete depths z
    iz = np.searchsorted(depth, z, side='right') - 1
    iz = np.clip(iz, 0, n_layers - 1)

    # Compute vertical travel-times vector
    vs_pro = np.where(vs == 0.0, 1e-9, vs)
    t0 = np.append(0.0, np.cumsum(thicknesses[:-1] / vs_pro[:-1]))
    t = t0[iz] + np.abs(z - depth[iz]) / vs_pro[iz]

    # Find quarter-wavelength depth and velocity
    t_pro = np.where(t == 0.0, 1e-9, t)
    vs_qwl = z / t_pro

    # Frequency depth matching
    min_i = np.zeros(n_freq, dtype=int)
    for i in range(n_freq):
        target = vs_qwl / (4.0 * f[i])
        min_i[i] = np.argmin(np.abs(z - target))

    z_f = z[min_i]
    vs_f = vs_qwl[min_i]
    iz_f = iz[min_i]

    # Find quarter-wavelength density
    rho_0 = np.append(0.0, np.cumsum(rho[:-1] * thicknesses[:-1]))
    rho_f = (rho_0[iz_f] + np.abs(z_f - depth[iz_f]) * rho[iz_f]) / z_f

    # Amplification factor
    denom_a = vs_f * rho_f
    denom_a_pro = np.where(denom_a == 0.0, 1e-9, denom_a)
    af = np.sqrt((vs_c * rho_c) / denom_a_pro)

    # Find quarter-wavelength Impedance contrast
    z2_f = np.zeros(n_freq)
    vs2_f = np.zeros(n_freq)
    for i in range(n_freq):
        idx_start = min_i[i]
        t_tmp = t[idx_start:] - t[idx_start]
        z_tmp = z[idx_start:] - z[idx_start]
        with np.errstate(divide='ignore', invalid='ignore'):
            vs_qwl_tmp = z_tmp / t_tmp

        # Match secondary QWL condition
        target_tmp = vs_qwl_tmp / (4.0 * f[i])
        min_c = np.nanargmin(np.abs(z_tmp - target_tmp))

        z2_f[i] = z[idx_start + min_c]
        vs2_f[i] = vs_qwl_tmp[min_c]

    # Impedance contrast
    ic_f = vs_f / vs2_f

    return z_f, vs_f, rho_f, af, z2_f, vs2_f, ic_f


# -----------------------------------------------------------------------------
def respsh(depth: np.ndarray, vs: np.ndarray, rho: np.ndarray,
           fox: np.ndarray, f: np.ndarray,
           th: float = 0.0) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Returns seismic response of damped soil layers on rock to incident SH-wave.

    Args:
        depth (np.ndarray): Top depth of layers [m].
        vs (np.ndarray): Shear wave velocity [m/s].
        rho (np.ndarray): Density [kg/m3].
        fox (np.ndarray): Layer damping ratio.
        f (np.ndarray): Frequency discretization [Hz].
        th (float): Wave incidence angle [deg] (0-60, beware higher angles).
    Returns:
        tuple[np.ndarray, np.ndarray, np.ndarray]: h1, h2, h3 transfer func.
              h1 - "surface to incident SH" transfer function.
              h2 - "surface to borehole" transfer function.
              h3 - "surface to rock outcrop" transfer function.
    """
    n_layers = len(depth)
    n_freq = len(f)

    # Calculate layer thicknesses
    thicknesses = np.zeros(n_layers)
    thicknesses[:-1] = depth[1:] - depth[:-1]
    thicknesses[-1] = 1e9

    # Complex shear wave velocity (1j is the imaginary unit)
    vs_s = vs * (1.0 + 1j * fox)

    # Horizontal and vertical incident shear wave slowness
    h_slow = np.sin(np.radians(th)) / vs[-1]
    v_slow = np.sqrt(((1.0 + 0j) / (vs_s ** 2)) - (h_slow ** 2))

    # Complex impedance ratio
    vs_s_imp = (vs_s ** 2) * v_slow
    alpha_s = (rho[:-1] * vs_s_imp[:-1]) / (rho[1:] * vs_s_imp[1:])

    # Allocate complex matrices for up-going (A) and down-going (B) waves
    a_mat = np.ones((n_layers, n_freq), dtype=np.complex128)
    b_mat = np.ones((n_layers, n_freq), dtype=np.complex128)

    # Recursive loop over layers
    for m in range(n_layers - 1):
        # Complex wave number times thickness
        ks_h = 2.0 * np.pi * f * thicknesses[m] * v_slow[m]

        # Exponential terms for up and down waves
        exp_p = np.exp(1j * ks_h)
        exp_m = np.exp(-1j * ks_h)

        # Matrix recursive update
        a_mat[m + 1, :] = (0.5 * a_mat[m, :] * (1.0 + alpha_s[m]) * exp_p +
                           0.5 * b_mat[m, :] * (1.0 - alpha_s[m]) * exp_m)
        b_mat[m + 1, :] = (0.5 * a_mat[m, :] * (1.0 - alpha_s[m]) * exp_p +
                           0.5 * b_mat[m, :] * (1.0 + alpha_s[m]) * exp_m)

    # Transfer functions
    # h1: "surface to incident SH"
    h1 = 2.0 / a_mat[-1, :]

    # h2: "surface to borehole"
    h2 = (a_mat[0, :] + b_mat[0, :]) / (a_mat[-1, :] + b_mat[-1, :])

    # h3: "surface to rock outcrop"
    h3 = 1.0 / a_mat[-1, :]

    return h1, h2, h3


# -----------------------------------------------------------------------------
# DUSK colormap for Probability Density Function (PDF)
def dusk() -> np.ndarray:
    """
    DUSK colormap for Probability Density Function (PDF).

    Returns:
        np.ndarray: color_matrix (RGB).
    """
    color_matrix = np.array([
        [1,  1,  0.99],
        [1,  1,  0.930929686420777],
        [1,  1,  0.861859372841553],
        [1,  1,  0.792789059262330],
        [1,  1,  0.723718745683105],
        [1,  1,  0.654648432103882],
        [1,  1,  0.585578118524658],
        [1,  1,  0.516507804945435],
        [1,  0.989791948919503,  0.5],
        [1,  0.976377952755906,  0.5],
        [1,  0.962963956592309,  0.5],
        [1,  0.949549960428711,  0.5],
        [1,  0.936135964265114,  0.5],
        [1,  0.922721968101517,  0.5],
        [1,  0.909307971937920,  0.5],
        [1,  0.895893975774323,  0.5],
        [1,  0.882479979610726,  0.5],
        [1,  0.869065983447129,  0.5],
        [1,  0.855651987283532,  0.5],
        [1,  0.842237991119935,  0.5],
        [1,  0.828823994956337,  0.5],
        [1,  0.815409998792740,  0.5],
        [1,  0.801996002629143,  0.5],
        [1,  0.788582006465546,  0.5],
        [1,  0.775168010301949,  0.5],
        [1,  0.761754014138352,  0.5],
        [1,  0.748340017974755,  0.5],
        [1,  0.734926021811158,  0.5],
        [1,  0.721512025647561,  0.5],
        [1,  0.708098029483964,  0.5],
        [1,  0.694684033320367,  0.5],
        [1,  0.681270037156770,  0.5],
        [1,  0.667856040993172,  0.5],
        [1,  0.654442044829575,  0.5],
        [1,  0.641028048665978,  0.5],
        [1,  0.627614052502381,  0.5],
        [1,  0.614200056338784,  0.5],
        [1,  0.600786060175187,  0.5],
        [1,  0.587372064011590,  0.5],
        [1,  0.573958067847993,  0.5],
        [1,  0.560544071684396,  0.5],
        [1,  0.547130075520799,  0.5],
        [1,  0.533716079357201,  0.5],
        [1,  0.520302083193604,  0.5],
        [1,  0.506888087030007,  0.5],
        [0.987188265346430,  0.5,  0.5],
        [0.960853764516894,  0.5,  0.5],
        [0.934519263687357,  0.5,  0.5],
        [0.908184762857820,  0.5,  0.5],
        [0.881850262028283,  0.5,  0.5],
        [0.855515761198747,  0.5,  0.5],
        [0.829181260369210,  0.5,  0.5],
        [0.802846759539673,  0.5,  0.5],
        [0.776512258710136,  0.5,  0.5],
        [0.750177757880599,  0.5,  0.5],
        [0.723843257051063,  0.5,  0.5],
        [0.697508756221526,  0.5,  0.5],
        [0.671174255391989,  0.5,  0.5],
        [0.644839754562453,  0.5,  0.5],
        [0.618505253732916,  0.5,  0.5],
        [0.592170752903379,  0.5,  0.5],
        [0.565836252073842,  0.5,  0.5],
        [0.539501751244305,  0.5,  0.5],
        [0.513167250414768,  0.5,  0.5],
        [0.486832749585232,  0.5,  0.5],
        [0.460498248755695,  0.5,  0.5],
        [0.434163747926158,  0.5,  0.5],
        [0.407829247096621,  0.5,  0.5],
        [0.381494746267085,  0.5,  0.5],
        [0.355160245437548,  0.5,  0.5],
        [0.328825744608011,  0.5,  0.5],
        [0.302491243778474,  0.5,  0.5],
        [0.276156742948937,  0.5,  0.5],
        [0.249822242119401,  0.5,  0.5],
        [0.223487741289864,  0.5,  0.5],
        [0.197153240460327,  0.5,  0.5],
        [0.170818739630790,  0.5,  0.5],
        [0.144484238801254,  0.5,  0.5],
        [0.118149737971717,  0.5,  0.5],
        [0.091815237142180,  0.5,  0.5],
        [0.065480736312643,  0.5,  0.5],
        [0.039146235483107,  0.5,  0.5],
        [0.012811734653570,  0.5,  0.5],
        [0,  0.493111912969993,  0.5],
        [0,  0.479697916806396,  0.5],
        [0,  0.466283920642799,  0.5],
        [0,  0.452869924479202,  0.5],
        [0,  0.439455928315605,  0.5],
        [0,  0.426041932152008,  0.5],
        [0,  0.412627935988410,  0.5],
        [0,  0.399213939824813,  0.5],
        [0,  0.385799943661216,  0.5],
        [0,  0.372385947497619,  0.5],
        [0,  0.358971951334022,  0.5],
        [0,  0.345557955170425,  0.5],
        [0,  0.332143959006828,  0.5],
        [0,  0.318729962843231,  0.5],
        [0,  0.305315966679634,  0.5],
        [0,  0.291901970516037,  0.5],
        [0,  0.278487974352439,  0.5],
        [0,  0.265073978188842,  0.5],
        [0,  0.251659982025245,  0.5],
        [0,  0.238245985861648,  0.5],
        [0,  0.224831989698051,  0.5],
        [0,  0.211417993534454,  0.5],
        [0,  0.198003997370857,  0.5],
        [0,  0.184590001207260,  0.5],
        [0,  0.171176005043663,  0.5],
        [0,  0.157762008880066,  0.5],
        [0,  0.144348012716468,  0.5],
        [0,  0.130934016552871,  0.5],
        [0,  0.117520020389274,  0.5],
        [0,  0.104106024225677,  0.5],
        [0,  0.0906920280620800,  0.5],
        [0,  0.0772780318984829,  0.5],
        [0,  0.0638640357348858,  0.5],
        [0,  0.0504500395712887,  0.5],
        [0,  0.0370360434076916,  0.5],
        [0,  0.0236220472440945,  0.5],
        [0,  0.0102080510804974,  0.5],
        [0,  0,  0.483492195054566],
        [0,  0,  0.414421881475342],
        [0,  0,  0.345351567896118],
        [0,  0,  0.276281254316895],
        [0,  0,  0.207210940737671],
        [0,  0,  0.138140627158447],
        [0,  0,  0.069070313579224],
        [0,  0,  0]
    ])

    return color_matrix


# -----------------------------------------------------------------------------
# Read POP headers (NEOPSY format)
def read_pop_headers(pop_prefix: Path) -> SimpleNamespace:
    """
    Reads inversion header metadata (NEOPSY format) and prepare for plotting.

    Args:
        pop_prefix (Path): Path prefix for population files.
    Returns:
        SimpleNamespace: Header metadata.
    """
    header_path = Path(f"{pop_prefix}_headers.txt")

    if not header_path.exists():
        raise FileNotFoundError(f"Missing header file: {header_path}")

    with open(header_path, 'r', encoding='utf-8') as f:
        next(f)  # Skip description line
        n_bins = int(next(f).strip())
        n_depth = int(next(f).strip())
        models_count = int(next(f).strip())
        next(f)  # Skip dummy line

        # Read all_bins matrix
        all_bins_list = []
        for _ in range(5):
            line = next(f).strip()
            all_bins_list.append(np.fromstring(line, sep=' '))
        # Transpose
        all_bins = np.array(all_bins_list).T

        # Read depth bins
        d_bins = np.fromstring(next(f).strip(), sep=' ')
        log_d_bins = np.fromstring(next(f).strip(), sep=' ')

    # Prepare P1 vectors
    d_step = d_bins[1] - d_bins[0]
    d_bins_p1 = np.append(d_bins, d_bins[-1] + d_step)

    # Prepare bin edges
    all_bins_p1 = np.zeros((n_bins + 1, 4))
    for i in range(4):
        step = all_bins[1, i] - all_bins[0, i]
        all_bins_p1[:, i] = np.append(all_bins[:, i], all_bins[-1, i] + step)

    return SimpleNamespace(
        n_bins=n_bins,
        n_depth=n_depth,
        models_count=models_count,
        all_bins=all_bins,
        d_bins=d_bins,
        log_d_bins=log_d_bins,
        d_bins_p1=d_bins_p1,
        all_bins_p1=all_bins_p1
    )


# -----------------------------------------------------------------------------
# Read DC, ELL, ELA curves (NEOPSY format)
def read_curves(file: Path) -> tuple[np.ndarray, np.ndarray]:
    """
    Read DC, ELL, ELA curves (NEOPSY one file format).

    Args:
        file (Path): Path to the file with curves.
    Returns:
        tuple: data (np.ndarray), f_n (np.ndarray).
    """
    with open(file, 'r', encoding='utf-8') as f:
        next(f)  # Skip description line
        f_n = np.fromstring(next(f).strip(), sep=' ', dtype=int)
        num_curves = len(f_n)
        max_fn = max(f_n)
        data = np.zeros((max_fn, 4, num_curves))
        for m in range(num_curves):
            n_points = f_n[m]
            if n_points <= 0:
                continue
            next(f)  # Skip curve header line
            for i in range(n_points):
                line = next(f).strip()
                data[i, :, m] = np.fromstring(line, sep=' ')
    return data, f_n


# -----------------------------------------------------------------------------
# Read input data (DC and ELL curves, NEOPSY format)
def read_input_data(datafile: Path) -> SimpleNamespace:
    """
    Reads input data (DC and ELL curves) and calculates sigmas.

    Args:
        datafile (Path): Path to file in_data.txt.
    Returns:
        SimpleNamespace: Input (measured) curves.
    """
    if not datafile.exists():
        raise FileNotFoundError(f"Input data file not found: {datafile}")

    # Read DC, ELL, ELA curves
    data, f_n = read_curves(datafile)

    # Calculate sigmas and masks
    max_fn = max(f_n)
    num_curves = len(f_n)
    data_s = np.zeros((max_fn, num_curves))
    okf = np.zeros((max_fn, num_curves), dtype=bool)
    for m in range(num_curves):
        n = f_n[m]
        if n > 0:
            sigma = data[:n, 2, m]
            valid_idx = sigma != 0
            data_s[:n, m][valid_idx] = 1.0 / sigma[valid_idx]
            okf[:n, m] = valid_idx

    return SimpleNamespace(
        f_n=f_n,
        data=data,
        data_s=data_s,
        okf=okf
    )


# -----------------------------------------------------------------------------
# Read model data (DC and ELL curves, NEOPSY format)
def read_model_data(datafile: Path) -> SimpleNamespace:
    """
    Reads model data (DC and ELL curves).

    Args:
        datafile (Path): Path to file with DC, ELL, ELA curves.
    Returns:
        SimpleNamespace: Model (resulting) curves.
    """
    if not datafile.exists():
        raise FileNotFoundError(f"Model data file not found: {datafile}")

    # Read DC, ELL, ELA curves
    data, f_n = read_curves(datafile)

    return SimpleNamespace(
        f_n=f_n,
        data=data
    )


# -----------------------------------------------------------------------------
# Read structural 1D velocity model (NEOPSY format)
def read_velocity_model(model_file: Path, depth_max: float) -> SimpleNamespace:
    """
    Read structural velocity model (NEOPSY) and prepares a staircase profile.

    Args:
        model_file (Path): Path to the model file (.txt).
        depth_max (float): Maximum depth for the last layer in the plot.
    Returns:
        SimpleNamespace: Velocity model fields with staircase plotting vectors.
    """
    if not model_file.exists():
        raise FileNotFoundError(f"Velocity model file not found: {model_file}")

    with open(model_file, 'r', encoding='utf-8') as f:
        n_layers = int(next(f).strip())
        vr = float(next(f).strip())
        tmp_mod_list = []
        for _ in range(5):
            line = next(f).strip()
            tmp_mod_list.append(np.fromstring(line, sep=' '))
        tmp_mod = np.array(tmp_mod_list).T  # Transpose

    # Prepare the velocity model for plots
    plot_mod = np.repeat(tmp_mod, 2, axis=0)
    thicknesses = tmp_mod[:, 0]
    cumulative_depths = np.cumsum(thicknesses)
    depths = np.zeros(2 * n_layers)
    # Set boundaries for intermediate layers
    for i in range(1, n_layers):
        depths[2 * i - 1] = cumulative_depths[i - 1]
        depths[2 * i] = cumulative_depths[i - 1]
    # Last point boundary
    if n_layers > 1:
        depths[-2] = cumulative_depths[-2]
    else:
        depths[-2] = 0.0
    depths[-1] = depth_max
    plot_mod[:, 0] = depths

    return SimpleNamespace(
        n_layers=n_layers,
        vr=vr,
        plot_mod=plot_mod,
        thicknesses=tmp_mod[:, 0],
        depth=plot_mod[:, 0],
        vs=plot_mod[:, 1],
        vp=plot_mod[:, 2],
        rho=plot_mod[:, 3],
        nu=plot_mod[:, 4]
    )


# -----------------------------------------------------------------------------
# Read synthetic target model (NEOPSY format)
def read_synthetic_model(synt_file: Path, depth_max: float) -> SimpleNamespace:
    """
    Reads the synthetic target model and prepares a staircase profile.

    Args:
        synt_file (Path): Path to the synthetic model file (.model).
        depth_max (float): Maximum depth for the last layer in the plot.
    Returns:
        SimpleNamespace: Synthetic model with staircase plotting vectors.
    """
    if not synt_file.exists():
        raise FileNotFoundError(f"Synthetic model file not found: {synt_file}")

    # Read model from a file
    with open(synt_file, 'r', encoding='utf-8') as f:
        next(f)  # Skip description line
        n_layers = int(next(f).strip())
        next(f)  # Skip description line
        # Read by line (thickness, Vp, Vs, rho)
        raw_layers = []
        for _ in range(n_layers-1):
            line = next(f).strip()
            raw_layers.append(np.fromstring(line, sep=' '))
        next(f)  # Skip description line
        line = next(f).strip()
        raw_layers.append(np.fromstring(line, sep=' '))
        tmp_mod = np.array(raw_layers)

    # Compute Poisson ratio nu
    vp = tmp_mod[:, 1]
    vs = tmp_mod[:, 2]
    denom = 2.0 * (vp**2 - vs**2)
    nu = np.zeros(n_layers)
    valid = denom != 0
    nu[valid] = (vp[valid]**2 - 2.0 * vs[valid]**2) / denom[valid]
    tmp_mod_full = np.zeros((n_layers, 5))
    tmp_mod_full[:, :4] = tmp_mod
    tmp_mod_full[:, 4] = nu

    # Prepare the velocity model for plots
    syn_mod = np.repeat(tmp_mod_full, 2, axis=0)
    thicknesses = tmp_mod_full[:, 0]
    cumulative_depths = np.cumsum(thicknesses)
    depths = np.zeros(2 * n_layers)
    # Set boundaries for intermediate layers
    for i in range(1, n_layers):
        depths[2 * i - 1] = cumulative_depths[i - 1]
        depths[2 * i] = cumulative_depths[i - 1]
    # Last point boundary
    if n_layers > 1:
        depths[-2] = cumulative_depths[-2]
    else:
        depths[-2] = 0.0
    depths[-1] = depth_max
    syn_mod[:, 0] = depths

    return SimpleNamespace(
        n_layers=n_layers,
        plot_mod=syn_mod,
        thicknesses=tmp_mod_full[:, 0],
        depth=syn_mod[:, 0],
        vp=syn_mod[:, 1],
        vs=syn_mod[:, 2],
        rho=syn_mod[:, 3],
        nu=syn_mod[:, 4]
    )


# -----------------------------------------------------------------------------
# Reads a Fortran binary ensemble (NEOPSY format)
def read_pop_binary_data(file: Path, f_n: np.ndarray, nm: int) -> np.ndarray:
    """
    Reads binary ensemble (NEOPSY format) and downsamples it to nm models.

    Args:
        file (Path): Path to the _data.bin file.
        f_n (np.ndarray): Vector with number of points per curve.
        nm (int): Target number of models to retain for plotting.
    Returns:
        np.ndarray: Downsampled data array of shape (nm, 2, sum_fn).
    """
    if not file.exists():
        raise FileNotFoundError(f"Ensemble file not found: {file}")

    sum_fn = int(np.sum(f_n))
    record_size = 2 * sum_fn

    # Read the entire file at once as float64 (Fortran double)
    raw_data = np.fromfile(file, dtype=np.float64)

    # Calculate how many full models are actually in the file
    count = len(raw_data) // record_size
    if count == 0:
        raise ValueError(f"Ensemble {file.name} contains no models.")

    bin_data = raw_data[:count * record_size].reshape(count, 2, sum_fn)
    nm = max(1, min(nm, count))
    step = max(1, count // nm)
    syn_data = bin_data[step - 1::step, :, :]
    syn_data = syn_data[:nm, :, :]

    return syn_data


# -----------------------------------------------------------------------------
# Read quarter-wavelength representation headers (NEOPSY format)
def read_qwl_headers(pop_prefix: Path) -> SimpleNamespace:
    """
    Reads quarter-wavelength (QWL) metadata headers and prepares for plotting.

    Args:
        pop_prefix (Path): Path prefix for population files.
    Returns:
        SimpleNamespace: QWL metadata and plotting bin arrays.
    """
    qwl_head_path = Path(f"{pop_prefix}_QWL_head.txt")

    if not qwl_head_path.exists():
        raise FileNotFoundError(f"Missing QWL header file: {qwl_head_path}")

    with open(qwl_head_path, 'r', encoding='utf-8') as f:
        next(f)  # Skip description line
        n_bins = int(next(f).strip())
        n_freq = int(next(f).strip())
        models_count = int(next(f).strip())
        freq = np.fromstring(next(f).strip(), sep=' ')
        vs_vec = np.fromstring(next(f).strip(), sep=' ')
        dep_vec = np.fromstring(next(f).strip(), sep=' ')
        imp_vec = np.fromstring(next(f).strip(), sep=' ')
        amp_vec = np.fromstring(next(f).strip(), sep=' ')
        n_x_bins = int(next(f).strip())
        vs30_bins = np.fromstring(next(f).strip(), sep=' ')

    # QWL bins for surf plots
    f_bins_w = np.log10(freq[1]) - np.log10(freq[0])
    f_bins_p1 = np.append(
        10.0 ** (np.log10(freq) - f_bins_w / 2.0),
        10.0 ** (np.log10(freq[-1]) + 0.5 * f_bins_w)
    )
    qwl_bins_p1 = np.zeros((n_bins + 1, 4))
    # VsVec
    step_vs = vs_vec[1] - vs_vec[0]
    qwl_bins_p1[:, 0] = np.append(vs_vec, vs_vec[-1] + step_vs)
    # DepVec
    s = np.log10(dep_vec[1]) - np.log10(dep_vec[0])
    qwl_bins_p1[:, 1] = np.append(dep_vec, 10.0 ** (np.log10(dep_vec[-1]) + s))
    # ImpVec
    step_imp = imp_vec[1] - imp_vec[0]
    qwl_bins_p1[:, 2] = np.append(imp_vec, imp_vec[-1] + step_imp)
    # AmpVec
    s = np.log10(amp_vec[1]) - np.log10(amp_vec[0])
    qwl_bins_p1[:, 3] = np.append(amp_vec, 10.0 ** (np.log10(amp_vec[-1]) + s))

    # QWL-frequency tics
    f_tick = np.array([
        0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0, 5.0,
        6.0, 7.0, 8.0, 9.0, 10.0, 20.0, 30.0, 40.0, 50.0
    ])
    f_tick_label = [
        '0.1', '', '0.3', '', '0.5', '', '', '', '', '1', '2', '3', '4', '5',
        '', '', '', '', '10', '20', '30', '40', '50'
    ]

    # QWL-depth tics
    d_tick = np.array([
        1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100,
        200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000,
        6000, 7000, 8000, 9000, 10000
    ])
    d_tick_label = [
        '1', '', '', '', '', '', '', '', '', '10', '', '', '', '', '', '',
        '', '', '100', '', '', '', '', '', '', '', '', '1000', '', '', '',
        '', '', '', '', '', '10000'
    ]

    # SH-amplification tics
    a_tick = np.array([
        0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,
        10.0, 20.0, 30.0, 40.0
    ])
    a_tick_label = [
        '0.5', '', '', '', '', '1', '2', '3', '4', '5', '', '', '', '', '10',
        '20', '30', '40'
    ]

    return SimpleNamespace(
        n_bins=n_bins,
        n_freq=n_freq,
        models_count=models_count,
        freq=freq,
        vs_vec=vs_vec,
        dep_vec=dep_vec,
        imp_vec=imp_vec,
        amp_vec=amp_vec,
        n_x_bins=n_x_bins,
        vs30_bins=vs30_bins,
        f_bins_p1=f_bins_p1,
        qwl_bins_p1=qwl_bins_p1,
        f_tick=f_tick,
        f_tick_label=f_tick_label,
        d_tick=d_tick,
        d_tick_label=d_tick_label,
        a_tick=a_tick,
        a_tick_label=a_tick_label
    )


# -----------------------------------------------------------------------------
# Read model QWL and SH amplification data (NEOPSY format)
def read_model_qwl(qwl_file: Path) -> SimpleNamespace:
    """
    Reads QWL model simulation data, VS30 values, and SH amplification.

    Args:
        qwl_file (Path): Path to the model QWL text file.
    Returns:
        SimpleNamespace: Object containing parsed vectors and matrices.
    """
    if not qwl_file.exists():
        raise FileNotFoundError(f"Model QWL file not found: {qwl_file}")

    with open(qwl_file, 'r', encoding='utf-8') as f:
        n_freq = int(next(f).strip())
        vs30 = np.fromstring(next(f).strip(), sep=' ')
        freq = np.fromstring(next(f).strip(), sep=' ')
        qwl_list = []
        for _ in range(4):
            line = next(f).strip()
            qwl_list.append(np.fromstring(line, sep=' '))
        qwl_data = np.array(qwl_list).T
        sh_list = []
        for _ in range(2):
            line = next(f).strip()
            sh_list.append(np.fromstring(line, sep=' '))
        # Transpose to shape (n_freq_m, 2)
        sh_data = np.array(sh_list).T

    return SimpleNamespace(
        n_freq=n_freq,
        vs30=vs30,
        freq=freq,
        qwl_data=qwl_data,
        sh_data=sh_data
    )


# -----------------------------------------------------------------------------
# Read 1D velocity model and compute QWL and SH amplification (NEOPSY format)
def compute_qwl(model_file: Path, ref_file: Path,
                freq: np.ndarray) -> SimpleNamespace:
    """
    Reads a velocity model, computes QWL, VS30, SH TF, and site amplification.

    Args:
        model_file (Path): Path to a velocity model (.model).
        ref_file (Path): Path to the reference rock velocity model.
        freq (np.ndarray): Frequencies for computation.
    Returns:
        SimpleNamespace: Derived QWL, Vs30, SH response, and ampl. factors.
    """
    if not model_file.exists():
        raise FileNotFoundError(f"Velocity model file not found: {model_file}")
    if not ref_file.exists():
        raise FileNotFoundError(f"Reference rock model not found: {ref_file}")

    # Read velocity model from a file
    with open(model_file, 'r', encoding='utf-8') as f:
        next(f)  # Skip description line
        n_layers = int(next(f).strip())
        next(f)  # Skip description line
        # Read by line (thickness, Vp, Vs, rho)
        raw_layers = []
        for _ in range(n_layers-1):
            line = next(f).strip()
            raw_layers.append(np.fromstring(line, sep=' '))
        next(f)  # Skip description line
        line = next(f).strip()
        raw_layers.append(np.fromstring(line, sep=' '))
        model = np.array(raw_layers)

    # Convert thickness to top depth
    thicknesses = model[:-1, 0]
    depths = np.append(0.0, np.cumsum(thicknesses))

    # Compute QWL
    z_f, vs_f, rho_f, af, z2_f, vs2_f, ic_f = respqwl(
        depths, model[:, 2], model[:, 3], freq
    )

    # Compute Vs30
    f30, vs30 = qwl30(depths, model[:, 2], 30.0)

    # Compute SH transfer functions with zero damping ratio (Fox)
    fox_zeros = np.zeros(n_layers)
    h1, h2, h3 = respsh(
        depths, model[:, 2], model[:, 3], fox_zeros, freq, 0.0
    )

    # Load reference rock profile (thickness and velocity)
    ref_mod = np.loadtxt(ref_file)

    # Convert reference thicknesses to top depth
    ref_thicknesses = ref_mod[:-1, 0]
    ref_depths = np.append(0.0, np.cumsum(ref_thicknesses))

    # Compute reference rock QWL (reference density = 2500 kg/m3)
    rho_ref = np.ones(len(ref_mod)) * 2500.0
    ref_zf, ref_vsf, ref_rhof, ref_a, ref_z2f, ref_vs2f, ref_icf = respqwl(
        ref_depths, ref_mod[:, 1], rho_ref, freq
    )

    # Amplification correction factor
    half_space_vs = model[-1, 2]
    denom_cor = ref_vsf * ref_rhof
    denom_cor_pro = np.where(denom_cor == 0.0, 1e-9, denom_cor)
    ampcor = np.sqrt((half_space_vs * 2500.0) / denom_cor_pro)

    return SimpleNamespace(
        model=model,
        depths=depths,
        z_f=z_f,
        vs_f=vs_f,
        rho_f=rho_f,
        af=af,
        z2_f=z2_f,
        vs2_f=vs2_f,
        ic_f=ic_f,
        f30=f30,
        vs30=vs30,
        h1=h1,
        h2=h2,
        h3=h3,
        ampcor=ampcor
    )
