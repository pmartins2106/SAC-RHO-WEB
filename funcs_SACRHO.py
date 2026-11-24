"""
SAC-RHO Core Analysis Module
Re-calculates RHO values using a fixed Seeding Activity Coefficient (SAC).
"""

from typing import Dict, List, Tuple
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import streamlit as st


def extract_curves(df_curve: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray]:
    """Drop empty cells and extract time (x) and signal (y) arrays."""
    df_clean = df_curve.dropna()
    xdata = np.array(df_clean.iloc[:, 0])
    ydata = np.array(df_clean.iloc[:, 1])
    return xdata, ydata


def plot_fit(
    p_x: np.ndarray,
    p_xfit: np.ndarray,
    p_y: np.ndarray,
    p_yfit: np.ndarray,
    curve_idx: int,
    yax_min: float,
    yax_max: float,
    colors: List[Tuple[float, float, float, float]],
    df: pd.DataFrame,
) -> None:
    """Plot experimental data points and fitted model curve."""
    plt.style.use("dark_background")
    label = df.columns[curve_idx + 1]
    plt.scatter(p_x, p_y, s=40, color=colors[curve_idx], label=label)
    plt.plot(p_xfit, p_yfit, color="g", linewidth=1)
    plt.xlabel("Time")
    plt.ylabel("Signal")
    plt.xlim(right=float(np.max(p_x)))
    plt.ylim([yax_min * 0.9, yax_max * 1.1])
    plt.legend(loc="lower center", bbox_to_anchor=(0.5, 1.05), ncol=4)


def segment_curves(
    time: np.ndarray,
    fluor: np.ndarray,
    f_min: float,
    f_max: float,
    segments: int,
    min_points: int,
) -> Dict[int, int]:
    """Divide normalized fluorescence curve into vertical segments."""
    idx_seg: Dict[int, int] = {}
    alpha = (fluor - f_min) / (f_max - f_min)
    alpha_seg = np.linspace(0, 1, segments + 1)

    for i, target_alpha in enumerate(alpha_seg):
        idx_temp = int(np.argmin(np.abs(alpha - target_alpha)))
        if i == 0:
            idx_seg[i] = 0
        elif idx_temp > idx_seg[i - 1] + min_points:
            idx_seg[i] = idx_temp
        else:
            idx_temp2 = idx_seg[i - 1] + min_points
            if idx_temp2 >= len(alpha):
                idx_seg[i] = len(alpha) - 1
            else:
                idx_seg[i] = int(np.argmin(np.abs(alpha[idx_temp2:] - target_alpha))) + idx_temp2

    return idx_seg


def clm_k0(
    x: np.ndarray,
    sac: float,
    log_a0k: float,
    baseline: float,
    endpoint: float,
    kb: float,
) -> np.ndarray:
    """Kinetic aggregation model with fixed baseline and endpoint."""
    a0k = 10**log_a0k
    return baseline + (endpoint - baseline) * (
        1 - 1 / ((kb + a0k) / (1 + a0k) * (np.exp(sac * (1 + a0k) * x) - 1) + 1)
    )


def handle_fit_error(stop_execution: bool = False) -> None:
    """Display standardized curve fitting error messages in Streamlit."""
    st.error("Error - curve fit failed.")
    st.error(
        "Possible causes: .ods file not properly filled, noisy data, "
        "or insufficient number of time points."
    )
    if stop_execution:
        st.stop()


def SACRHO2(
    df: pd.DataFrame,
    baseline: List[float],
    endpoint: List[float],
    segments: int,
    kac0: List[float],
    yfit0: List[float],
    n_curves: int,
    baseline0: float,
    kb: float,
    colors: List[Tuple[float, float, float, float]],
    backc: str,
    sac_mean: float,
    min_points: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, List[List[float]], np.ndarray, List[List[np.ndarray]]]:
    """Re-calculate RHO values for each curve using a fixed SAC value."""
    fig3 = plt.figure(facecolor=backc)
    results = np.array([], dtype=np.int64).reshape(0, 1)

    xseed0_all: List[List[float]] = []
    results_all: List[List[float]] = []
    t0_all: List[List[float]] = []
    sac_all: List[List[float]] = []
    f0_all: List[List[object]] = []
    mmt_all: List[List[object]] = []
    sigma_all: List[List[np.ndarray]] = []
    y_max_all: List[float] = []
    y_min_all: List[float] = []

    for curve in range(n_curves):
        f0_k: List[object] = []
        t0_k: List[float] = []
        sac_k: List[float] = []
        xseed0_k: List[float] = []
        results_k: List[object] = []
        mmt_k: List[object] = []
        sigma_k: List[np.ndarray] = []

        df_curve = df.iloc[:, [0, curve + 1]].astype(float)
        xdata, ydata = extract_curves(df_curve)
        y_max = float(np.max(ydata))
        y_min = float(np.min(ydata))

        idx_max = int(np.argmin(np.abs(y_max - ydata)))
        idx_max = idx_max if idx_max != 0 else min_points

        if y_max < y_min * 1.15:
            sac_k.append(0.0)
            sac_all.append(sac_k)
            results_k.append(0)
            f0_k.append("")
            t0_k.append(0.0)
            mmt_k.append("")
            xseed0_k.append(0.0)
        else:
            ydata_trim = ydata[:idx_max]
            xdata_trim = xdata[:idx_max]
            idx_seg = segment_curves(
                xdata_trim, ydata_trim, baseline[curve], endpoint[curve], segments, min_points
            )

            if all(np.diff(list(idx_seg.values())) < min_points):
                sac_k.append(0.0)
                sac_all.append(sac_k)
                results_k.append(0)
                f0_k.append("")
                t0_k.append(0.0)
                mmt_k.append("")
                xseed0_k.append(0.0)
                results_all.append(results_k)
            else:
                try:
                    t0_k.append(float(xdata_trim[0]))
                    xdata_k = xdata_trim[0 : idx_seg[1]] - t0_k[0]
                    ydata_k = ydata_trim[0 : idx_seg[1]]

                    ig = np.array([sac_mean, np.log10(yfit0[curve])])
                    bounds = ([sac_mean * 0.99, -50], [sac_mean * 1.01, 0])

                    parameters, covariance = curve_fit(
                        lambda time, a, b: clm_k0(
                            time, a, b, baseline[curve], endpoint[curve], kb
                        ),
                        xdata_k,
                        ydata_k,
                        p0=ig,
                        bounds=bounds,
                        maxfev=10000,
                        method="dogbox",
                    )

                    x00 = 10 ** parameters[1]
                    sigma_ab = np.sqrt(np.diagonal(covariance))

                    sigma_k.append(sigma_ab)
                    sac_k.append(float(parameters[0]))
                    f0_k.append(baseline[curve])
                    xseed0_k.append(float(x00))
                    mmt_k.append(float(x00))
                    results_k.append(float(xdata_trim[idx_seg[1]]))

                    x_ode = np.linspace(0, float(np.max(xdata_k)), 100)
                    y_max_all.append(y_max)
                    y_min_all.append(y_min)
                    y_ode = clm_k0(
                        x_ode, parameters[0], parameters[1], baseline[curve], endpoint[curve], kb
                    )
                    plot_fit(
                        xdata_trim,
                        x_ode + t0_k[0],
                        ydata_trim,
                        y_ode,
                        curve,
                        min(y_min_all),
                        max(y_max_all),
                        colors,
                        df,
                    )

                except RuntimeError:
                    results_k.append("")
                    t0_k.append(0.0)
                    sac_k.append(0.0)
                    f0_k.append("")
                    xseed0_k.append(0.0)
                    mmt_k.append("")
                    handle_fit_error(stop_execution=False)

        sac_all.append(sac_k)
        results_all.append(results_k)
        f0_all.append(f0_k)
        t0_all.append(t0_k)
        results = np.vstack((results, results_all[curve]))
        mmt_all.append(mmt_k)
        xseed0_all.append(xseed0_k)
        sigma_all.append(sigma_k)

    st.pyplot(fig3)

    return (
        results,
        np.array(sac_all),
        np.array(f0_all),
        np.array(xseed0_all),
        t0_all,
        np.array(mmt_all),
        sigma_all,
    )
