"""
Analysis page for the SAC-RHO Streamlit web application.
Determines Seeding Activity Coefficient (SAC) and RHO values from uploaded protein aggregation curves.
"""

from typing import Any, Dict, List, Optional, Tuple, Union
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import streamlit as st
from uncertainties import ufloat

from funcs_SACRHO import SACRHO2

MAX_COLS = 100
MAX_ROWS = 1000
DEFAULT_SEGMENTS = 4
TEMPLATE_PATH = "datasets/SAC-RHO_template.ods"


# ---------------------------------------------------------------------------
# Kinetic Models & Mathematical Helpers
# ---------------------------------------------------------------------------

def clm_model(
    x: np.ndarray, sac: float, log_y0: float, f0: float, fmax: float, kb: float = 0.0
) -> np.ndarray:
    """Kinetic model for protein aggregation curves."""
    y0 = 10**log_y0
    return f0 + (fmax - f0) * (
        1 - 1 / ((kb + y0) / (1 + y0) * (np.exp(sac * (1 + y0) * x) - 1) + 1)
    )


def clm_k0_model(
    x: np.ndarray, sac: float, log_a0k: float, baseline: float, endpoint: float, kb: float = 0.0
) -> np.ndarray:
    """Kinetic model with fixed baseline and endpoint."""
    a0k = 10**log_a0k
    return baseline + (endpoint - baseline) * (
        1 - 1 / ((kb + a0k) / (1 + a0k) * (np.exp(sac * (1 + a0k) * x) - 1) + 1)
    )


def clm_k0_a0k_model(
    x: np.ndarray, sac: float, baseline: float, endpoint: float, a0k: float, kb: float = 0.0
) -> np.ndarray:
    """Kinetic model with fixed baseline, endpoint, and seed ratio (a0k)."""
    return baseline + (endpoint - baseline) * (
        1 - 1 / ((kb + a0k) / (1 + a0k) * (np.exp(sac * (1 + a0k) * x) - 1) + 1)
    )


def format_scientific(val: Any, threshold: float = 0.1) -> str:
    """Format small non-zero numbers in scientific notation and others to 4 decimal places."""
    if isinstance(val, (int, float, np.number)) and abs(val) < threshold and val != 0:
        return f"{val:.2e}"
    if isinstance(val, str):
        return val
    return f"{val:.4f}"


# ---------------------------------------------------------------------------
# Data Loading & Curve Processing
# ---------------------------------------------------------------------------

def load_odf(file_source: Any) -> Tuple[Optional[pd.DataFrame], int, List[str], int]:
    """Load an .ods file into a cleaned Pandas DataFrame."""
    try:
        df_load = pd.read_excel(file_source, skiprows=4, engine="odf", header=None)
        df_load.columns = df_load.iloc[0]
        df_load = df_load.drop(0)
        df_load.dropna(how="all", axis=1, inplace=True)

        n_curves = df_load.shape[1] - 1
        n_points = df_load.shape[0]

        if n_curves > MAX_COLS or n_points > MAX_ROWS:
            st.warning(
                f"Warning: Your file exceeds the limit of {MAX_COLS} curves and "
                f"{MAX_ROWS - 5} time points. Only data within these limits will be analyzed."
            )
            df_load = pd.read_excel(
                file_source,
                skiprows=4,
                engine="odf",
                header=None,
                usecols=list(range(min(MAX_COLS, n_curves) + 1)),
                nrows=min(MAX_ROWS, n_points),
            )
            df_load.columns = df_load.iloc[0]
            df_load = df_load.drop(0)
            df_load.dropna(how="all", axis=1, inplace=True)
            n_curves = df_load.shape[1] - 1

        col_list = [f"Curve {i}" for i in range(1, n_curves + 1)]
        df_load.columns = np.concatenate((["Time"], col_list), axis=0)
        return df_load, n_curves, col_list, 0

    except Exception as exc:
        st.error(f"Failed to load file: {exc}")
        return None, 0, [], 1


def extract_curves(df_curve: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray]:
    """Drop NaN values and return time (x) and signal (y) numpy arrays."""
    df_clean = df_curve.dropna()
    xdata = np.array(df_clean.iloc[:, 0])
    ydata = np.array(df_clean.iloc[:, 1])
    return xdata, ydata


def segment_curves(
    fluor: np.ndarray, f_min: float, f_max: float, segments: int, min_points: int
) -> Dict[int, int]:
    """Segment normalized fluorescence curves into equal vertical intervals."""
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


def handle_fit_error(stop_execution: bool = True) -> None:
    """Display standardized curve fitting error messages in Streamlit."""
    st.error("Error - curve fit failed.")
    st.error(
        "Possible causes: .ods file not properly filled, noisy data, "
        "or insufficient number of time points."
    )
    if stop_execution:
        st.stop()


def fit_initial_parameters(
    df: pd.DataFrame, n_curves: int, kb: float, endpoint0: float, min_points: int
) -> Tuple[List[float], List[float], List[float], List[float]]:
    """Estimate baseline, endpoint, initial ka, and initial seed ratio for each curve."""
    baseline: List[float] = []
    endpoint: List[float] = []
    kac0: List[float] = []
    yfit0: List[float] = []
    error_count = 0

    for curve in range(n_curves):
        df_curve = df.iloc[:, [0, curve + 1]].astype(float)
        xdata, ydata = extract_curves(df_curve)
        y_max, y_min = float(np.max(ydata)), float(np.min(ydata))

        if y_max < y_min * 1.15:
            baseline.append(y_min)
            endpoint.append(y_max)
            kac0.append(0.0)
            yfit0.append(1e-24)
            continue

        idx_50 = int(np.argmin(np.abs(ydata - y_min - (y_max - y_min) * 0.5)))
        idx_50 = idx_50 if idx_50 != 0 else min_points

        ydata50 = ydata[:idx_50]
        xdata50 = xdata[:idx_50] - xdata[0]

        tmax = float(np.max(xdata))
        ymin_lim = min(y_min, 0.0)
        initial_guess = np.asarray([1.0 / tmax, -2.0, y_min, y_max])
        bounds = ([0, -50, ymin_lim, ymin_lim], [np.inf, 0, np.inf, np.inf])

        try:
            # Fit lower half for baseline and initial kinetic parameters
            params_low, _ = curve_fit(
                lambda x, p0, p1, p2, p3: clm_model(x, p0, p1, p2, p3, kb),
                xdata50,
                ydata50,
                p0=initial_guess,
                bounds=bounds,
                maxfev=10000,
            )
            kac0.append(params_low[0] / (10 ** params_low[1]))
            yfit0.append(10 ** params_low[1])
            baseline.append(params_low[2])

            # Fit upper half for endpoint
            ydata100 = ydata[idx_50:]
            xdata100 = xdata[idx_50:] - xdata[idx_50]
            params_high, _ = curve_fit(
                lambda x, p0, p1, p2, p3: clm_model(x, p0, p1, p2, p3, kb),
                xdata100,
                ydata100,
                p0=initial_guess,
                bounds=bounds,
                maxfev=10000,
            )
            if not endpoint0 or params_high[3] < 1.5 * endpoint0:
                endpoint.append(params_high[3])
            else:
                endpoint.append(endpoint0)

        except RuntimeError:
            error_count += 1

    if error_count > 0:
        handle_fit_error(stop_execution=True)

    return baseline, endpoint, kac0, yfit0


def run_sacrho_fitting(
    df: pd.DataFrame,
    baseline: List[float],
    endpoint: List[float],
    kac0: List[float],
    yfit0: List[float],
    n_curves: int,
    baseline0: float,
    kb: float,
    min_points: int,
) -> Tuple[np.ndarray, np.ndarray, List[List[np.ndarray]]]:
    """Determine individual ka and seed ratio values for each curve."""
    sac_all: List[List[float]] = []
    xseed0_all: List[List[float]] = []
    sigma_all: List[List[np.ndarray]] = []

    for curve in range(n_curves):
        sac_k: List[float] = []
        xseed0_k: List[float] = []
        sigma_k: List[np.ndarray] = []

        df_curve = df.iloc[:, [0, curve + 1]].astype(float)
        xdata, ydata = extract_curves(df_curve)
        y_max, y_min = float(np.max(ydata)), float(np.min(ydata))

        idx_max = int(np.argmin(np.abs(y_max - ydata)))
        idx_max = idx_max if idx_max != 0 else min_points

        if y_max < y_min * 1.15:
            sac_all.append([0.0])
            xseed0_all.append([0.0])
            sigma_all.append([])
            continue

        ydata_trim = ydata[:idx_max]
        xdata_trim = xdata[:idx_max]
        idx_seg = segment_curves(
            ydata_trim, baseline[curve], endpoint[curve], DEFAULT_SEGMENTS, min_points
        )

        if all(np.diff(list(idx_seg.values())) < min_points):
            sac_all.append([0.0])
            xseed0_all.append([0.0])
            sigma_all.append([])
            continue

        try:
            t0 = xdata_trim[0]
            xdata_k = xdata_trim[0 : idx_seg[1]] - t0
            ydata_k = ydata_trim[0 : idx_seg[1]]

            if not baseline0 or baseline0 > baseline[curve] or yfit0[curve] < 0.03:
                ig = np.array([kac0[curve] * yfit0[curve], np.log10(yfit0[curve])])
                bds = ([0, -60], [np.inf, 0])
                params, cov = curve_fit(
                    lambda time, a, b: clm_k0_model(
                        time, a, b, baseline[curve], endpoint[curve], kb
                    ),
                    xdata_k,
                    ydata_k,
                    p0=ig,
                    bounds=bds,
                    maxfev=10000,
                    method="dogbox",
                )
                x00 = 10 ** params[1]
            else:
                x00 = (baseline[curve] - baseline0) / (endpoint[curve] - baseline[curve])
                ig = np.array([kac0[curve] * yfit0[curve]])
                bds = ([0], [np.inf])
                params, cov = curve_fit(
                    lambda time, a: clm_k0_a0k_model(
                        time, a, baseline[curve], endpoint[curve], x00, kb
                    ),
                    xdata_k,
                    ydata_k,
                    p0=ig,
                    bounds=bds,
                    maxfev=10000,
                    method="dogbox",
                )

            sigma_ab = np.sqrt(np.diagonal(cov))
            sigma_k.append(sigma_ab)
            sac_k.append(params[0])
            xseed0_k.append(x00)

        except RuntimeError:
            sac_k.append(0.0)
            xseed0_k.append(0.0)
            handle_fit_error(stop_execution=False)

        sac_all.append(sac_k)
        xseed0_all.append(xseed0_k)
        sigma_all.append(sigma_k)

    st.success("Individual Curve Analysis")
    return np.array(sac_all), np.array(xseed0_all), sigma_all


# ---------------------------------------------------------------------------
# Streamlit UI Entrypoint
# ---------------------------------------------------------------------------

def analyse() -> None:
    """Render the SAC-RHO analysis page in Streamlit."""
    st.markdown(
        "<h2 style='text-align: center;'>🧪 SAC-RHO Web Tool</h2>",
        unsafe_allow_html=True,
    )

    st.info(
        "- Upload a dataset of protein aggregation curves (use the template)\n"
        "- Find the seeding activity coefficient (SAC) that characterizes the whole dataset\n"
        "- Find the RHO values for each progress curve\n"
        "- RHO should increase linearly with seed concentration in a log-log scale"
    )

    # Sidebar configuration
    with st.sidebar.expander("**Primary Nucleation Options**"):
        kb_mode = st.radio("Options", ("$k_b$ = 0", "Set a $k_b$ value:"), key="kb_radio")
        kb = 0.0 if kb_mode == "$k_b$ = 0" else st.number_input("$k_b$ = ", format="%e", value=1e-20)

    with st.sidebar.expander("**SAC Normalization**"):
        sac_mode = st.radio(
            "Options", ("Use the average $k_a$ value", "Set a SAC value:"), key="sac_radio"
        )
        sac_norm = (
            0.0
            if sac_mode == "Use the average $k_a$ value"
            else st.number_input("SAC = ", format="%e", value=0.1)
        )

    with st.sidebar.expander("**Reference Endpoint Values**"):
        ep_mode = st.radio(
            "Options",
            ("Use Fitted $F_\\infty^0$ Values", "Define Reference Endpoint:"),
            key="ep_radio",
        )
        endpoint0 = (
            0.0
            if ep_mode == "Use Fitted $F_\\infty^0$ Values"
            else st.number_input("Reference Endpoint = ", format="%2.1f", value=10.0)
        )

    with st.sidebar.expander("**Minimum Number of Datapoints**"):
        pts_mode = st.radio(
            "Options", ("Default ($N_{min}$ = 15)", "Configure:"), key="pts_radio"
        )
        min_points = (
            14
            if pts_mode == "Default ($N_{min}$ = 15)"
            else int(st.number_input("$N_{min}$ = ", format="%d", value=15, min_value=8, max_value=200)) - 1
        )

    # File Upload & Example Controls
    uploaded_file = st.file_uploader("Upload your dataset (.ods)")
    st.write("---------- OR ----------")
    run_example = st.checkbox("Run a prefilled example")

    df: Optional[pd.DataFrame] = None
    n_curves = 0
    col_list: List[str] = []
    flag_return = 1

    if uploaded_file is None and not run_example:
        if st.checkbox("Download template"):
            try:
                with open(TEMPLATE_PATH, "rb") as fp:
                    st.download_button(
                        label="Download",
                        data=fp,
                        file_name="SAC-RHO_template.ods",
                        mime="application/vnd.oasis.opendocument.spreadsheet",
                    )
                st.markdown(
                    "**Fill the template, save the file (keeping the .ods format), "
                    "and upload it into SAC-RHO.**"
                )
            except FileNotFoundError:
                st.error(f"Template file not found at `{TEMPLATE_PATH}`.")

    if run_example:
        df, n_curves, col_list, flag_return = load_odf(TEMPLATE_PATH)
    elif uploaded_file is not None:
        with st.spinner("Loading data..."):
            df, n_curves, col_list, flag_return = load_odf(uploaded_file)

    # Run Analysis & Display Results
    if df is not None and flag_return == 0:
        colormap = plt.cm.coolwarm
        colors = [colormap(i) for i in np.linspace(0, 1, n_curves)]
        backc = "black"
        baseline0 = 0.0

        with st.spinner("SAC-RHO running... Please wait"):
            baseline, endpoint, kac0, yfit0 = fit_initial_parameters(
                df, n_curves, kb, endpoint0, min_points
            )
            sac_all, xseed0_all, sigma_all = run_sacrho_fitting(
                df, baseline, endpoint, kac0, yfit0, n_curves, baseline0, kb, min_points
            )

        baseline_df = pd.DataFrame(baseline, columns=["Estimated Initial Signal"], index=col_list)
        endpoint_df = pd.DataFrame(endpoint, columns=["Estimated Endpoint Signal"], index=col_list)

        sac0_all: List[List[Union[float, str]]] = []
        current_baseline0 = baseline0

        for curve in range(n_curves):
            if not current_baseline0 or current_baseline0 > baseline[curve] or yfit0[curve] < 0.03:
                current_baseline0 = baseline[curve]

            unc = sigma_all[curve][0].tolist() if len(sigma_all[curve]) > 0 else []
            a_val: Any = float(sac_all[curve, 0])
            a_val = a_val if sum(unc) == 0 else ufloat(a_val, unc[0])

            if endpoint0 and len(unc) > 0:
                correction = (endpoint[curve] - current_baseline0) / (endpoint0 - current_baseline0)
                unc[0] = unc[0] * correction

            a_std = "" if isinstance(a_val, (float, np.floating)) else a_val.std_dev
            a_mean = a_val if isinstance(a_val, (float, np.floating)) else a_val.n
            sac0_all.append([a_mean, a_std])

        sac0_arr = np.array(sac0_all, dtype=object)
        sac0_df = pd.DataFrame(sac0_arr, columns=["ka", "Std Error"], index=col_list)

        st.dataframe(baseline_df.apply(lambda col: col.map(format_scientific)))
        st.dataframe(endpoint_df.apply(lambda col: col.map(format_scientific)))
        st.write(sac0_df.apply(lambda col: col.map(format_scientific)))

        st.success("SAC-RHO Results")

        if sac_norm == 0.0:
            arr = np.array(sac0_arr[:, 0], dtype=float)
            valid_values = arr[(arr != 0) & (~np.isnan(arr))]
            kac_norm = float(np.mean(valid_values)) if len(valid_values) > 0 else 0.0
            st.write(f"SAC value = {format_scientific(kac_norm)} (average ka value)")
        else:
            kac_norm = sac_norm
            st.write(f"SAC value = {format_scientific(kac_norm)} (user-defined)")

        _, _, _, xseed0_all_2, _, _, sigma_all_2 = SACRHO2(
            df,
            baseline,
            endpoint,
            DEFAULT_SEGMENTS,
            kac0,
            yfit0,
            n_curves,
            current_baseline0,
            kb,
            colors,
            backc,
            kac_norm,
            min_points,
        )

        rho_all2: List[List[Union[float, str]]] = []
        dfs = np.array(endpoint) - np.array(baseline)
        df_max = float(np.max(dfs))

        for curve in range(n_curves):
            unc2 = sigma_all_2[curve][0].tolist() if len(sigma_all_2[curve]) > 0 else []
            b2: Any = float(xseed0_all_2[curve, 0])
            correction2 = dfs[curve] / df_max if df_max != 0 else 1.0
            b2 = b2 * correction2

            if len(unc2) > 1:
                unc2[1] = unc2[1] * np.log(10) * xseed0_all_2[curve, 0] * correction2
                unc2[1] = unc2[1] * abs(1 - kb) / (1 + xseed0_all_2[curve, 0] ** 2)
                b2 = ufloat(b2, unc2[1])

            b_std2 = "" if isinstance(b2, (float, np.floating)) else b2.std_dev
            b_mean2 = b2 if isinstance(b2, (float, np.floating)) else b2.n
            rho_all2.append([b_mean2, b_std2])

        rho0_df2 = pd.DataFrame(rho_all2, columns=["RHO", "Std Error"], index=col_list)
        st.info("RHO values")
        st.dataframe(rho0_df2.apply(lambda col: col.map(format_scientific)))
