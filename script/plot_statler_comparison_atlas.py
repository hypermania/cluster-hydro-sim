#!/usr/bin/env python3
"""Build a figure-by-figure HydroSim comparison to Statler et al. (1987)."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from textwrap import fill

import fitz
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np


HISTORY_FILES = {
    "time": "history_time_trh.dat",
    "rho0": "history_rho0.dat",
    "vms2": "history_vms2_0.dat",
    "rc": "history_rc.dat",
    "rh": "history_rh.dat",
    "nb": "history_nb.dat",
    "ratio": "history_ratio.dat",
    "age": "history_mean_binary_age_trh.dat",
    "capture_rate": "history_capture_number_rate_per_year.dat",
    "threebody_rate": "history_threebody_number_rate_per_year.dat",
    "capture_energy": "history_capture_energy_rate.dat",
    "direct_s_bs": "history_direct_s_bs.dat",
    "direct_b_bs": "history_direct_b_bs.dat",
    "direct_b_bb": "history_direct_b_bb.dat",
    "mass": "history_mass.dat",
    "energy": "history_energy.dat",
}

SNAPSHOT_FILES = {
    "radius": "snapshot_radius.dat",
    "rho_s": "snapshot_rho_s.dat",
    "rho_b": "snapshot_rho_b.dat",
    "u_s": "snapshot_u_s.dat",
    "u_b": "snapshot_u_b.dat",
    "luminosity": "snapshot_luminosity.dat",
}


@dataclass(frozen=True)
class RunData:
    history: dict[str, np.ndarray]
    snapshot_time: np.ndarray
    snapshots: dict[str, np.ndarray]


def _read_float64(path: Path) -> np.ndarray:
    if not path.is_file() or path.stat().st_size % 8:
        raise ValueError(f"missing or malformed float64 input: {path}")
    result = np.fromfile(path, dtype="<f8")
    if not result.size or not np.all(np.isfinite(result)):
        raise ValueError(f"empty or non-finite input: {path}")
    return result


def load_run(directory: str | Path) -> RunData:
    directory = Path(directory).expanduser().resolve()
    history = {
        key: _read_float64(directory / filename)
        for key, filename in HISTORY_FILES.items()
    }
    lengths = {values.size for values in history.values()}
    if len(lengths) != 1 or np.any(np.diff(history["time"]) < 0.0):
        raise ValueError("history arrays must have one length and ordered time")
    zone_path = directory / "snapshot_zone_count.dat"
    if not zone_path.is_file() or zone_path.stat().st_size != 8:
        raise ValueError(f"missing or malformed zone count: {zone_path}")
    zones = int(np.fromfile(zone_path, dtype="<i8")[0])
    times = _read_float64(directory / "snapshot_time_trh.dat")
    snapshots = {}
    for key, filename in SNAPSHOT_FILES.items():
        flat = _read_float64(directory / filename)
        if flat.size != times.size * zones:
            raise ValueError(f"snapshot shape mismatch: {directory / filename}")
        snapshots[key] = flat.reshape(times.size, zones)
    if np.any(np.diff(times) < 0.0):
        raise ValueError("snapshot times are not ordered")
    return RunData(history, times, snapshots)


def projected_density(radius: np.ndarray, density: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Abel-project a positive spherical density profile using a z integral."""
    mask = (radius > 0.0) & (density > 0.0) & np.isfinite(radius * density)
    radius = radius[mask]
    density = density[mask]
    if radius.size < 4 or np.any(np.diff(radius) <= 0.0):
        raise ValueError("projection needs at least four ordered positive samples")
    projected_radius = radius[:-2]
    log_radius = np.log(radius)
    log_density = np.log(density)
    sigma = np.empty_like(projected_radius)
    for index, projected in enumerate(projected_radius):
        maximum_z = np.sqrt(max(radius[-1] ** 2 - projected**2, 0.0))
        positive_z = np.geomspace(max(projected * 1.0e-7, 1.0e-14), maximum_z, 900)
        z = np.concatenate(([0.0], positive_z))
        spherical_radius = np.sqrt(projected**2 + z**2)
        interpolated = np.exp(np.interp(np.log(spherical_radius), log_radius, log_density))
        sigma[index] = 2.0 * np.trapezoid(interpolated, z)
    return projected_radius, sigma


# Each entry is (PDF page, normalized left, top, right, bottom, quarter-turns CCW).
REFERENCE_CROPS: dict[int, list[tuple[int, float, float, float, float, int]]] = {
    1: [(11, 0.02, 0.00, 0.80, 1.00, 0)],
    2: [(12, 0.46, 0.00, 1.00, 0.49, 0)],
    3: [(12, 0.46, 0.46, 1.00, 1.00, 0)],
    4: [(13, 0.16, 0.055, 0.77, 0.445, 0),
        (14, 0.17, 0.095, 0.79, 0.47, 0),
        (13, 0.16, 0.475, 0.77, 0.845, 0)],
    5: [(15, 0.03, 0.04, 0.95, 0.95, 3)],
    6: [(16, 0.00, 0.00, 0.58, 0.62, 0)],
    7: [(17, 0.04, 0.05, 0.95, 0.94, 3)],
    8: [(18, 0.04, 0.05, 0.95, 0.94, 3)],
    9: [(19, 0.00, 0.00, 0.72, 0.57, 0)],
    10: [(19, 0.00, 0.48, 0.72, 1.00, 0)],
    11: [(20, 0.00, 0.00, 0.74, 1.00, 0)],
    12: [(21, 0.25, 0.08, 0.74, 0.75, 0)],
    13: [(22, 0.03, 0.04, 0.94, 0.94, 3)],
    14: [(23, 0.00, 0.00, 0.69, 1.00, 0)],
    15: [(24, 0.00, 0.00, 0.69, 1.00, 0)],
    16: [(25, 0.00, 0.00, 0.74, 1.00, 0)],
    17: [(26, 0.05, 0.09, 0.94, 0.36, 0), (26, 0.05, 0.38, 0.49, 0.65, 0)],
    18: [(26, 0.49, 0.43, 1.00, 1.00, 0)],
    19: [(27, 0.18, 0.04, 0.79, 0.86, 0)],
    20: [(28, 0.03, 0.04, 0.95, 0.95, 3)],
    21: [(29, 0.03, 0.04, 0.95, 0.95, 3)],
    22: [(30, 0.00, 0.00, 0.69, 1.00, 0)],
    23: [(31, 0.00, 0.00, 0.74, 0.58, 0)],
}


# Limits read from the *frames*, not just the outermost labelled ticks.
# Logs of dimensional quantities printed in the paper are represented with
# log-scaled positive quantities here. Each pair of panels uses identical limits.
@dataclass(frozen=True)
class PanelAxes:
    xlim: tuple[float, float]
    ylim: tuple[float, float]
    xscale: str = "log"
    yscale: str = "log"


GLOBAL_AXES = [PanelAxes((1, 1e4), y) for y in [(1e-1, 1e4), (1e-3, 10), (1e-3, 100)]]
COUNT_AXES = PanelAxes((1, 1e4), (1, 1e4))
RATIO_AXES = PanelAxes((1, 1e4), (1e-3, 100))
AGE_AXES = PanelAxes((1, 1e4), (1e-1, 100))
HEATING_AXES = [PanelAxes(x, (1e-9, 1e-1)) for x in [(20, 1e-2), (1e-2, 1e4)]]
DENSITY_AXES = PanelAxes((1e-7, 1e5), (1e-25, 1e10))
VELOCITY_AXES = PanelAxes((1e-7, 1e5), (1e-5, 10))
PANEL_AXES = {
    1: GLOBAL_AXES, 2: [COUNT_AXES], 3: [RATIO_AXES],
    4: [PanelAxes((20, 1e-2), (1e-12, 1)),
        PanelAxes((1e-2, 1e4), (1e-12, 1)),
        PanelAxes((1e-2, 1e4), (1e-3, 10))],
    5: HEATING_AXES, 6: [AGE_AXES], 7: [DENSITY_AXES] * 4,
    8: [VELOCITY_AXES, VELOCITY_AXES,
        PanelAxes((1e-7, 1e5), (1e-3, 10)), VELOCITY_AXES],
    9: [PanelAxes((-7, 5), (-2e-3, 0.65e-3), "linear", "linear")],
    10: [PanelAxes((1e-2, 1e4), (1e-2, 1))],
    11: GLOBAL_AXES, 12: [COUNT_AXES, RATIO_AXES], 13: HEATING_AXES,
    14: [DENSITY_AXES] * 2, 15: [VELOCITY_AXES] * 2,
    16: GLOBAL_AXES, 17: [COUNT_AXES, RATIO_AXES, AGE_AXES],
    18: [PanelAxes((1, 1e4), (1e-15, 1))], 19: HEATING_AXES,
    20: [PanelAxes(x, (-1.1, 1.1), yscale="linear")
         for x in [(20, 1e-2), (1e-2, 1e4)]],
    21: [DENSITY_AXES] * 2, 22: [VELOCITY_AXES] * 2, 23: [DENSITY_AXES],
}

# (PDF page, frame rectangle in PDF points, quarter-turns CCW).
# Calibrated against the 612 x 792 point ADS scan, excluding tick labels.
# These are raster registrations, not digitized or fitted reference curves.
OVERLAY_FRAMES = {
    11: [(20, (63, 73, 266, 236), 0), (20, (67, 300, 270, 464), 0),
         (20, (63, 531, 265, 695), 0)],
    12: [(21, (199, 95, 412, 268), 0), (21, (199, 375, 414, 552), 0)],
    13: [(22, (170, 403, 391, 678), 3), (22, (169, 66, 390, 342), 3)],
    14: [(23, (175, 94, 415, 290), 0), (23, (176, 363, 417, 559), 0)],
    15: [(24, (178, 95, 421, 292), 0), (24, (179, 370, 422, 567), 0)],
}

# Paper primary rate unit: initial N_star / initial half-mass relaxation time.
STATLER_RATE_PER_YEAR = 3e5 / 225e6


def apply_panel_axes(axes: list[plt.Axes], figure_number: int) -> None:
    for axis, spec in zip(axes, PANEL_AXES[figure_number], strict=True):
        axis.set_xscale(spec.xscale)
        axis.set_yscale(spec.yscale)
        axis.set_xlim(spec.xlim)
        axis.set_ylim(spec.ylim)


def overlay_reference(axis: plt.Axes, document: fitz.Document,
                      frame: tuple) -> None:
    page_number, rectangle, turns = frame
    page = document.load_page(page_number - 1)
    # Scale the calibration if the same scan has been uniformly resized.
    rect = fitz.Rect(*(v * (page.rect.width / 612 if i % 2 == 0 else page.rect.height / 792)
                       for i, v in enumerate(rectangle)))
    pixmap = page.get_pixmap(matrix=fitz.Matrix(3, 3), clip=rect, alpha=False)
    image = np.frombuffer(pixmap.samples, dtype=np.uint8).reshape(
        pixmap.height, pixmap.width, pixmap.n)[..., :3]
    image = np.rot90(image, turns)
    # Axes coordinates preserve linear pixel spacing in log10(x), log10(y),
    # including the reversed precollapse x axis. A data-space imshow extent
    # would incorrectly warp the raster on logarithmic axes.
    axis.imshow(image, extent=(0, 1, 0, 1), transform=axis.transAxes,
                origin="upper", aspect="auto", alpha=0.55, zorder=0)
    # Fix the physical axes rectangle, not the data-unit aspect (which would
    # be incorrect for unequal logarithmic ranges). Account for scan rotation.
    axis.set_box_aspect(image.shape[0] / image.shape[1])


def reference_images(document: fitz.Document, figure_number: int) -> list[np.ndarray]:
    images = []
    for page_number, left, top, right, bottom, turns in REFERENCE_CROPS[figure_number]:
        page = document.load_page(page_number - 1)
        pixmap = page.get_pixmap(matrix=fitz.Matrix(1.6, 1.6), alpha=False)
        image = np.frombuffer(pixmap.samples, dtype=np.uint8).reshape(
            pixmap.height, pixmap.width, pixmap.n
        )[..., :3]
        x0, x1 = int(left * image.shape[1]), int(right * image.shape[1])
        y0, y1 = int(top * image.shape[0]), int(bottom * image.shape[0])
        image = image[y0:y1, x0:x1]
        if turns:
            image = np.rot90(image, turns)
        images.append(image)
    return images


def configure_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Latin Modern Roman", "Computer Modern Roman"],
            "mathtext.fontset": "cm",
            "font.size": 8.5,
            "axes.labelsize": 9,
            "axes.titlesize": 9,
            "legend.fontsize": 7,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "lines.linewidth": 1.35,
        }
    )


def _history_mask(data: RunData) -> np.ndarray:
    return data.history["time"] > 0.0


def _base_axis(axis: plt.Axes, xlabel: str, ylabel: str) -> None:
    axis.set_xlabel(xlabel)
    axis.set_ylabel(ylabel)
    axis.grid(which="major", color="0.88", linewidth=0.45)


def plot_evolution(axes: list[plt.Axes], data: RunData) -> None:
    h = data.history
    mask = _history_mask(data)
    axes[0].loglog(h["time"][mask], h["rho0"][mask], color="#1f77b4")
    _base_axis(axes[0], r"$t/t_{rh}$", r"$\rho_0/(M/r_0^3)$")
    axes[1].loglog(h["time"][mask], h["vms2"][mask], color="#d62728")
    _base_axis(axes[1], r"$t/t_{rh}$", r"$v_{m,s}^2/(GM/r_0)$")
    axes[2].loglog(h["time"][mask], h["rc"][mask], label=r"$r_c$", color="#9467bd")
    axes[2].loglog(h["time"][mask], h["rh"][mask], label=r"$r_h$", color="#2ca02c")
    axes[2].legend(frameon=False)
    _base_axis(axes[2], r"$t/t_{rh}$", r"$r/r_0$")
    for axis in axes:
        axis.set_xlim(1.0, 1.0e4)


def plot_binary_number(axis: plt.Axes, data: RunData) -> None:
    mask = _history_mask(data)
    # No primordial binaries: zero counts cannot be displayed on a log axis.
    mask &= data.history["nb"] > 0.0
    axis.loglog(data.history["time"][mask], data.history["nb"][mask], color="#1f77b4")
    axis.set_xlim(1.0, 1.0e4)
    _base_axis(axis, r"$t/t_{rh}$", r"$N_b$")


def plot_ratio(axis: plt.Axes, data: RunData) -> None:
    mask = _history_mask(data)
    axis.loglog(data.history["time"][mask], data.history["ratio"][mask], color="#d62728")
    axis.set_xlim(1.0, 1.0e4)
    _base_axis(axis, r"$t/t_{rh}$", r"$\rho_{bc}/\rho_{sc}$")


def _peak_time(data: RunData) -> float:
    return float(data.history["time"][np.argmax(data.history["rho0"])])


def _phase(data: RunData, before: bool) -> tuple[np.ndarray, np.ndarray]:
    time = data.history["time"]
    separation = (_peak_time(data) - time) if before else (time - _peak_time(data))
    mask = (separation > 1.0e-3) & (time > 0.0)
    return mask, separation


def plot_interaction_rates(axes: list[plt.Axes], data: RunData) -> None:
    h = data.history
    for axis, before in zip(axes[:2], [True, False], strict=True):
        mask, separation = _phase(data, before)
        axis.loglog(separation[mask], h["capture_rate"][mask] / STATLER_RATE_PER_YEAR, label="tidal capture")
        axis.loglog(separation[mask], h["threebody_rate"][mask] / STATLER_RATE_PER_YEAR, linestyle="--", label="3-body")
        axis.invert_xaxis() if before else None
        axis.legend(frameon=False)
        _base_axis(axis, r"$|t-t_c|/t_{rh}$", r"rate ($N_*/t_{rh}$)")
    mask, separation = _phase(data, False)
    if np.count_nonzero(mask) > 3:
        time_year = h["time"] * 225.0e6
        net_rate = np.gradient(h["nb"], time_year)
        ratio = net_rate / np.maximum(h["capture_rate"], np.finfo(float).tiny)
        good = mask & (ratio > 0.0)
        axes[2].loglog(separation[good], ratio[good], color="#9467bd")
    else:
        axes[2].text(0.5, 0.5, "No postcollapse state", ha="center", va="center", transform=axes[2].transAxes)
    _base_axis(axes[2], r"$(t-t_c)/t_{rh}$", r"$\dot N_b/\dot N_{b,f}$")


def _heating_series(data: RunData) -> list[tuple[str, np.ndarray, str]]:
    h = data.history
    return [
        ("$-$ tidal-capture cooling", -h["capture_energy"], "-"),
        ("single heating, binary--single", h["direct_s_bs"], "--"),
        ("binary heating, binary--single", h["direct_b_bs"], "-."),
        ("binary--binary heating", h["direct_b_bb"], ":"),
    ]


def plot_heating(axes: list[plt.Axes], data: RunData) -> None:
    for axis, before in zip(axes, [True, False], strict=True):
        mask, separation = _phase(data, before)
        for label, values, linestyle in _heating_series(data):
            good = mask & (values > 0.0)
            if np.any(good):
                axis.loglog(separation[good], values[good], label=label, linestyle=linestyle)
        axis.invert_xaxis() if before else None
        if axis.lines:
            axis.legend(frameon=False)
        else:
            axis.text(0.5, 0.5, "No postcollapse state", ha="center", va="center", transform=axis.transAxes)
        _base_axis(axis, r"$|t-t_c|/t_{rh}$", r"$|\dot E|/[GM^2/(r_0t_{rh})]$")


def plot_age(axis: plt.Axes, data: RunData) -> None:
    mask = _history_mask(data) & (data.history["age"] > 0.0)
    axis.loglog(data.history["time"][mask], data.history["age"][mask], color="#9467bd")
    axis.set_xlim(1.0, 1.0e4)
    _base_axis(axis, r"$t/t_{rh}$", r"estimated $\langle\tau_b\rangle/t_{rh}$")


def _profile_indices(data: RunData, postcollapse: bool, count: int = 4) -> np.ndarray:
    peak = _peak_time(data)
    if postcollapse:
        available = np.flatnonzero(data.snapshot_time > peak * 1.02)
    else:
        available = np.flatnonzero((data.snapshot_time > 0.0) & (data.snapshot_time < peak))
    if available.size <= count:
        return available
    return available[np.linspace(0, available.size - 1, count).round().astype(int)]


def _plot_profiles(axis: plt.Axes, data: RunData, quantity: str,
                   component: str, postcollapse: bool) -> None:
    ylabel = r"$\rho/(M/r_0^3)$" if quantity.startswith("density") else r"$v_m^2/(GM/r_0)$"
    indices = _profile_indices(data, postcollapse)
    if not indices.size:
        axis.text(0.5, 0.5, "No postcollapse state", ha="center", va="center", transform=axis.transAxes)
        _base_axis(axis, r"$r/r_0$", ylabel)
        return
    colors = plt.colormaps["turbo"](np.linspace(0.05, 0.92, indices.size))
    for color, index in zip(colors, indices, strict=True):
        radius = data.snapshots["radius"][index]
        if quantity.startswith("density"):
            values = (3.0 / (4.0 * np.pi)) * data.snapshots[f"rho_{component}"][index]
        else:
            values = 6.0 * data.snapshots[f"u_{component}"][index]
        axis.loglog(radius, values, color=color, label=rf"$t/t_{{rh}}={data.snapshot_time[index]:.2g}$")
    axis.legend(frameon=False, ncol=2)
    _base_axis(axis, r"$r/r_0$", ylabel)
    axis.set_xlim(1.0e-6, 1.0e4)


def plot_four_profiles(axes: list[plt.Axes], data: RunData, quantity: str) -> None:
    combinations = [("s", False), ("s", True), ("b", False), ("b", True)]
    for axis, (component, postcollapse) in zip(axes, combinations, strict=True):
        _plot_profiles(axis, data, quantity, component, postcollapse)
        axis.set_title(("postcollapse " if postcollapse else "precollapse ") + ("singles" if component == "s" else "binaries"))


def plot_two_profiles(axes: list[plt.Axes], data: RunData, quantity: str) -> None:
    for axis, component in zip(axes, ["s", "b"], strict=True):
        _plot_profiles(axis, data, quantity, component, True)
        axis.set_title("postcollapse " + ("singles" if component == "s" else "binaries"))


def plot_luminosity(axis: plt.Axes, data: RunData) -> None:
    indices = _profile_indices(data, True, 3)
    multipliers = [1.0, 10.0, 100.0]
    selected_multipliers = multipliers[-len(indices) :] if len(indices) else []
    for index, multiplier in zip(indices[-3:], selected_multipliers, strict=True):
        radius = data.snapshots["radius"][index]
        luminosity = multiplier * data.snapshots["luminosity"][index]
        axis.plot(np.log10(radius), luminosity, label=rf"$t/t_{{rh}}={data.snapshot_time[index]:.2g}$, $\times{multiplier:g}$")
    if indices.size:
        axis.legend(frameon=False)
    else:
        axis.text(0.5, 0.5, "No postcollapse state", ha="center", va="center", transform=axis.transAxes)
    _base_axis(axis, r"$\log_{10}(r/r_0)$", r"$L_c/[GM^2/(r_0t_{rh})]$")


def plot_mass_slope(axis: plt.Axes, data: RunData) -> None:
    mask, separation = _phase(data, False)
    h = data.history
    positive_time = h["time"] > 0.0
    slope = np.full_like(h["time"], np.nan)
    slope[positive_time] = -np.gradient(np.log(h["mass"][positive_time]), np.log(h["time"][positive_time]))
    good = mask & (slope > 0.0)
    axis.loglog(separation[good], slope[good], color="#2ca02c")
    _base_axis(axis, r"$(t-t_c)/t_{rh}$", r"$-d\ln M/d\ln t$")


def plot_fractional_heating(axes: list[plt.Axes], data: RunData) -> None:
    h = data.history
    components = {
        "formation": h["capture_energy"],
        "direct heating": h["direct_s_bs"] + h["direct_b_bs"] + h["direct_b_bb"],
    }
    total = sum(components.values())
    absolute_sum = sum(np.abs(values) for values in components.values())
    for axis, before in zip(axes, [True, False], strict=True):
        mask, separation = _phase(data, before)
        good_total = np.abs(total) > 0.02 * absolute_sum
        for label, values in components.items():
            good = mask & good_total
            axis.semilogx(separation[good], values[good] / total[good], label=label)
        axis.invert_xaxis() if before else None
        axis.axhline(0.0, color="0.5", linewidth=0.6)
        axis.set_ylim(-1.1, 1.1)
        axis.legend(frameon=False)
        _base_axis(axis, r"$|t-t_c|/t_{rh}$", "fractional resolved heating")


def plot_projected(axis: plt.Axes, data: RunData) -> None:
    target = 13000.0 / 225.0
    index = int(np.argmin(np.abs(data.snapshot_time - target)))
    radius = data.snapshots["radius"][index]
    projected = []
    for component, label, linestyle in [("s", "singles", "-"), ("b", "binaries", "--")]:
        projected_radius, sigma = projected_density(radius, data.snapshots[f"rho_{component}"][index])
        sigma *= 3.0 / (4.0 * np.pi)
        projected.append(sigma)
        axis.loglog(projected_radius, sigma, label=label, linestyle=linestyle)
    axis.loglog(projected_radius, projected[0] + projected[1], label="total", linestyle=":")
    axis.legend(frameon=False)
    axis.set_xlim(1.0e-6, 1.0e4)
    _base_axis(axis, r"$R/r_0$", r"$\Sigma/(M/r_0^2)$")
    axis.set_title(rf"closest state: $t/t_{{rh}}={data.snapshot_time[index]:.3g}$")


FIGURE_TITLES = {
    1: "Ejection-only global evolution", 2: "Ejection-only binary count",
    3: "Ejection-only central binary ratio", 4: "Ejection-only interaction rates",
    5: "Ejection-only heating rates", 6: "Ejection-only mean binary age",
    7: "Ejection-only density profiles", 8: "Ejection-only velocity profiles",
    9: "Ejection-only conductive luminosity", 10: "Ejection-only mass-loss exponent",
    11: "Direct-heating-only global evolution", 12: "Direct-heating-only binary population",
    13: "Direct-heating-only heating rates", 14: "Direct-heating-only density profiles",
    15: "Direct-heating-only velocity profiles", 16: "Full-model global evolution",
    17: "Full-model binary population", 18: "Full-model single-star source rates",
    19: "Full-model heating rates", 20: "Full-model fractional heating",
    21: "Full-model density profiles", 22: "Full-model velocity profiles",
    23: "Full-model projected density at 13 Gyr",
}

EJECTION_NOTE = "Not apples-to-apples: HydroSim has no encounter-driven ejection. The blue result is a tidal-capture/cooling control with direct heating also disabled."
DIRECT_NOTE = "Closest available matched experiment: identical Plummer scaling and Press--Teukolsky capture, no ejection; HydroSim uses a local conducting-fluid heating closure rather than the paper's orbit-averaged encounter kernels."
FULL_NOTE = "Missing-physics baseline: the source figure includes ejection plus direct heating, while HydroSim includes direct heating only. Initial model, physical scaling, capture law, observable, and axes are matched."


def _make_right_axes(figure: plt.Figure, right_spec, count: int) -> list[plt.Axes]:
    if count == 1:
        return [figure.add_subplot(right_spec)]
    if count == 2:
        grid = right_spec.subgridspec(2, 1, hspace=0.38)
    elif count == 3:
        grid = right_spec.subgridspec(3, 1, hspace=0.48)
    elif count == 4:
        grid = right_spec.subgridspec(2, 2, hspace=0.38, wspace=0.35)
    else:
        raise ValueError("unsupported panel count")
    return [figure.add_subplot(grid[index]) for index in range(count)]


def _source_panel(figure: plt.Figure, source_spec, images: list[np.ndarray]) -> None:
    grid = source_spec.subgridspec(len(images), 1, hspace=0.04)
    for index, image in enumerate(images):
        axis = figure.add_subplot(grid[index])
        axis.imshow(image)
        axis.axis("off")
        if index == 0:
            axis.set_title("Statler, Ostriker & Cohn (1987) scan", pad=4)


def make_page(document: fitz.Document, figure_number: int,
              control: RunData, direct: RunData, *, overlay: bool = False) -> plt.Figure:
    if overlay and figure_number not in OVERLAY_FRAMES:
        raise ValueError("reference overlays are calibrated for Figures 11--15 only")
    figure = plt.figure(figsize=(8.27, 11.69) if overlay else (11.69, 8.27))
    outer = figure.add_gridspec(1, 1 if overlay else 2,
                                left=0.09 if overlay else 0.035, right=0.96,
                                top=0.88, bottom=0.15 if overlay else 0.105,
                                wspace=0.22)
    if not overlay:
        _source_panel(figure, outer[0], reference_images(document, figure_number))
    right = outer[-1]
    data = control if figure_number <= 10 else direct
    note = EJECTION_NOTE if figure_number <= 10 else (DIRECT_NOTE if figure_number <= 15 else FULL_NOTE)

    if figure_number in {1, 11, 16}:
        axes = _make_right_axes(figure, right, 3); plot_evolution(axes, data)
    elif figure_number in {2}:
        axes = _make_right_axes(figure, right, 1); plot_binary_number(axes[0], data)
    elif figure_number in {3}:
        axes = _make_right_axes(figure, right, 1); plot_ratio(axes[0], data)
    elif figure_number == 4:
        axes = _make_right_axes(figure, right, 3); plot_interaction_rates(axes, data)
    elif figure_number in {5, 13, 19}:
        axes = _make_right_axes(figure, right, 2); plot_heating(axes, data)
    elif figure_number == 6:
        axes = _make_right_axes(figure, right, 1); plot_age(axes[0], data)
    elif figure_number == 7:
        axes = _make_right_axes(figure, right, 4); plot_four_profiles(axes, data, "density")
    elif figure_number == 8:
        axes = _make_right_axes(figure, right, 4); plot_four_profiles(axes, data, "velocity")
    elif figure_number == 9:
        axes = _make_right_axes(figure, right, 1); plot_luminosity(axes[0], data)
    elif figure_number == 10:
        axes = _make_right_axes(figure, right, 1); plot_mass_slope(axes[0], data)
    elif figure_number in {12}:
        axes = _make_right_axes(figure, right, 2)
        plot_binary_number(axes[0], data); plot_ratio(axes[1], data)
    elif figure_number in {14, 21}:
        axes = _make_right_axes(figure, right, 2); plot_two_profiles(axes, data, "density")
    elif figure_number in {15, 22}:
        axes = _make_right_axes(figure, right, 2); plot_two_profiles(axes, data, "velocity")
    elif figure_number == 17:
        axes = _make_right_axes(figure, right, 3)
        plot_binary_number(axes[0], data); plot_ratio(axes[1], data); plot_age(axes[2], data)
    elif figure_number == 18:
        axes = _make_right_axes(figure, right, 1)
        mask = _history_mask(data)
        axes[0].loglog(data.history["time"][mask], 2.0 * data.history["capture_rate"][mask] / STATLER_RATE_PER_YEAR, label="singles consumed by capture")
        axes[0].loglog(data.history["time"][mask], data.history["threebody_rate"][mask] / STATLER_RATE_PER_YEAR, linestyle="--", label="3-body binary events")
        axes[0].legend(frameon=False); _base_axis(axes[0], r"$t/t_{rh}$", r"rate ($N_*/t_{rh}$)")
    elif figure_number == 20:
        axes = _make_right_axes(figure, right, 2); plot_fractional_heating(axes, data)
    elif figure_number == 23:
        axes = _make_right_axes(figure, right, 1); plot_projected(axes[0], data)
    else:
        raise AssertionError(f"unhandled figure {figure_number}")

    apply_panel_axes(axes, figure_number)
    if overlay:
        for axis, frame in zip(axes, OVERLAY_FRAMES[figure_number], strict=True):
            overlay_reference(axis, document, frame)
        # imshow may change limits even with an axes transform: restore the
        # same explicit contract used for the side-by-side plots.
        apply_panel_axes(axes, figure_number)
        figure.text(0.09, 0.925, "Colored curves: HydroSim    |    Gray scan: Statler et al. (1987)", fontsize=10)
        figure.text(0.09, 0.085,
                    "Scan-registered overlay, not digitized data; reference annotations remain gray.\n"
                    "Profiles show available saved HydroSim epochs (legend), not matched reference epochs.",
                    fontsize=8)

    figure.suptitle(f"Statler et al. Figure {figure_number}: {FIGURE_TITLES[figure_number]}", fontsize=12 if overlay else 13, y=0.965)
    figure.text(0.035, 0.035, fill(note, 112 if overlay else 165),
                ha="left", va="bottom", fontsize=8.2)
    figure.text(0.98, 0.018, "HydroSim comparison atlas", ha="right", va="bottom", fontsize=7, color="0.35")
    return figure


def build_atlas(reference_pdf: Path, control: RunData, direct: RunData,
                output: Path) -> None:
    configure_style()
    output.parent.mkdir(parents=True, exist_ok=True)
    document = fitz.open(reference_pdf)
    try:
        with PdfPages(output) as pdf:
            for figure_number in range(1, 24):
                figure = make_page(document, figure_number, control, direct)
                pdf.savefig(figure)
                plt.close(figure)
        overlay_output = output.with_name(output.stem + "_overlays.pdf")
        with PdfPages(overlay_output) as pdf:
            for figure_number in OVERLAY_FRAMES:
                figure = make_page(document, figure_number, control, direct, overlay=True)
                pdf.savefig(figure)
                plt.close(figure)
    finally:
        document.close()


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference_pdf", type=Path)
    parser.add_argument("control_directory", type=Path)
    parser.add_argument("direct_directory", type=Path)
    parser.add_argument("--output", type=Path, default=Path("output/statler_atlas/plots/statler_hydrosim_comparison_atlas.pdf"))
    return parser


def main(argv: list[str] | None = None) -> int:
    arguments = build_argument_parser().parse_args(argv)
    control = load_run(arguments.control_directory)
    direct = load_run(arguments.direct_directory)
    output = arguments.output.expanduser().resolve()
    build_atlas(arguments.reference_pdf.expanduser().resolve(), control, direct, output)
    print(f"Wrote {output}")
    print(f"Wrote {output.with_name(output.stem + '_overlays.pdf')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
