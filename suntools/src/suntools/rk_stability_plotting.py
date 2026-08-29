#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2025-2026, Lawrence Livermore National Security,
# University of Maryland Baltimore County, and the SUNDIALS contributors.
# Copyright (c) 2013-2025, Lawrence Livermore National Security
# and Southern Methodist University.
# Copyright (c) 2002-2013, Lawrence Livermore National Security.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------------------

"""Matplotlib-facing plotting utilities for Runge-Kutta stability regions."""

from __future__ import annotations

from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from .rk_butcher_table import ButcherTable
from .rk_stability import stability_magnitude

# ---------------------------------------------------------------------------
# Output canvas
# ---------------------------------------------------------------------------

# Each figure uses a fixed canvas with the axes box pinned to fixed margins, so all plots have
# the same shape with the axes in the same place. Margins are in inches.
_FIG_SIZE = (5.5, 6.5)
_FIG_DPI = 150

_MARGIN_LEFT = 0.78
_MARGIN_RIGHT = 0.22
_MARGIN_BOTTOM = 0.55
_MARGIN_TOP = 0.45

_COLORBAR_STRIP = 1.30
_COLORBAR_WIDTH = 0.20
_COLORBAR_PAD = 0.12


# ---------------------------------------------------------------------------
# Options
# ---------------------------------------------------------------------------


@dataclass
class PlotOptions:
    """Everything that varies between stability-region plots."""

    bounds: tuple[float, float, float, float] | None = None
    grid_points: int = 500
    # Output canvas, in inches and dots per inch.
    figure_size: tuple[float, float] = _FIG_SIZE
    dpi: int = _FIG_DPI
    # Zeros/poles are informative annotations only; they do not change the shaded region.
    show_zeros: bool = False
    show_poles: bool = False
    show_embedding: bool = True
    shade: bool = False
    shade_embedding: bool = False
    suppress_tiny_islands: bool = True
    font_size: float = 14

    def __post_init__(self):
        if self.shade and self.shade_embedding:
            raise ValueError("shade and shade_embedding are mutually exclusive")
        if self.bounds is not None:
            self.bounds = tuple(float(v) for v in self.bounds)
        self.figure_size = tuple(float(v) for v in self.figure_size)
        # Fail here rather than with a negative axes later.
        _axes_geometry(self.figure_size, colorbar=True)


# ---------------------------------------------------------------------------
# Styling
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class RegionStyle:
    """Color and line properties for one method's region."""

    fill_color: str
    fill_alpha: float
    boundary_color: str
    boundary_lw: float
    boundary_ls: str
    shade_label: str
    legend_fill_alpha: float | None = None
    # Shade mode only: the boundary sits on the colormap rather than a flat fill, so the main
    # region overrides its to stay readable. Unset means keep boundary_color.
    shade_boundary_color: str | None = None

    @property
    def shade_color(self) -> str:
        return self.shade_boundary_color or self.boundary_color


MAIN_REGION_STYLE = RegionStyle(
    fill_color="#cfe3ff",
    fill_alpha=0.7,
    boundary_color="#1f5fb0",
    boundary_lw=1.8,
    boundary_ls="-",
    shade_label=r"$|R(z)|$ within the main stable region",
    shade_boundary_color="k",
)

EMBEDDED_REGION_STYLE = RegionStyle(
    fill_color="#ffcccc",
    fill_alpha=0.45,
    boundary_color="#cc0000",
    boundary_lw=1.6,
    boundary_ls="--",
    shade_label=r"$|\tilde{R}(z)|$ within the embedded stable region",
    legend_fill_alpha=0.7,
)

_SHADE_LEVELS = np.linspace(0.0, 1.0, 21)

# Auto-framing: scan a coarse grid wide enough to hold every bounded feature, then crop.
_BOUNDS_SCAN_POINTS = 420
_BOUNDS_MIN_SPAN = 10.0
_BOUNDS_STAGE_PAD = 4.0
_BOUNDS_GROWTH = 2.0
_BOUNDS_MAX_SPAN = 80.0
_BOUNDS_PAD_FRAC = 0.20
_BOUNDS_PAD_MIN = 0.5
_BOUNDS_FALLBACK_MIN = 6.0
_BOUNDS_FALLBACK_STAGE_PAD = 2.0

# Any value above 1 is unstable; used to paint over suppressed tiny islands.
_UNSTABLE_FILL = 2.0
_ISLAND_MIN_RADIUS = 0.06


# ---------------------------------------------------------------------------
# Region rendering
# ---------------------------------------------------------------------------


def _region_magnitude(phi, Z, X, Y, suppress_tiny_islands: bool):
    """|phi| on the plot grid, optionally with tiny stable islands masked out."""
    R = stability_magnitude(phi, Z)
    if suppress_tiny_islands:
        R = _suppress_tiny_islands(R, X, Y, phi.p.r)
    return R


def _render_region(
    ax,
    fig,
    X,
    Y,
    R,
    *,
    style: RegionStyle,
    mode: str,
    label: str | None,
    cax=None,
    font_size: float,
):
    """Draw one method's region and return its legend artist (or None).

    *mode* is "shade" (colormap under the boundary, with a colorbar), "fill" (flat tint under the
    boundary), or "boundary" (the curve only).
    """
    # Default to method color for the boundary color.
    color = style.boundary_color

    if mode == "shade":
        # Gradient fill for |phi| <= 1.
        color = style.shade_color
        contours = ax.contourf(X, Y, R, levels=_SHADE_LEVELS, cmap="viridis")
        # cax is the pre-placed colorbar axes; without one the bar is carved out of ax.
        bar = fig.colorbar(contours, cax=cax) if cax is not None else fig.colorbar(contours, ax=ax)
        bar.set_label(style.shade_label, fontsize=font_size)
        bar.ax.tick_params(labelsize=font_size)
    elif mode == "fill":
        # Solid fill for |phi| <= 1.
        ax.contourf(X, Y, R, levels=[0.0, 1.0], colors=[style.fill_color], alpha=style.fill_alpha)

    # |phi| = 1 boundary
    ax.contour(
        X,
        Y,
        R,
        levels=[1.0],
        colors=[color],
        linewidths=style.boundary_lw,
        linestyles=style.boundary_ls,
    )

    if label is None:
        return None
    if mode == "fill":
        alpha = {} if style.legend_fill_alpha is None else {"alpha": style.legend_fill_alpha}
        return Patch(
            facecolor=style.fill_color, edgecolor=style.boundary_color, label=label, **alpha
        )
    return Line2D([0], [0], color=color, lw=style.boundary_lw, ls=style.boundary_ls, label=label)


def _suppress_tiny_islands(R, X, Y, zeros, min_radius: float = _ISLAND_MIN_RADIUS):
    """Mask out tiny stable islands around zeros of phi."""
    if zeros is None or len(zeros) == 0:
        return R
    ny, nx = R.shape
    # Grid spacing and origin, for turning an (x, y) position into row/column indices.
    dx = X[0, 1] - X[0, 0]
    dy = Y[1, 0] - Y[0, 0]
    x0, y0 = X[0, 0], Y[0, 0]
    # Max cell count for an island considered "tiny": area of a disk of min_radius in cells. Floor
    # of 8 keeps very coarse grids workable.
    cap = max(8, int(np.pi * min_radius**2 / abs(dx * dy)))
    stable = R <= 1.0
    R = R.copy()
    for z in zeros:
        j = int(round((z.real - x0) / dx))
        i = int(round((z.imag - y0) / dy))
        # Skip zeros off-grid or not in a stable cell.
        if not (0 <= i < ny and 0 <= j < nx) or not stable[i, j]:
            continue
        # Flood fill outward (up/down/left/right) from this zero through stable cells
        island = {(i, j)}  # every cell found so far
        stack = [(i, j)]  # found cells whose neighbours have not been checked yet
        overflow = False
        while stack:
            ci, cj = stack.pop()
            # Exit early if the region is not tiny.
            if len(island) > cap:
                overflow = True
                break
            for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                ni, nj = ci + di, cj + dj
                # Record a cell when it's queued, so the same cell is never queued twice.
                if 0 <= ni < ny and 0 <= nj < nx and (ni, nj) not in island and stable[ni, nj]:
                    island.add((ni, nj))
                    stack.append((ni, nj))
        # Overwrite the island as unstable so it is not drawn.
        if not overflow:
            for ci, cj in island:
                R[ci, cj] = _UNSTABLE_FILL
    return R


# ---------------------------------------------------------------------------
# Framing and figure setup
# ---------------------------------------------------------------------------


def _auto_bounds(table: ButcherTable, include_embedding: bool = True):
    """Choose a plotting box that frames the interesting bounded features."""
    methods = [table.stability_function]
    emb = table.embedded_stability_function if include_embedding else None
    if emb is not None:
        methods.append(emb)

    def touches_edge(mask) -> bool:
        """True if any flagged cell is on the outer row or column of the grid (not fully closed)."""
        return bool(mask[0, :].any() or mask[-1, :].any() or mask[:, 0].any() or mask[:, -1].any())

    def feature_mask(phi, Z):
        """Return whichever side of |phi| = 1 closes inside the grid, or None if neither does."""
        R = stability_magnitude(phi, Z)
        for mask in (R <= 1.0, R > 1.0):
            if mask.any() and not touches_edge(mask):
                return mask
        return None

    # Start with a window sized off the stage count, then grow it until any features fit.
    span = min(max(_BOUNDS_MIN_SPAN, table.stages + _BOUNDS_STAGE_PAD), _BOUNDS_MAX_SPAN)
    while True:
        axis = np.linspace(-span, span, _BOUNDS_SCAN_POINTS)
        X, Y = np.meshgrid(axis, axis)
        masks = [m for m in (feature_mask(phi, X + 1j * Y) for phi in methods) if m is not None]
        # Stop once every method has a closed shape or the window reaches the max size.
        if len(masks) == len(methods) or span >= _BOUNDS_MAX_SPAN:
            break
        span = min(_BOUNDS_GROWTH * span, _BOUNDS_MAX_SPAN)

    # Union of the method and embedding feature regions
    features = np.logical_or.reduce(masks) if masks else None

    if features is not None and features.any():
        # Tight box around the flagged cells, including the origin
        xv, yv = X[features], Y[features]
        xmin, xmax = min(xv.min(), 0.0), max(xv.max(), 0.0)
        ymin, ymax = min(yv.min(), 0.0), max(yv.max(), 0.0)
        # Padding proportional to the box, with a floor for very flat shapes.
        padx = max(_BOUNDS_PAD_MIN, _BOUNDS_PAD_FRAC * (xmax - xmin))
        pady = max(_BOUNDS_PAD_MIN, _BOUNDS_PAD_FRAC * (ymax - ymin))
        return (xmin - padx, xmax + padx, ymin - pady, ymax + pady)

    # Nothing closed inside the largest window, fall back to a plain square.
    w = max(_BOUNDS_FALLBACK_MIN, table.stages + _BOUNDS_FALLBACK_STAGE_PAD)
    return (-w, w, -w, w)


def _axes_geometry(figure_size, colorbar: bool):
    """Axes rect in figure fractions, and the data aspect ratio that fills it exactly."""
    fig_w, fig_h = figure_size
    # The colorbar needs extra width on the right.
    right = _MARGIN_RIGHT + (_COLORBAR_STRIP if colorbar else 0.0)
    axes_w = fig_w - _MARGIN_LEFT - right
    axes_h = fig_h - _MARGIN_BOTTOM - _MARGIN_TOP
    if axes_w <= 0.0 or axes_h <= 0.0:
        # Margins exceed the figure; a rect here would be zero or negative.
        raise ValueError(
            f"figure size {tuple(figure_size)} leaves no room for the axes; "
            f"it must exceed {_MARGIN_LEFT + right:.2f} x "
            f"{_MARGIN_BOTTOM + _MARGIN_TOP:.2f} inches"
        )
    # add_axes wants (left, bottom, width, height) as fractions of the figure.
    rect = (_MARGIN_LEFT / fig_w, _MARGIN_BOTTOM / fig_h, axes_w / fig_w, axes_h / fig_h)
    return rect, axes_h / axes_w


def _colorbar_rect(figure_size):
    """Rect for the colorbar axes inside the strip reserved on the right."""
    fig_w, fig_h = figure_size
    # Step in from the right edge past the outer margin and the reserved strip, then back out by
    # the gap that separates the bar from the plot.
    left = fig_w - _MARGIN_RIGHT - _COLORBAR_STRIP + _COLORBAR_PAD
    height = fig_h - _MARGIN_BOTTOM - _MARGIN_TOP  # matches the plot axes exactly
    return (left / fig_w, _MARGIN_BOTTOM / fig_h, _COLORBAR_WIDTH / fig_w, height / fig_h)


def _fit_bounds_to_aspect(bounds, aspect: float):
    """Grow *bounds* about its centre until its height/width equals *aspect*."""
    xmin, xmax, ymin, ymax = bounds
    width, height = xmax - xmin, ymax - ymin
    if height < aspect * width:
        # Too short: pad the y range, split evenly so the centre doesn't move.
        grow = 0.5 * (aspect * width - height)
        return (xmin, xmax, ymin - grow, ymax + grow)
    # Otherwise too narrow: pad the x range instead. Only ever enlarging means the framing from
    # _auto_bounds is preserved and no feature gets cropped.
    grow = 0.5 * (height / aspect - width)
    return (xmin - grow, xmax + grow, ymin, ymax)


# ---------------------------------------------------------------------------
# Axis annotations
# ---------------------------------------------------------------------------


def _draw_zeros_poles(ax, phi, options: PlotOptions):
    """Optionally mark the zeros and poles of phi; returns their legend handles."""
    handles = []
    p, q = phi.p, phi.q
    if options.show_zeros and p.r.size:
        handles.append(
            ax.plot(p.r.real, p.r.imag, "kx", ms=7, mew=1.5, label=r"zeros of $R(z)$")[0]
        )
    if options.show_poles and q.r.size:
        handles.append(
            ax.plot(
                q.r.real,
                q.r.imag,
                "o",
                mfc="none",
                mec="#b03030",
                ms=8,
                mew=1.5,
                label=r"poles of $R(z)$",
            )[0]
        )
    return handles


def _decorate_axes(ax, table: ButcherTable, bounds, legend_handles, font_size: float):
    """Axis limits, labels, title, legend, and grid."""
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(bounds[0], bounds[1])
    ax.set_ylim(bounds[2], bounds[3])
    ax.set_xlabel(r"Re$(z)$", fontsize=font_size)
    ax.set_ylabel(r"Im$(z)$", fontsize=font_size)
    # Always use two lines as long names overrun one line and a constant title height keeps the top
    # margin and axes box fixed.
    ax.set_title(f"Linear stability region for\n{table.name}", fontsize=font_size)
    ax.tick_params(axis="both", labelsize=font_size)
    if legend_handles:
        ax.legend(handles=legend_handles, loc="upper right", fontsize=0.8 * font_size)
    # Put the grid under the contourf fills so the fill alpha fades it; the default ordering
    # ('line', zorder 1.5) draws the grid over the fill.
    ax.set_axisbelow(True)
    ax.grid(True, alpha=0.25)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def plot_stability_region(
    table: ButcherTable,
    options: PlotOptions | None = None,
    *,
    filename: str | None = None,
    ax=None,
    close: bool = False,
):
    """Plot { z : |phi(z)| <= 1 } with the boundary curve.

    Returns (fig, ax). When *ax* is supplied the caller owns the figure, so *filename* and *close*
    are ignored.
    """
    options = options or PlotOptions()

    phi = table.stability_function
    emb = table.embedded_stability_function if options.show_embedding else None
    shade_emb = options.shade_embedding and emb is not None

    bounds = options.bounds or _auto_bounds(table, include_embedding=options.show_embedding)

    # The canvas is fixed, so the framing is adjusted as needed. This has to happen before the grid
    # is built: padding the box afterwards would leave the added strip outside the evaluated grid,
    # and a region that runs off the scan would be cut short mid-axes.
    created = ax is None
    cax = None
    if created:
        colorbar = options.shade or shade_emb
        rect, axes_aspect = _axes_geometry(options.figure_size, colorbar)
        fig = plt.figure(figsize=options.figure_size)
        ax = fig.add_axes(rect)
        if colorbar:
            # Do not use colorbar(ax=ax), it will shrink the plot axes and break the fixed geometry
            cax = fig.add_axes(_colorbar_rect(options.figure_size))
        bounds = _fit_bounds_to_aspect(bounds, axes_aspect)
    else:
        # The caller owns the geometry, so its axes box is the only thing that can decide the
        # aspect; leave the framing tight and let equal aspect letterbox it.
        fig = ax.figure

    # Evaluate |phi| over the final framing, so the grid and the axes limits agree.
    axis_x = np.linspace(bounds[0], bounds[1], options.grid_points)
    axis_y = np.linspace(bounds[2], bounds[3], options.grid_points)
    X, Y = np.meshgrid(axis_x, axis_y)
    Z = X + 1j * Y

    # Whichever region is shaded gets the colormap; the other is drawn as a bare contour so the two
    # do not fight over the same visual channel.
    region_handles = []
    main_handle = _render_region(
        ax,
        fig,
        X,
        Y,
        _region_magnitude(phi, Z, X, Y, options.suppress_tiny_islands),
        style=MAIN_REGION_STYLE,
        mode="shade" if options.shade else ("boundary" if shade_emb else "fill"),
        label=None if options.shade else "main method",
        cax=cax,
        font_size=options.font_size,
    )
    if main_handle is not None:
        region_handles.append(main_handle)

    if emb is not None:
        emb_handle = _render_region(
            ax,
            fig,
            X,
            Y,
            _region_magnitude(emb, Z, X, Y, options.suppress_tiny_islands),
            style=EMBEDDED_REGION_STYLE,
            mode="shade" if shade_emb else ("boundary" if options.shade else "fill"),
            label="embedded method",
            cax=cax,
            font_size=options.font_size,
        )
        if emb_handle is not None:
            region_handles.append(emb_handle)

    marker_handles = _draw_zeros_poles(ax, phi, options)

    # Real and imaginary axes
    ax.axhline(0.0, color="0.4", lw=0.8)
    ax.axvline(0.0, color="0.4", lw=0.8)

    show_region = options.shade or shade_emb or len(region_handles) > 1
    _decorate_axes(
        ax,
        table,
        bounds,
        (region_handles if show_region else []) + marker_handles,
        options.font_size,
    )

    # Only touch the file and the figure lifetime when we made the figure.
    if created and filename:
        fig.savefig(filename, dpi=options.dpi)
    if created and close:
        plt.close(fig)
    return fig, ax
