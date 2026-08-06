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
from matplotlib import patheffects
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.transforms import ScaledTranslation

from .rk_butcher_table import ButcherTable
from .rk_stability import (
    TOL,
    axis_extent,
    imag_stable_intervals,
    max_axis_crossing,
    stability_magnitude,
)

# ---------------------------------------------------------------------------
# Output canvas
# ---------------------------------------------------------------------------

# Every figure is rendered on one fixed canvas with the axes box pinned to fixed
# margins, so a directory of plots is a single image shape with the axes in the same
# place in each one -- what documentation needs to lay them out without per-figure
# scaling.  Uniformity is paid for in framing, not in shape: the data box is *grown* to
# the aspect of the axes box (see :func:`_fit_bounds_to_aspect`) so `set_aspect("equal")`
# is satisfied exactly.  Nothing is stretched and nothing is cropped.
#
# The default axes aspect, 5.5 / 4.5 = 1.22, is the median data-box aspect over the
# shipped ERK and DIRK tables; stability regions run tall because the real extent is
# one-sided while the imaginary extent is symmetric.  Sizes are inches.
_FIG_SIZE = (5.5, 6.5)
_FIG_DPI = 150

# Margins must hold the tick labels, axis labels, and title at every framing.
_MARGIN_LEFT = 0.78
_MARGIN_RIGHT = 0.22
_MARGIN_BOTTOM = 0.55
_MARGIN_TOP = 0.45

# Reserved on the right when a region is shaded; the bar itself is narrower than the
# strip, which also has to hold the tick labels and the bar's own axis label.
_COLORBAR_STRIP = 1.30
_COLORBAR_WIDTH = 0.20
_COLORBAR_PAD = 0.12


# ---------------------------------------------------------------------------
# Options
# ---------------------------------------------------------------------------


@dataclass
class PlotOptions:
    """Everything that varies between stability-region plots.

    Field names match the CLI's argparse destinations so the two stay in step without
    the option list being written out a third time.
    """

    bounds: tuple[float, float, float, float] | None = None
    grid_points: int = 500
    # Output canvas, in inches and dots per inch.  Fixed rather than fitted, so
    # figure_size * dpi is the pixel size of every image this module writes.
    figure_size: tuple[float, float] = _FIG_SIZE
    dpi: int = _FIG_DPI
    # Zeros/poles are informative annotations only; they do not change the shaded region.
    show_zeros: bool = False
    show_poles: bool = False
    show_embedding: bool = True
    shade: bool = False
    shade_embedding: bool = False
    suppress_pinholes: bool = True
    show_axis_extent: bool = True
    highlight_imag_axis_intervals: bool = False
    label_detached_imag_axis_intervals: bool = False

    def __post_init__(self):
        if self.shade and self.shade_embedding:
            raise ValueError("shade and shade_embedding are mutually exclusive")
        if self.bounds is not None:
            self.bounds = tuple(float(v) for v in self.bounds)
        self.figure_size = tuple(float(v) for v in self.figure_size)
        # Fail here rather than with a negative axes rect deep inside the draw.
        _axes_geometry(self.figure_size, colorbar=True)


# ---------------------------------------------------------------------------
# Styling
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class RegionStyle:
    """Colours, line properties, and annotation side for one method's region.

    ``side`` is +1 for the main method and -1 for the embedding; every mirrored label
    offset is derived from it, so the two annotation passes share one code path.
    """

    fill_color: str
    fill_alpha: float
    boundary_color: str
    boundary_lw: float
    boundary_ls: str
    shade_label: str
    side: int
    legend_fill_alpha: float | None = None
    # Only the main region needs a different boundary colour on top of the viridis fill.
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
    shade_label=r"$|\phi(z)|$  within the main stable region",
    side=+1,
    shade_boundary_color="k",
)

EMBEDDED_REGION_STYLE = RegionStyle(
    fill_color="#ffcccc",
    fill_alpha=0.45,
    boundary_color="#cc0000",
    boundary_lw=1.6,
    boundary_ls="--",
    shade_label=r"$|\hat\phi(z)|$  within the embedded stable region",
    side=-1,
    legend_fill_alpha=0.7,
)

_SHADE_LEVELS = np.linspace(0.0, 1.0, 21)

# An origin-attached imaginary interval shorter than this is drawn with a compact single
# label instead of two: below this length the +/- endpoints are visually indistinguishable
# from a mere tangency at the origin.
_MIN_IMAG_EXTENT_LABEL = 1e-2

# Annotation geometry, in points; all mirrored by RegionStyle.side.
_REAL_LABEL_DY = -11
_IMAG_LABEL_DX = 7
_SMALL_IMAG_LABEL_DY = 10
_IMAG_RAIL_DX = 4.0
_IMAG_RAIL_LW = 4.0
_AXIS_LABEL_FONTSIZE = 8
_AXIS_MARKER_SIZE = 5

# Auto-framing: scan a coarse grid over a span wide enough to hold every bounded feature,
# then crop to what was found.  The scan costs one extra evaluation of phi at lower
# resolution than the plot grid itself.
_BOUNDS_SCAN_POINTS = 420
_BOUNDS_MIN_SPAN = 10.0
_BOUNDS_STAGE_PAD = 4.0
_BOUNDS_CROSSING_PAD = 1.3
_BOUNDS_MAX_SPAN = 80.0
_BOUNDS_PAD_FRAC = 0.20
_BOUNDS_PAD_MIN = 0.5
_BOUNDS_FALLBACK_MIN = 6.0
_BOUNDS_FALLBACK_STAGE_PAD = 2.0

# Any value above 1 reads as "unstable"; used to paint over suppressed pinhole islands.
_UNSTABLE_FILL = 2.0
_ISLAND_MIN_RADIUS = 0.06

_AXIS_LABEL_BOX = dict(boxstyle="round,pad=0.15", fc="white", ec="0.6", lw=0.5, alpha=0.85)


# ---------------------------------------------------------------------------
# Region rendering
# ---------------------------------------------------------------------------


def _region_magnitude(p, q, Z, X, Y, suppress_pinholes: bool):
    """|phi| on the plot grid, optionally with microscopic stable islands masked out."""
    R = stability_magnitude(p, q, Z)
    if suppress_pinholes:
        R = _suppress_tiny_islands(R, X, Y, p.r)
    return R


def _render_region(ax, fig, X, Y, R, *, style: RegionStyle, mode: str, label: str | None, cax=None):
    """Draw one method's region and return its legend artist (or None).

    ``mode`` is 'shade' (|phi| colormap plus colorbar), 'fill' (flat stable fill), or
    'boundary' (contour only).  All three draw the |phi| = 1 contour.  *cax* is the
    pre-placed colorbar axes; without one the bar is carved out of *ax*, which is only
    appropriate when the caller owns the figure geometry.
    """
    if mode == "shade":
        color = style.shade_color
        contours = ax.contourf(X, Y, R, levels=_SHADE_LEVELS, cmap="viridis")
        bar = fig.colorbar(contours, cax=cax) if cax is not None else fig.colorbar(contours, ax=ax)
        bar.set_label(style.shade_label)
    else:
        color = style.boundary_color
        if mode == "fill":
            ax.contourf(
                X, Y, R, levels=[0.0, 1.0], colors=[style.fill_color], alpha=style.fill_alpha
            )

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
    """Mask out microscopic stable islands around zeros of phi.

    They render as stray dots rather than meaningful features.  Each zero seeds a
    flood fill over the stable mask; components that stay below the area a disc of
    ``min_radius`` would cover are painted unstable.
    """
    if zeros is None or len(zeros) == 0:
        return R
    ny, nx = R.shape
    dx = X[0, 1] - X[0, 0]
    dy = Y[1, 0] - Y[0, 0]
    x0, y0 = X[0, 0], Y[0, 0]
    cap = max(8, int(np.pi * min_radius**2 / abs(dx * dy)))
    stable = R <= 1.0
    R = R.copy()
    for z in zeros:
        j = int(round((z.real - x0) / dx))
        i = int(round((z.imag - y0) / dy))
        if not (0 <= i < ny and 0 <= j < nx) or not stable[i, j]:
            continue
        island = {(i, j)}
        stack = [(i, j)]
        overflow = False
        while stack:
            ci, cj = stack.pop()
            if len(island) > cap:
                overflow = True
                break
            for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                ni, nj = ci + di, cj + dj
                if 0 <= ni < ny and 0 <= nj < nx and (ni, nj) not in island and stable[ni, nj]:
                    island.add((ni, nj))
                    stack.append((ni, nj))
        if not overflow:
            for ci, cj in island:
                R[ci, cj] = _UNSTABLE_FILL
    return R


# ---------------------------------------------------------------------------
# Framing and figure setup
# ---------------------------------------------------------------------------


def _auto_bounds(table: ButcherTable, include_embedding: bool = True):
    """Choose a plotting box that frames the interesting bounded features.

    The scan span is widened to hold every axis crossing, then capped at
    ``_BOUNDS_MAX_SPAN``. The cap is a legibility limit, not an estimate: several ARKODE
    embeddings cross the real axis past t = 100 while their main method turns over inside
    t = 20, and honouring both would give a box in which the main region is a few pixels
    wide.

    A method whose region does not close inside the capped scan is *deliberately not
    framed* -- its mask reaches the scan edge, so it is dropped and the box is fitted to
    whichever methods do close. Its contour is still drawn and will simply leave the axes.
    Five of the shipped tables hit this, all of them embeddings. If every method overflows
    there is nothing to fit and the stage-count fallback box is used instead.
    """
    methods = [table.stability_function]
    emb = table.embedded_stability_function if include_embedding else None
    if emb is not None:
        methods.append(emb)

    crossings = [max_axis_crossing(phi) for phi in methods]
    span = max(
        [_BOUNDS_MIN_SPAN, table.stages + _BOUNDS_STAGE_PAD]
        + [_BOUNDS_CROSSING_PAD * c for c in crossings if c > 0]
    )
    span = min(span, _BOUNDS_MAX_SPAN)

    axis = np.linspace(-span, span, _BOUNDS_SCAN_POINTS)
    X, Y = np.meshgrid(axis, axis)
    Z = X + 1j * Y

    def touches_edge(mask) -> bool:
        return bool(mask[0, :].any() or mask[-1, :].any() or mask[:, 0].any() or mask[:, -1].any())

    def feature_mask(p, q):
        """The bounded side of the boundary, or None if neither side closes inside the scan.

        None is the overflow signal described above: the region runs past the capped span,
        so this method contributes nothing to the framing.
        """
        R = stability_magnitude(p, q, Z)
        for mask in (R <= 1.0, R > 1.0):
            if mask.any() and not touches_edge(mask):
                return mask
        return None

    masks = [m for m in (feature_mask(phi.p, phi.q) for phi in methods) if m is not None]
    features = np.logical_or.reduce(masks) if masks else None

    if features is not None and features.any():
        xv, yv = X[features], Y[features]
        xmin, xmax = min(xv.min(), 0.0), max(xv.max(), 0.0)
        ymin, ymax = min(yv.min(), 0.0), max(yv.max(), 0.0)
        padx = max(_BOUNDS_PAD_MIN, _BOUNDS_PAD_FRAC * (xmax - xmin))
        pady = max(_BOUNDS_PAD_MIN, _BOUNDS_PAD_FRAC * (ymax - ymin))
        return (xmin - padx, xmax + padx, ymin - pady, ymax + pady)

    w = max(_BOUNDS_FALLBACK_MIN, table.stages + _BOUNDS_FALLBACK_STAGE_PAD)
    return (-w, w, -w, w)


def _axes_geometry(figure_size, colorbar: bool):
    """Axes rect in figure fractions, and the data aspect that fills it exactly.

    Returns ``(rect, aspect)`` where *aspect* is the height/width of the axes box in
    inches.  A data box of that same aspect satisfies ``set_aspect("equal")`` with no
    room left over, so the axes are never shrunk away from *rect*.
    """
    fig_w, fig_h = figure_size
    right = _MARGIN_RIGHT + (_COLORBAR_STRIP if colorbar else 0.0)
    axes_w = fig_w - _MARGIN_LEFT - right
    axes_h = fig_h - _MARGIN_BOTTOM - _MARGIN_TOP
    if axes_w <= 0.0 or axes_h <= 0.0:
        raise ValueError(
            f"figure size {tuple(figure_size)} leaves no room for the axes; "
            f"it must exceed {_MARGIN_LEFT + right:.2f} x "
            f"{_MARGIN_BOTTOM + _MARGIN_TOP:.2f} inches"
        )
    rect = (_MARGIN_LEFT / fig_w, _MARGIN_BOTTOM / fig_h, axes_w / fig_w, axes_h / fig_h)
    return rect, axes_h / axes_w


def _colorbar_rect(figure_size):
    """Rect for the colorbar's own axes, inside the strip reserved on the right.

    The bar is placed rather than carved out of the plot axes with ``colorbar(ax=ax)``,
    which would shrink the axes and undo the fixed geometry.
    """
    fig_w, fig_h = figure_size
    left = fig_w - _MARGIN_RIGHT - _COLORBAR_STRIP + _COLORBAR_PAD
    height = fig_h - _MARGIN_BOTTOM - _MARGIN_TOP
    return (left / fig_w, _MARGIN_BOTTOM / fig_h, _COLORBAR_WIDTH / fig_w, height / fig_h)


def _fit_bounds_to_aspect(bounds, aspect: float):
    """Grow *bounds* about its centre until its height/width equals *aspect*.

    The framing this widens has already been chosen to hold the interesting features, so
    the fit only ever enlarges: no feature is cropped and no axis is rescaled relative to
    the other.  The surplus is filled with more of the complex plane, which is the same
    picture the letterboxing alternative would show, minus the blank margin beside it.
    """
    xmin, xmax, ymin, ymax = bounds
    width, height = xmax - xmin, ymax - ymin
    if height < aspect * width:
        grow = 0.5 * (aspect * width - height)
        return (xmin, xmax, ymin - grow, ymax + grow)
    grow = 0.5 * (height / aspect - width)
    return (xmin - grow, xmax + grow, ymin, ymax)


# ---------------------------------------------------------------------------
# Axis annotations
# ---------------------------------------------------------------------------


def _plot_axis_markers(ax, xs, ys, color):
    ax.plot(xs, ys, marker="D", color=color, ms=_AXIS_MARKER_SIZE, ls="none", zorder=6)


def _annotate_axis_label(ax, text, xy, offset, *, ha, va, color):
    ax.annotate(
        text,
        xy=xy,
        xytext=offset,
        textcoords="offset points",
        ha=ha,
        va=va,
        fontsize=_AXIS_LABEL_FONTSIZE,
        color=color,
        bbox=_AXIS_LABEL_BOX,
        zorder=6,
    )


def _imag_endpoint_label_specs(t: float):
    """(y, text, dy, va) for the +/- labels of an imaginary-axis endpoint at height t."""
    return ((t, f"${t:.2f}\\,i$", 4, "bottom"), (-t, f"$-{t:.2f}\\,i$", -4, "top"))


def _small_imag_extent_label(value: float, precision: int = 2) -> str:
    mantissa, exponent = f"{value:.{precision}e}".split("e")
    return f"$\\pm{mantissa} \\times 10^{{{int(exponent)}}}\\,i$"


def _highlight_imag_axis_intervals(ax, intervals, style: RegionStyle):
    """Draw offset vertical rails marking finite imaginary-axis stability intervals."""
    dx_points = _IMAG_RAIL_DX * style.side
    transform = ax.transData + ScaledTranslation(dx_points / 72.0, 0.0, ax.figure.dpi_scale_trans)
    effects = [patheffects.withStroke(linewidth=6.0, foreground="white", alpha=0.95)]

    for t_left, t_right in intervals:
        if t_left <= TOL.root_zero:
            rails = [(-t_right, t_right)]
        else:
            rails = [(t_left, t_right), (-t_right, -t_left)]
        for y0, y1 in rails:
            (line,) = ax.plot(
                [0.0, 0.0],
                [y0, y1],
                transform=transform,
                color=style.boundary_color,
                lw=_IMAG_RAIL_LW,
                ls=style.boundary_ls,
                solid_capstyle="round",
                zorder=5.5,
            )
            line.set_path_effects(effects)


def _mark_axis_extents(ax, phi, style: RegionStyle, options: PlotOptions):
    """Mark and label the real/imaginary stability interval endpoints for one method."""
    side = style.side
    color = style.boundary_color
    imag_dx = _IMAG_LABEL_DX * side
    imag_ha = "left" if side > 0 else "right"

    intervals = imag_stable_intervals(phi)
    if not options.label_detached_imag_axis_intervals:
        intervals = [iv for iv in intervals if iv[0] <= TOL.root_zero]

    real_extent = axis_extent(phi, "real")
    if real_extent:
        real_dy = _REAL_LABEL_DY * side
        _plot_axis_markers(ax, [-real_extent], [0.0], color)
        _annotate_axis_label(
            ax,
            f"$-{real_extent:.2f}$",
            (-real_extent, 0.0),
            (0, real_dy),
            ha="center",
            va="top" if real_dy < 0 else "bottom",
            color=color,
        )

    if options.highlight_imag_axis_intervals:
        _highlight_imag_axis_intervals(ax, intervals, style)

    for t_left, t_right in intervals:
        attached_to_origin = t_left <= TOL.root_zero

        if attached_to_origin and t_right <= _MIN_IMAG_EXTENT_LABEL:
            # Too short to label at its endpoints: one compact label at the origin.
            small_dy = _SMALL_IMAG_LABEL_DY * side
            _plot_axis_markers(ax, [0.0], [0.0], color)
            _annotate_axis_label(
                ax,
                _small_imag_extent_label(t_right),
                (0.0, 0.0),
                (imag_dx, small_dy),
                ha=imag_ha,
                va="bottom" if small_dy > 0 else "top",
                color=color,
            )
            continue

        endpoints = (t_right,) if attached_to_origin else (t_left, t_right)
        _plot_axis_markers(
            ax, [0.0] * (2 * len(endpoints)), [s * t for t in endpoints for s in (1, -1)], color
        )
        for t in endpoints:
            for y, text, dy, va in _imag_endpoint_label_specs(t):
                _annotate_axis_label(
                    ax, text, (0.0, y), (imag_dx, dy), ha=imag_ha, va=va, color=color
                )


def _draw_zeros_poles(ax, p, q, options: PlotOptions):
    """Optionally mark the zeros and poles of phi; returns their legend handles."""
    handles = []
    if options.show_zeros and p.r.size:
        handles.append(
            ax.plot(p.r.real, p.r.imag, "kx", ms=7, mew=1.5, label=r"zeros of $\phi$")[0]
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
                label=r"poles of $\phi$",
            )[0]
        )
    return handles


def _decorate_axes(ax, table: ButcherTable, bounds, legend_handles):
    """Axis limits, labels, title, legend, and grid."""
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(bounds[0], bounds[1])
    ax.set_ylim(bounds[2], bounds[3])
    ax.set_xlabel(r"Re$(z)$")
    ax.set_ylabel(r"Im$(z)$")
    # Two lines, always: the longest ARKODE names overrun a fixed canvas on one line, and
    # a title of constant height keeps the top margin -- and so the axes box -- constant.
    ax.set_title(f"Linear stability region for\n{table.name}", fontsize=10)
    if legend_handles:
        ax.legend(handles=legend_handles, loc="upper right", fontsize=8)
    # Matplotlib's default axisbelow of 'line' puts the grid at zorder 1.5, above the
    # contourf fills at zorder 1, so the grid is composited over a region at full
    # strength while outside it only has white to work against.  Dropping it to 0.5 puts
    # it under the fills, which then attenuate it by their own alpha -- the grid reads as
    # background everywhere, and a region reads as something laid over it.  The axis
    # lines below stay above the fills; they are annotation, not background.
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

    Returns (fig, ax).  When *ax* is supplied the caller owns the figure, so *filename*
    and *close* are ignored.
    """
    options = options or PlotOptions()

    phi = table.stability_function
    p, q = phi.p, phi.q
    emb = table.embedded_stability_function if options.show_embedding else None
    shade_emb = options.shade_embedding and emb is not None

    bounds = options.bounds or _auto_bounds(table, include_embedding=options.show_embedding)

    # The canvas is fixed, so the framing is what adapts.  This has to happen before the
    # grid is built: padding the box afterwards would leave the added strip outside the
    # evaluated grid, and a region that runs off the scan would be cut short mid-axes.
    created = ax is None
    cax = None
    if created:
        colorbar = options.shade or shade_emb
        rect, axes_aspect = _axes_geometry(options.figure_size, colorbar)
        fig = plt.figure(figsize=options.figure_size)
        ax = fig.add_axes(rect)
        if colorbar:
            cax = fig.add_axes(_colorbar_rect(options.figure_size))
        bounds = _fit_bounds_to_aspect(bounds, axes_aspect)
    else:
        # The caller owns the geometry, so its axes box is the only thing that can
        # decide the aspect; leave the framing tight and let equal aspect letterbox it.
        fig = ax.figure

    axis_x = np.linspace(bounds[0], bounds[1], options.grid_points)
    axis_y = np.linspace(bounds[2], bounds[3], options.grid_points)
    X, Y = np.meshgrid(axis_x, axis_y)
    Z = X + 1j * Y

    # Whichever region is shaded gets the colormap; the other is drawn as a bare contour
    # so the two do not fight over the same visual channel.
    region_handles = []
    main_handle = _render_region(
        ax,
        fig,
        X,
        Y,
        _region_magnitude(p, q, Z, X, Y, options.suppress_pinholes),
        style=MAIN_REGION_STYLE,
        mode="shade" if options.shade else ("boundary" if shade_emb else "fill"),
        label=None if options.shade else "main method",
        cax=cax,
    )
    if main_handle is not None:
        region_handles.append(main_handle)

    if emb is not None:
        pe, qe = emb.p, emb.q
        label = "embedded method"
        if table.embedding_order is not None:
            label += f" (order {table.embedding_order})"
        emb_handle = _render_region(
            ax,
            fig,
            X,
            Y,
            _region_magnitude(pe, qe, Z, X, Y, options.suppress_pinholes),
            style=EMBEDDED_REGION_STYLE,
            mode="shade" if shade_emb else ("boundary" if options.shade else "fill"),
            label=label,
            cax=cax,
        )
        if emb_handle is not None:
            region_handles.append(emb_handle)

    marker_handles = _draw_zeros_poles(ax, p, q, options)

    ax.axhline(0.0, color="0.4", lw=0.8)
    ax.axvline(0.0, color="0.4", lw=0.8)

    if options.show_axis_extent:
        _mark_axis_extents(ax, phi, MAIN_REGION_STYLE, options)
        if emb is not None:
            _mark_axis_extents(ax, emb, EMBEDDED_REGION_STYLE, options)

    # A single unshaded region needs no legend entry of its own; its label is the title.
    show_region = options.shade or shade_emb or len(region_handles) > 1
    _decorate_axes(ax, table, bounds, (region_handles if show_region else []) + marker_handles)

    if created and filename:
        fig.savefig(filename, dpi=options.dpi)
    if created and close:
        plt.close(fig)
    return fig, ax
