import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

import util
from util import *
from itertools import combinations
from extq import projection

cm_seq = sns.cubehelix_palette(
    start=0, rot=-0.70, gamma=0.40, light=0.9, dark=0.1, as_cmap=True, reverse=True
)
colors = mpl.colors.to_rgba_array(
    [
        "#364B9A",
        "#4A7BB7",
        "#6EA6CD",
        "#98CAE1",
        "#C2E4EF",
        "#EAECCC",
        "#FEDA8B",
        "#FDB366",
        "#F67E4B",
        "#DD3D2D",
        "#A50026",
    ]
)
cm_div = mpl.colors.LinearSegmentedColormap.from_list("", colors)


def format_cvs(ax, centroids=True, **kwargs):
    # sets the axis labels and includes the locations of the 4 states for a CV plot
    # ax: matplotlib.pyplot Axis object
    if not centroids:
        ax.plot(-8.919, -109.9, "s", mec="k", mfc="w", **kwargs)
        ax.plot(-4.174, -50.8, "o", mec="k", mfc="w", **kwargs)
        ax.plot(0, 0, "^", mec="k", mfc="w", **kwargs)
        ax.plot(4.535, 43.7, "D", mec="k", mfc="w", **kwargs)
    elif centroids:
        ax.plot(-7.863, -108.6, "s", mec="k", mfc="w", **kwargs)
        ax.plot(-4.240, -56.95, "o", mec="k", mfc="w", **kwargs)
        ax.plot(-0.506, 3.940, "^", mec="k", mfc="w", **kwargs)
        ax.plot(6.447, 60.713, "D", mec="k", mfc="w", **kwargs)

    ax.set(
        xlabel="Translocation (Å)",
        ylabel="Rotation ($^\circ$)",
        xlim=[-10, 10],
        ylim=[-150, 100],
    )


def plot_models(ax, points, **kwargs):
    """Plots four formatted points corresponding to the models in some projection
    ax: matplotlib.pyplot Axis object
    points: np.ndarray(4, 2)
    """
    if points.shape != (4, 2):
        raise ValueError
    ax.plot(points[0, 0], points[0, 1], "s", mec="k", mfc="w", label="Down-", **kwargs)
    ax.plot(points[1, 0], points[1, 1], "o", mec="k", mfc="w", label="Down", **kwargs)
    ax.plot(points[2, 0], points[2, 1], "^", mec="k", mfc="w", label="Up", **kwargs)
    ax.plot(points[3, 0], points[3, 1], "D", mec="k", mfc="w", label="Up+", **kwargs)


def plot_single_model(ax, point, label, **kwargs):
    """
    ax : matplotlib.pyplot Axis
    point : array-like of length 2
    label : str, "dd", "d", "u", or "uu"
    """
    if label not in ("dd", "d", "u", "uu"):
        raise ValueError("Label must be one of 'dd', 'd', 'u', or 'uu'")
    if len(point) != 2:
        raise ValueError("Point must be of length 2")
    if label == "dd":
        ax.plot(point[0], point[1], "s", mec="k", mfc="w", label="Down-", **kwargs)
    elif label == "d":
        ax.plot(point[0], point[1], "o", mec="k", mfc="w", label="Down", **kwargs)
    elif label == "u":
        ax.plot(point[0], point[1], "^", mec="k", mfc="w", label="Up", **kwargs)
    elif label == "uu":
        ax.plot(point[0], point[1], "D", mec="k", mfc="w", label="Up+", **kwargs)


def get_sb_models(centroids=True):
    if centroids:
        feat2_models = np.load(
            "/project/dinner/scguo/ci-vsd/data/models_centroids_feat2.npy"
        )
    else:
        feat2_models = np.load("/project/dinner/scguo/ci-vsd/data/models_feat2.npy")
    return feat2_models


def plot_sb_models(axes, mode="du", centroids=True, **kwargs):
    # plot models in salt bridge distances
    # use 6 different combinations
    # mode specifies for which transition to plot, chooses
    # appropriate salt bridges
    # centroids plots MD-centroids or crystal structures
    feat2_models = get_sb_models(centroids=centroids)
    if mode == "dd":
        plot_sb_models_helper(axes, [36, 42, 41, 47], feat2_models, **kwargs)
    elif mode == "du":
        plot_sb_models_helper(axes, [42, 48, 47, 53], feat2_models, **kwargs)
    elif mode == "uu":
        plot_sb_models_helper(axes, [48, 54, 53, 59], feat2_models, **kwargs)
    else:
        raise ValueError("Invalid mode")


def plot_sb_models_helper(axes, idxs, models, **kwargs):
    # axes must be of length 6
    # idxs must be of length 4
    # plots all combinations
    # models must be of shape[4, ...]
    if len(axes) != 6:
        raise ValueError("Number of axes must be 6")
    if len(idxs) != 4:
        raise ValueError("Number of salt bridge ids must be 4")
    for (ax, (id1, id2)) in zip(axes, combinations(idxs, 2)):
        skip = id2 - id1
        if skip >= 0:
            final = id1 + skip + 1
        elif skip < 0:
            final = id1 + skip - 1
        plot_models(ax, models[:, id1:final:skip], **kwargs)


def plot_pmf(ax, pmf, xlim, ylim, units="kT", clines=None, cmap=None):
    if clines is None:
        clines = np.linspace(0, 6, 7)
    if cmap is None:
        cmap = cm_div
    if units != "kT" and units != "kcal":
        raise ValueError("Units must be in kT or kcal!")

    # compute grid
    centerx = (xlim[1:] + xlim[:-1]) / 2
    centery = (ylim[1:] + ylim[:-1]) / 2

    # calculate energy difference against minimum
    min_energy = np.min(-np.log(pmf[np.nonzero(pmf)]))
    diff = -np.log(pmf.T) - min_energy
    if units == "kcal":
        diff *= 0.593
    h = ax.pcolor(xlim, ylim, diff, cmap=cmap)
    ax.contour(
        centerx,
        centery,
        diff,
        levels=clines,
        colors="black",
        linestyles="solid",
    )

    return ax, h


def plot_sb_pmfs(
    sb_ids,
    pmfs_sb,
    sb1_lim,
    sb2_lim,
    clines=None,
    cmap=None,
    units="kT",
    fig_kwargs=None,
    **kwargs,
):
    """Plot PMFs for a series of salt bridge variables.

    Parameters
    ---------
    sb_ids : list of 4 int
    pmfs_sb : list of histogrammed weights
    sb1_lim :
    sb2_lim :
    clines : optional, contour lines
    cmap : optional, colormap
    units : str, 'kT' or 'kcal'
    fig_kwargs : dictionary of keyword arguments to create the figure

    Returns
    -------
    f : Figure object
    axes : axes
    """
    if clines is None:
        clines = np.linspace(0, 6, 7)
    if cmap is None:
        cmap = cm_div
    if units != "kT" and units != "kcal":
        raise ValueError("Units must be in kT or kcal!")
    if fig_kwargs is None:
        fig_kwargs = {"sharex": True, "sharey": True, "constrained_layout": True}

    f, axes = plt.subplots(2, 3, **fig_kwargs)

    for ((sb1, sb2), ax, p, xlim, ylim) in zip(
        combinations(sb_ids, 2), axes.flat, pmfs_sb, sb1_lim, sb2_lim
    ):
        ax, h = plot_pmf(ax, p, xlim, ylim, units=units, clines=clines, cmap=cmap)
        ax.set_xlabel(f"{sb_label(sb1)}")
        ax.set_ylabel(f"{sb_label(sb2)}")
        ax.set_xlim([0.3, 2.0])
        ax.set_ylim([0.3, 2.0])

    cb = plt.colorbar(h, ax=axes[:, -1])
    if units == "kT":
        cb.set_label("kT", rotation=-90, labelpad=10)
    elif units == "kcal":
        cb.set_label("kcal / mol", rotation=-90, labelpad=10)

    return f, axes


def plot_sb_q(
    sb_ids,
    sb_trajs,
    q,
    weights,
    sb1_lim=(3, 20),
    sb2_lim=(3, 20),
    bins=100,
    cmap=cm_div,
    centroids=True,
    fig_kwargs=None,
    **kwargs,
):
    """Plot histogrammed committor for a series of salt bridge variables.

    Parameters
    ---------
    sb_ids : list of 4 int
    sb_trajs : list of ndarray (nframes, n_dim)
        Salt bridge trajectories
    q : list of ndarray (nframes,)
        Committors for each frame
    sb1_lim : optional, array-like of length 2 or shape (6, 2)
        Range for first axis. Can be the same for each plot or
        specified individually.
    sb2_lim : optional, array-like of length 2 or shape (6, 2)
    bins : int
    cmap : optional, colormap
    centroids : bool
        Whether to plot centroids or crystal structure models
    fig_kwargs : dictionary of keyword arguments to create the figure
    **kwargs :
        Keyword arguments to plot models.

    Returns
    -------
    f : Figure object
    axes : axes
    """
    if cmap is None:
        cmap = cm_div
    if fig_kwargs is None:
        fig_kwargs = {"sharex": True, "sharey": True, "constrained_layout": True}
    if len(sb1_lim) == 2 and len(sb2_lim) == 2:
        sb1_lim = [sb1_lim] * 6
        sb2_lim = [sb2_lim] * 6

    f, axes = plt.subplots(2, 3, **fig_kwargs)

    sb_models = get_sb_models(centroids=centroids)
    for ((sb1, sb2), ax, lim1, lim2) in zip(
        combinations(sb_ids, 2), axes.flat, sb1_lim, sb2_lim
    ):
        traj_1 = [traj[:, sb1] * 10 for traj in sb_trajs]
        traj_2 = [traj[:, sb2] * 10 for traj in sb_trajs]
        xlim = np.linspace(*lim1, bins + 1)
        ylim = np.linspace(*lim2, bins + 1)
        qbin = projection.average2d(traj_1, traj_2, q, weights, xlim, ylim)
        h = ax.pcolor(xlim, ylim, qbin.T, cmap=cmap, vmin=0, vmax=1)
        ax.set_xlabel(f"{util.sb_label(sb1)} distance / $\AA$")
        ax.set_ylabel(f"{util.sb_label(sb2)} distance / $\AA$")
        ax.set_xlim(lim1)
        ax.set_ylim(lim2)
        plot_models(ax, sb_models[:, [sb1, sb2]] * 10, **kwargs)

    cb = f.colorbar(h, ax=axes[:, -1])
    cb.set_label("$q_+$", rotation=-90, labelpad=10)

    return f, axes


def multiline(xs, ys, c, ax=None, **kwargs):
    """Plot lines with different colorings

    Parameters
    ----------
    xs : iterable container of x coordinates
    ys : iterable container of y coordinates
    c : iterable container of numbers mapped to colormap
    ax (optional): Axes to plot on.
    kwargs (optional): passed to LineCollection

    Notes:
        len(xs) == len(ys) == len(c) is the number of line segments
        len(xs[i]) == len(ys[i]) is the number of points for each line (indexed by i)

    Returns
    -------
    lc : LineCollection instance.
    """

    # find axes
    ax = plt.gca() if ax is None else ax

    # create LineCollection
    segments = [np.column_stack([x, y]) for x, y in zip(xs, ys)]
    lc = mpl.collections.LineCollection(segments, **kwargs)

    # set coloring of line segments
    #    Note: I get an error if I pass c as a list here... not sure why.
    lc.set_array(np.asarray(c))

    # add lines to axes and rescale
    #    Note: adding a collection doesn't autoscalee xlim/ylim
    ax.add_collection(lc)
    ax.autoscale()
    return lc
