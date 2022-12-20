# Useful functions for loading files, plotting, etc.

import glob
import numpy as np
import scipy
from scipy import sparse
from scipy import signal

# def load_ds_trajs():
#     # load trajectories for distance and spin CVs
#     # returns a list of numpy arrays, one for each trajectory
#     cv_files_amber = sorted(
#         glob.glob(
#             "/project2/roux/scguo/ci-vsd/amber-gpu/\
#                 msm-projection-300ns/TRANSFORM_2D_npy/*"
#         )
#     )
#     cv_trajs = []
#     for f in cv_files_amber:
#         data = np.load(f, allow_pickle=True)
#         cv_trajs.append(data)
#
#     ds_anton_files = sorted(
#         glob.glob("/project2/roux/scguo/ci-vsd/data/projection_anton/*")
#     )
#     for f in ds_anton_files:
#         data = np.loadtxt(f)
#         cv_trajs.append(data[1:])
#
#     return cv_trajs


def frame(idx):
    # Find the trajectory number and corresponding frame number to the data
    # set index (between 0 and 6 580 000)

    if idx < 3000000:
        # in amber-gpu trajectories (1000 x 3000 frames)
        traj = idx // 3000
        frame = idx % 3000
        return (traj, frame)
    elif 3000000 <= idx < 4790000:
        # first group of Anton trajectories (179 x 10 000 frames)
        idx -= 3000000
        traj = idx // 10000
        frame = idx % 10000
        return (traj + 1000, frame)
    elif 4790000 <= idx < 5490000:
        # 10 us Anton trajectories (7 x 100 000 frames)
        idx -= 4790000
        traj = idx // 100000
        frame = idx % 100000
        return (traj + 1179, frame)
    else:
        # last group of Anton trajectories (109 x 10 000 frames)
        idx -= 5490000
        traj = idx // 10000
        frame = idx % 10000
        return (traj + 1186, frame)


def svd_whiten(X, n, frac_retain=1):
    # From John Strahan - performs whitening of matrix by SVD
    U, L, V = np.linalg.svd(X, full_matrices=False)
    s = np.cumsum(L) / np.sum(L)
    n = np.where(s >= frac_retain)[0][0]
    return (U @ V)[:, 0:n], L


def orthogonalize(basis, pi):
    """Orthgonalize basis functions with respect to a sampling measure.

    Parameters
    ---------
    basis : array-like of ndarray
    pi : array-like of ndarray
    """
    numer = 0
    denom = 0
    for b, p in zip(basis, pi):
        numer += np.einsum("m,mi,mj->ij", p, b, b)
        denom += np.sum(p)
    evals, evecs = scipy.linalg.eigh(numer / denom)
    coeffs = evecs / np.sqrt(evals)[None, :]
    return [b @ coeffs for b in basis]


def delay_embed(tlist, n_embed, lag=1):
    # make delay embedded data (stolen from Erik Thiede's repo)
    # tlist is list of trajectories
    embed_traj_list = []
    for i, traj_i in enumerate(tlist):
        N_i = len(traj_i)
        if N_i - (lag * n_embed) <= 0:  # Must be longer than max embedding
            continue
        embed_traj_i = []
        for n in range(n_embed + 1):
            start_ndx = lag * (n_embed - n)
            stop_ndx = N_i - (lag * n)
            embed_traj_i.append(traj_i[start_ndx:stop_ndx])
        embed_traj_i = np.concatenate(embed_traj_i, axis=1)
        embed_traj_list.append(embed_traj_i)
    return embed_traj_list


def lift_function(tlist, n_embed, lag=1):
    # evaluates a function on the delay-embedded space
    lifted_fxn = []
    for i, fxn_i in enumerate(tlist):
        N_i = len(fxn_i)
        if N_i - (lag * n_embed) <= 0:  # Must be longer than max embedding
            continue
        sub_fxn = fxn_i[int(n_embed * lag / 2) : int(N_i - (n_embed * lag / 2))]
        lifted_fxn.append(sub_fxn)
    return lifted_fxn


def get_index(traj):
    # returns the index of the correct Anton trajectory
    # with the removed duplicates
    duplicates = [
        1282,
        1283,
        1284,
        1285,
        1286,
        1288,
        1289,
        1290,
        1187,
        1188,
        1189,
        1190,
        1191,
        1197,
        1198,
        1199,
        1203,
        1205,
        1206,
        1207,
        1211,
        1212,
        1213,
        1214,
        1215,
        1225,
        1226,
        1227,
        1228,
        1231,
        1232,
        1233,
        1236,
        1237,
        1238,
        1242,
        1245,
        1246,
        1252,
        1253,
        1260,
        1261,
        1262,
        1263,
        1266,
        1267,
        1268,
        1269,
        1270,
        1271,
        1272,
        1273,
        1274,
        1275,
        1276,
        1277,
        1278,
        1279,
    ]

    remove = sorted([i - 1000 for i in duplicates], reverse=True)
    all_ids = list(range(0, 295))
    for i in remove:
        del all_ids[i]
    return all_ids[traj]


def anton_frame(idx):
    if idx < 1790000:
        # first group of Anton trajectories (179 x 10 000 frames)
        traj = idx // 10000
        frame = idx % 10000
        return (traj, frame)
    elif 1790000 <= idx < 2490000:
        # 10 us Anton trajectories (7 x 100 000 frames)
        idx -= 1790000
        traj = idx // 100000
        frame = idx % 100000
        return (traj + 179, frame)
    else:
        # last group of Anton trajectories (109 x 10 000 frames)
        idx -= 2490000
        traj = get_index(idx // 10000 + 186)
        frame = idx % 10000
        return (traj, frame)


def find_closest_points(point, data, n=1):
    """Find ID of closest point(s) in data to point
    point must be of the same dimensionality as data, i.e.
    point.shape = d, data.shape = (N, d) where N in the number
    of data points

    Parameters
    ----------
    point : array-like of length d
        point to which to find closest values
    data : array-like of shape (N, d)
        data to search
    n : int, optional
        number of indices to return

    Returns
    -------
    indices : array-like of length n
    """
    assert len(point) == data.shape[1]
    distances = np.sum((data - point) ** 2, axis=1)
    return np.argsort(distances)[:n]


def make_sparse_basis(dtrajs):
    """Converts a discretized trajectory (e.g. from k-means clustering)
    into a sparse basis of indicator functions.

    Parameters
    ----------
    dtrajs : ndarray
        discretized trajectories

    Return
    ------
    basis : scipy.sparse.csr_matrix
    """
    nclusters = len(np.unique(dtrajs))
    rows, cols = [], []
    for i in range(nclusters):
        pts = np.argwhere(dtrajs == i)
        # indices of which frames are in the cluster i
        rows.append(pts.squeeze())
        # all assigned as 1 in the basis
        cols.append(np.repeat(i, len(pts)))
    rows = np.hstack(rows)
    cols = np.hstack(cols)
    data = np.ones(len(rows), dtype=float)
    basis = sparse.csr_matrix((data, (rows, cols)), shape=(len(dtrajs), nclusters))
    return basis


def make_dense_basis(dtrajs):
    """Converts a discretized trajectory (e.g. from k-means clustering)
    into a sparse basis of indicator functions.

    Parameters
    ----------
    dtrajs : np.ndarray
        discretized trajectories

    Return
    ------
    basis : np.ndarray
    """

    n_basis = len(np.unique(dtrajs))
    basis = np.zeros((len(dtrajs), n_basis))
    basis[np.arange(len(dtrajs)), dtrajs] += 1.0
    return basis


def split_indices(arrays):
    """Gets the indices for np.split from a
    list of arrays.

    Parameters
    ----------
    arrays : ndarray or list/tuple of ndarray
        Arrays from which to get indices

    Returns
    -------
    traj_inds : list of int
        Frame separators to use in np.split
    """
    traj_lens = [len(traj) for traj in arrays]
    traj_inds = []
    subtot = 0
    for length in traj_lens[:-1]:
        subtot += length
        traj_inds.append(subtot)
    return traj_inds


def kdesum2d(
    x,
    y,
    w,
    *,
    xmin=None,
    xmax=None,
    ymin=None,
    ymax=None,
    xstd=None,
    ystd=None,
    nx=100,
    ny=100,
    cut=4.0,
):
    """Compute a 2D kernel density estimate.

    This function histograms the data, then uses a Gaussian filter to
    approximate a kernel density estimate with a Gaussian kernel.

    Credit to Chatipat Lorpaiboon for this code.

    Parameters
    ----------
    x, y : ndarray or list/tuple of ndarray
        Coordinates of each frame.
    w : ndarray or list/tuple of ndarray
        Weight or value of each frame. The output is the sum of these
        values in each bin, after smoothing.
    xmin, xmax, ymin, ymax : float, optional
        Limits of kernel density estimate. If None, takes the min/max
        of the data along the coordinate.
    xstd, ystd : float, optional
        Standard deviation of the Gaussian filter. If None, these are
        set to (xmax - xmin) / nx and (ymax - ymin) / ny, respectively.
        Increase this to smooth the results more.
    nx, ny : int, optional
        Number of bins in each dimension. This should be set as high as
        reasonable, since xstd/ystd takes care of the smoothing.
    cut : float, optional
        Number of standard deviations at which to truncate the Gaussian
        filter. The default, 4, usually doesn't need to be changed.

    Returns
    -------
    kde : (nx, ny) ndarray
        Kernel density estimate, given as bins.
    xedges : (nx+1,) ndarray
        Bin edges along the x dimension.
    yedges : (ny+1,) ndarray
        Bin edges along the y dimension.

    """

    # flatten input to 1D arrays
    x = _flatten(x)
    y = _flatten(y)
    w = _flatten(w)

    # limits
    _xmin = np.min(x)
    _xmax = np.max(x)
    _ymin = np.min(y)
    _ymax = np.max(y)
    if xmin is None:
        xmin = _xmin
    if xmax is None:
        xmax = _xmax
    if ymin is None:
        ymin = _ymin
    if ymax is None:
        ymax = _ymax

    # separation between grid points
    xsep = (xmax - xmin) / nx
    ysep = (ymax - ymin) / ny

    # number of grid points to pad the boundaries,
    # since the Gaussian filter extends beyond the boundaries
    # usually overestimates the padding, but whatever
    ax = max(0, int(np.ceil((xmin - _xmin) / xsep + 1e-6)))
    bx = max(0, int(np.ceil((_xmax - xmax) / xsep + 1e-6)))
    ay = max(0, int(np.ceil((ymin - _ymin) / ysep + 1e-6)))
    by = max(0, int(np.ceil((_ymax - ymax) / ysep + 1e-6)))

    # output bin edges
    xedges = np.linspace(xmin, xmax, nx + 1)
    yedges = np.linspace(ymin, ymax, ny + 1)

    # bin edges, with the added padding
    xedges_padded = np.concatenate(
        [
            xmin + xsep * np.arange(-ax, 0),
            xedges,
            xmax + xsep * np.arange(1, bx + 1),
        ]
    )
    yedges_padded = np.concatenate(
        [
            ymin + ysep * np.arange(-ay, 0),
            yedges,
            ymax + ysep * np.arange(1, by + 1),
        ]
    )
    assert np.allclose(xedges_padded[1:] - xedges_padded[:-1], xsep)
    assert np.allclose(yedges_padded[1:] - yedges_padded[:-1], ysep)
    assert xedges_padded[0] <= _xmin and _xmax <= xedges_padded[-1]
    assert yedges_padded[0] <= _ymin and _ymax <= yedges_padded[-1]

    # construct 2D histogram on padded edges
    hist_padded, _, _ = np.histogram2d(
        x, y, weights=w, bins=(xedges_padded, yedges_padded)
    )

    # Gaussian kernel parameters
    if xstd is None:
        xstd = xsep
    if ystd is None:
        ystd = ysep

    # apply Gaussian filter to histogram
    kde_padded = scipy.ndimage.gaussian_filter(
        hist_padded,
        sigma=(xstd / xsep, ystd / ysep),  # in units of grid points
        mode="constant",
        truncate=cut,
    )

    # remove the padding
    assert ax + nx + bx == kde_padded.shape[0]
    assert ay + ny + by == kde_padded.shape[1]
    kde = kde_padded[ax : ax + nx, ay : ay + ny]

    return kde, xedges, yedges


def _flatten(a):
    if isinstance(a, np.ndarray):
        # avoid creating a new array (and using twice the memory)
        return np.ravel(a)
    else:
        return np.ravel(np.concatenate(a))


def moving_average(x, w):
    """Computes a moving average for an array.

    Parameters
    ----------
    x : np.ndarray, dimension (n_frames,)
    w : int
        Window size
    """
    return signal.convolve(x, np.ones(w) / w, mode="same")


def smooth_moving_average(trajs, w, n=2):
    """Smooths trajectory by computing repeated moving averages
    (i.e. multiple convolutions with a constant array).

    Parameters
    ----------
    trajs : np.ndarray or array-like of np.ndarray
        Trajectory or list of trajectories
    w : int
        Window size
    n : int, optional
        Number of times to perform convolution, which controls
        the amount of smoothing.

    Returns
    -------
    ans : np.ndarray or array-like of np.ndarray
        Smoothed trajectory(s)
    """
    if isinstance(trajs, np.ndarray):
        ans = trajs
        for i in range(n):
            ans = moving_average(ans, w)
        return ans
    else:
        ans = []
        for arr in trajs:
            assert arr.ndim == 1
            for i in range(n):
                arr = moving_average(arr, w)
            ans.append(arr)
        return ans


def zero_h_occupancy(pdb_file, new_file_name=None, savelines=3):
    import Bio.PDB
    from Bio.PDB import PDBParser, StructureBuilder

    if new_file_name is None:
        new_file_name = pdb_file
    parser = PDBParser()
    pdb = parser.get_structure("", pdb_file)
    for atom in pdb.get_atoms():
        if atom.element == "H" or atom.element == "HE":
            atom.set_bfactor(0.0)
        else:
            atom.set_bfactor(1.0)
    io = Bio.PDB.PDBIO()
    io.set_structure(pdb)
    with open(new_file_name, mode="w+") as f:
        f.writelines(open(pdb_file).readlines()[:savelines])
        io.save(f)


def sb_label(i):
    sb_labels = []
    for r in ("R217", "R223", "R226", "R229", "R232"):
        for n in ("D129", "D136", "D151", "D164", "E183", "D186"):
            sb_labels.append(f"{r} - {n}")
    return sb_labels[i - 30]
