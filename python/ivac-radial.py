import numpy as np
import extq
import ivac
import pyemma
import glob
import sys
import matplotlib.pyplot as plt
from matplotlib import ticker

sys.path.insert(1, "./")
sys.path.insert(1, "../../")
import util
import kde
import seaborn as sns

plt.style.use("seaborn-ticks")
sns.set_palette("colorblind")

# define radial basis functions
def gauss_short(r, r0, d0):
    return np.where(r < d0, 1, np.exp(-((r - d0) ** 2) / (2 * r0 ** 2)))


def gauss_long(r, r0, d0):
    return 1 - gauss_short(r, r0, d0)


def transform_data(data_sb_arr, r0, d0):
    data_short = gauss_short(data_sb_arr, r0, d0)
    data_long = gauss_long(data_sb_arr, r0, d0)
    return (data_short, data_long)


def plot_its(its):
    f, ax = plt.subplots(figsize=(8, 6))
    plt.plot(its[:10] * 0.0001, ".")  # in units of us
    plt.ylabel("Implied timescales / $\mu s$")
    plt.xlabel("Eigenvalue index")
    plt.ylim([0, 10])
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    plt.savefig("../fig/ivac-radial-its.png", dpi=300)


def plot_evals(evals):
    f, ax = plt.subplots(figsize=(8, 6))
    plt.plot(evals[:10], ".")
    plt.ylabel("Eigenvalue")
    plt.xlabel("Eigenvalue index")
    plt.ylim([1000, 3100])
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    plt.savefig("../fig/ivac-radial-evals.png", dpi=300)


def plot_evecs(evecs):
    f, axs = plt.subplots(1, 4, figsize=(20, 4))
    for i, ax in enumerate(axs):
        ax.plot(evecs[i + 1], ".-")
        ax.set_xlim([0, 63])
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax.set_title(f"iTIC {i + 1}")
    plt.savefig("../fig/ivac-radial-evecs.png", dpi=300)


def plot_proj(livac_arr, models_livac):
    f, axs = plt.subplots(1, 3, figsize=(18, 6))
    for i, ax in enumerate(axs):
        h = ax.scatter(livac_arr[::10, 1], livac_arr[::10, i + 2], alpha=0.1)
        ax.set_xlabel("iTIC 1")
        ax.set_ylabel(f"iTIC {i + 2}")

    util.plot_models(axs[0], models_livac[:, 1:3])
    util.plot_models(axs[1], models_livac[:, 1:4:2])
    util.plot_models(axs[2], models_livac[:, 1:5:3])
    axs[0].legend(["Down-", "Down", "Up", "Up+"])
    f.tight_layout()
    plt.savefig("../fig/ivac-radial-dataproj.png", dpi=300)


def plot_cv_proj(cv_arr, livac_arr):
    f, axs = plt.subplots(1, 3, figsize=(21, 6))
    for i, ax in enumerate(axs):
        h = ax.scatter(
            cv_arr[::20, 0],
            cv_arr[::20, 1],
            c=livac_arr[::20, i + 1],
            cmap="rocket",
            alpha=0.2,
        )
        cbar = plt.colorbar(h, ax=ax)
        cbar.ax.set_ylabel(f"iTIC {i+1}")
        util.format_cvs(ax)
    axs[0].legend(["Down-", "Down", "Up", "Up+"])
    f.tight_layout()
    plt.savefig("../fig/ivac-radial-ds-iTIC-proj.png", dpi=300)


def main():
    # load data
    with np.load("../data/raw_feat/feat2_raw.npz", allow_pickle=True) as f:
        data = f["arr_0"]
    data_arr = np.concatenate(data)
    with np.load("../data/raw_feat/cv_dist_spin.npz", allow_pickle=True) as f:
        cv_arr = f["arr_0"]
    traj_lens = [len(traj) for traj in data]
    traj_inds = []
    subtot = 0
    for length in traj_lens[:-1]:
        subtot += length
        traj_inds.append(subtot)
    cv_trajs = np.split(cv_arr, traj_inds)
    data_cz = [traj[:, 30:] for traj in data]
    feat2_models = np.load("../data/models_feat2.npy")
    models_cz = feat2_models[:, 30:]

    sb_ids = [36, 42, 48, 46, 52, 58]
    data_sb = [traj[:, sb_ids] for traj in data]
    data_sb_arr = np.concatenate(data_sb)

    r0 = float(sys.argv[1])
    d0 = float(sys.argv[2])
    data_short, data_long = transform_data(data_sb_arr, r0, d0)
    data_basis = [np.ones(data_arr.shape[0])]
    for i in range(len(sb_ids)):
        new_basis = []
        for temp in data_basis:
            new_basis.append(np.multiply(temp, data_short[:, i]))
            new_basis.append(np.multiply(temp, data_long[:, i]))
            data_basis = new_basis

    data_basis = np.array(data_basis).T
    data_basis_trajs = np.split(data_basis, traj_inds)

    # IVAC calculations
    livac = ivac.LinearIVAC(
        minlag=1, maxlag=3000, nevecs=5, reweight=False, adjust=True, method="fft"
    )
    livac.fit(data_basis_trajs)

    # plot implied timescales
    plot_its(livac.its)
    # plot eigenvalues
    plot_evals(livac.evals)
    # plot eigenvectors
    plot_evecs(livac.evecs)

    livac_trajs = livac.transform(data_basis_trajs)
    livac_arr = np.concatenate(livac_trajs)

    # transform models
    models_sb = feat2_models[:, sb_ids]
    models_short = gauss_short(models_sb, r0, d0)
    models_long = gauss_long(models_sb, r0, d0)
    models_basis = [np.ones(4)]
    for i in range(len(sb_ids)):
        new_basis = []
        for temp in models_basis:
            new_basis.append(np.multiply(temp, models_short[:, i]))
            new_basis.append(np.multiply(temp, models_long[:, i]))
        models_basis = new_basis
    models_basis = np.array(models_basis).T
    models_livac = livac.transform([models_basis])[0]

    # plot data projections
    plot_proj(livac_arr, models_livac)

    # plot projection onto CVs
    plot_cv_proj(cv_arr, livac_arr)


if __name__ == "__main__":
    main()
