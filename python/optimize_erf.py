import sys
import os
import numpy as np
import scipy
from scipy import optimize

from jax.config import config
import jax
import jax.numpy as jnp
from jax import jacfwd, jacrev, jit

import util

config.update("jax_enable_x64", True)

DATA_DIR = "/project/dinner/scguo/ci-vsd/data"


def newton_method(fn, jac_fn, hess_fn, x0, n_its, args, kwargs):
    x = x0
    for i in range(n_its):
        x = x - jnp.linalg.inv(hess_fn(*args, **kwargs)) @ jac_fn(*args, **kwargs)


def jac_hess(fn):
    jac_fn = jacrev(fn)
    hess_fn = jacfwd(jac_fn)
    return jac_fn, hess_fn


def corr_q(e, z_traj, z0, lag=500, sigma=1):
    """Compute committor correlation function
    using erf ansatz

    Returns
    -------
    C : correlation function (negative for minimization)
    """
    C = 0
    for z in z_traj:
        q_0 = q(e, z[:-lag], z0, sigma=sigma)
        q_tau = q(e, z[lag:], z0, sigma=sigma)
        C += jnp.mean((q_tau - q_0) ** 2)
    return -0.5 * C  # want to maximize, so flip sign


def Lagrangian(x, z_traj, z0, lag=500, sigma=1):
    e = x[:-1]
    L = corr_q(e, z_traj, z0, lag=lag, sigma=sigma)
    L -= x[-1] * (jnp.linalg.norm(e) - 1)
    return L


@jit
def q(e, z, z0, sigma=1):
    """Compute the value of the committor for a trajectory z given
    a unit vector e

    z : ndarray of shape (n_frames, n_features)
    z0 : ndarray of shape (n_features,)
    """
    arg = jnp.dot((z - z0), e) / sigma
    return 0.5 * (1 + jax.scipy.special.erf(arg))


def load_data():
    # S4 translocation/rotation data
    cv_trajs = list(
        np.load(f"{DATA_DIR}/raw_feat/cv_dist_spin_anton.npy", allow_pickle=True)
    )
    cv_trajs.extend(np.load(f"{DATA_DIR}/raw_feat/cv_dist_spin_anton2.npy"))
    cv_arr = np.concatenate(cv_trajs)
    # salt bridge distances for states
    sb_trajs = list(
        np.load(f"{DATA_DIR}/raw_feat/feat2_raw_anton.npy", allow_pickle=True)
    )
    sb_trajs.extend(np.load(f"{DATA_DIR}/raw_feat/feat2_raw_anton2.npy"))
    sb_arr = np.concatenate(sb_trajs)
    # sb_models = np.load(f"{DATA_DIR}/models_centroids_feat2.npy")
    rf161 = list(np.load(f"{DATA_DIR}/raw_feat/rf161.npy", allow_pickle=True))
    rf161.extend(np.load(f"{DATA_DIR}/raw_feat/rf161_anton2.npy"))
    rf161_arr = np.concatenate(rf161)
    # committors
    qp_du = np.load(
        f"{DATA_DIR}/feat2_dist_du_anton2/qp_downup_3.npy", allow_pickle=True
    )[
        8
    ]  # 50 ns lag time
    # weights
    weights = np.load(
        f"{DATA_DIR}/feat2_dist_du_anton2/weights_3_feat5ivac.npy", allow_pickle=True
    )[
        0
    ]  # 0.1 ns lag time
    X = np.hstack((cv_arr, sb_arr[:, 30:], rf161_arr))
    y = np.concatenate(qp_du)
    return X, y, weights


def run_optimizer(e0, z_trajs, z0):
    def cons_f(x):
        return jnp.sum(x * x)

    def cons_J(x):
        return 2 * x

    def cons_H(x, v):
        return jnp.eye(len(x)) * 2 * v

    constraint = optimize.NonlinearConstraint(
        cons_f, 1, 1, jac=cons_J, hess=cons_H
    )

    grad_cq, hess_cq = jac_hess(corr_q)

    jit_grad = jit(grad_cq, static_argnames=("lag", "sigma"))
    jit_hess = jit(hess_cq, static_argnames=("lag", "sigma"))
    solve_ans = optimize.minimize(
        corr_q,
        e0,
        args=(z_trajs, z0),
        method="trust-constr",
        jac=jit_grad,
        hess=jit_hess,
        constraints=constraint,
        options={"verbose": 2},
    )
    return solve_ans


def main():
    X, q_arr, weights = load_data()
    offset = 0.05
    ts_ids = (q_arr > (0.5 - offset)) & (q_arr < (0.5 + offset))
    traj_inds = util.split_indices(weights)
    z_arr = scipy.sparse.diags(np.concatenate(weights)) @ X
    z_trajs = np.split(z_arr, traj_inds, axis=0)
    z0 = np.mean(z_arr[ts_ids], axis=0)
    n_feat = z_arr.shape[-1]
    vec = np.random.normal(size=n_feat)
    e0 = vec / np.linalg.norm(vec)

    print(f"z0: {z0}")
    print(f"Initial guess: {e0}")
    print(f"Data size: {X.shape}")
    ans = run_optimizer(e0, z_trajs, z0)
    print(ans)


if __name__ == "__main__":
    main()
