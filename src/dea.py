"""dea.py -- DEA library: BoD, AR-BoD, super-efficiency, GR-DEA (QP), common weights, AUC.

Part of SpatialDEA-MPM (Nazarpour, submitted to Mathematical Geosciences).
Run from the repository root (see README.md). Set SDEA_BASE to override the base directory.
"""
"""Spatial DEA models for catchment-based mineral prospectivity.

Models
------
bod        : Benefit-of-the-Doubt DEA (CRS, unit input)  -- per-DMU LP
bod_ar     : BoD with virtual-weight (assurance-region) restrictions
super_eff  : Andersen-Petersen super-efficiency version of bod / bod_ar
grdea      : Geographically Regularised DEA (joint QP over all DMUs)
"""
import numpy as np
from scipy.optimize import linprog
import cvxpy as cp


def normalise(Y, eps=1e-3):
    Y = np.asarray(Y, float)
    Yn = Y / Y.max(0, keepdims=True)
    return np.maximum(Yn, eps)


def _lp(Y, o, alpha=0.0, exclude_self=False):
    n, s = Y.shape
    yo = Y[o]
    rows = np.arange(n) != o if exclude_self else np.ones(n, bool)
    A = [Y[rows]]
    b = [np.ones(rows.sum())]
    if alpha > 0:
        # alpha * (u.yo) - u_r yo_r <= 0   for each r
        Aar = alpha * np.tile(yo, (s, 1)) - np.diag(yo)
        A.append(Aar); b.append(np.zeros(s))
    res = linprog(-yo, A_ub=np.vstack(A), b_ub=np.concatenate(b),
                  bounds=[(0, None)] * s, method="highs")
    if res.status != 0:
        return np.nan, np.full(s, np.nan)
    return -res.fun, res.x


def bod(Y, alpha=0.0, super_eff=False):
    n, s = Y.shape
    th = np.empty(n); U = np.empty((n, s))
    for o in range(n):
        th[o], U[o] = _lp(Y, o, alpha, exclude_self=super_eff)
    return th, U


def knn_weights(xy, k=6):
    from scipy.spatial import cKDTree
    t = cKDTree(xy)
    d, idx = t.query(xy, k + 1)
    b = np.median(d[:, 1:])
    n = len(xy)
    W = np.zeros((n, n))
    for i in range(n):
        for dd, j in zip(d[i, 1:], idx[i, 1:]):
            w = np.exp(-(dd / b) ** 2)
            W[i, j] = max(W[i, j], w); W[j, i] = max(W[j, i], w)
    W *= n / W.sum()          # total edge mass = n  (makes lambda dimensionless)
    return W, b


def grdea(Y, W, lam, alpha=0.0, solver="CLARABEL", ridge=1e-6):
    """max sum_j u_j.y_j - lam * sum_{j<k} w_jk ||D(u_j-u_k)||^2
       s.t. u_j.y_l <= 1 for all j,l ; u_j >= 0 ; virtual-weight shares >= alpha."""
    n, s = Y.shape
    D = Y.mean(0)                                # scale weights to 'virtual' units
    U = cp.Variable((n, s), nonneg=True)
    V = cp.multiply(U, np.tile(D, (n, 1)))
    iu, ju = np.nonzero(np.triu(W, 1))
    wv = W[iu, ju]
    diff = V[iu, :] - V[ju, :]
    energy = cp.sum(cp.multiply(wv[:, None], cp.square(diff)))
    theta = cp.sum(cp.multiply(U, Y), axis=1)
    obj = cp.Maximize(cp.sum(theta) - lam * energy - ridge * cp.sum_squares(V))
    cons = [U @ Y.T <= 1]
    if alpha > 0:
        cons.append(cp.multiply(U, Y) >= alpha * cp.reshape(theta, (n, 1), order="F") @ np.ones((1, s)))
    prob = cp.Problem(obj, cons)
    prob.solve(solver=solver)
    Uv = U.value
    th = (Uv * Y).sum(1)
    Vv = Uv * D
    E = float((wv * ((Vv[iu] - Vv[ju]) ** 2).sum(1)).sum())
    return th, Uv, E, prob.status


def common_weights(Y, alpha=0.0):
    """lambda -> infinity limit: one weight vector maximising total score."""
    n, s = Y.shape
    c = -Y.sum(0)
    A = [Y]; b = [np.ones(n)]
    res = linprog(c, A_ub=np.vstack(A), b_ub=np.concatenate(b), bounds=[(0, None)] * s, method="highs")
    u = res.x
    return Y @ u, u


def auc(score, pos, neg=None):
    score = np.asarray(score); pos = np.asarray(pos, bool)
    neg = ~pos if neg is None else np.asarray(neg, bool)
    sp, sn = score[pos], score[neg]
    gt = (sp[:, None] > sn[None, :]).sum() + 0.5 * (sp[:, None] == sn[None, :]).sum()
    return gt / (len(sp) * len(sn))


def auc_perm_p(score, pos, neg=None, nperm=20000, seed=0):
    rng = np.random.default_rng(seed)
    pos = np.asarray(pos, bool)
    neg = ~pos if neg is None else np.asarray(neg, bool)
    pool = np.where(pos | neg)[0]; npos = pos.sum()
    a0 = auc(score, pos, neg)
    cnt = 0
    for _ in range(nperm):
        sel = rng.choice(pool, npos, replace=False)
        p = np.zeros_like(pos); p[sel] = True
        q = np.zeros_like(pos); q[pool] = True; q[sel] = False
        cnt += auc(score, p, q) >= a0
    return a0, (cnt + 1) / (nperm + 1)
