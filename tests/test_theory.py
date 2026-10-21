"""Numerical checks of the propositions and theorem in the manuscript (runs in < 1 min).

    python3 tests/test_theory.py
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src'))
import numpy as np
from dea import normalise, bod, grdea, knn_weights, common_weights

rng = np.random.default_rng(42)
n, s = 40, 5
Y = normalise(rng.lognormal(0, 1, (n, s)))
xy = rng.uniform(0, 10000, (n, 2)); W, _ = knn_weights(xy, 5)
tol = 1e-6

# Proposition 2(a,b,c): feasibility alpha <= 1/s, weakest-link bound, monotone in alpha
prev = None
for a in [0.0, 0.05, 0.1, 0.2]:
    th, U = bod(Y, alpha=a)
    if a > 0:
        assert np.all(th <= Y.min(1) / a + tol), 'weakest-link bound violated'
    if prev is not None:
        assert np.all(th <= prev + tol), 'score not non-increasing in alpha'
    prev = th
print('Proposition 2: feasibility, weakest-link bound and monotonicity in alpha  OK')

# Proposition 3(a): invariance to column rescaling;  3(b): multiplicative stability bound
th0, _ = bod(Y)
th1, _ = bod(Y * np.array([3.0, 0.2, 1.0, 7.0, 0.5]))
assert np.allclose(th0, th1, atol=1e-7), 'units invariance violated'
d = 0.1
Yp = Y * rng.uniform(1 - d, 1 + d, Y.shape)
thp, _ = bod(Yp)
r = thp / th0
assert np.all(r >= (1 - d) / (1 + d) - tol) and np.all(r <= (1 + d) / (1 - d) + tol), 'stability bound violated'
print('Proposition 3: units invariance and multiplicative stability bound  OK')

# Theorem 1(ii): lambda = 0 decouples into AR-BoD
th_ar, _ = bod(Y, alpha=0.05)
th_g0, _, _, _ = grdea(Y, W, 0.0, alpha=0.05)
assert np.allclose(th_g0, th_ar, atol=1e-3), 'decoupling at lambda=0 violated'
# Theorem 1(iii): total efficiency and Dirichlet energy non-increasing in lambda
S_prev, E_prev = np.inf, np.inf
for lam in [0.0, 0.1, 1, 10, 100, 1000]:
    th, U, E, st = grdea(Y, W, lam, alpha=0.05)
    assert st == 'optimal'
    assert th.sum() <= S_prev + 1e-4 and E <= E_prev + 1e-6, 'monotonicity in lambda violated'
    S_prev, E_prev = th.sum(), E
# Theorem 1(iv): large lambda -> (near) common weights
th_big, U_big, E_big, _ = grdea(Y, W, 1e5, alpha=0.0)
assert E_big < 1e-4, 'weights not equalised at large lambda'
print('Theorem 1: decoupling, monotonicity in lambda and common-weight limit  OK')
print('All checks passed.')
