# tests/test_operators.py
import numpy as np
from hamsim.operators import kron_at

X = np.array([[0, 1], [1, 0]])
I = np.eye(2)

def test_identity_at_any_site_is_identity_matrix():
    assert np.allclose(kron_at(I, 0, 2), np.eye(4))
    assert np.allclose(kron_at(I, 1, 2), np.eye(4))

def test_X_squares_to_identity():
    Xi = kron_at(X, 0, 2)
    assert np.allclose(Xi @ Xi, np.eye(4))

def test_X_on_qubit0_flips_qubit0_only():
    # |11> (index 3) -> |01> (index 1)
    M = kron_at(X, 0, 2)
    e3 = np.zeros(4); e3[3] = 1
    result = M @ e3
    expected = np.zeros(4); expected[1] = 1
    assert np.allclose(result, expected)

def test_commutation():
    # X on qubit 0 and X on qubit 1 should act independently
    X0 = kron_at(X, 0, 2)
    X1 = kron_at(X, 1, 2)
    assert np.allclose(X0 @ X1, X1 @ X0)
