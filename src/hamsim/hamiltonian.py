# hamiltonian
from numpy import array, pi
from hamsim.operators import kron_at


def ising(J=pi, h=2*pi):
    """A transverse-field Ising Hamiltonian for 2 qubits with
    interaction constant J and degenerate energy h

    """
    n = 2
    X = array([[0, 1], [1, 0]])
    Z = array([[1, 0], [0, -1]])
    Xn = [kron_at(X, i, n) for i in range(n)]
    Zn = [kron_at(Z, i, n) for i in range(n)]

    H = J*(Zn[0] @ Zn[1]) + h*(Xn[0] + Xn[1])

    return H
