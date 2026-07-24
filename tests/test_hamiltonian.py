# test hamiltonian
from numpy import allclose
from hamsim import hamiltonian


def test_hermiticity():
    H = hamiltonian.ising()
    H_dagger = H.conj().T
    assert allclose(H, H_dagger), (
        f"the Hamiltonian is not Hermitian"
    )

    return
