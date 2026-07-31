# test hamiltonian
from numpy import allclose
from hamsim import hamiltonian


def test_hermiticity():
    H = hamiltonian.ising(n=10)
    Hk = hamiltonian.ising_trotter(n=10)
    Hk_dagger = [Hk[i].conj().T for i in range(len(Hk))]
    H_dagger = H.conj().T
    assert allclose(H, H_dagger), (
        f"the Hamiltonian is not Hermitian"
    )
    assert allclose(Hk, Hk_dagger), (
        f"the subspace Hamiltonian H_k is not Hermitian"
    )

    return

def test_sum_local_H():
    H = hamiltonian.ising(n=10)
    sum_Hk = sum(hamiltonian.ising_trotter(n=10))

    assert allclose(H, sum_Hk), ( f"the summation of set of local "
    f"Hamiltonian is not equal to the the full Hamiltonian"
    )
