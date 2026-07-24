# evolution
from scipy import linalg
from numpy import allclose, ndarray, zeros

def exact_evolve(H, psi0, t):
    """unitary evolution of an initial state psi0 and Hamiltonian H to a
    state psi_t at time t

    """
    U = linalg.expm(-1j * H * t)
    spectral_norm = linalg.norm(U, ord=2)
    assert allclose(1.0, spectral_norm), (
        f"Spectral Norm: {spectral_norm:.4f} should be 1"
    )

    return U @ psi0

def exact_evolve_series(H, psi0, times: ndarray):
    psi_t = zeros((len(times), len(psi0)), dtype=complex)
    for i, t in enumerate(times):
        psi_t[i, :] = exact_evolve(H, psi0, t)
    return psi_t
        
