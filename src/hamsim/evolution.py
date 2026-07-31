# evolution
from scipy import linalg
from numpy import allclose, ndarray, zeros, log2, eye

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
        
def trotter_step(H, dt):
    """Single first-order trotter step: product of exp(-i*Hi*dt) for
    each local term Hi in set H, applied in sequence
    
    """
    dim = H[0].shape[0]
    U = eye(dim, dtype=complex)
    for Hi in H:
        U = linalg.expm(-1j * Hi * dt) @ U

    return U

def trotter_evolve(H, psi0, t, n_steps):
    """Evolve psi0 for local time t using n_steps Trotter steps."""
    dt = t / n_steps
    U_step = trotter_step(H, dt)
    psi = psi0.copy()
    for _ in range(n_steps): psi = U_step @ psi
        
    return psi

def trotter_true_error(H, Hk, psi0, t, r):
    psi_exact = exact_evolve(H, psi0, t)
    psi_trotter = trotter_evolve(Hk, psi0, t, r)

    return linalg.norm(psi_exact - psi_trotter)

def trotter_error_bound(H, t, r):
    commutator_sum = 0
    for i in range(len(H)):
        for j in range(i+1, len(H)):
            comm = H[i] @ H[j] - H[j] @ H[i]
            commutator_sum += linalg.norm(comm, ord=2)

    return (t**2 / (2*r)) * commutator_sum
