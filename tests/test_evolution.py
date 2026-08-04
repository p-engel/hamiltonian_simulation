# test_evolution
from numpy import eye, allclose, linalg, zeros
from hamsim import hamiltonian, evolution

def test_trotter_error_bound():
    n = 4; t = 0.5; n_steps = 200
    H = hamiltonian.ising(n=n)
    Hl = hamiltonian.ising_trotter(n=n)
    psi0 = zeros(2**n, dtype=complex); psi0[1] = 1

    assert (
        evolution.trotter_true_error(H, Hl, psi0, t, n_steps)
            < evolution.trotter_error_bound(Hl, t, n_steps)
    )

def test_trotter_error_scales_linearly_with_dt():
    n = 4; t = 0.5
    H = hamiltonian.ising(n=n)
    Hl = hamiltonian.ising_trotter(n=n)
    psi0 = zeros(2**n, dtype=complex); psi0[1] = 1

    errors = []
    step_counts = [10, 20, 40, 80]
    for r in step_counts:
        errors.append(
            evolution.trotter_true_error(H, Hl, psi0, t, r)
        )

    # error is proportional to 1/r
    ratios = [errors[i]/errors[i+1] for i in range(len(errors)-1)]
    assert all(1.5 < ratio < 2.5 for ratio in ratios)

# def test_trotter_step_is_unitary():
#     H = hamiltonian.ising_trotter(n=4)
#     U = evolution.trotter_step(H, dt=0.01)
#     assert allclose( U @ U.conj().T, eye(U.shape[0]) )
