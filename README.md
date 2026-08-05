# Hamiltonian Simulation

Exact and Trotter-Suzuki simulation of the transverse-field Ising
model, with error and runtime analysis benchmarked against theory.

## Objective

The aim is to build toward simulating open quantum systems. As a
concrete first step, this repository implements and benchmarks
Hamiltonian simulation using Trotterization. The local interaction
implementation from Trotter product formula scales far better than the
dense matrix exponential simulation. The result highlights the
connection to the polynomial-time advantage Trotterization provides on
quantum hardware.

## Result

**Trotter error converges at the expected first-order rate, and
   within the theoretical commutator bound.**

For a first-order product formula, Lloyd (1996)	bounds the error by

```math
\vert\vert U - U_{\rm trott} \vert\vert \leq (t^2 / 2n)
\sum_{i \lt j} \vert\vert[H_i, H_j ] \vert\vert
```

where $U = e^{-iHt}$ is the exact evolution operator, and $U_{\rm
trott} = (e^{-iH_1t/n} \cdots e^{-iH_lt/n})^n$ is the $n$-step Trotter
product. Each local operator $e^{-iH_it/n}$ acts on a Hilbert space of
only $m=4$ dimensions (of nearest neighbors), rather than the full
$2^N$ dimensional Hilbert space in $H$, of the whole system of N
qubits. This locality is what the efficient implementation exploits
directly, via tensor contraction rather than dense matrix
exponentiation.

The measured error tracks the bound's scaling closely; both
curves are parallel on a log-log plot, confirming the correct
first-order convergence rate, and the true error sits consistently
below the bound across the full range of step counts tested.

<!-- ![Trotter error vs. theoretical
bound](docs/figures/error_scaling.png) -->
<div align="center">
  <img src="docs/figures/error_scaling.png" width="500">
</div>

**The efficient (tensor-contraction) Trotter implementation is
dramatically faster.**

At $N=12$, exact matrix exponentiation of the systems evolution takes
$\approx 27 {\rm s}$ while the evolution via Trotterization takes
$\approx 5 {\rm ms}$, roughly a 5400x speedup. The exact method
grows approximately as $e^{(1.53)N}$, which is somewhat less than the
expected cubic scaling $e^{ln(8)N}$, with cost $\mathcal{O}(8^N)$, for
dense matrix exponentiation. (The discrepancy plausibly reflects
optimized linear algebra routines in scipy's module: `expm`, rather
than a departure from the theoretical asymptotic.) The implementation
of Trotter product formula costs $\mathcal{O}(N \cdot 2^N)$ per
trotter step, since a classical computer must store the full
$\mathcal{O}(2^N)$ dimensional state. Both algorithms are exponential
in $N$, but trotterization has a much smaller effective exponent and
thus yields the dramatic speedup observed in the Figure on runtime
scaling cost.

<div align="center">
  <img src="docs/figures/runtime_scaling.png" width="500">
</div>
<!-- ![Runtime scaling](docs/figures/runtime_scaling.png) -->

The genuine polynomial-time advantage Trotterization is known for is
specific to real quantum hardware, where a single local gate costs
$\mathcal{O}(1)$, independent of $N$, and a full Trotter step needs
only $\mathcal{O}(N)$ such gates, polynomial overall.

## Setup
pip install -e .
pip install jupyter matplotlib ipykernel

## Repository structure

```
hamiltonian-simulation/
├── README.md
├── pyproject.toml
├── docs/
│   └── figures/
├── notebooks/
│   ├── exact_time_evolution.ipynb
│   ├── first_order_trotter.ipynb
│   └── error_scaling.ipynb
├── src/
│   └── hamsim/
│       ├── operators.py
│       ├── hamiltonian.py
│       └── evolution.py
└── tests/
    ├── test_operators.py
    ├── test_hamiltonian.py
    └── test_evolution.py
```

## Notebooks
- `exact_time_evolution.ipynb` — exact diagonalization and time evolution
- `first_order_trotter.ipynb` — first-order Trotter (dense vs. exact)
- `error_scaling.ipynb` — error scaling and runtime benchmark
(Second-order Suzuki: follow-up)