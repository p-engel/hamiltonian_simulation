# operators
from numpy import kron, eye
from functools import reduce

def kron_at(op, site, n_qubits):
    """Place `op` (a 2x2 array) at position `site` in an n_qubits-qubit
    tensor product, identity (2x2) everywhere else. Returns a
    (2**n_qubits, 2**n_qubits) array.

    """
    matrices = []
    for i in range(n_qubits):
        matrices.append(eye(2))

    matrices[site] = op

    return reduce(kron, matrices)
