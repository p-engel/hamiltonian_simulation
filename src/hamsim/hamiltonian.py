# hamiltonian
from numpy import array, pi, zeros
from hamsim.operators import kron_at

def X(): return array([[0, 1], [1, 0]])

def Z(): return array([[1, 0], [0, -1]])

def ops(A, n):
    """define set of operators of at site n in hilbert space 2**n by
    2**n

    """
    return [kron_at(A, i, n) for i in range(n)]
    
def ising(J=pi, h=2*pi, n=2):
    """A transverse-field Ising Hamiltonian for n qubits with
    interaction constant J and degenerate energy h

    """
    Zn = ops(Z(), n); Xn = ops(X(), n)
    H = h * sum(Xn)
    
    for i in range(n - 1):
        H += J * ( Zn[i] @ Zn[i+1] )

    return H

# def ising_trotter(J=pi, h=2*pi, n=10):
#    """An n-qubit transverse field Ising Hamiltonian expressed in terms
#    of set of local interactions in subspace m by m where m = 2
# 
#    """
#    Zn = ops(Z(), n); Xn = ops(X(), n)
# 
#    l = n - 1  # number of local interactions or gates
#    A = [ J * Zn[i] @ Zn[i+1] for i in range(l) ]
#    B = [ h * Xn[i]  for i in range(n) ]
#    H = []  # set of local interactions
#    for i in range(l):
#        Hi = A[i] + (B[i] + B[i+1])/2 
#        H.append(Hi)
# 
#    # correction terms to edge qubits
#    H[0] += (h/2) * Xn[0]; H[-1] += (h/2) * Xn[n-1]
# 
#    return H

def ising_trotter_local(
        J, h, is_left_end=False, is_right_end=False
):
    """Local interaction Hamiltonian in subspace (m by m) of n-qubit
    transverse field Ising model

    """
    m = 2
    Xn = ops(X(), m); Zn = ops(Z(), m)
    Hl = J * (Zn[0] @ Zn[1]) + h * (Xn[0] + Xn[1]) / 2

    if is_left_end: Hl += (h/2) * Xn[0]
    if is_right_end: Hl += (h/2) * Xn[-1]
    
    return Hl

def ising_trotter(n, J=pi, h=2*pi):
    l = n - 1
    H = []
    for i in range(l):
        H.append(
            ising_trotter_local(
                J, h, is_left_end=(i==0), is_right_end=(i==l-1)
            )
        )

    return H
