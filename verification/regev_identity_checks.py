"""Independent exact/numerical checks for the Regev order-polynomial protocol.

The general identity is inherited from Procesi (2021), A note on the Formanek
Weingarten function, Theorem 7, Eq.(2.85), DOI 10.1285/i15900932v41n1p69.
It also appears in arXiv:2102.10657v2, Theorem 4.13.
These checks do not substitute for that theorem.
"""
from collections import defaultdict
from functools import reduce
from itertools import permutations
from math import factorial
from fractions import Fraction
from pathlib import Path
import argparse
import json
import numpy as np


def regev_constant(d):
    C = Fraction(1)
    for k in range(1, d + 1):
        C *= factorial(2*k - 1)
    for k in range(1, d):
        C /= factorial(k)
    return (-1)**(d-1) * C*C / (factorial(d)**2 * (2*d-1))


def basis_regev(d):
    """Sum signed basis-matrix products by sparse exact dynamic programming.

    X and Y each enumerate E_00,E_01,... in row-major order.
    Regev word groups have lengths X1,Y1,X3,Y3,...,X(2d-1),Y(2d-1).
    A matrix-unit product is a directed path, so only compatible edges survive.
    """
    pattern = []
    for k in range(1, d+1):
        pattern += [0]*(2*k-1) + [1]*(2*k-1)
    states = {(0, 0, start, start): 1 for start in range(d)}
    peak = 0
    for group in pattern:
        new = defaultdict(int)
        for (mx, my, start, current), val in states.items():
            mask = mx if group == 0 else my
            for nxt in range(d):
                idx = current*d+nxt
                if mask & (1 << idx):
                    continue
                sign = (-1)**((mask >> (idx+1)).bit_count())
                nmx = mx | (1 << idx) if group == 0 else mx
                nmy = my | (1 << idx) if group == 1 else my
                new[nmx, nmy, start, nxt] += sign*val
        states = {k: v for k, v in new.items() if v}
        peak = max(peak, len(states))
    out = np.zeros((d, d), dtype=object)
    for (_, _, a, b), v in states.items():
        out[a,b] += v
    return out, peak


def parity(p):
    return (-1)**sum(p[i] > p[j] for i in range(len(p)) for j in range(i+1,len(p)))


def regev2(X, Y):
    out = np.zeros((2,2), dtype=complex)
    for px in permutations(range(4)):
        for py in permutations(range(4)):
            word = [X[px[0]],Y[py[0]],X[px[1]],X[px[2]],X[px[3]],Y[py[1]],Y[py[2]],Y[py[3]]]
            out += parity(px)*parity(py)*reduce(np.matmul, word)
    return out


def verify():
    exact=[]
    for d in (2,3):
        out, peak = basis_regev(d)
        c = regev_constant(d)
        assert np.array_equal(out, int(c)*np.eye(d,dtype=int)), (d,out,c)
        n = d*d
        p_dep = c*c / (factorial(n)**2 * d**(2*n))
        exact.append(dict(d=d, basis_value=str(c), exact_depolarizing_success=str(p_dep),
                          success_float=float(p_dep), dp_peak_states=peak))
    rng = np.random.default_rng(20261028)
    max_residual = 0.
    for _ in range(6):
        X = rng.normal(size=(4,2,2)) + 1j*rng.normal(size=(4,2,2))
        Y = rng.normal(size=(4,2,2)) + 1j*rng.normal(size=(4,2,2))
        rhs = -3*np.linalg.det(X.reshape(4,4).T)*np.linalg.det(Y.reshape(4,4).T)*np.eye(2)
        lhs = regev2(X,Y)
        residual = np.linalg.norm(lhs-rhs)/(1+np.linalg.norm(rhs))
        max_residual = max(max_residual,residual)
        assert residual < 1e-11
    # Cauchy--Binet normalization for a generic full-rank qubit channel.
    Z = rng.normal(size=(12,2)) + 1j*rng.normal(size=(12,2))
    Q, _ = np.linalg.qr(Z)
    K = Q.reshape(6,2,2)
    V = K.reshape(6,4).T
    J = V@V.conj().T
    det_sum = sum(abs(np.linalg.det(V[:,p]))**2 for p in permutations(range(6),4))
    assert np.isclose(det_sum, factorial(4)*np.linalg.det(J).real,rtol=1e-12)
    p_from_sum = 9 * det_sum**2 / factorial(4)**4
    p_from_choi = 9 * np.linalg.det(J).real**2 / factorial(4)**2
    return dict(exact_matrix_unit_checks=exact,
                d2_random_tuple_checks=6,maximum_relative_residual=max_residual,
                choi_determinant=float(np.linalg.det(J).real),
                ordered_tuple_determinant_sum=float(det_sum),
                success_from_tuple_sum=float(p_from_sum),success_from_choi=float(p_from_choi))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,help='Optional JSON report path')
    args=parser.parse_args()
    text=json.dumps(verify(),indent=2)+'\n'
    print(text,end='')
    if args.output:
        args.output.write_text(text)


if __name__ == '__main__':
    main()
