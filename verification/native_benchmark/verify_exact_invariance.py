#!/usr/bin/env python3
"""Independent exact audit of the full three-call controller converse.

Standard-library only: python3 verify_exact_invariance.py --output-dir DIR
Writes exact_invariance_audit.json in DIR (default: beside this script).
Words use matrix-product order 012,021,102,120,201,210.
This checks finite algebra, not literature novelty or physical hardware costs.
"""
from fractions import Fraction as F
from pathlib import Path
import itertools as it
import json
import argparse
import math

N = 6
ORDERS = list(it.permutations(range(3)))
PARITY = [(-1) ** sum(o[i] > o[j] for i in range(3) for j in range(i + 1, 3)) for o in ORDERS]


def mm(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def transpose(a):
    return list(map(list, zip(*a)))


def zero(a):
    return all(v == 0 for row in a for v in row)


def scale(c, a):
    return [[c * v for v in row] for row in a]


def add(a, b):
    return [[x + y for x, y in zip(ar, br)] for ar, br in zip(a, b)]


I = [[F(i == j) for j in range(N)] for i in range(N)]
REV = [[F(a == tuple(reversed(b))) for b in ORDERS] for a in ORDERS]
T = [[F(1, 6) for _ in range(N)] for _ in range(N)]
S = [[F(PARITY[i] * PARITY[j], 6) for j in range(N)] for i in range(N)]
O = add(scale(F(1, 2), add(I, scale(-1, REV))), scale(-1, S))
E = add(scale(F(1, 2), add(I, REV)), scale(-1, T))
STD = add(O, E)
PROJECTORS = {'T': (T, 1), 'S': (S, 1), 'O': (O, 2), 'E': (E, 2)}


def integer_word_vectors():
    # K_a=M_a/sqrt(2); three-letter words have common factor 1/(2 sqrt(2)).
    matrices = [[[1, 0], [1, 0]], [[0, 1], [0, -1]]]
    out = []
    for history in it.product(range(2), repeat=3):
        row = []
        for order in ORDERS:
            a = [[1, 0], [0, 1]]
            for j in order:
                a = mm(a, matrices[history[j]])
            row.append([v for ar in a for v in ar])
        out.append(row)
    return out


def native_flags(b):
    return [('O', (4 - b) / 12, (1 + 2 * b) / 12),
            ('E', (4 - b) / 12, (1 + b) / 12),
            ('T', (1 + 2 * b) / 6, (1 + 5 * b) / 24),
            ('S', (1 - b) / 6, (1 - b) / 24)]


def vertex_useful(a, b):
    """Direct rational comparison of every native cumulative cost vertex."""
    flags = sorted((z for z in native_flags(b) if z[1] > 0), key=lambda z: z[2] / z[1], reverse=True)
    p = F(0)
    f = F(0)
    floor2 = 3 * (1 - a) / (2 * (2 - a)) if a > F(1, 2) else F(1, 2)
    for _, prob, overlap in flags:
        p += prob
        f += overlap
        error = 1 - f / p
        cost = 3 / p
        if error < floor2:
            return True
        if a > F(1, 2) and floor2 <= error < F(1, 2):
            baseline = 1 + 2 * (6 + a) * (F(1, 2) - error) / (2 * a - 1)
            if cost < baseline:
                return True
    return False


def closed_useful(a, b):
    if a <= F(1, 2):
        return b > F(2, 5)
    return b > min(2 * a / (3 - a), (57 * a + 4) / (13 * a + 52))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parent,help='directory for the regenerated exact audit JSON')
    args=parser.parse_args()
    output_dir=args.output_dir.resolve();output_dir.mkdir(parents=True,exist_ok=True)
    for name, (p, rank) in PROJECTORS.items():
        assert mm(p, p) == p
        assert sum(p[i][i] for i in range(N)) == rank
        assert all(p[i][i] == F(rank, 6) for i in range(N))
    assert add(add(T, S), add(O, E)) == I
    for (name_a, (a, _)), (name_b, (b, _)) in it.combinations(PROJECTORS.items(), 2):
        assert zero(mm(a, b)), (name_a, name_b)

    ivecs = integer_word_vectors()
    left_actions = [
        [[F(a == tuple(g[j] for j in b)) for b in ORDERS] for a in ORDERS]
        for g in ORDERS
    ]
    for a, b in it.product(range(4), repeat=2):
        kernel = [[sum(v[i][a] * v[j][b] for v in ivecs) for j in range(N)] for i in range(N)]
        for action in left_actions:
            assert mm(mm(transpose(action), kernel), action) == kernel
        assert zero(mm(mm(T, kernel), STD))
        assert zero(mm(mm(STD, kernel), T))
        assert zero(mm(S, kernel)) and zero(mm(kernel, S))

    # Bell vectors for I,U,Y,V; omitted scalar/phase factors do not affect zeros.
    integer_bells = [[1, 0, 0, 1], [-1, 1, 1, 1], [0, -1, 1, 0], [1, 1, 1, -1]]
    scores = [[[sum(q[a] * v[i][a] for a in range(4)) for i in range(N)] for v in ivecs]
              for q in integer_bells]

    def bell_kernel(a, b):
        return [[sum(scores[a][h][i] * scores[b][h][j] for h in range(8))
                 for j in range(N)] for i in range(N)]

    wants = [scale(F(3, 2), T), scale(F(3, 4), O),
             scale(F(1, 2), add(T, E)), add(T, scale(F(1, 4), E))]
    for a, denom, want in zip(range(4), [32, 64, 32, 64], wants):
        assert scale(F(1, denom), bell_kernel(a, a)) == want
    for a, b in it.combinations(range(4), 2):
        if (a, b) != (1, 3):
            assert zero(bell_kernel(a, b))
    uv = bell_kernel(1, 3)
    assert mm(mm(O, uv), E) == uv
    assert mm(mm(E, transpose(uv)), O) == transpose(uv)

    # Test the full unequal-noise strict-usefulness predicate independently.
    # Rational comparisons cover the grid AND exact boundaries plus nearby points.
    cases = {(F(i, 50), F(j, 50)) for i in range(51) for j in range(51)}
    for i in range(51):
        a = F(i, 50)
        bound = F(2, 5) if a <= F(1, 2) else min(2 * a / (3 - a), (57 * a + 4) / (13 * a + 52))
        for shift in [F(0), F(-1, 10**10), F(1, 10**10)]:
            if 0 <= bound + shift <= 1:
                cases.add((a, bound + shift))
    for a, b in sorted(cases):
        assert vertex_useful(a, b) == closed_useful(a, b), (a, b)

    # Optional duration-law boundary, using only standard-library bisection.
    polynomial = lambda t: 13 * t**5 + 52 * t**3 - 57 * t**2 - 4
    lo, hi = 0.9, 1.0
    assert polynomial(lo) < 0 < polynomial(hi)
    for _ in range(70):
        mid = (lo + hi) / 2
        if polynomial(mid) > 0:
            hi = mid
        else:
            lo = mid
    t = (lo + hi) / 2

    result = {
        'status': 'PASS',
        'arithmetic': 'fractions.Fraction; all algebraic invariance and phase-predicate comparisons exact',
        'orders': [''.join(map(str, order)) for order in ORDERS],
        'target_coordinate_kernels': 16,
        'left_actions_per_kernel': 6,
        'left_invariance_equalities': 96,
        'trivial_standard_zero_block_equalities': 32,
        'sign_zero_side_equalities': 32,
        'bell_diagonal_kernel_identities': ['M_I=3T/2', 'M_U=3O/4', 'M_Y=(T+E)/2', 'M_V=T+E/4'],
        'zero_cross_bell_pairs': ['I,U', 'I,Y', 'I,V', 'U,Y', 'Y,V'],
        'sole_cross_bell_kernel': 'C_UV=O C_UV E; C_VU=E C_VU O (transpose/conjugate orientation)',
        'integer_UV_kernel': uv,
        'projector_ranks_and_constant_diagonals': {name: {'rank': rank, 'diagonal': str(F(rank, 6))}
                                                 for name, (_, rank) in PROJECTORS.items()},
        'complex_gram_scope': 'The exact kernels cover arbitrary complex PSD G. U,V coherence may be complex; the PSD two-dimensional block trace bounds its largest eigenvalue.',
        'unequal_noise_predicate': 'a<=1/2: b>2/5; a>1/2: b>min(2a/(3-a),(57a+4)/(13a+52))',
        'unequal_noise_exact_test_cases': len(cases),
        'unequal_noise_mismatches': 0,
        'duration_law': {
            'assumption': 'eta_m=exp(-m gamma), with no additional architecture-dependent duration overhead',
            'threshold_a': t*t,
            'threshold_b': t*t*t,
            'threshold_gamma': -math.log(t),
            'equation': '13t^5+52t^3-57t^2-4=0',
        },
        'physical_scope': [
            'Independent compulsory channel slots, inaccessible environments, initially independent target/reference and controller.',
            'Global unflagged logical order-basis dephasing after routing and before ANY final controller processing.',
            'Architecture-specific retained coherences a=eta_2,b=eta_3 are independent model parameters unless a duration law is explicitly imposed.',
            'Fine controller outcomes may have input-dependent success; aggregate accepted channel must have input-independent p.',
            'The corrected accepted TARGET-ONLY channel averages/discards the fine controller/correction label.',
            'No hardware, equal-time, device-independent causal, capacity, or novelty conclusion is certified by this audit.',
        ],
    }
    path = output_dir/'exact_invariance_audit.json'
    path.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ['status', 'left_invariance_equalities', 'unequal_noise_exact_test_cases', 'unequal_noise_mismatches', 'duration_law']}, indent=2))
    print(path)


if __name__ == '__main__':
    main()
