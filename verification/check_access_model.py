"""Independent finite checks of the order-only access formulas.

These checks validate conventions, gauge invariance, and boundary examples.
They do not replace the analytic proofs.
Only NumPy is required.
"""

from itertools import permutations, product
from pathlib import Path
import argparse
import json
import numpy as np


RNG = np.random.default_rng(202610280314)
I2 = np.eye(2, dtype=complex)
PAULI = [
    I2,
    np.array([[0, 1], [1, 0]], dtype=complex),
    np.array([[0, -1j], [1j, 0]], dtype=complex),
    np.diag([1, -1]).astype(complex),
]


def vec(a):
    return a.reshape(-1, order="F")


def random_channel(d=2, rank=4):
    z = RNG.normal(size=(d * rank, d)) + 1j * RNG.normal(size=(d * rank, d))
    q, _ = np.linalg.qr(z)
    return [q[i * d : (i + 1) * d, :] for i in range(rank)]


def choi(kraus):
    return sum(np.outer(vec(k), vec(k).conj()) for k in kraus)


def apply_channel(kraus, rho):
    return sum(k @ rho @ k.conj().T for k in kraus)


def compose(after, before):
    return [a @ b for a in after for b in before]


def fidelity(kraus):
    d = kraus[0].shape[0]
    return float(sum(abs(np.trace(k)) ** 2 for k in kraus) / d**2)


def polynomial4(a, b, c, d):
    x = a @ b - b @ a
    y = c @ d - d @ c
    return x @ y + y @ x


def success4(slots):
    """Eight signed words, normalized by their coefficient l1 norm, 8."""
    bs = [polynomial4(*ks) / 8 for ks in product(*slots)]
    j = choi(bs)
    p = float(np.trace(j).real / 2)
    ideal = p * np.outer(vec(I2), vec(I2).conj())
    return p, j, float(np.linalg.norm(j - ideal))


def gauge_expand(kraus, new_rank):
    old_rank = len(kraus)
    z = RNG.normal(size=(new_rank, old_rank)) + 1j * RNG.normal(
        size=(new_rank, old_rank)
    )
    v, _ = np.linalg.qr(z)
    return [
        sum(v[b, a] * kraus[a] for a in range(old_rank))
        for b in range(new_rank)
    ]


def partial_trace_control(x, d, control_dimension):
    return np.trace(
        x.reshape(d, control_dimension, d, control_dimension),
        axis1=1,
        axis2=3,
    )


def mixed_ancillary_control_check():
    """Compare a direct switch with the polynomial Gram decomposition."""
    d, nc, na = 2, 2, 2
    dim = nc * na
    slots = [random_channel(), random_channel()]
    orders = list(permutations(range(2)))
    z = RNG.normal(size=(dim, dim)) + 1j * RNG.normal(size=(dim, dim))
    sigma = z @ z.conj().T
    sigma /= np.trace(sigma)
    z = RNG.normal(size=(dim, dim)) + 1j * RNG.normal(size=(dim, dim))
    effect = z @ z.conj().T
    effect /= 1.1 * np.linalg.eigvalsh(effect)[-1]
    ew, eu = np.linalg.eigh(effect)
    effect_root = (eu * np.sqrt(ew)) @ eu.conj().T
    projectors = [
        np.kron(np.diag([int(j == p) for j in range(nc)]), np.eye(na))
        for p in range(nc)
    ]
    gram = np.array(
        [
            [np.trace(effect @ pp @ sigma @ pq) for pq in projectors]
            for pp in projectors
        ]
    )
    assert np.linalg.eigvalsh(gram)[0] > -1e-12
    residuals = []
    tp = np.zeros((d * dim, d * dim), dtype=complex)
    histories = []
    for ks in product(*slots):
        words = [ks[o[1]] @ ks[o[0]] for o in orders]
        w = sum(np.kron(word, proj) for word, proj in zip(words, projectors))
        tp += w.conj().T @ w
        histories.append((words, w))
    for x in range(d):
        for y in range(d):
            rho = np.zeros((d, d), dtype=complex)
            rho[x, y] = 1
            direct = np.zeros_like(rho)
            decomposed = np.zeros_like(rho)
            for words, w in histories:
                joint = w @ np.kron(rho, sigma) @ w.conj().T
                f = np.kron(np.eye(d), effect_root)
                direct += partial_trace_control(f @ joint @ f, d, dim)
                for p, q in product(range(nc), repeat=2):
                    decomposed += gram[p, q] * (
                        words[p] @ rho @ words[q].conj().T
                    )
            residuals.append(float(np.linalg.norm(direct - decomposed)))
    return {
        "general_effect_ancillary_control_max_residual": max(residuals),
        "switch_trace_preservation_residual": float(
            np.linalg.norm(tp - np.eye(d * dim))
        ),
        "coefficient_gram_min_eigenvalue": float(np.linalg.eigvalsh(gram)[0]),
    }


def classical_bound_check():
    slots = [random_channel() for _ in range(3)]
    eta = [2 * np.linalg.eigvalsh(choi(s))[0] for s in slots]
    w = 1 - np.prod([1 - x for x in eta])
    lower_bound = float(w * (1 - 1 / 4))
    values = []
    for order in permutations(range(3)):
        ks = [I2]
        for j in order:
            ks = compose(slots[j], ks)
        infidelity = 1 - fidelity(ks)
        assert infidelity >= lower_bound - 1e-12
        values.append(float(infidelity))
    dep = [p / 2 for p in PAULI]
    assert abs((1 - fidelity(dep)) - 0.75) < 1e-14
    return {
        "slot_depolarizing_fractions": [float(x) for x in eta],
        "classical_infidelity_lower_bound": lower_bound,
        "minimum_infidelity_over_six_orders": min(values),
        "depolarizer_infidelity": float(1 - fidelity(dep)),
    }


def correlated_noise_check():
    n = 4
    theta = np.pi / (2 * n)
    unitaries = [I2] + [
        np.cos(theta) * I2 + 1j * np.sin(theta) * p
        for p in PAULI[1:]
    ]
    marginal = [u / 2 for u in unitaries]
    shared = [np.linalg.matrix_power(u, n) / 2 for u in unitaries]
    independent_probability, _, residual = success4([marginal] * n)
    correlated_polynomial_max = max(
        np.linalg.norm(polynomial4(u, u, u, u)) for u in unitaries
    )
    return {
        "single_slot_choi_min_eigenvalue": float(
            np.linalg.eigvalsh(choi(marginal))[0]
        ),
        "shared_noise_output_choi_distance_from_depolarizer": float(
            np.linalg.norm(choi(shared) - np.eye(4) / 2)
        ),
        "independent_slot_success_probability": independent_probability,
        "independent_slot_identity_residual": residual,
        "shared_noise_central_polynomial_max_norm": float(
            correlated_polynomial_max
        ),
    }


def verify():
    global RNG
    RNG = np.random.default_rng(202610280314)
    slots = [random_channel() for _ in range(4)]
    p, j, residual = success4(slots)
    expanded = [gauge_expand(k, 5) for k in slots]
    pg, jg, residual_gauge = success4(expanded)
    dep = [pauli / 2 for pauli in PAULI]
    pd, _, rd = success4([dep] * 4)
    near_identity = []
    for eps in [0.1, 0.3, 0.8]:
        channel = [np.sqrt(1 - eps) * I2] + [
            np.sqrt(eps) * k for k in dep
        ]
        pe, _, re = success4([channel] * 4)
        near_identity.append(
            {
                "epsilon": eps,
                "probability": pe,
                "expected_epsilon_four_scaling": eps**4 * pd,
                "identity_residual": re,
            }
        )
        assert abs(pe - eps**4 * pd) < 1e-13
    assert p > 0
    assert residual < 1e-12
    assert residual_gauge < 1e-12
    assert np.linalg.norm(j - jg) < 1e-12
    assert abs(pd - 3 / 64) < 1e-14
    result = {
        "random_heterogeneous_fullrank_success": p,
        "identity_choi_residual": residual,
        "expanded_kraus_gauge_success": pg,
        "kraus_gauge_choi_difference": float(np.linalg.norm(j - jg)),
        "depolarizing_success": pd,
        "depolarizing_identity_residual": rd,
        "near_identity_scaling": near_identity,
        **mixed_ancillary_control_check(),
        **classical_bound_check(),
        **correlated_noise_check(),
    }
    assert result["general_effect_ancillary_control_max_residual"] < 1e-12
    assert result["switch_trace_preservation_residual"] < 1e-12
    assert result["single_slot_choi_min_eigenvalue"] > 0
    assert result["shared_noise_output_choi_distance_from_depolarizer"] < 1e-12
    assert result["independent_slot_success_probability"] > 0
    assert result["shared_noise_central_polynomial_max_norm"] < 1e-12
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Optional JSON report path")
    args = parser.parse_args()
    text = json.dumps(verify(), indent=2) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
