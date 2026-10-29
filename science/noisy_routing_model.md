# Noisy routing of raw channel slots

The ideal [central-polynomial construction](polynomial_realization.md) cancels unwanted Kraus products exactly. Controller dephasing weakens the interference responsible for that cancellation. The recovery problem is then to determine the best conditional Bell fidelity at a charged acceptance probability, allowing every controller preparation and final measurement compatible with the supplied router.

This page fixes the common physical model for the recovery results here and the [retained-versus-accepted entanglement comparison](https://github.com/GoGoKo699/Quantum-Order-Entanglement/blob/main/science/noisy_separation.md). The two projects keep different outputs: this repository optimizes an accepted reference–target pair; the sister project studies a retained controller–target state.

<a id="shared-access"></a>
## 1. Access, independent slots, and accepted operations

The target is $`T\simeq\mathbb C^2`$, with an untouched reference $`R`$. The routing controller $`C`$ and any finite controller ancilla $`A`$ begin independently of the input:

```math
\rho_{RT}\otimes\omega_{CA}.
```

Each supplied slot is an independent ordinary quantum channel with an inaccessible environment. Every occupied order traverses every slot exactly once. No target encoding, target-carrying bypass, intermediate target gate, per-slot unitary wrapper, or accessible noise flag is supplied. Arbitrary controller processing and an instrument are allowed after the noisy switch. Each accepted classical outcome $`x`$ may select one final target unitary $`U_x`$. The controller can be mixed and entangled with $`A`$; collective final effects on $`CA`$ are included.

These conditions specialize the [general routing model](model.md). They distinguish the problem from arbitrary target circuits or a reusable noisy device with correlated memory. The supplied router is otherwise ideal; the result is not a gate-noise or fault-tolerance theorem.

For a Kraus tuple $`\boldsymbol a=(a_1,\ldots,a_m)`$ and occupied word $`\pi`$, define

```math
A_{\boldsymbol a,\pi}
=K^{(\pi_1)}_{a_{\pi_1}}\cdots K^{(\pi_m)}_{a_{\pi_m}},
\qquad
W_{\boldsymbol a}
=\sum_{\pi\in\Pi}A_{\boldsymbol a,\pi}\otimes|\pi\rangle\langle\pi|.
```

Words specify **matrix-product order**: their rightmost slot acts first chronologically. A slot's environmental index is the same across all order branches and independent of the other slots' indices. Independent Kraus changes leave the supplied channel unchanged, as proved for the normalized switch in [polynomial realization](polynomial_realization.md#1-the-normalized-switch).

<a id="physical-codebook"></a>
## 2. The physical eight-order controller

The four-call results use precisely these occupied words and codes:

| Matrix-product word | Physical control code |
|---|---|
| 1234 | 000 |
| 1243 | 001 |
| 2134 | 010 |
| 2143 | 011 |
| 3412 | 100 |
| 3421 | 110 |
| 4312 | 101 |
| 4321 | 111 |

The code bits are $`(c,x,y)`$. The outer bit exchanges the two blocks; $`x`$ and $`y`$ reverse the fixed channel pairs $`12`$ and $`34`$. The inner bits stay attached to those pairs when the blocks exchange. This is the nested-switch architecture used in [Kechrimparis et al. (2024)](https://doi.org/10.1088/1367-2630/ada05e) and [Kechrimparis, Moran, and Kwon (2025)](https://doi.org/10.1103/mvrg-yj3l).

In displayed-word order, the integer codes are $`[0,1,2,3,4,6,5,7]`$. Assigning successive binary integers $`[0,1,\ldots,7]`$ to this word list changes the physical local-noise model. All local-noise conclusions below use the displayed codebook. The optimum is over controllers on these eight words, not over all 24 permutations or all encodings.

## 3. Two unobserved controller-noise geometries

A Schur controller channel multiplies order coherences by a positive semidefinite correlation matrix $`\Gamma`$ with diagonal entries one:

```math
|\pi\rangle\langle\sigma|
\longmapsto\Gamma_{\pi\sigma}|\pi\rangle\langle\sigma|.
```

The noise acts immediately before or immediately after the switch, with no additional coherence-restoration stage before routing. These placements commute because every $`W_{\boldsymbol a}`$ is diagonal in the actual routing labels. Its event flag and environment are inaccessible. Final controller processing after the noisy switch is already included in the accepted effects.

Let $`\mathcal P_C`$ pinch the order basis. Global dephasing is

```math
\mathcal D_\eta=\eta\,\mathrm{id}+(1-\eta)\mathcal P_C,
\qquad
\Gamma^{\mathrm{glob}}_{\pi\sigma}
=\eta+(1-\eta)\delta_{\pi\sigma},
\qquad 0\le\eta\le1.
```

Independent phase flips of probability $`z`$ on each physical control qubit instead give

```math
\gamma=1-2z\in[0,1],
\qquad
\Gamma^{\mathrm{loc}}_{\pi\sigma}
=\gamma^{d_H(b_\pi,b_\sigma)}.
```

Here $`b_\pi`$ is the physical three-bit code and $`d_H`$ its Hamming distance. Thus $`\gamma=4/5`$ means a phase-flip probability $`1/10`$ per control qubit. It is not the global model with $`\eta=4/5`$. Their different coherence geometries can [reverse the preferred recovery herald](recovery_frontier.md#preferred-herald).

<a id="channel-family"></a>
## 4. A full-rank entanglement-breaking channel family

Let $`D_Z`$ measure in the computational basis and reprepare the measured basis state. With $`\Delta(X)=I\mathrm{Tr}(X)/2`$, define

```math
\mathcal M_\theta
=\mathrm{Ad}_{R_y(\theta)}\circ D_Z,
\qquad
\mathcal N_{\epsilon,\theta}
=(1-\epsilon)\mathcal M_\theta+\epsilon\Delta,
\qquad 0\le\epsilon\le1.
```

The two measure–prepare Kraus operators are

```math
|n_{\theta,+}\rangle\langle0|,
\qquad |n_{\theta,-}\rangle\langle1|,
\qquad
|n_{\theta,+}\rangle=
\begin{pmatrix}\cos(\theta/2)\\\sin(\theta/2)\end{pmatrix},
\qquad
|n_{\theta,-}\rangle=
\begin{pmatrix}-\sin(\theta/2)\\\cos(\theta/2)\end{pmatrix}.
```

The depolarizing component has Kraus operators $`\sqrt{\epsilon/2}\,E_{ab}`$, and the two displayed pure Kraus operators carry a factor $`\sqrt{1-\epsilon}`$. Every slot is entanglement breaking. Its normalized Choi eigenvalues are $`(2-\epsilon)/4`$ twice and $`\epsilon/4`$ twice, so every $`\epsilon>0`$ has full Choi rank.

At $`\theta=\pi/2`$, write $`s=1-\epsilon`$. The Bloch action is $`r\mapsto s\,\hat x r_z`$, so

```math
\mathcal N_{\epsilon,\pi/2}^{\,2}=\Delta,
\qquad
\mathcal N_{\epsilon,\pi/2}^{\,m}=\Delta
\quad(m\ge2).
```

This is the orthogonal measure–prepare point. The tilted family sets $`\theta=\pi/2+\delta`$. An output wrapper $`R_y(-\theta)`$ after **each** supplied slot would change the raw channel into a Pauli channel. Those wrappers are excluded. One final correction after the complete routed process does not generally implement gates between its uses. This access distinction is essential to the comparison with earlier Pauli-channel switch formulas.

## 5. Recovery objectives and resource accounting

Use the normalized Bell vector and the unnormalized accepted Choi output

```math
|\Phi_2\rangle=(|00\rangle+|11\rangle)/\sqrt2,
\qquad
\tau=(\mathrm{id}_R\otimes\mathcal E)
(|\Phi_2\rangle\langle\Phi_2|).
```

After incorporating each accepted herald's final correction, define

```math
p=\mathrm{Tr}\tau,
\qquad u=\langle\Phi_2|\tau|\Phi_2\rangle,
\qquad F=u/p\quad(p>0).
```

The recovery task requires input-independent acceptance, $`\mathcal E^\dagger(I)=pI`$. An arbitrary complex controller can have input-dependent success. The [controller reduction](controller_optimization.md#real-controller-reduction) proves that every achievable Bell pair $`(p,u)`$ in this real, $`Y`$-covariant family has a real, constant-success representative. It makes the unrestricted Bell-input upper bound tight for the required recovery task without assigning constant success to every allowed controller.

The upper frontier is $`u_*(p)=pF_*(p)`$. A four-call trial consumes four raw channel slots on success and failure. Repeating with a fresh input pair costs $`4/p`$ expected channel calls per accepted pair. It does not transmit an already lost unknown input again. Bell fidelity above $`1/2`$ witnesses reference–target entanglement; a fidelity at or below $`1/2`$ alone does not classify every possible output state's separability.

## 6. Proof route

| Question | Authoritative result |
|---|---|
| Which coefficient matrices can any controller realize? | [Controller feasibility, symmetry, and correction reduction](controller_optimization.md) |
| Can noisy routing still yield an exact nonzero unitary branch? | [Strict Schur-noise obstruction](noisy_recovery_obstruction.md) |
| What is the complete orthogonal-family recovery frontier? | [Physical Walsh readout and the attaining allocation](recovery_frontier.md) |
| Does simple readout remain optimal away from simultaneous diagonalization? | [Tilted diagonal-envelope theorem](tilted_readout.md) |
| How does accepted reference–target entanglement compare with a retained controller–target state? | [Protocol-specific robust separation](https://github.com/GoGoKo699/Quantum-Order-Entanglement/blob/main/science/noisy_separation.md) |

The finite coefficient identities and rational certificates used by the last three results are checked by [the independent verifier](../verification/audit_verify.py). The [verification guide](../verification/README.md) distinguishes exact reconstruction, analytical arguments, and numerical evaluations.
