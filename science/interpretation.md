# What exact heralded recovery means

The [model](model.md) supplies an ideal coherent ordering resource for independent ordinary channels. The [construction](polynomial_realization.md) gives an accepted operation $`p\,\mathrm{id}_d`$, and the [query theorem](query_threshold.md) characterizes the minimum number of compulsory channel traversals. This page states the operational consequence and the boundaries of those claims.

## Heralded entanglement generation

Prepare a fresh maximally entangled state

```math
|\Phi_d\rangle
=\frac1{\sqrt d}\sum_{i=1}^d|i\rangle|i\rangle
```

and send one half through the protocol. Every accepted trial returns exactly the same maximally entangled state across the retained reference and the target output. The classical success flag identifies that event and must be made available to the users of the pair. Repeating on fresh pairs gives a positive heralded entanglement yield whenever $`p>0`$.

A failed trial need not preserve the target. An arbitrary unknown message therefore cannot generally be retried after failure. The fresh-pair procedure also does not by itself establish positive unassisted quantum capacity: a flagged channel may have a perfect successful branch while its overall failure probability precludes unassisted transmission.

The coherent routing and its control remain charged resources of this statement. An abstract switch of channel boxes is not automatically a protocol on an ordinary fixed one-way communication network. [Probability and resources](probability.md) separates noisy-slot counts, accepted probabilities, and controller costs.

## Classical-order comparison

The classical comparison uses the same compulsory channel slots and an initially independent controller, but chooses a complete order through classical randomness. Its success flag may depend on the classical controller record. It cannot inspect the target or the noise environments.

Conditioned on any nonzero accepted event, the normalized target channel consequently has the form

```math
\mathcal C
=\sum_{\pi\in S_m}r_\pi\,
 \mathcal N_{\pi_1}\circ\cdots\circ\mathcal N_{\pi_m},
\qquad
r_\pi\ge0,\qquad \sum_\pi r_\pi=1.
\qquad\text{(1)}
```

The rightmost channel acts first, consistently with the [word convention](model.md#2-compulsory-traversals-and-word-order). For identical slots, every term is $`\mathcal N^m`$. Postselecting the independent classical record cannot change that channel.

For a full-Choi-rank slot, put

```math
\eta_j=d\,\lambda_{\min}(J_{\mathcal N_j})\in(0,1].
```

The unnormalized Choi convention gives a channel decomposition

```math
\mathcal N_j
=\eta_j\Delta+(1-\eta_j)\mathcal M_j.
\qquad\text{(2)}
```

When $`\eta_j\lt 1`$, positivity follows from
$`J_{\mathcal N_j}-\eta_jI_{d^2}/d\ge0`$, and its channel normalization follows by taking the appropriate partial trace. The case $`\eta_j=1`$ is $`\mathcal N_j=\Delta`$.

### A quantitative gap from the identity

Define

```math
w=1-\prod_{j=1}^m(1-\eta_j).
```

Expanding any complete composition in (1) using (2), the terms containing at least one depolarizer have total weight $`w`$. Every such term is a replacer channel

```math
\mathcal R_\sigma(\rho)=\sigma\mathrm{Tr}\rho.
```

The first depolarizer removes dependence on the input; subsequent trace-preserving maps only change its fixed output. Combining those terms across all selected orders yields

```math
\mathcal C=(1-w)\mathcal T+w\mathcal R_\sigma
\qquad\text{(3)}
```

for a channel $`\mathcal T`$ and a state $`\sigma`$.

Let the entanglement fidelity relative to identity be

```math
F_e(\mathcal C)=
\langle\Phi_d|
(\mathrm{id}\otimes\mathcal C)
(|\Phi_d\rangle\langle\Phi_d|)
|\Phi_d\rangle.
```

Every replacer has $`F_e(\mathcal R_\sigma)=1/d^2`$, and $`F_e(\mathcal T)\le1`$. Hence

```math
\boxed{
1-F_e(\mathcal C)
\ge
\left[1-\prod_{j=1}^m(1-\eta_j)\right]
\left(1-\frac1{d^2}\right).
}
\qquad\text{(4)}
```

The same right side lower-bounds
$`\tfrac12\|\mathcal C-\mathrm{id}\|_\diamond`$: use the maximally entangled input and measure its projector at the output. Known unitary corrections depending on a classical outcome do not evade the bound because they take replacers to replacers.

The gap is strictly positive for full-rank slots, so an accepted classical-order channel in this model cannot be exactly identity. For completely depolarizing slots, it is exactly $`\Delta`$ and equality holds in (4). Coherent order instead permits a positive-probability accepted branch with unit fidelity. This comparison concerns conditional exactness and does not assert an advantage in throughput.

An unrestricted causal circuit with protected target-carrying memory is outside both sides of this comparison. Such a memory could simply retain the input while the required noisy calls act elsewhere.

## Correlated instances

Full-rank single-slot marginals do not specify the joint action of a device whose errors are correlated across uses. The following qubit example gives an exact failure of the independent-instance conclusion.

Let

```math
U_0=I_2,\qquad
U_j=\exp\!\left(\frac{i\pi}{8}\sigma_j\right),
\quad j=1,2,3.
```

Choose $`\mu\in\{0,1,2,3\}`$ uniformly once per trial, and let all four slots apply the same $`U_\mu`$. Each one-slot marginal is

```math
\mathcal N(\rho)
=\frac14\sum_{\mu=0}^3U_\mu\rho U_\mu^\dagger.
\qquad\text{(5)}
```

The four matrices span $`M_2(\mathbb C)`$: subtracting the appropriate multiple of $`U_0`$ from $`U_j`$ gives a nonzero multiple of each Pauli matrix. Thus (5) has full Choi rank.

Nevertheless every order of the four correlated slots applies $`U_\mu^4`$. Since

```math
U_0^4=I_2,\qquad U_j^4=i\sigma_j,
```

the complete target channel is

```math
\frac14\left(\rho+\sigma_x\rho\sigma_x+
\sigma_y\rho\sigma_y+\sigma_z\rho\sigma_z\right)
=\Delta(\rho).
\qquad\text{(6)}
```

For each shared label, all order words are identical. The target therefore stays independent of the controller, and every nonzero control-only accepted branch is proportional to $`\Delta`$. No herald can restore identity.

Four independent copies of the same marginal (5) satisfy the qubit theorem. These four correlated instances do not. Independence of the supplied channels, rather than full rank of their separate Choi matrices, is the required hypothesis.

Some restricted correlations can preserve scalarity while making success vanish. More general shared memory can make a path amplitude cease to be a permutation of a fixed tuple of $`d\times d`$ matrices. The theorem does not classify either situation or establish a general memory-channel result.

## Rank-deficient channels

Full Choi rank guarantees positivity for every nonidentity central-polynomial controller and enables the spanning converse. It is not a necessary condition for every possible instance of exact recovery.

For example, a unitary channel generated by $`U`$ with $`U^m`$ proportional to identity is exactly recovered by $`m`$ ordinary consecutive uses. This is a rank-one Choi channel, and it does not obey the full-rank minimum-use conclusion.

Other rank-deficient channels cannot be recovered by order control. For independent copies of a completely dephasing channel in one fixed basis, use its diagonal projector Kraus operators. They commute, so every permutation gives the same product for each history. The switch output factors into the ordinary target channel and the initial controller. A control-only accepted outcome can only multiply the dephasing channel by a scalar.

This also explains a boundary of the central-polynomial construction. If the arguments commute, every permutation word equals the same product, so the polynomial reduces to

```math
P(X_1,\ldots,X_m)
=\left(\sum_\pi c_\pi\right)X_1\cdots X_m.
```

For a positive-degree multilinear central polynomial in dimension $`d>1`$, setting all but one variable to identity forces $`\sum_\pi c_\pi=0`$. Such a polynomial therefore vanishes on every commuting tuple.

The [probability account](probability.md) gives the full-rank noise-floor bound and the precise suppression near the identity-channel boundary. Its probability and expected-call laws concern the stated full-rank, unaugmented model; they are not classifications of rank-deficient channels.

## Randomizing wrappers

If known output postprocessing may be attached inside each coherently routed box, a finite unitary 1-design implements the completely depolarizing channel. For instance, uniform conjugation by the $`d^2`$ Weyl operators gives $`\Delta`$. Each original slot can then be replaced by

```math
\widetilde{\mathcal N}_j
=\Delta\circ\mathcal N_j
=\Delta.
\qquad\text{(7)}
```

Applying a central-polynomial controller to these composite slots gives an exact positive-probability branch for arbitrary original CPTP channels. Every order still uses each original channel and its attached randomizer once. This introduces no protected target-carrying bypass.

It does, however, change the supplied resource. The original channel becomes informationally irrelevant, and the randomizing devices together with coherent order provide the capability. This corollary is subordinate to the unmodified full-rank theorem. In particular, the unaugmented small-noise lower bounds in [probability and resources](probability.md) do not apply after these wrappers are allowed.

Each randomizing environment belongs to its physical slot, consistently across all coherent order branches. Choosing unrelated randomizations for different branches can reveal the order and would not implement the switch of the composite channels in (7).

## Ideal control and conditional error

The central-polynomial identity remains exact when the independent channel maps inside the declared model change: it holds for every matrix tuple. Such a change can alter the success probability, or make it zero at a rank-deficient boundary.

Routing mistakes, incorrect control amplitudes, and additional target interactions are different errors. They can invalidate the prescribed polynomial. A simple conditional-error bound displays the role of a small success probability.

Suppose the implemented accepted operation $`\widetilde{\mathcal E}`$ is completely positive and trace nonincreasing, with

```math
\|\widetilde{\mathcal E}
-p\,\mathrm{id}\|_\diamond\le\delta.
```

For any normalized input including a reference, let $`\widetilde\sigma`$ be the unnormalized accepted output and $`\widetilde p=\mathrm{Tr}\widetilde\sigma`$. Then

```math
\|\widetilde\sigma-p\rho\|_1\le\delta,
\qquad |\widetilde p-p|\le\delta.
```

If $`p>\delta`$, normalization and the triangle inequality give

```math
\left\|\frac{\widetilde\sigma}{\widetilde p}-\rho\right\|_1
\le
\frac{\|\widetilde\sigma-p\rho\|_1+|p-\widetilde p|}
{\widetilde p}
\le\frac{2\delta}{p-\delta}.
\qquad\text{(8)}
```

The normalized successful transformation can depend nonlinearly on its input through its success probability, so (8) is a bound for each input, including references, rather than a diamond distance between normalized channels.

This is a conditioning estimate, not a hardware fault-tolerance threshold. It shows why a small positive accepted probability needs a separate control-accuracy account. The [resource analysis](probability.md) also distinguishes noisy-call counts from the synthesis and coherence costs of many order states.

## Noisy exactness and useful recovery

The ideal central-polynomial theorem and a noisy accepted Bell state describe different endpoints. Under the [strict Schur assumptions](noisy_recovery_obstruction.md#strict-schur-obstruction), every nonzero accepted branch contains an unavoidable full-Choi-rank diagonal component, so no final target unitary can make it exact. The [frontier](recovery_frontier.md) nevertheless quantifies useful approximate recovery, charging the accepted probability instead of retaining only its best conditional fidelity.

The physical geometry of controller noise also matters. A single scalar controller coherence or fidelity does not determine the best herald: the [global/local comparison](recovery_frontier.md) exhibits a branch-ranking reversal even at matched scalar coherence. The actual three-bit codebook is therefore a load-bearing part of the stated local-noise theorem.

For the same raw channel family and per-bit phase noise, [Quantum Order Entanglement proves](https://github.com/GoGoKo699/Quantum-Order-Entanglement/blob/main/science/noisy_separation.md#uniform-separation) a robust specified-protocol distinction: three-call retained $`C:T`$ outputs are separable while the fixed four-call accepted $`R:T`$ pair remains NPT. This compares different bipartitions, conditioning, controller sizes and call counts. It gives neither an optimized equal-resource separation nor a no-go theorem for all retained-controller protocols.

Bell fidelity above one half is a sufficient two-qubit entanglement witness. A lower optimized Bell fidelity does not show that every complex-controller output is separable. Likewise, a retained $`C:T`$ witness on a known target input does not certify recovery of an unknown input. These inference boundaries join the two tracks without conflating their tasks.

## Relation to established communication results

Coherent-order communication activation is an established phenomenon. In particular, cyclic-order protocols already provide higher-dimensional heralded transmission with an error that decreases as more channels are used. Central-polynomial quantum resetting also predates this channel model.

The result here is finite exact identity heralding through ordinary inaccessible-environment channel slots, together with the minimum-use characterization in that restricted architecture. The precise inherited results, source versions, and implication boundaries are given in [Related work](../references/RELATED_WORK.md), with bibliographic data in [references.bib](../references/references.bib).
