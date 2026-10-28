# Channels, order control, and the recovery task

The target is a quantum system of known dimension $`d\ge2`$. The available resource is a generalized switch of independent ordinary quantum channels. Every order uses every supplied channel once. An independent controller selects a coherent superposition of those orders and is measured at the end.

The [polynomial realization](polynomial_realization.md) proves how this resource produces an exact successful identity branch. The [query characterization](query_threshold.md) proves the converse and the minimum number of channel uses. This page fixes the common model and conventions.

## 1. Independent ordinary channel slots

Let $`T`$ denote the target. There are $`m\ge1`$ channel slots

```math
\mathcal N_j:\mathcal B(T)\longrightarrow\mathcal B(T),
\qquad
\mathcal N_j(\rho)
=\sum_{a_j}K^{(j)}_{a_j}\rho K^{(j)\dagger}_{a_j},
\qquad j=1,\ldots,m.
```

Each $`\mathcal N_j`$ is completely positive and trace preserving. The channels may differ. When all slots contain the same channel $`\mathcal N`$, they are independent fresh instances of that map. A reusable device is covered only when its uses implement these independent channel instances, for example through an appropriate environment reset.

The Kraus labels describe inaccessible noise environments. The controller receives neither those labels nor a coherent Stinespring implementation. It cannot measure, select, or reverse a channel's environmental interaction. Independence is a property of the joint supplied resource, not merely of its individual marginal descriptions. The [correlated-instance counterexample](interpretation.md#correlated-instances) shows that full-rank marginal channels alone are insufficient.

## 2. Compulsory traversals and word order

The order register $`C`$ has an orthonormal basis $`\{|\pi\rangle:\pi\in S_m\}`$. A word
$`\pi=(\pi_1,\ldots,\pi_m)`$ lists **mathematical matrix-product order**. Its rightmost slot acts first in physical time. Thus

```math
A_{\boldsymbol a,\pi}
=K^{(\pi_1)}_{a_{\pi_1}}\cdots K^{(\pi_m)}_{a_{\pi_m}},
\qquad
W_{\boldsymbol a}
=\sum_{\pi\in S_m}
A_{\boldsymbol a,\pi}\otimes|\pi\rangle\langle\pi|.
\tag{1}
```

The supplied switch is the channel

```math
\mathfrak S_{\mathcal N_1,\ldots,\mathcal N_m}(X)
=\sum_{\boldsymbol a}W_{\boldsymbol a}XW_{\boldsymbol a}^{\dagger}.
\tag{2}
```

The tuple $`\boldsymbol a=(a_1,\ldots,a_m)`$ uses the same environmental label $`a_j`$ for slot $`j`$ in every order branch. All $`m`$ slots occur exactly once in every word. Only a subset of the $`m!`$ order states need be populated.

Equations (1) and (2) are the standard generalized-switch resource used in [Kechrimparis, Moran, and Kwon, published Appendix A](../references/RELATED_WORK.md#kechrimparis2025). Its trace preservation and invariance under independent changes of Kraus representation are proved in [polynomial realization, Section 1](polynomial_realization.md#1-the-normalized-switch). A chosen Kraus decomposition is a mathematical representation, not additional operational access.

The count $`m`$ measures compulsory noisy-slot uses. A zero-slot identity wire is excluded from this conversion task. An additional compulsory slot cannot be ignored. Consequently a theorem about the minimum feasible positive $`m`$ does not by itself establish feasibility at every larger integer.

## 3. Allowed controller and successful instrument

The controller may include an auxiliary register $`A`$ of arbitrary finite dimension. Its initial density operator $`\omega_{CA}`$ may be mixed or entangled between $`C`$ and $`A`$. It must be independent of the target and every reference: for an arbitrary input $`\rho_{RT}`$, the initial state is

```math
\rho_{RT}\otimes\omega_{CA}.
```

The switch acts on $`T`$ and $`C`$, and as identity on $`R`$ and $`A`$. After the switch, arbitrary processing and a quantum instrument on $`CA`$ are allowed. If the controller's quantum output is discarded, an accepted outcome or collection of outcomes is specified by an effect $`0\le F_{CA}\le I_{CA}`$. The corresponding target operation is

```math
\mathcal E_F(\rho)
=\operatorname{Tr}_{CA}\!\left[
(I_T\otimes F_{CA}^{1/2})
(\mathfrak S_{\mathcal N_1,\ldots,\mathcal N_m}
 \otimes\operatorname{id}_A)(\rho\otimes\omega_{CA})
(I_T\otimes F_{CA}^{1/2})
\right].
\tag{3}
```

This expression includes pure or mixed preparations, rank-one postselection, general positive effects, and coarse-grained successful events. The form of the final instrument matters only through its accepted effect when its quantum output is discarded. Its outcomes may supply a classical success flag.

The controller may depend on the complete descriptions of the channels when proving a lower bound for a fixed instance. The universal construction instead uses a controller depending only on the target dimension and a chosen central polynomial. These are distinct quantifiers in the [query theorem](query_threshold.md).

## 4. Target access and explicit exclusions

The only target action between the initial input and final output is the supplied routing through the complete channel slots. In particular, the core model permits no:

- Encoding of the input into a protected ancillary system.
- Target-carrying bypass around a channel slot.
- Known target gates or target-control interactions inserted between slots.
- Access to a channel's environment, Kraus label, or coherent dilation.
- Shared channel memory in place of the independent resource in Section 1.

An initially independent controller ancilla is allowed because it does not carry an encoded input around the noisy slots. Arbitrary target-processing quantum memories are different resources. Without this distinction, an unrestricted fixed-order circuit could store the input in a noiseless memory while using the channels on another system.

An optional extension allows a fixed known target unitary correction chosen from the classical successful herald. [The query characterization](query_threshold.md#4-final-unitary-corrections-do-not-lower-the-count) proves that these corrections do not lower the minimum for full-rank channels. That result does not authorize target gates between channel slots.

Known randomizing postprocessing attached inside each routed box is a separate augmented access model. Its [all-channel corollary](interpretation.md#randomizing-wrappers) does not change the assumptions of the unaugmented minimum-use or resource bounds.

Ideal coherent order control is supplied. The slot count is not a gate-level implementation cost and does not include a theorem that an ordinary fixed-order circuit can synthesize the switch with one call to each unknown channel.

## 5. Exact successful recovery

An exact identity branch means equality of completely positive maps,

```math
\mathcal E_F=p\,\operatorname{id}_d,\qquad p>0.
\tag{4}
```

Its success probability is $`p`$ on every normalized input. For every untouched reference $`R`$ and joint input,

```math
(\operatorname{id}_R\otimes\mathcal E_F)(\rho_{RT})
=p\,\rho_{RT}.
```

Thus the task preserves arbitrary unknown states and their entanglement on the accepted branch. It does not merely preserve a selected alphabet or an average fidelity. The failure branch is unrestricted apart from belonging to the physical instrument. In particular, failure need not preserve the input for another attempt.

The [operational interpretation](interpretation.md) uses fresh entangled-pair trials. Exact heralded recovery does not, on its own, imply a positive unassisted quantum capacity.

## 6. Choi and full-rank conventions

Choose

```math
|\Omega_d\rangle=\sum_{i=1}^d|i\rangle|i\rangle,
\qquad
|K\rangle\!\rangle=(K\otimes I)|\Omega_d\rangle.
```

The Choi matrix is **unnormalized**, with the output factor written first:

```math
J_{\mathcal N}
 =(\mathcal N\otimes\operatorname{id})
   (|\Omega_d\rangle\langle\Omega_d|)
 =\sum_a|K_a\rangle\!\rangle\langle\!\langle K_a|,
\qquad
\operatorname{Tr}J_{\mathcal N}=d.
\tag{5}
```

Trace preservation is $`\operatorname{Tr}_{\mathrm{output}}J_{\mathcal N}=I_d`$.

For the completely depolarizing channel
$`\Delta(\rho)=I_d\operatorname{Tr}(\rho)/d`$,

```math
J_\Delta=\frac{I_{d^2}}d,
\qquad
J_{p\,\operatorname{id}}=
p|\Omega_d\rangle\langle\Omega_d|.
```

Full Choi rank means $`\operatorname{rank}J_{\mathcal N}=d^2`$. Equivalently, a Kraus family spans $`M_d(\mathbb C)`$. This equivalence and its use in the positivity proof are established in [Theorem 1](polynomial_realization.md#theorem-1-universal-exact-branch). Full rank is sufficient for the universal construction to succeed, and it is the hypothesis of the fixed-channel converse. It is not a claim that every rank-deficient channel is impossible to distill.

## 7. Result map

Each proof has one main location.

| Result | Authoritative account |
|---|---|
| Switch normalization, Kraus-representation invariance, and realization of order polynomials | [Polynomial realization](polynomial_realization.md) |
| Universal exact branch and full-rank positivity | [Theorem 1](polynomial_realization.md#theorem-1-universal-exact-branch) |
| Fixed-channel converse, general control instruments, final unitary corrections, and minimum $`c(d)`$ | [Query characterization](query_threshold.md) |
| Explicit qubit construction and finite lower-degree certificate | [Qubit account](qubit.md) |
| Explicit eight-call qutrit construction and inherited degree minimum | [Qutrit account](qutrit.md) |
| Success formulas, probability laws, and expected-call accounting | [Probability and resources](probability.md) |
| Classical-order gap, correlated-noise counterexample, boundary cases, and augmented wrappers | [Interpretation](interpretation.md) |
| Antecedents and the distinction between channel access and coherent Kraus access | [Related work](../references/RELATED_WORK.md) |

Bibliographic records are collected in [references.bib](../references/references.bib).
