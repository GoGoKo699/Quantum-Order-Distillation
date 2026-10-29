# From quantum channels to exact recovery by coherent order

Suppose a noisy box replaces every input qubit with the maximally mixed state $`I_2/2`$. Sending the same qubit through four independent copies in any fixed order still gives $`I_2/2`$. The target alone carries no information about its initial state.

This project asks what changes when a supplied quantum controller coherently chooses the order of those four complete channel uses. A suitable final control measurement has a successful outcome whose target operation is exactly

```math
\mathcal E_{\mathrm{succ}}(\rho)=\frac3{64}\rho.
```

The coefficient is the probability of that outcome. Conditioned on success, every input state is preserved, including its entanglement with another system. The [worked example below](#worked-example) calculates the probability and specifies the preparation, routing labels, and measurement.

The ideal-control result concerns independent fresh channels with full Choi rank in a precise [order-only model](../science/model.md). It connects the minimum positive number of compulsory channel calls to the minimum degree of a central polynomial. The noisy-control results then identify when exactness is obstructed and determine complete attainable acceptance–fidelity frontiers for a specified eight-order qubit architecture. The two parts share the channel and instrument foundations, while requiring additional polynomial and positive-matrix arguments.

The selected educational anchor is **John Watrous, *The Theory of Quantum Information* (Cambridge University Press, 2018)**. The [author's page](https://jhwatrous.github.io/books-and-courses.html) provides the [standard-layout manuscript](https://jhwatrous.github.io/TQI.pdf) and [errata](https://jhwatrous.github.io/TQI.errata.pdf). The manuscript is offered for personal use and must not be redistributed. All page numbers below are its **printed page numbers**, rather than PDF viewer indexes. The errata are dated 23 July 2025.

This bridge supplies the connections needed for the project; it does not require reading the entire book. The [30-milestone map](MILESTONES.md) organizes the complete scientific story into 13 background milestones and 17 narrative milestones. The authoritative statements and proofs remain in [science/](../science/README.md).

<a id="source-map"></a>

## 1. A focused route through the anchor

Use the source entries when their concepts become relevant. The topic labels in the first column describe their role here.

| Concept needed for this project | Watrous location | Where it enters the project |
|---|---|---|
| Operators, tensor products, vectorization, and positive matrices | §§1.1.2–1.1.3, starting pp. 8 and 24; vectorization on pp. 23–24 | Choi matrices, Kraus spans, and positive sums |
| States, retained reference systems, and discarding a subsystem | §§2.1.2–2.1.3, starting pp. 61 and 67 | Recovery of unknown inputs and their entanglement |
| Channel definitions and their equivalent representations | §§2.2.1–2.2.3, starting pp. 72, 77, and 91 | The supplied noisy boxes and their Choi rank |
| Measurement effects and selected quantum operations | §§2.3.1–2.3.2, starting pp. 100 and 105; instruments on pp. 111–113 | The successful control outcome and its probability |
| State overlap and errors including a reference | §3.2, p. 139; §3.3.2, starting p. 166 | Fidelity conventions and conditional control errors |
| Random unitary descriptions of depolarization | §§4.1.1–4.1.2, starting pp. 202 and 212 | The Pauli worked example and the Weyl-channel comparison |
| Entangled pairs and classical communication | §§6.1.2, 6.2.2, and 6.3.1, starting pp. 324, 345, and 359 | Fresh-pair attempts and the use of an accepted pair |
| Positive-matrix optimization | §1.2.3, semidefinite programming | Exact feasible controller effects and operator upper bounds |
| Schur maps | §4.1.3 | Physical controller dephasing as an entrywise multiplier |
| Separability and partial transpose | §§6.1.1–6.1.2 and §6.2.3 | The sister comparison’s different entanglement cuts |
| Probability and expected resource use | §1.2.2, starting p. 47 | Independent attempts and charged channel calls |
| Capacity as a separate asymptotic task | §8.2.1, starting p. 512, when needed | Why one heralded branch does not establish unassisted quantum capacity |

The main representation results are Theorem 2.22 on page 82, Corollary 2.23 on page 84, and Theorem 2.26 on page 87. They respectively supply complete-positivity tests, freedom to change a Kraus representation, and trace-preservation tests. Their roles are explained below. They do not assume or construct a coherent permutation of unknown channel boxes.

<a id="notation"></a>

## 2. Translating notation and normalization

### Basic dictionary

The target register is $`T`$, with state space $`\mathbb C^d`$. A state $`\rho`$ is a positive matrix of trace one. An untouched register $`R`$ can hold a reference, so the input may be a joint state $`\rho_{RT}`$ rather than a state of $`T`$ alone.

| Watrous notation or convention | Repository notation or convention |
|---|---|
| $`\mathrm L(\mathcal X)`$ and $`\mathrm D(\mathcal X)`$ | Operators $`\mathcal B(T)=M_d(\mathbb C)`$ and density matrices $`\rho`$ |
| Adjoint $`A^*`$ | Adjoint $`A^\dagger`$ |
| Basis vector $`e_i`$, outer product $`uv^*`$ | $`\lvert i\rangle`$, $`\lvert u\rangle\langle v\rvert`$ |
| Operator identity $`1_{\mathcal X}`$, map identity $`1_{\mathrm L(\mathcal X)}`$ | $`I_d`$, $`\mathrm{id}_d`$ |
| $`\mathrm{vec}(K)`$ | $`\lvert K\rangle\!\rangle=(K\otimes I)\lvert\Omega_d\rangle`$ |
| $`J(\mathcal N)`$ | $`J_{\mathcal N}`$, with the same unnormalized, output-first convention |
| Completely depolarizing channel $`\Omega`$ | Completely depolarizing channel $`\Delta`$ |
| Completely dephasing channel $`\Delta`$ | Dephasing is written out explicitly to distinguish it from the repository's $`\Delta`$ |
| Completely bounded trace norm, written with triple bars | Diamond norm $`\lVert\mathcal E\rVert_\diamond`$ |

The distinction between an operator and a map matters: $`I_d`$ multiplies state vectors, whereas $`\mathrm{id}_d`$ sends an input matrix to itself. Likewise, the project's $`|\Omega_d\rangle`$ is a vector, not Watrous's depolarizing map $`\Omega`$.

### Choi matrices are not normalized states here

Watrous defines the Choi representation in Eq. (2.64), page 78. With this repository's symbols,

```math
|\Omega_d\rangle=\sum_{i=1}^d|i\rangle|i\rangle,
\qquad
J_{\mathcal N}
=(\mathcal N\otimes\mathrm{id})
  (|\Omega_d\rangle\langle\Omega_d|).
```

The vector has squared norm $`d`$. A trace-preserving channel therefore has $`\mathrm{Tr}J_{\mathcal N}=d`$, not one. Its normalized Choi state is $`J_{\mathcal N}/d`$. Two useful checks are

```math
J_{\mathrm{id}_d}=|\Omega_d\rangle\langle\Omega_d|,
\qquad
J_\Delta=\frac{I_{d^2}}d,
\qquad
J_{p\,\mathrm{id}_d}=p|\Omega_d\rangle\langle\Omega_d|.
```

Thus a successful identity branch has probability $`p=\mathrm{Tr}J_{\mathcal E}/d`$. Omitting the divisor would give the wrong probability. The complete [Choi conventions](../science/model.md#6-choi-and-full-rank-conventions) are fixed in the model page.

Full Choi rank also has a specific meaning. It is **not invertibility of the channel as a linear transformation**. The depolarizer sends all normalized inputs to the same state and nevertheless has Choi rank $`d^2`$. The role of full Choi rank is to make the Kraus matrices span the entire matrix algebra.

### Fidelity, distance, and selected outcomes

Watrous's fidelity in Definition 3.11, page 139, is the root convention,

```math
F(\rho,\sigma)=\|\sqrt\rho\sqrt\sigma\|_1.
```

The repository's entanglement fidelity relative to identity is an overlap. For the normalized vector $`|\Phi_d\rangle=|\Omega_d\rangle/\sqrt d`$ and output $`\tau=(\mathrm{id}\otimes\mathcal C)(|\Phi_d\rangle\langle\Phi_d|)`$,

```math
F_e(\mathcal C)
=\langle\Phi_d|\tau|\Phi_d\rangle
=F(|\Phi_d\rangle\langle\Phi_d|,\tau)^2.
```

The square is necessary when moving between the two conventions. The diamond norm is Watrous's completely bounded trace norm, Definition 3.43 on page 170. It measures a map's effect while allowing an untouched reference. A factor $`1/2`$ is included only when explicitly written, as in the [classical-order distance bound](../science/interpretation.md#a-quantitative-gap-from-the-identity).

An accepted operation is kept **unnormalized** until its probability is taken. If $`\mathcal E(\rho)`$ is the accepted output, its trace is the acceptance probability and its conditional state is $`\mathcal E(\rho)/\mathrm{Tr}\mathcal E(\rho)`$. This is the instrument description in §2.3.2. For $`\mathcal E=p\,\mathrm{id}`$, the probability is input independent and the conditional output is exactly $`\rho`$. For an imperfect implementation, that probability can depend on the input; conditioning need not define a linear channel.

<a id="conceptual-gaps"></a>

## 3. The missing steps from ordinary channels to the theorem

### A channel representation does not grant environment access

A channel can be written

```math
\mathcal N(\rho)=\sum_a K_a\rho K_a^\dagger,
\qquad
\sum_a K_a^\dagger K_a=I.
```

One way to interpret this is to let a larger interaction correlate the target with an environment label $`a`$ and then discard the environment. The sum describes what remains when that label is unavailable. It does not grant the experimenter a measurement of $`a`$, a controlled use of $`K_a`$, or the ability to reverse the larger interaction.

The same channel admits different Kraus lists. Corollary 2.23 relates equal-length lists by a unitary mixing; lists can first be padded with zero operators when necessary. Consequently a proposed operation on ordinary channels must not acquire a different physical meaning merely because one uses a different Kraus list to calculate it.

The project supplies coherent **order control** as an additional resource. Each slot has its own inaccessible environment, and each branch uses that slot once. The same environment label for a slot is used consistently across all its order branches. These assumptions are why the complete switch is invariant under independent Kraus changes. The [switch normalization and invariance proof](../science/polynomial_realization.md#1-the-normalized-switch) establishes this explicitly.

### Coherence combines amplitudes within a fixed noise history

For two slots, let $`|12\rangle`$ label the matrix product $`K^{(1)}_aK^{(2)}_b`$ and $`|21\rangle`$ the reversed product. Matrix multiplication acts on the input from right to left, so word $`12`$ means that slot 2 is traversed first.

Preparing the sum of the two order labels and postselecting their difference gives

```math
|\alpha\rangle=\frac{|12\rangle+|21\rangle}{\sqrt2},
\qquad
|\beta\rangle=\frac{|12\rangle-|21\rangle}{\sqrt2},
\qquad
B_{ab}=\frac{K^{(1)}_aK^{(2)}_b-K^{(2)}_bK^{(1)}_a}{2}.
```

The subtraction occurs between two orders of the **same** pair of noise labels. The accepted channel still sums $`B_{ab}\rho B_{ab}^\dagger`$ over the unobserved pairs. It does not take a coherent sum over $`a,b`$.

A classical random choice of order adds the corresponding output density matrices with nonnegative weights. It has no cross terms with which to make this subtraction. This is the physical distinction developed in the [classical-order comparison](../science/interpretation.md#classical-order-comparison). A commutator alone need not preserve the target, so further algebra is required.

### Final control measurements produce order polynomials

With $`m`$ slots, each order is a word containing every slot label once. A controller can therefore produce a polynomial

```math
P(X_1,\ldots,X_m)=\sum_{\pi\in S_m}c_\pi
X_{\pi_1}\cdots X_{\pi_m}.
```

The variables represent separate slots. Even when two slots contain the same channel, their environment labels are independent. The polynomial is *multilinear*: it is linear in any one matrix argument when the others are held fixed.

The coefficients arise from a preparation amplitude multiplied by the conjugate of a measurement amplitude. They need physical normalization. If $`L=\sum_\pi|c_\pi|`$, the construction uses

```math
\alpha_\pi=\sqrt{|c_\pi|/L},
\qquad
\beta_\pi=e^{-i\arg(c_\pi)}\sqrt{|c_\pi|/L},
\qquad
B_{\boldsymbol a}=\frac{P(K^{(1)}_{a_1},\ldots,K^{(m)}_{a_m})}{L}.
```

Zero coefficients correspond to unused order labels. The [realization proof](../science/polynomial_realization.md#2-every-multilinear-order-polynomial-is-realizable) verifies the valid preparation and measurement. The [normalization bound](../science/probability.md#optimality-for-a-fixed-coefficient-ray) concerns this fixed coefficient representation and one rank-one accepted outcome; it is not a global optimization of success probability.

### Centrality makes every accepted noise history harmless

A central polynomial has a scalar-matrix value for every matrix tuple:

```math
P(X_1,\ldots,X_m)=\lambda(X_1,\ldots,X_m)I_d.
```

Here the polynomial has positive degree, and it must be nonzero on at least one tuple. This last condition distinguishes a useful central polynomial from a polynomial identity that vanishes on all tuples. The scalar $`\lambda`$ may vary with the tuple and may vanish on many tuples.

Each accepted history now has $`B_{\boldsymbol a}=b_{\boldsymbol a}I_d`$. The environment may retain its history, but all surviving histories act identically on the target up to a scalar. Summing their contributions gives

```math
\mathcal E_{\mathrm{succ}}(\rho)
=\sum_{\boldsymbol a}|b_{\boldsymbol a}|^2\rho
=p\rho.
```

This is equality of maps, so it preserves correlations with every reference as well. No environmental record has to be measured or erased by the controller.

### Full Choi rank guarantees a surviving history and enables the converse

The Kraus and Choi representations give

```math
J_{\mathcal N}
=\sum_a|K_a\rangle\!\rangle\langle\!\langle K_a|.
```

This positive matrix has range equal to the span of the Kraus vectors. Full rank therefore means that every matrix can be expressed as a linear combination of the Kraus matrices. If a multilinear polynomial vanished on every independent Kraus tuple, expanding arbitrary matrix arguments in those spanning families would make it vanish everywhere. A nonidentity central polynomial must consequently have at least one nonzero successful tuple. This is the positivity mechanism in [Theorem 1](../science/polynomial_realization.md#theorem-1-universal-exact-branch).

For the opposite direction, an exact identity branch has a rank-one Choi matrix. A positive sum of Kraus-vector outer products can have that one-dimensional support only if every contributing vector lies in it. Thus every accepted Kraus matrix is scalar. The same spanning argument extends scalarity from the channel's Kraus tuples to all matrix tuples, forcing a central polynomial.

This yields [Theorem 2](../science/query_threshold.md#theorem-2-exact-query-characterization):

```math
m_\star(\mathcal N)=c(d)
\qquad\text{for every full-Choi-rank }\mathcal N.
```

Here $`m_\star`$ is the minimum positive number of compulsory calls, and $`c(d)`$ is the minimum positive multilinear central degree. The full proof includes mixed controls, ancillary control registers, general accepted effects, and final herald-dependent unitary corrections. It allows the lower-bound controller to know the particular channel. These extensions belong to the canonical proof, rather than being assumptions that a pure-control example can silently replace.

The conclusion is architecture specific. The target cannot be encoded into protected memory, bypass a slot, or undergo inserted target gates. The controller starts independently of the target. Correlated device uses also need a different joint description. These conditions are part of the [model](../science/model.md#4-target-access-and-explicit-exclusions), not consequences of the textbook's representation theorems.

<a id="worked-example"></a>

## 4. Complete worked example: four depolarizing qubit slots

### The channel and the supplied resource

Write the Pauli matrices as

```math
\sigma_0=I_2,\qquad
\sigma_x=\begin{pmatrix}0&1\\1&0\end{pmatrix},\qquad
\sigma_y=\begin{pmatrix}0&-i\\i&0\end{pmatrix},\qquad
\sigma_z=\begin{pmatrix}1&0\\0&-1\end{pmatrix}.
```

For every normalized qubit state,

```math
\Delta_2(\rho)=\frac14\sum_{\mu\in\{0,x,y,z\}}
\sigma_\mu\rho\sigma_\mu=\frac{I_2}{2}.
```

Use four independent slots, each with Kraus list $`K_\mu=\sigma_\mu/2`$. The supplied generalized switch coherently routes the target through each slot exactly once. The controller starts independently of the target and any retained reference. Nothing measures the four noise labels.

### Preparation, orders, and measurement

Use the polynomial

```math
P_2(X_1,X_2,X_3,X_4)
=[X_1,X_2][X_3,X_4]+[X_3,X_4][X_1,X_2],
\qquad [A,B]=AB-BA.
```

Its eight words divide into four positive and four negative coefficients:

| Coefficient | Matrix-product words |
|---|---|
| $`+1`$ | $`1234,\ 2143,\ 3412,\ 4321`$ |
| $`-1`$ | $`1243,\ 2134,\ 3421,\ 4312`$ |

Call this set $`S`$ and its signs $`c_\pi`$. Prepare and measure

```math
|\alpha\rangle=\frac1{\sqrt8}\sum_{\pi\in S}|\pi\rangle,
\qquad
|\beta\rangle=\frac1{\sqrt8}\sum_{\pi\in S}c_\pi|\pi\rangle,
\qquad
F_{\mathrm{succ}}=|\beta\rangle\langle\beta|.
```

The other measurement effect is $`I_C-F_{\mathrm{succ}}`$. Both effects are positive and sum to the control identity. For each tuple $`\boldsymbol\mu=(\mu_1,\mu_2,\mu_3,\mu_4)`$, the accepted target Kraus operator is

```math
B_{\boldsymbol\mu}
=\frac{P_2(K_{\mu_1},K_{\mu_2},K_{\mu_3},K_{\mu_4})}{8}.
```

There are eight occupied order states, which fit in a three-qubit register. This is a description of the ideal controller and routing resource, not a gate count for constructing the router.

### One history and then the complete sum

Take $`\boldsymbol\mu=(x,y,x,y)`$. Since $`[\sigma_x,\sigma_y]=2i\sigma_z`$ and $`\sigma_z^2=I`$,

```math
P_2(\sigma_x,\sigma_y,\sigma_x,\sigma_y)=-8I,
\qquad
P_2(K_x,K_y,K_x,K_y)=-\frac12I,
\qquad
B_{xyxy}=-\frac1{16}I.
```

The middle factor follows from multilinearity: four factors of $`1/2`$ contribute $`1/16`$. The final factor of $`1/8`$ comes from the controller, not the channel noise. This one accepted Kraus contribution has weight $`1/256`$.

To sum all histories, a Pauli commutator is zero if either entry is identity or if its two nonidentity entries are equal. Otherwise it is $`\pm2i`$ times the remaining Pauli axis. The anticommutator of two different Pauli axes is zero, while that of an axis with itself is $`2I`$. Therefore a tuple survives exactly when its two ordered pairs produce the same commutator axis.

There are three choices of that axis, two orders of the other two axes for the first pair, and two orders for the second pair. Exactly $`3\cdot2\cdot2=12`$ of the $`4^4=256`$ tuples survive. Each gives $`B_{\boldsymbol\mu}=\pm I/16`$. Consequently

```math
\mathcal E_{\mathrm{succ}}(\rho)
=12\frac{\rho}{256}
=\frac3{64}\rho.
```

This counts contributions in a Kraus representation; it does not assume that the controller observes a random Pauli label. The result is a channel statement and is unchanged by changing the Kraus representation.

### Output and costs

For a fresh entangled input pair, acceptance returns the same pair with probability $`3/64`$. Independent fresh-pair attempts therefore use an average of $`64/3`$ attempts and $`256/3`$ channel calls per accepted pair, before routing and control costs. A failed trial has probability $`61/64`$ and need not preserve the input. An unknown message cannot simply be reused after failure.

Every classical order of these same slots gives $`\Delta_2^4=\Delta_2`$. Conditioning on an independent classical order record does not restore the input. The successful quantum branch thus has a different exact action, even though it is rare.

The general centrality proof, the inherited superswitch attribution, and the exact certificate ruling out one, two, or three calls are in [the qubit account](../science/qubit.md). The example computes one explicit protocol's success probability. It does not establish that this probability is the largest possible.

<a id="resource-scaling"></a>

## 5. Minimum calls and expected calls answer different questions

The query theorem minimizes calls **within one successful-capable attempt**. It does not minimize the expected number of charged calls after failures. The two explicit minimum-call constructions illustrate the difference:

| Target | Calls per attempt | Occupied order labels | Success on complete depolarization | Mean calls per accepted fresh pair |
|---|---:|---:|---:|---:|
| Qubit | 4 | 8 | $`3/64`$ | $`256/3`$ |
| Qutrit | 8 | 2880 | $`1/43740`$ | 349920 |

These are the probabilities of the specified controllers. The [qutrit construction](../science/qutrit.md) gives the compact polynomial and its centrality proof; the [probability calculation](../science/probability.md#eight-call-qutrit-probability) supplies two exact contractions yielding $`1/43740`$.

### Why less noise can make a fixed exact filter succeed less often

Consider independent slots containing

```math
\mathcal N_\epsilon=(1-\epsilon)\mathrm{id}+\epsilon\Delta,
\qquad 0<\epsilon\le1.
```

At the minimum degree $`c(d)`$, substituting identity for any one variable of a central polynomial must make it vanish. A nonzero remaining polynomial would otherwise be a central polynomial of lower degree. In the identity-mixture channel, the minimum-degree filter therefore accepts only the all-depolarizing part of the expansion. Its success scales as $`\epsilon^{c(d)}`$.

For the worked qubit controller this gives $`3\epsilon^4/64`$. The exact qutrit controller gives $`\epsilon^8/43740`$. These laws describe a filter that demands **exact** recovery. They do not say that ordinary transmission through a nearly noiseless channel is poor. At $`\epsilon=0`$, the channel itself is identity and the full-rank hypothesis no longer applies. The [small-noise proof](../science/probability.md#noise-approaching-the-identity) establishes the law for every minimum-call exact controller and its optimized supremum.

### Longer contiguous blocks improve the charged cost

The same order controller can treat a fixed sequence of $`k`$ supplied channels as one block. Every order visits all channels within each block in the same internal order and permutes the complete blocks. No additional operation acts on the target. For identity-mixture noise,

```math
\mathcal N_\epsilon^k=\mathcal N_\delta,
\qquad
\delta=1-(1-\epsilon)^k.
```

For example, the qubit controller with four such blocks has $`4k`$ calls, the same eight active block orders, and success $`3\delta^4/64`$. Choosing $`k=\lceil1/\epsilon\rceil`$ makes $`\delta\ge1-e^{-1}`$. Success is then bounded below by a positive constant, while the number of charged calls grows as $`1/\epsilon`$.

The [fixed-batch resource theorem](../science/probability.md#expected-charged-calls) extends this reasoning to each fixed dimension and proves a matching lower scaling:

```math
\mathcal Q_d(\epsilon)
=\inf_{m,\mathcal C:\ p_{\mathcal C}>0}\frac{m}{p_{\mathcal C}}
=\Theta_d(\epsilon^{-1}).
```

The infimum is over exact identity controllers in the unaugmented order-only model, with a fixed number $`m`$ of calls per attempt. The batch size can depend on the known noise parameter. The constants can be large and depend on dimension; this does not optimize their prefactors or charge the hardware cost of coherent paths. Its lower bound also does not extend to the [augmented model with randomizing wrappers](../science/interpretation.md#randomizing-wrappers).

<a id="proof-route"></a>

## 6. Where the complete arguments live

The bridge explains the conceptual transitions. Each scientific claim has an authoritative proof or source account elsewhere in the repository.

| Question | Canonical route |
|---|---|
| Exactly what may the controller and target do? | [Model, Sections 1–5](../science/model.md) |
| Why is the switch a well-defined channel independent of Kraus coordinates? | [Normalization and invariance](../science/polynomial_realization.md#1-the-normalized-switch) |
| How do coefficients become a physical preparation and accepted effect? | [Order-polynomial realization](../science/polynomial_realization.md#2-every-multilinear-order-polynomial-is-realizable) |
| Why does full rank ensure nonzero exact success? | [Theorem 1](../science/polynomial_realization.md#theorem-1-universal-exact-branch) |
| Why can no allowed controller succeed at a lower degree? | [Theorem 2 and its general-instrument proof](../science/query_threshold.md) |
| What establishes the four- and eight-call values? | [Qubit construction and finite certificate](../science/qubit.md); [qutrit construction and imported minimality](../science/qutrit.md) |
| How do probability, noise strength, and total charged calls interact? | [Probability and resources](../science/probability.md) |
| What changes for classical orders, correlated noise, or rank-deficient channels? | [Operational interpretation and counterexamples](../science/interpretation.md) |
| What if the controller is imperfect? | [General conditional-error estimate](../science/interpretation.md#ideal-control-and-conditional-error); [physical noisy-routing model](../science/noisy_routing_model.md#shared-access) |
| What is the exact domain of all accepted controller instruments? | [Controller feasibility and attaining real representatives](../science/controller_optimization.md#controller-feasibility) |
| Which controller noise excludes finite exact recovery? | [Strict Schur-floor obstruction](../science/noisy_recovery_obstruction.md#strict-schur-obstruction) |
| Which readout attains the complete noisy recovery frontier? | [Orthogonal frontier](../science/recovery_frontier.md#frontier-theorem); [certified tilted interval](../science/tilted_readout.md#tilted-optimality) |
| How does this connect to retained entanglement? | [Specified-protocol uniform separation](https://github.com/GoGoKo699/Quantum-Order-Entanglement/blob/main/science/noisy_separation.md#uniform-separation) |
| Which computations can be reproduced, and what do they establish? | [Verification guide](../verification/README.md) |

The control-error distinction is especially relevant when success is small. A small error in the unnormalized accepted operation can become larger after conditioning. The bound on the conditional trace norm error in the linked proof is $`2\delta/(p-\delta)`$ when the operation error is at most $`\delta<p`$ in diamond norm. It is an input-by-input bound including references, rather than a diamond distance between generally nonlinear normalized transformations.

<a id="source-boundaries"></a>

## 7. Textbook foundations and research dependencies

Watrous is the single educational anchor. The primary research sources remain necessary because several ingredients are outside the ordinary-channel foundations used here.

| Layer | What supports it |
|---|---|
| States, positive operators, channels, Kraus and Choi representations, instruments, and channel distances | The selected Watrous sections above |
| The generalized switch as a supplied higher-order resource | The [model](../science/model.md) and [published generalized-switch reference](../references/RELATED_WORK.md#kechrimparis2025) |
| Central-polynomial resetting as an antecedent and the inherited qubit branch | The [Navascués and superswitch comparisons](../references/RELATED_WORK.md#navascues2018) |
| Existence and minimum degrees of central polynomials | The [algebraic source account](../references/RELATED_WORK.md#central-polynomials) |
| Explicit qutrit centrality | Drensky–Kasparian's construction, with the [complete local proof](../science/qutrit.md#2-centrality-and-its-scalar-value) |
| The qutrit lower-degree exclusion | The historical theorem reported in Drensky–Kostadinov; its [source-access qualification](../references/RELATED_WORK.md#central-polynomials) remains explicit |
| General realization, fixed-channel converse, and fixed-batch charged-cost scaling | The canonical project proofs linked above |
| Regev's determinant evaluation and its probability | The [published algebraic formula](../references/RELATED_WORK.md#regev-procesi) and [channel probability derivation](../science/probability.md#regev-construction) |

In particular, $`c(2)=4`$ and $`c(3)=8`$ are exact algebraic degrees, whereas $`c(4)\le13`$ is an upper bound. The dimension-four equality is not established here. The query characterization transfers the algebraic minimum to the channel model; it does not prove new central-polynomial identities or imply feasibility at every larger compulsory call count.

Verification scripts check coefficients, exact witnesses, contractions, normalization, and numerical channel implementations. The qubit lower-degree certificate is a finite algebraic proof certificate. Sampled matrix evaluations are not proofs of universal centrality, and the qutrit computations do not reproduce the historical minimality proof. The [verification guide](../verification/README.md) distinguishes these roles.

The ideal-control route goes from channel representations to controlled order amplitudes, from those amplitudes to central polynomials, and from rank-one Choi support back to the same algebraic condition. The noisy-control route below adds a complete feasible instrument domain and tight operator bounds. The access model and resource account determine both routes’ physical conclusions.

<a id="noisy-control-bridge"></a>

## 8. From instruments and Schur maps to a noisy recovery frontier

**Milestones M22–M26.** The textbook’s instruments and semidefinite-programming language meet in a concrete question: how much accepted Bell fidelity can *any* permitted controller obtain at a chosen acceptance probability? A selected herald gives a lower bound. To establish an optimum, the calculation must include arbitrary mixed preparations, ancillary controllers, collective effects and outcome-dependent final target unitaries.

### Fix the physical register before specifying its noise

Let the raw qubit slot measure the target’s Z axis and prepare the two opposite states of a tilted axis, with a small depolarizing admixture:

```math
\mathcal N_{\epsilon,\theta}
=(1-\epsilon)\operatorname{Ad}_{R_y(\theta)}\circ\mathcal D_Z
+\epsilon\Delta,
\qquad \theta=\pi/2+\delta.
```

Here $`\mathcal D_Z`$ deletes target off-diagonal entries, whereas $`\Delta(X)=\operatorname{Tr}(X)I_2/2`$ replaces the target by the maximally mixed state. The slot is entanglement breaking; positive $`\epsilon`$ gives full Choi rank. A fixed target unitary after the complete experiment is a correction. A unitary attached to every supplied slot would change the resource.

The eight matrix words in lexical order are $`1234,1243,2134,2143,3412,3421,4312,4321`$. Their physical three-bit codes are respectively **000, 001, 010, 011, 100, 110, 101, 111**. The inner reversal bits stay attached to the same channel pairs when the outer bit exchanges the two blocks. Replacing this assignment by lexical binary counting changes the physical local-noise model.

Watrous §4.1.3 supplies the Schur-map viewpoint: a matrix element between physical order labels is multiplied entrywise by a positive kernel. Global dephasing has off-diagonal factor $`\eta`$. Independent phase flips of probability $`z`$ on each physical bit instead give

```math
\Gamma_{\pi\sigma}=\gamma^{d_H(b_\pi,b_\sigma)},
\qquad \gamma=1-2z.
```

Thus $`\gamma=4/5`$ means a 10% phase-flip probability on each bit. It is a different noise process from global dephasing with $`\eta=4/5`$. The [physical codebook and common family](../science/noisy_routing_model.md#physical-codebook) fix these definitions fully.

### An instrument has a shared coefficient budget

An accepted outcome is represented by a positive matrix $`H_x`$ on the eight order labels. All outcomes must fit one initial preparation:

```math
H_x\succeq0,\qquad
\sum_x H_x\preceq\operatorname{diag}(t),
\qquad t_\pi\ge0,\quad\sum_\pi t_\pi=1.
```

The vector $`t`$ gives preparation populations. The bound accounts for mutually compatible effects of one instrument, rather than allowing every outcome its own independent preparation budget. The [feasibility theorem](../science/controller_optimization.md#controller-feasibility) proves both directions: every allowed controller gives such matrices, and every feasible family can be implemented. Mixed preparations and ancillary controllers do not enlarge this set. Symmetry of these identical slots and their physical noise permits uniform preparation and the budget $`I_8/8`$ without losing an optimum.

Acceptance and the unnormalized corrected Bell overlap are linear functionals of $`H_x`$. Their Gram matrices encode all raw Kraus-history contractions. For a general complex Schur kernel, the appropriate entrywise contraction uses its transpose; this convention is fixed in the scientific derivation. The physical kernels here are real and symmetric.

The Bell input is an optimization probe, not a reason to assume every branch has input-independent success. Arbitrary complex $`H_x`$ can produce input-dependent acceptance. For this real, Y-covariant family, the [real-controller reduction](../science/controller_optimization.md#real-controller-reduction) proves that the optimum has a real attaining representative, with a scalar success effect. Its accepted normalization therefore defines a recovery channel on arbitrary inputs. This is an attainment theorem, not a classification of every feasible controller.

### A readout frontier becomes fractional knapsack

At the orthogonal point $`\delta=0`$, the acceptance and relevant corrected-reward matrices are diagonal in the physical Walsh basis. This is the basis measured by three physical X measurements. Let $`b_j`$ be a mode’s acceptance coefficient and $`a_j`$ its best corrected Bell-reward coefficient. If $`x_j`$ is its accepted coefficient weight, the complete frontier is

```math
\max\ pF=\sum_j a_jx_j,
\qquad
0\le x_j\le\frac18,
\qquad \sum_j b_jx_j=p.
```

Each mode has acceptance capacity $`b_j/8`$ and fidelity $`a_j/b_j`$ when $`b_j>0`$. Accept modes in descending fidelity, randomizing acceptance of the final boundary mode when needed; zero-capacity modes contribute no acceptance. An operator upper bound covers every allowed controller, and this simple readout attains it. The [frontier theorem](../science/recovery_frontier.md#frontier-theorem) supplies the coefficient tables and proof.

For example, at $`\epsilon=.01`$ and physical $`\gamma=.8`$, the highest-fidelity mode accepts with probability approximately **.06196376409** and corrected Bell fidelity **.758264362729**. Increasing acceptance to **.1245025** gives optimal fidelity approximately **.688117778572**. These are two points on one tradeoff curve, not conflicting controller recommendations. The deterministic optimum is approximately **.47329633218**. A Bell-fidelity optimum below one half does **not** imply that every possible complex-controller Choi output is separable.

The inherited P2 branch and the alternative Q branch illustrate why the physical kernel matters. Their fidelity ranking reverses between the stated global and local noise experiments, even under a matched scalar-coherence comparison. The [exact branch formulas](../science/recovery_frontier.md) retain both fidelity and acceptance, so a claim about a favorable herald is not mistaken for the full frontier or for a cost-free accepted output.

<a id="noise-obstruction-bridge"></a>

## 9. Positive noise floors, tilted optimality and the sister connection

**Milestones M24 and M27–M30.** The same positivity foundation supports three distinct arguments.

First, if the Schur kernel has a strictly positive identity floor, every nonzero accepted effect leaves a nonzero diagonal contribution from ordinary composed slot channels. For full-Choi-rank slots that contribution has full Choi rank. A positive map containing it cannot have the rank-one Choi operator of an exact unitary branch. The [strict Schur obstruction](../science/noisy_recovery_obstruction.md#strict-schur-obstruction) gives the precise finite-call theorem. Global dephasing with $`\eta<1`$ and physical local dephasing with $`\gamma<1`$ have such floors. This proof concerns unobserved phase noise on the actual order register, with no noise flag or coherence restoration; it is not a claim about every type of controller noise.

Second, the optimum survives a finite target-axis tilt even after the reward matrices stop commuting. A diagonal matrix made from the best individual mode rewards can still dominate every corrected reward matrix. The [tilted proof](../science/tilted_readout.md#tilted-optimality) reduces this operator-envelope condition to two coupled 2-by-2 Walsh blocks. Exact polynomial identities and positive Bernstein coefficients establish the required margin. At $`\epsilon=1/100`$ and physical $`\gamma=4/5`$, the complete three-X readout frontier remains optimal throughout

```math
t=\tan(\delta/2),\qquad
|t|\le\frac15,
\qquad |\delta|\le2\arctan(1/5)\simeq .39479112.
```

The angle $`\delta`$ is in radians. The conclusion optimizes the final correction separately for every outcome and angle. Holding the orthogonal-point corrections fixed would not test the stated optimal readout.

Third, the sister project keeps the *same raw channel family and per-bit noise* but compares different outputs. Its specified three-call experiment retains C and T unconditionally, while the four-call P2 experiment measures C and accepts a Bell pair across R:T. Exact rational grid certificates plus a Stinespring continuity bound prove a uniform separation on the smaller angular interval $`|\delta|\le .2`$ radians. This interval is stated in $`\delta`$, not in $`t`$. The proof and its full resource comparison live in [Entanglement’s uniform theorem](https://github.com/GoGoKo699/Quantum-Order-Entanglement/blob/main/science/noisy_separation.md#uniform-separation).

The three-call occupied 3-by-2 output is PPT and therefore separable for all three physical encodings. The four-call accepted R:T output is NPT, with $`.035<p<.095`$ and $`F>101/190`$. Calls, controller bits, bipartitions and conditional endpoints differ. The theorem concerns those experiments, not every retained protocol or an optimized equal-resource comparison. It establishes a connection between the tracks without making retained entanglement a recovery certificate.

Watrous supplies positive operators, instruments, Schur maps, norms and separability language. The controller-feasibility theorem, tight frontiers, envelope certificate and rational-grid interpolation are the missing project arguments. The [verification guide](../verification/README.md) distinguishes independent exact reconstruction, interval certificates and numerical evaluations from the analytical proofs that connect them.
