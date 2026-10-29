# Scientific overview

A noisy channel leaves information in an inaccessible environment. Choosing a classical order for independent uses generally cannot remove that noise exactly. A coherent order control changes the available amplitudes: for each fixed environment history, its successful effect combines different products of Kraus operators. This raises a precise question: **when can every successful history act as a scalar on an unknown target?**

A central polynomial supplies exactly that condition. Its value is a scalar matrix for every tuple of matrices, even when the individual matrices do not commute. The same identity can therefore apply separately to every inaccessible noise history. The environment labels are still summed incoherently; within each fixed label tuple, the orders interfere.

This account establishes the resulting channel capability and its exact minimum-use characterization. It also records the probability cost, an explicit qutrit construction, and the assumptions under which the characterization applies. A second connected argument determines noisy-controller recovery frontiers and their relation to retained entanglement.

## The argument in order

1. **Fix the physical access.** The [model](model.md) supplies a generalized switch of independent ordinary CPTP channel slots. Every branch traverses every slot once. The target has no protected bypass, and only the initially independent controller is prepared and measured.
2. **Realize an order polynomial.** The [construction](polynomial_realization.md) turns arbitrary multilinear coefficients into normalized control amplitudes. Centrality makes every successful Kraus operator scalar. Full Choi rank supplies the spanning property that makes at least one scalar nonzero.
3. **Prove the converse.** The [query threshold](query_threshold.md) uses positivity of a rank-one successful Choi operator. Each contributing Kraus vector must lie in the identity line. Spanning and multilinearity then force centrality on all matrix tuples, even for a controller tailored to one known channel.
4. **Make the first nontrivial cases explicit.** The [qubit account](qubit.md) gives the four-call commutator polynomial and an exact finite lower-degree certificate. The [qutrit account](qutrit.md) gives all 2,880 words of the eight-call Drensky–Kasparian polynomial through two independent coefficient constructions, a centrality proof, and a nonzero witness.
5. **Charge success and failure.** The [probability account](probability.md) computes the chosen branches' success probabilities and separates calls per attempt from expected calls per heralded pair. It proves the small-noise suppression law and matching fixed-dimension expected-call scaling.
6. **Interpret the boundary.** The [physical interpretation](interpretation.md) gives the operational use, a classical-order comparison, and examples showing why independence and the target-access restrictions matter.

## Main claims and proof dependencies

| Claim | Assumptions doing the work | Proof or construction |
|---|---|---|
| A dimension-dependent controller heralds an exact identity branch for arbitrary full-rank slots | Independent channels; multilinear central polynomial; full Kraus spans | [Universal exact branch](polynomial_realization.md#3-centrality-gives-exact-recovery) |
| For a fixed full-rank channel, the minimum positive number of calls is $`c(d)`$ | Control-only processing; compulsory traversals; rank-one successful Choi operator | [Exact query characterization](query_threshold.md#theorem-2-exact-query-characterization) |
| The minima are four calls for qubits and eight for qutrits | Characterization plus inherited central-degree facts | [Known degrees](query_threshold.md#5-known-degrees-and-their-quantum-meaning) |
| General call counts obey $`2d-1\le c(d)\le(d-1)^2+4`$ for $`d>2`$ | Matrix-unit lower bound; inherited central-polynomial construction | [Dimension bounds](query_threshold.md#7-a-dimension-dependent-lower-bound) |
| Minimum-call success on $`(1-\epsilon)\mathrm{id}+\epsilon\Delta`$ is suppressed by $`\epsilon^{c(d)}`$ | Minimal central degree; independent mixture expansion | [Noise approaching the identity](probability.md#noise-approaching-the-identity) |
| Optimized fixed-batch expected calls scale as $`\Theta_d(1/\epsilon)`$ | Fresh entangled-pair attempts; known or comparably bounded noise; no extra target wrappers | [Expected charged calls](probability.md#expected-charged-calls) |

The known algebraic degrees are not new polynomial-identity results. The characterization explains their operational meaning in this channel model. The qutrit lower-degree exclusion is an imported theorem; the repository separately proves and verifies the explicit degree-eight construction.

## Why the success probability remains visible

Exactness describes the normalized successful output, not how often that output occurs. The explicit qutrit protocol succeeds on completely depolarizing slots with probability $`1/43{,}740`$. Independent repetition on fresh entangled pairs therefore costs an average of 349,920 original channel calls per successful pair before routing and controller costs are charged.

Approaching the identity creates a second resource distinction. A minimum-call central polynomial vanishes whenever one argument is the identity. Its exact branch consequently becomes rare as the noise strength tends to zero. Larger batches can first combine contiguous noisy uses into a channel with a constant amount of noise, then apply the central-polynomial controller to the blocks. This yields the matching expected-call scaling in the probability account. At the exactly noiseless endpoint, the channel itself is an identity, and the full-rank converse no longer applies.

These statements explain capability and cost within one model. They do not estimate the engineering overhead of a generalized switch or establish a rate for an ordinary directed communication link.

## Relation to the literature

Central-polynomial resetting with returned quantum probes was already established by Navascués. Coherent cyclic orders already activate communication through complete depolarization and approach perfect heralded transmission with increasing uses. The qubit superswitch distillation branch is also prior work. The comparison here concerns **finite exact identity heralding from ordinary inaccessible-environment channel slots and its minimum compulsory traversal count**.

The [related-work account](../references/RELATED_WORK.md) compares the physical access and exact claims of these sources, including current versions. It distinguishes the order-polynomial construction from a coherent Kraus-oracle compiler, whose environment access is stronger.

## Noisy control: a second connected argument

The ideal exact branch and its minimum-call characterization do not settle useful recovery after controller coherence is lost. The noisy extension keeps the channel slots raw and fixes the physical controller encoding before optimization.

1. [Specify the shared family and physical noise](noisy_routing_model.md). The four-slot nested-switch architecture is inherited; its eight physical order codes determine local dephasing distances. Final target corrections are permitted, while per-slot wrappers and target interleavings are excluded.
2. [Characterize all physical controllers](controller_optimization.md#controller-feasibility). Positive accepted coefficient matrices obey one diagonal budget. Symmetry gives uniform preparation; correction-sector splitting and realification produce input-independent attaining representatives without falsely assigning that property to every controller.
3. [Locate the exactness boundary](noisy_recovery_obstruction.md#strict-schur-obstruction). Strict Schur noise leaves unavoidable positive diagonal acceptance. Independent full-rank slots make that component full Choi rank, excluding an exact unitary branch.
4. [Solve the orthogonal recovery frontier](recovery_frontier.md#frontier-theorem). Coefficient-wise Walsh identities and an optimized correction envelope reduce the whole acceptance–fidelity problem to a fractional knapsack, including boundary modes and deterministic acceptance.
5. [Prove tilted optimality](tilted_readout.md#tilted-optimality). A diagonal matrix envelope survives noncommuting overlap matrices. An exact degree-16 Bernstein certificate establishes its strict block margin throughout the stated interval.
6. [Compare the operational outputs](https://github.com/GoGoKo699/Quantum-Order-Entanglement/blob/main/science/noisy_separation.md#uniform-separation). The sister repository proves that specified retained-$`C:T`$ outputs can be separable while an accepted four-call $`R:T`$ output is entangled under the same per-bit noise. The protocols' distinct resources and conditioning are part of the comparison.

The shared optimization and certificate bundle is one common contribution across the two repositories. This repository owns recovery and its optimization proofs; the sister repository owns the retained-entanglement criterion and the full uniform separation proof. Neither the original central-polynomial threshold nor the order-support criterion alone implies this noisy comparison.

## The rank-deficient native benchmark

At the fixed measure-$`Z`$/prepare-$`X`$ endpoint, the channel has Choi rank two. The full-rank central-polynomial minimum consequently does not control its call count. The [native three-call frontier](native_three_call_frontier.md) instead uses every permutation of three raw slots and a global order-dephasing mask. Its four controller projectors provide a finite attaining protocol and an all-controller converse for accepted target-only half-diamond error.

The [randomized call-cost account](randomized_call_cost.md) supplies the complete two-call baseline proof and then convexifies all one-/two-/three-call strategies. It distinguishes raw attempt-selection probabilities from fractions of accepted outputs, proves the unequal logical-noise activation condition and checks an explicitly conditional duration law. The target-only output averages fine labels, and only aggregate acceptance is required to be input independent.

These are additional scoped benchmark results. A Bell deficit equals half-diamond error here because the attaining corrected mixtures have error Choi support orthogonal to the identity Bell state. That equality is proved for these maps; it is not transferred automatically to the four-call tilted Bell frontier or to a retained controller output.

## Verification route

Read the [verification guide](../verification/README.md) before interpreting a numerical result. The evidence has three distinct roles:

- The analytical proofs establish channel validity, centrality, positivity, the converse, and resource bounds.
- Exact finite calculations verify the coefficient transcription, nonidentity witness, qubit lower-degree certificate, and rational success probabilities.
- Numerical channel calculations check general nonunital instances, independent Kraus-isometry changes, reference preservation, and ancillary control decompositions.

The [root overview](../README.md) supplies the shortest entry route. For the prerequisites and their connection to this argument, use the [Watrous bridge](../learning/WATROUS_BRIDGE.md) with the [30-milestone map](../learning/MILESTONES.md). The [bibliography](../references/references.bib) supplies structured source metadata.
