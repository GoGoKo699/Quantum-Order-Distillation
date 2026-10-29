# Thirty milestones for quantum order distillation

This map follows one operational story: **how coherent order can herald an exact identity channel, what that capability costs, and which recovery tradeoffs remain attainable when the controller is noisy.** It connects the central-polynomial characterization to a complete controller optimization and a proved comparison with retained entanglement.

The sole educational anchor, shared with the sister project, is **John Watrous, *The Theory of Quantum Information* (2018)**. The [Watrous bridge](WATROUS_BRIDGE.md) identifies its version and relevant sections, translates notation, supplies worked explanations, and marks the steps developed in the project. Primary papers in [Related work](../references/RELATED_WORK.md) provide research attribution rather than additional educational anchors.

## Phase structure

| Phase | Milestones | Scientific purpose | Closing synthesis |
|---|---|---|---|
| Background/preliminaries | [M01](#m01)–[M13](#m13) | Join states, channel representations, instruments, independent ordered slots and matrix polynomials into the recovery problem. | [M13: The complete recovery specification](#m13) |
| Narrative/innovation | [M14](#m14)–[M30](#m30) | Build exact recovery and its costs, then controller-noise obstructions, attaining frontiers and the sister comparison. | [M30: The complete capability, frontier and resource statement](#m30) |

The **13 + 17 allocation** reflects the meeting of channel mathematics, polynomial identities and constrained optimization. The background proceeds from local objects to the complete model. The narrative combines proof components into physical capabilities and limits; its later milestones build on the ideal-control story rather than replacing it. The map is common to every reader and targets global understanding, while specialist calculations remain in their authoritative scientific pages.

## Navigation

| Milestones | Scientific connection |
|---|---|
| [M01](#m01)–[M06](#m06) | States, references, channels, positive Choi matrices and accepted operations |
| [M07](#m07)–[M13](#m13) | Independent ordered slots, polynomials, exactness and the complete task |
| [M14](#m14)–[M18](#m18) | Physical construction, full-rank converse, qubit and qutrit realizations |
| [M19](#m19)–[M21](#m21) | Success cost, batching and model boundaries |
| [M22](#m22)–[M27](#m27) | Physical noise, complete controller feasibility, no-go and attaining frontiers |
| [M28](#m28)–[M30](#m30) | Robust sister comparison, evidence and whole-project synthesis |

The [proof route](WATROUS_BRIDGE.md#proof-route) and [scientific overview](../science/README.md) supply direct specialist routes.

## Phase I: Background/preliminaries

<a id="m01"></a>

### M01. States and target matrices

**Phase:** Background/preliminaries.

**Central question:** Which matrices describe physical states, and which matrices must a channel identity handle?

**Understanding target:** Distinguish a density operator from a general linear operator on a $`d`$-dimensional target. Connect adjoints, positivity, trace, matrix units, and linear combinations. In particular, Kraus operators and polynomial arguments need not be density operators, Hermitian matrices, or unitary matrices.

**Dependencies and story role:** No earlier milestone. This supplies the local objects from which channels, noise histories, and polynomial evaluations will be built.

**Conceptual connection:** A matrix basis describes all possible operator amplitudes, while positivity and normalization select physical states. Confusing those roles would weaken the later centrality condition.

**Source pointers:** [Watrous Sections 1.1.2, 1.1.3, and 2.1.2](WATROUS_BRIDGE.md#source-map); [project target and channel spaces](../science/model.md#1-independent-ordinary-channel-slots).

<a id="m02"></a>

### M02. Entangled reference systems

**Phase:** Background/preliminaries.

**Central question:** Why must recovery preserve a target's correlations with a system that the protocol never touches?

**Understanding target:** Explain tensor-product systems, reduced states, and the role of an untouched reference. Relate the normalized maximally entangled state $`|\Phi_d\rangle`$ to preservation of an unknown quantum input. Distinguish the target's dimension $`d`$ from the dimension $`d^2`$ of its joint space with an equal-sized reference.

**Dependencies and story role:** [M01](#m01). Enlarges the local target description to the correlated inputs needed for a complete quantum task.

**Conceptual connection:** A protocol that preserves selected target states can still damage entanglement. Channel-level recovery includes those correlations from the start.

**Source pointers:** [Watrous Section 2.1.3](WATROUS_BRIDGE.md#source-map); [exact successful recovery](../science/model.md#5-exact-successful-recovery); [heralded entanglement generation](../science/interpretation.md#heralded-entanglement-generation).

<a id="m03"></a>

### M03. Ordinary channels and inaccessible environments

**Phase:** Background/preliminaries.

**Central question:** What information does an ordinary noisy channel expose to its user?

**Understanding target:** Read a CPTP map through its Kraus representation and the normalization $`\sum_a K_a^\dagger K_a=I`$. Explain complete positivity using an untouched reference. Distinguish the mathematical noise label $`a`$ from an observed outcome or a controllable environment register.

**Dependencies and story role:** [M01](#m01), [M02](#m02). Introduces the physical boxes whose order will be controlled.

**Conceptual connection:** Tracing out an environment gives a sum of Kraus operations. Having a formula for that sum does not supply coherent access to its terms. That distinction will separate the present protocol from stronger environment-access models.

**Source pointers:** [Watrous Sections 2.2.1–2.2.2, including Theorem 2.22](WATROUS_BRIDGE.md#source-map); [ordinary channel slots](../science/model.md#1-independent-ordinary-channel-slots); [access distinctions](WATROUS_BRIDGE.md#conceptual-gaps).

<a id="m04"></a>

### M04. A channel as a Choi matrix

**Phase:** Background/preliminaries.

**Central question:** How can a linear map be represented by one positive operator without losing its action on entangled inputs?

**Understanding target:** Connect Kraus and Choi descriptions and use the project's unnormalized, output-first convention:

```math
J_{\mathcal N}
=(\mathcal N\otimes\mathrm{id})(|\Omega_d\rangle\langle\Omega_d|),
\qquad
|\Omega_d\rangle=\sum_i|i\,i\rangle.
```

For a channel, $`\mathrm{Tr}J_{\mathcal N}=d`$ and the output partial trace is $`I_d`$. Distinguish this matrix from the normalized Choi state $`J_{\mathcal N}/d`$.

**Dependencies and story role:** [M02](#m02), [M03](#m03). Converts the channel problem into positive-matrix geometry.

**Conceptual connection:** The identity channel has a rank-one Choi matrix even though it transmits an arbitrary target state. The depolarizer has $`J_\Delta=I_{d^2}/d`$.

**Source pointers:** [Watrous Section 2.2.2, Eq. (2.64) and Theorem 2.26](WATROUS_BRIDGE.md#source-map); [notation translation](WATROUS_BRIDGE.md#notation); [project Choi conventions](../science/model.md#6-choi-and-full-rank-conventions).

<a id="m05"></a>

### M05. Positivity, rank, and spanning

**Phase:** Background/preliminaries.

**Central question:** What does the rank of a positive sum reveal about all of its contributing vectors?

**Understanding target:** Explain why the range of $`\sum_a|v_a\rangle\langle v_a|`$ is the span of the $`v_a`$. Apply this to full Choi rank and to a rank-one Choi operator. In the latter case, positivity forces every contributing vector into the same line; contributions outside it cannot cancel.

**Dependencies and story role:** [M01](#m01), [M04](#m04). Supplies the common geometric ingredient of the later achievability and converse arguments.

**Conceptual connection:** Full Choi rank means that Kraus operators span the entire matrix algebra. Rank-one successful output means that every successful Kraus amplitude has one prescribed operator direction.

**Source pointers:** [Watrous Sections 1.1.3 and 2.2.2](WATROUS_BRIDGE.md#source-map); [full-rank spanning in Theorem 1](../science/polynomial_realization.md#theorem-1-universal-exact-branch); [rank-one Choi argument](../science/query_threshold.md#2-the-rank-one-choi-argument).

<a id="m06"></a>

### M06. Effects, instruments, and accepted operations

**Phase:** Background/preliminaries.

**Central question:** What mathematical operation describes accepting some outcomes of a controller measurement?

**Understanding target:** Distinguish a trace-preserving channel from one completely positive, trace-nonincreasing branch of an instrument. Explain why a final controller instrument, with its quantum output discarded, is represented by an accepted effect $`0\le F\le I`$. Include mixed control preparations, auxiliary control registers, and unions of accepted outcomes.

**Dependencies and story role:** [M03](#m03), [M04](#m04). Defines the success mechanism before any particular rank-one construction is chosen.

**Conceptual connection:** Rank-one postselection is a convenient realization. The later lower bound must also cover more general preparations and accepted effects.

**Source pointers:** [Watrous Sections 2.3.1–2.3.2, especially instruments on pp. 111–113](WATROUS_BRIDGE.md#source-map); [allowed controller and instrument](../science/model.md#3-allowed-controller-and-successful-instrument).

<a id="m07"></a>

### M07. Independent uses and matrix-product order

**Phase:** Background/preliminaries.

**Central question:** What joint resource is specified by several independent channel slots?

**Understanding target:** Translate a mathematical word $`K^{(\pi_1)}_{a_{\pi_1}}\cdots K^{(\pi_m)}_{a_{\pi_m}}`$ into physical time, with its rightmost factor acting first. Explain why a separate noise label belongs to each slot and why repeated use of one device requires independent instances for this model. Equal marginal channels do not determine their joint noise.

**Dependencies and story role:** [M03](#m03), [M06](#m06). Builds multi-call processes from individual channels while fixing the independence assumption.

**Conceptual connection:** The once-per-slot word convention makes degree count channel traversals. Independent environments make the complete noise label a tuple, which is essential for both Kraus invariance and the spanning proof.

**Source pointers:** [Watrous channel and representation foundations, Sections 2.2.1–2.2.2](WATROUS_BRIDGE.md#source-map); [compulsory traversals and word order](../science/model.md#2-compulsory-traversals-and-word-order); [independence requirement](../science/model.md#1-independent-ordinary-channel-slots).

<a id="m08"></a>

### M08. Classical choice and coherent order

**Phase:** Background/preliminaries.

**Central question:** Which alternatives can interfere when the environment remains inaccessible?

**Understanding target:** Contrast a classical mixture of complete orders with a superposition stored in an order register and read by a final effect. For one fixed noise tuple, order amplitudes can add and cancel. Different inaccessible noise tuples remain separate Kraus contributions and are summed incoherently.

**Dependencies and story role:** [M06](#m06), [M07](#m07). Identifies the physical interference mechanism that the polynomial language will describe.

**Conceptual connection:** The controller combines different products of the same slot operators. It neither selects a favorable Kraus outcome nor turns a classical distribution over orders into the same resource.

**Source pointers:** [Watrous Sections 2.1.2, 2.2.2, and 2.3.2](WATROUS_BRIDGE.md#source-map); [conceptual bridge to order control](WATROUS_BRIDGE.md#conceptual-gaps); [switch definition](../science/polynomial_realization.md#1-the-normalized-switch).

<a id="m09"></a>

### M09. Multilinear words and compulsory slots

**Phase:** Background/preliminaries.

**Central question:** Why is the relevant polynomial a sum of permutation words?

**Understanding target:** Interpret

```math
P(X_1,\ldots,X_m)
=\sum_{\pi\in S_m}c_\pi X_{\pi_1}\cdots X_{\pi_m}
```

as a polynomial linear in each variable separately, with every variable occurring once in every word. Connect its degree to the number of compulsory slots and its nonzero coefficients to occupied order labels. Explain how multilinearity extends statements from a spanning set to arbitrary operator tuples.

**Dependencies and story role:** [M01](#m01), [M05](#m05), [M07](#m07), [M08](#m08). Supplies the algebraic language for the physical interference.

**Conceptual connection:** Separate channel instances correspond to separate variables even when their channel maps, or particular matrix values, coincide.

**Source pointers:** [Watrous operator algebra background, Section 1.1.2](WATROUS_BRIDGE.md#source-map); [the polynomial step beyond the textbook](WATROUS_BRIDGE.md#conceptual-gaps); [realizable order polynomials](../science/polynomial_realization.md#2-every-multilinear-order-polynomial-is-realizable).

<a id="m10"></a>

### M10. Central polynomials and polynomial identities

**Phase:** Background/preliminaries.

**Central question:** What is the difference between always producing a scalar and always producing zero?

**Understanding target:** Distinguish a formal nonzero polynomial, a polynomial identity that vanishes on every matrix tuple, and a nonidentity central polynomial that is always scalar and sometimes nonzero. Recognize the separate roles of a centrality proof and a nonzero evaluation. Understand why characteristic-zero polarization converts repeated formal variables into multilinear variables while preserving degree and nonidentity.

**Dependencies and story role:** [M09](#m09). Identifies the precise algebraic resource needed for a useful successful channel.

**Conceptual connection:** Scalarity will give exact recovery; nonidentity together with spanning will give positive success. Neither condition alone supplies both conclusions.

**Source pointers:** [Scope of the Watrous anchor](WATROUS_BRIDGE.md#source-boundaries); [centrality definition](../science/polynomial_realization.md#3-centrality-gives-exact-recovery); [linearization](../science/polynomial_realization.md#4-why-finite-dimensional-constructions-exist); [primary algebraic sources](../references/RELATED_WORK.md#central-polynomials).

<a id="m11"></a>

### M11. Exact recovery, fidelity, and distance

**Phase:** Background/preliminaries.

**Central question:** How does an exact successful channel differ from a high-fidelity output?

**Understanding target:** Read $`\mathcal E=p\,\mathrm{id}_d`$ as equality of completely positive maps, with input-independent success probability and preservation of every reference. Distinguish that equality from entanglement fidelity and from a norm bound on an imperfect accepted operation. The diamond norm includes reference-assisted inputs; conditioning on success can introduce input-dependent normalization.

**Dependencies and story role:** [M02](#m02), [M04](#m04), [M06](#m06). Combines the objects into an exact task with meaningful approximate comparisons.

**Conceptual connection:** The project's $`F_e`$ is an overlap with a maximally entangled state. Watrous's fidelity uses the root convention, so that overlap is the square of the corresponding textbook fidelity.

**Source pointers:** [Watrous Sections 3.2.1 and 3.3.2](WATROUS_BRIDGE.md#source-map); [notation and fidelity conventions](WATROUS_BRIDGE.md#notation); [exact recovery](../science/model.md#5-exact-successful-recovery).

<a id="m12"></a>

### M12. Heralding and distinct resource counts

**Phase:** Background/preliminaries.

**Central question:** What must be counted when a perfect accepted output occurs only sometimes?

**Understanding target:** Distinguish calls per attempt $`m`$, success probability $`p`$, occupied order labels, their register encoding, and physical routing cost. Explain $`m/p`$ as mean noisy-channel calls under independent repetition on fresh entangled pairs. A failed unknown input cannot generally be reused, and a successful identity branch alone is not an unassisted-capacity theorem.

**Dependencies and story role:** [M06](#m06), [M08](#m08), [M11](#m11). Establishes the operational quantities against which the algebraic result will be assessed.

**Conceptual connection:** Encoding many order labels in a small number of qubits does not provide a small routing circuit. Minimum positive use count excludes a zero-call identity wire.

**Source pointers:** [Watrous instrument and capacity context](WATROUS_BRIDGE.md#source-map); [resource bridge](WATROUS_BRIDGE.md#resource-scaling); [resource accounting](../science/probability.md#resource-accounting); [fresh-pair interpretation](../science/interpretation.md#heralded-entanglement-generation).

<a id="m13"></a>

### M13. The complete recovery specification

**Phase:** Background/preliminaries. **Phase synthesis.**

**Central question:** What complete physical and mathematical problem do the preceding ideas define?

**Understanding target:** Assemble an arbitrary reference-target input, an independent controller and ancilla, independent ordinary CPTP slots, compulsory once-per-branch traversals, inaccessible environments, and a final control-only instrument. State the goal $`p\,\mathrm{id}_d`$ with $`p>0`$. Distinguish this architecture from protected input encodings, intermediate target gates, target-carrying bypass memory, and coherent environment access.

**Dependencies and story role:** [M01](#m01) through [M12](#m12), especially the joint description in M07, the exact criterion in M11, and the costs in M12. Integrates the background into one resource-conversion question.

**Conceptual connection:** A dimension-dependent universal construction and a lower bound for one completely known fixed channel have different quantifiers. The supplied generalized switch is itself a resource, not a consequence of the textbook channel axioms.

**Source pointers:** [Watrous-to-project conceptual gaps](WATROUS_BRIDGE.md#conceptual-gaps); [complete model](../science/model.md); [scientific question and proof dependencies](../science/README.md).

## Phase II: Narrative/innovation

<a id="m14"></a>

### M14. Physical order amplitudes become a matrix polynomial

**Phase:** Narrative/innovation.

**Central question:** How does a normalized controller turn a formal word sum into an accepted quantum operation?

**Understanding target:** Explain Kraus independence of the supplied switch and the preparation/effect construction for arbitrary multilinear coefficients. Keep coherent addition of orders within one history separate from the incoherent sum over independent histories.

**Dependencies and story role:** M06–M10 and M13. Starts the narrative with the physical realization needed by every later algebraic claim.

**Conceptual connection:** Coefficient normalization is part of a protocol, not a free rescaling of a mathematical identity.

**Source pointers:** Watrous §§2.2.2 and 2.3.2 via the [bridge](WATROUS_BRIDGE.md#conceptual-gaps); [normalized switch and realization](../science/polynomial_realization.md#1-the-normalized-switch).

<a id="m15"></a>

### M15. Centrality gives an exact branch with positive probability

**Phase:** Narrative/innovation.

**Central question:** Why does a nonidentity central polynomial recover every input rather than merely a selected state?

**Understanding target:** Combine scalar successful Kraus operators with complete positivity to obtain a multiple of the identity channel. Use the full Kraus span and multilinearity to prove that at least one scalar is nonzero. Distinguish centrality from nonidentity.

**Dependencies and story role:** M05, M10, M11 and M14. Joins the algebraic mechanism to the complete quantum recovery goal.

**Conceptual connection:** A scalar on every inaccessible history protects an untouched reference; spanning supplies positive success.

**Source pointers:** Watrous §§1.1.3 and 2.2.2; [universal exact branch](../science/polynomial_realization.md#theorem-1-universal-exact-branch).

<a id="m16"></a>

### M16. The converse and the minimum call threshold

**Phase:** Narrative/innovation.

**Central question:** Why must every exact full-rank protocol obey the same central-polynomial degree limit?

**Understanding target:** Explain the rank-one Choi support argument for general accepted effects, then extend scalarity from the Kraus span to arbitrary matrices. Interpret the minimum positive call count as c(d), including the imported algebraic degree facts and their limits.

**Dependencies and story role:** M05–M07, M09–M11 and M15. Reverses the construction to produce a characterization rather than only an example.

**Conceptual connection:** Positivity prevents cancellation between unwanted successful histories and converts exact recovery into an algebraic necessity.

**Source pointers:** Watrous §§1.1.3, 2.2.2 and 2.3.2; [query characterization](../science/query_threshold.md#theorem-2-exact-query-characterization); [known degrees and dimension bounds](../science/query_threshold.md#5-known-degrees-and-their-quantum-meaning).

<a id="m17"></a>

### M17. The four-call qubit branch

**Phase:** Narrative/innovation.

**Central question:** What makes the eight-word qubit polynomial both central and physically useful?

**Understanding target:** Relate the anticommutator of two commutators to scalar matrices in dimension two, verify one nonzero evaluation, and specify the eight occupied orders and normalized signed herald. Explain what the exact lower-degree certificate establishes.

**Dependencies and story role:** M14–M16. Makes construction and minimality concrete in the smallest target dimension.

**Conceptual connection:** Four channel calls and eight order alternatives are different resource counts; the construction does not prove eight alternatives minimal.

**Source pointers:** Watrous operator and Choi foundations; [qubit polynomial and certificate](../science/qubit.md); [worked depolarizing example](WATROUS_BRIDGE.md#worked-example).

<a id="m18"></a>

### M18. Qutrit and all-dimensional realizations

**Phase:** Narrative/innovation.

**Central question:** How do inherited central polynomials become finite exact protocols beyond a qubit?

**Understanding target:** Explain the degree-eight qutrit construction, its 2,880-word support and nonidentity witness, and distinguish the local centrality proof from imported degree minimality. Place the general existence and Regev constructions alongside their separate degree and probability costs.

**Dependencies and story role:** M14–M17. Extends the same mechanism while retaining the provenance of the algebraic ingredients.

**Conceptual connection:** A universal realization theorem turns algebraic identities into channel capabilities without making all dimensions equally economical.

**Source pointers:** Watrous §§1.1.2–1.1.3 and §2.2.2; [qutrit construction](../science/qutrit.md); [general existence](../science/polynomial_realization.md#4-why-finite-dimensional-constructions-exist); [Regev probability](../science/probability.md#regev-construction).

<a id="m19"></a>

### M19. Probability is part of the exact result

**Phase:** Narrative/innovation.

**Central question:** What does a perfect accepted channel cost when most attempts fail?

**Understanding target:** Use coefficient normalization and exact depolarizing contractions to interpret qubit success 3/64 and qutrit success 1/43,740. Separate p, calls per attempt, occupied orders, register size, and m/p for independent fresh-pair trials.

**Dependencies and story role:** M12, M17 and M18. Converts existence into an explicitly charged operational statement.

**Conceptual connection:** A perfect normalized branch can coexist with a large expected cost; a failed unknown input is not assumed reusable.

**Source pointers:** Watrous §§1.2.2 and 2.3.2; [success probability](../science/probability.md#qubit-success-probability); [resource accounting](../science/probability.md#resource-accounting).

<a id="m20"></a>

### M20. Noise floors, batching and the near-identity boundary

**Phase:** Narrative/innovation.

**Central question:** Why can weaker slot noise make a fixed exact filter rarer, and how can a larger batch help?

**Understanding target:** Distinguish a CP full-rank noise floor from a controller-noise floor. Explain the minimum-call suppression law near identity and the matching fixed-dimension expected-call scaling obtained by contiguous batching. Retain its fixed-batch, fresh-input and known-noise assumptions.

**Dependencies and story role:** M15, M16 and M19. Integrates probability with call-count optimization rather than treating minimum calls as the only cost.

**Conceptual connection:** Batching changes the effective slot strength while still charging every original traversal; it does not optimize routing overhead or establish communication capacity.

**Source pointers:** Watrous §1.2.2 and §2.2.2; [CP noise floor](../science/probability.md#cp-noise-floor); [near-identity suppression](../science/probability.md#noise-approaching-the-identity); [expected charged calls](../science/probability.md#expected-charged-calls).

<a id="m21"></a>

### M21. Which changes leave the theorem’s model?

**Phase:** Narrative/innovation.

**Central question:** Which physical resources would invalidate a transfer of the exact recovery theorem?

**Understanding target:** Explain the correlated-noise and rank-deficient boundary examples, classical-order comparison, and randomizing wrappers. Keep an unknown input, inaccessible environments and compulsory raw traversals visible when comparing protocols.

**Dependencies and story role:** M13, M16 and M20. Closes the ideal-control story by testing the assumptions on which its characterization depends.

**Conceptual connection:** A valid theorem transfer requires matching access and resources; equal marginal channels or a convenient Kraus formula do not suffice.

**Source pointers:** Watrous §§2.2.2 and 3.3.2; [interpretation and boundary cases](../science/interpretation.md); [raw-slot access](../science/noisy_routing_model.md#shared-access).

<a id="m22"></a>

### M22. Physical controller noise and the shared channel family

**Phase:** Narrative/innovation.

**Central question:** Why do order labels alone fail to specify a locally noisy controller?

**Understanding target:** Specify the four raw qubit slots, the eight matrix words, the physical code assignment [0,1,2,3,4,6,5,7], and the measure-Z/prepare-tilted-axis channel with depolarizing admixture. Distinguish global dephasing eta from independent physical-bit dephasing gamma and explain why per-slot wrappers alter the task.

**Dependencies and story role:** M07–M08, M13 and M21. Introduces the matched noisy experiment with an operational encoding rather than a mere relabeling.

**Conceptual connection:** Schur noise weights interference by physical Hamming distance, so scalar controller coherence does not identify all useful order correlations.

**Source pointers:** Watrous §§2.2.2 and 4.1.3; [channel family](../science/noisy_routing_model.md#channel-family); [physical codebook](../science/noisy_routing_model.md#physical-codebook).

<a id="m23"></a>

### M23. All allowed controllers have one positive-matrix budget

**Phase:** Narrative/innovation.

**Central question:** How can mixed preparations, ancillary registers and final instruments be optimized without guessing a controller?

**Understanding target:** Interpret the exact feasible set H_x positive semidefinite with their sum bounded by diag(t), where t is a probability vector. Explain the converse realization, symmetry reduction to I/8, and the Gram formulas for acceptance and Bell reward. Distinguish arbitrary complex branches from real attaining branches with input-independent success.

**Dependencies and story role:** M05–M06, M11 and M22. Replaces a chosen postselection with a complete optimization domain.

**Conceptual connection:** Positive-matrix feasibility charges the entire outcome instrument; ancillas do not enlarge this domain. Real and Y-symmetry reductions justify an attaining recovery channel rather than assuming all branches normalize linearly.

**Source pointers:** Watrous §§1.2.3, 2.3.2 and 3.2; [controller feasibility](../science/controller_optimization.md#controller-feasibility); [real-controller reduction](../science/controller_optimization.md#real-controller-reduction).

<a id="m24"></a>

### M24. A strictly positive Schur floor excludes exact recovery

**Phase:** Narrative/innovation.

**Central question:** Why does any nonzero diagonal acceptance obstruct an exact unitary branch under the specified controller noise?

**Understanding target:** Use the CP lower bound from a strictly positive Schur kernel and full Choi rank of composed raw slots to rule out a nonzero rank-one accepted Choi operator. State the finite-call scope and the excluded access to noise flags or coherence restoration.

**Dependencies and story role:** M05, M16 and M22–M23. Separates the noiseless exact mechanism from its mathematically precise noisy obstruction.

**Conceptual connection:** The crucial fact is unavoidable diagonal contribution, not a claim that every imaginable controller noise destroys exactness.

**Source pointers:** Watrous §§1.1.3, 2.2.2 and 4.1.3; [strict Schur obstruction](../science/noisy_recovery_obstruction.md#strict-schur-obstruction).

<a id="m25"></a>

### M25. The complete orthogonal recovery frontier

**Phase:** Narrative/innovation.

**Central question:** When does a simple physical readout attain the optimum over every allowed controller?

**Understanding target:** Explain physical Walsh diagonalization of the orthogonal-axis acceptance and corrected reward matrices. Reduce the optimization to fractional knapsack with mode capacities b_j/8, rewards a_j/8, and modes sorted by a_j/b_j; include partial acceptance of the boundary mode and zero-capacity cases.

**Dependencies and story role:** M22–M23. Turns the controller feasible set into an exact acceptance–fidelity tradeoff with an explicit attaining measurement.

**Conceptual connection:** An upper bound and a physical three-X readout attain the same curve, so the result is a full frontier for this architecture rather than only a good herald.

**Source pointers:** Watrous §§1.2.3 and 2.3.1–2.3.2; [frontier theorem](../science/recovery_frontier.md#frontier-theorem); [controller reduction](../science/controller_optimization.md#real-controller-reduction).

<a id="m26"></a>

### M26. Herald choice depends on noise geometry

**Phase:** Narrative/innovation.

**Central question:** Why can matching a controller’s scalar coherence fail to match its best recovery outcome?

**Understanding target:** Compare the inherited P2 branch and the alternative Q branch under global and physical local dephasing. Keep success probability and corrected fidelity together, distinguish a selected-branch comparison from the full frontier, and explain the reversal of the preferred branch at the stated parameters.

**Dependencies and story role:** M19, M22 and M25. Exposes which physical correlations the optimized frontier uses.

**Conceptual connection:** A scalar summary can agree while Schur kernels weight the same order pairs differently; fidelity alone also omits acceptance cost.

**Source pointers:** Watrous §§3.2 and 4.1.3; [exact frontier and branch evaluations](../science/recovery_frontier.md); [physical noise conventions](../science/noisy_routing_model.md#physical-codebook).

<a id="m27"></a>

### M27. Tilted readout optimality without commutation

**Phase:** Narrative/innovation.

**Central question:** How can the same readout remain globally optimal after the reward matrices cease to commute?

**Understanding target:** For epsilon=1/100 and physical gamma=4/5, explain the diagonal-envelope upper bound, the two coupled Walsh blocks, and the exact polynomial margin. Distinguish the interval in t=tan(delta/2), abs(t) at most 1/5, from an interval stated directly in the angle delta. Explain why Bernstein positivity certifies the entire interval.

**Dependencies and story role:** M23 and M25–M26. Generalizes the optimality argument through an operator majorant rather than an unjustified continuity claim.

**Conceptual connection:** Noncommutation alone does not give a collective-measurement advantage; a common dominating diagonal reward can still be tight.

**Source pointers:** Watrous §§1.1.3 and 1.2.3; [tilted optimality theorem and certificate](../science/tilted_readout.md#tilted-optimality); [verification guide](../verification/README.md).

<a id="m28"></a>

### M28. Recovery and retained entanglement separate in a robust experiment

**Phase:** Narrative/innovation.

**Central question:** What can a matched channel-and-noise comparison actually conclude about C:T entanglement and R:T recovery?

**Understanding target:** Explain the sister result: for abs(delta) at most .2 radians, the specified unconditional three-call C:T outputs are separable on occupied 3-by-2 support, while the fixed four-call P2 branch yields NPT R:T output with positive acceptance bounds. Preserve the different calls, controller bits, bipartitions and conditional endpoints.

**Dependencies and story role:** M11–M12, M22, M24 and M27. Connects the recovery track to its sister through a proved comparison of specified protocols.

**Conceptual connection:** Retained controller–target entanglement is not a recovery certificate. This experiment does not optimize equal resources or classify all retained protocols.

**Source pointers:** Watrous §§3.2, 6.1.1 and 6.2.3; [uniform separation](https://github.com/GoGoKo699/Quantum-Order-Entanglement/blob/main/science/noisy_separation.md#uniform-separation); [operational meaning](https://github.com/GoGoKo699/Quantum-Order-Entanglement/blob/main/science/noisy_separation.md#operational-meaning).

<a id="m29"></a>

### M29. Proof, certificate, evaluation and attribution

**Phase:** Narrative/innovation.

**Central question:** Which parts of the scientific account are analytical theorems, exact finite certificates or numerical crosschecks?

**Understanding target:** Trace the degree imports, inherited qubit branch, controller reduction, exact coefficient reconstruction, Bernstein positivity and rational CT/RT certificates to their proper sources. Explain what the independent rebuild verifies and why Bell fidelity at most one half does not prove separability of every possible complex-controller Choi output.

**Dependencies and story role:** M16–M18 and M23–M28. Integrates evidence and attribution before the final scientific synthesis.

**Conceptual connection:** Finite exact obligations can support a uniform theorem when combined with its analytical argument; numerical agreement and bounded literature search retain narrower meanings.

**Source pointers:** Watrous positivity and fidelity foundations; [source boundaries](WATROUS_BRIDGE.md#source-boundaries); [related work](../references/RELATED_WORK.md); [verification commands](../verification/README.md).

<a id="m30"></a>

### M30. The complete capability, frontier and resource statement

**Phase:** Narrative/innovation. **Phase synthesis.**

**Central question:** What does coherent order enable, what does controller noise forbid, and what is the charged operational conclusion?

**Understanding target:** Connect ideal exact recovery and its central-degree threshold to the strict-noise no-go, the attainable noisy recovery frontier, tilted robustness and the specified entanglement comparison. State the independent raw-slot model, final-only corrections and fresh-input accounting in one coherent explanation, without inferring unassisted capacity or unrestricted controller–target impossibility.

**Dependencies and story role:** M13–M29. Phase synthesis: combines construction, converse, noisy optimization, sister connection and evidence boundaries.

**Conceptual connection:** The common mathematical foundation supports distinct scientific tasks; it is the access model, endpoint, error and success convention together that determine what transfers.

**Source pointers:** Watrous channel, instrument, fidelity and separability sections through the [bridge](WATROUS_BRIDGE.md); [scientific route](../science/README.md); [shared access](../science/noisy_routing_model.md#shared-access); [verification](../verification/README.md).
