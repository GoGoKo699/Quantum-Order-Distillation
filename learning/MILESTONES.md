# Thirty milestones for quantum order distillation

This map follows the scientific question: **when can control of channel order alone herald an exact identity channel, and what does that capability cost?** It connects the physical access model, the central-polynomial mechanism, the minimum-use theorem, explicit protocols, and their operational limits.

The single educational anchor is **John Watrous, *The Theory of Quantum Information* (2018)**. The [Watrous bridge](WATROUS_BRIDGE.md) identifies the relevant sections, translates notation, supplies a worked example, and explains which research ingredients lie beyond the textbook. Primary research sources remain necessary for the inherited polynomial identities and their degree bounds; they are distinguished in [Related work](../references/RELATED_WORK.md).

## Phase structure

| Phase | Milestones | Scientific purpose | Closing synthesis |
|---|---|---|---|
| Background/preliminaries | [M01](#m01) through [M13](#m13) | Connect states, channel representations, conditional operations, independent ordered uses, and matrix polynomials into one precise recovery problem. | [M13: The complete recovery specification](#m13) |
| Narrative/innovation | [M14](#m14) through [M30](#m30) | Build the construction and converse, make their dimension and probability consequences concrete, and determine their physical scope. | [M30: The integrated operational theorem](#m30) |

The **13 + 17 allocation** reflects this project's meeting of two subjects. Quantum-channel mathematics alone does not explain centrality, while a matrix identity alone does not establish a physical channel protocol or its cost. The first phase builds the common language and ends with the complete task. The second needs room for both directions of the characterization, the qubit and qutrit realizations, and the distinction between minimum calls, expected calls, and routing resources.

Within each phase, the route proceeds from individual objects and arguments to their combined meaning. Narrative milestones include inherited mathematical tools and known special cases as well as the project consequences. The map targets global understanding; specialist derivations remain in their authoritative scientific pages.

## Navigation

| ID | Background/preliminaries | ID | Narrative/innovation |
|---|---|---|---|
| [M01](#m01) | States and target matrices | [M14](#m14) | A valid switch independent of Kraus representation |
| [M02](#m02) | Entangled reference systems | [M15](#m15) | Turning coefficients into control amplitudes |
| [M03](#m03) | Ordinary channels and inaccessible environments | [M16](#m16) | Centrality gives exact positive success |
| [M04](#m04) | A channel as a Choi matrix | [M17](#m17) | Why exact full-rank recovery forces centrality |
| [M05](#m05) | Positivity, rank, and spanning | [M18](#m18) | The minimum degree becomes a call threshold |
| [M06](#m06) | Effects, instruments, and accepted operations | [M19](#m19) | Four calls for a qubit |
| [M07](#m07) | Independent uses and matrix-product order | [M20](#m20) | Eight calls for a qutrit |
| [M08](#m08) | Classical choice and coherent order | [M21](#m21) | Evaluating specified controllers |
| [M09](#m09) | Multilinear words and compulsory slots | [M22](#m22) | An all-dimensional construction and its cost |
| [M10](#m10) | Central polynomials and polynomial identities | [M23](#m23) | Two consequences of a CP noise floor |
| [M11](#m11) | Exact recovery, fidelity, and distance | [M24](#m24) | Exact recovery near the identity boundary |
| [M12](#m12) | Heralding and distinct resource counts | [M25](#m25) | Optimizing expected original-channel calls |
| [M13](#m13) | The complete recovery specification | [M26](#m26) | Correlated and rank-deficient boundaries |
| | | [M27](#m27) | Randomizing wrappers change the resource |
| | | [M28](#m28) | Exact algebra with imperfect control |
| | | [M29](#m29) | Finite exact recovery in the communication landscape |
| | | [M30](#m30) | The integrated operational theorem |

The dependency links in each record identify earlier milestones that supply its main ideas. The [proof route](WATROUS_BRIDGE.md#proof-route) connects this map to the authoritative statements, and the [scientific overview](../science/README.md) provides a shorter research route.

## Phase I: Background/preliminaries

<a id="m01"></a>

### M01. States and target matrices

**Phase:** Background/preliminaries.

**Central question:** Which matrices describe physical states, and which matrices must a channel identity handle?

**Understanding target:** Distinguish a density operator from a general linear operator on a $`d`$-dimensional target. Connect adjoints, positivity, trace, matrix units, and linear combinations. In particular, Kraus operators and polynomial arguments need not be density operators, Hermitian matrices, or unitary matrices.

**Dependencies and story role:** No earlier milestone. This supplies the local objects from which channels, noise histories, and polynomial evaluations will be built.

**Conceptual connection:** A matrix basis describes all possible operator amplitudes, while positivity and normalization select physical states. Confusing those roles would weaken the later centrality condition.

**Sources:** [Watrous Sections 1.1.2, 1.1.3, and 2.1.2](WATROUS_BRIDGE.md#source-map); [project target and channel spaces](../science/model.md#1-independent-ordinary-channel-slots).

<a id="m02"></a>

### M02. Entangled reference systems

**Phase:** Background/preliminaries.

**Central question:** Why must recovery preserve a target's correlations with a system that the protocol never touches?

**Understanding target:** Explain tensor-product systems, reduced states, and the role of an untouched reference. Relate the normalized maximally entangled state $`|\Phi_d\rangle`$ to preservation of an unknown quantum input. Distinguish the target's dimension $`d`$ from the dimension $`d^2`$ of its joint space with an equal-sized reference.

**Dependencies and story role:** [M01](#m01). Enlarges the local target description to the correlated inputs needed for a complete quantum task.

**Conceptual connection:** A protocol that preserves selected target states can still damage entanglement. Channel-level recovery includes those correlations from the start.

**Sources:** [Watrous Section 2.1.3](WATROUS_BRIDGE.md#source-map); [exact successful recovery](../science/model.md#5-exact-successful-recovery); [heralded entanglement generation](../science/interpretation.md#heralded-entanglement-generation).

<a id="m03"></a>

### M03. Ordinary channels and inaccessible environments

**Phase:** Background/preliminaries.

**Central question:** What information does an ordinary noisy channel expose to its user?

**Understanding target:** Read a CPTP map through its Kraus representation and the normalization $`\sum_a K_a^\dagger K_a=I`$. Explain complete positivity using an untouched reference. Distinguish the mathematical noise label $`a`$ from an observed outcome or a controllable environment register.

**Dependencies and story role:** [M01](#m01), [M02](#m02). Introduces the physical boxes whose order will be controlled.

**Conceptual connection:** Tracing out an environment gives a sum of Kraus operations. Having a formula for that sum does not supply coherent access to its terms. That distinction will separate the present protocol from stronger environment-access models.

**Sources:** [Watrous Sections 2.2.1–2.2.2, including Theorem 2.22](WATROUS_BRIDGE.md#source-map); [ordinary channel slots](../science/model.md#1-independent-ordinary-channel-slots); [access distinctions](WATROUS_BRIDGE.md#conceptual-gaps).

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

**Sources:** [Watrous Section 2.2.2, Eq. (2.64) and Theorem 2.26](WATROUS_BRIDGE.md#source-map); [notation translation](WATROUS_BRIDGE.md#notation); [project Choi conventions](../science/model.md#6-choi-and-full-rank-conventions).

<a id="m05"></a>

### M05. Positivity, rank, and spanning

**Phase:** Background/preliminaries.

**Central question:** What does the rank of a positive sum reveal about all of its contributing vectors?

**Understanding target:** Explain why the range of $`\sum_a|v_a\rangle\langle v_a|`$ is the span of the $`v_a`$. Apply this to full Choi rank and to a rank-one Choi operator. In the latter case, positivity forces every contributing vector into the same line; contributions outside it cannot cancel.

**Dependencies and story role:** [M01](#m01), [M04](#m04). Supplies the common geometric ingredient of the later achievability and converse arguments.

**Conceptual connection:** Full Choi rank means that Kraus operators span the entire matrix algebra. Rank-one successful output means that every successful Kraus amplitude has one prescribed operator direction.

**Sources:** [Watrous Sections 1.1.3 and 2.2.2](WATROUS_BRIDGE.md#source-map); [full-rank spanning in Theorem 1](../science/polynomial_realization.md#theorem-1-universal-exact-branch); [rank-one Choi argument](../science/query_threshold.md#2-the-rank-one-choi-argument).

<a id="m06"></a>

### M06. Effects, instruments, and accepted operations

**Phase:** Background/preliminaries.

**Central question:** What mathematical operation describes accepting some outcomes of a controller measurement?

**Understanding target:** Distinguish a trace-preserving channel from one completely positive, trace-nonincreasing branch of an instrument. Explain why a final controller instrument, with its quantum output discarded, is represented by an accepted effect $`0\le F\le I`$. Include mixed control preparations, auxiliary control registers, and unions of accepted outcomes.

**Dependencies and story role:** [M03](#m03), [M04](#m04). Defines the success mechanism before any particular rank-one construction is chosen.

**Conceptual connection:** Rank-one postselection is a convenient realization. The later lower bound must also cover more general preparations and accepted effects.

**Sources:** [Watrous Sections 2.3.1–2.3.2, especially instruments on pp. 111–113](WATROUS_BRIDGE.md#source-map); [allowed controller and instrument](../science/model.md#3-allowed-controller-and-successful-instrument).

<a id="m07"></a>

### M07. Independent uses and matrix-product order

**Phase:** Background/preliminaries.

**Central question:** What joint resource is specified by several independent channel slots?

**Understanding target:** Translate a mathematical word $`K^{(\pi_1)}_{a_{\pi_1}}\cdots K^{(\pi_m)}_{a_{\pi_m}}`$ into physical time, with its rightmost factor acting first. Explain why a separate noise label belongs to each slot and why repeated use of one device requires independent instances for this model. Equal marginal channels do not determine their joint noise.

**Dependencies and story role:** [M03](#m03), [M06](#m06). Builds multi-call processes from individual channels while fixing the independence assumption.

**Conceptual connection:** The once-per-slot word convention makes degree count channel traversals. Independent environments make the complete noise label a tuple, which is essential for both Kraus invariance and the spanning proof.

**Sources:** [Watrous channel and representation foundations, Sections 2.2.1–2.2.2](WATROUS_BRIDGE.md#source-map); [compulsory traversals and word order](../science/model.md#2-compulsory-traversals-and-word-order); [independence requirement](../science/model.md#1-independent-ordinary-channel-slots).

<a id="m08"></a>

### M08. Classical choice and coherent order

**Phase:** Background/preliminaries.

**Central question:** Which alternatives can interfere when the environment remains inaccessible?

**Understanding target:** Contrast a classical mixture of complete orders with a superposition stored in an order register and read by a final effect. For one fixed noise tuple, order amplitudes can add and cancel. Different inaccessible noise tuples remain separate Kraus contributions and are summed incoherently.

**Dependencies and story role:** [M06](#m06), [M07](#m07). Identifies the physical interference mechanism that the polynomial language will describe.

**Conceptual connection:** The controller combines different products of the same slot operators. It neither selects a favorable Kraus outcome nor turns a classical distribution over orders into the same resource.

**Sources:** [Watrous Sections 2.1.2, 2.2.2, and 2.3.2](WATROUS_BRIDGE.md#source-map); [conceptual bridge to order control](WATROUS_BRIDGE.md#conceptual-gaps); [switch definition](../science/polynomial_realization.md#1-the-normalized-switch).

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

**Sources:** [Watrous operator algebra background, Section 1.1.2](WATROUS_BRIDGE.md#source-map); [the polynomial step beyond the textbook](WATROUS_BRIDGE.md#conceptual-gaps); [realizable order polynomials](../science/polynomial_realization.md#2-every-multilinear-order-polynomial-is-realizable).

<a id="m10"></a>

### M10. Central polynomials and polynomial identities

**Phase:** Background/preliminaries.

**Central question:** What is the difference between always producing a scalar and always producing zero?

**Understanding target:** Distinguish a formal nonzero polynomial, a polynomial identity that vanishes on every matrix tuple, and a nonidentity central polynomial that is always scalar and sometimes nonzero. Recognize the separate roles of a centrality proof and a nonzero evaluation. Understand why characteristic-zero polarization converts repeated formal variables into multilinear variables while preserving degree and nonidentity.

**Dependencies and story role:** [M09](#m09). Identifies the precise algebraic resource needed for a useful successful channel.

**Conceptual connection:** Scalarity will give exact recovery; nonidentity together with spanning will give positive success. Neither condition alone supplies both conclusions.

**Sources:** [Scope of the Watrous anchor](WATROUS_BRIDGE.md#source-boundaries); [centrality definition](../science/polynomial_realization.md#3-centrality-gives-exact-recovery); [linearization](../science/polynomial_realization.md#4-why-finite-dimensional-constructions-exist); [primary algebraic sources](../references/RELATED_WORK.md#central-polynomials).

<a id="m11"></a>

### M11. Exact recovery, fidelity, and distance

**Phase:** Background/preliminaries.

**Central question:** How does an exact successful channel differ from a high-fidelity output?

**Understanding target:** Read $`\mathcal E=p\,\mathrm{id}_d`$ as equality of completely positive maps, with input-independent success probability and preservation of every reference. Distinguish that equality from entanglement fidelity and from a norm bound on an imperfect accepted operation. The diamond norm includes reference-assisted inputs; conditioning on success can introduce input-dependent normalization.

**Dependencies and story role:** [M02](#m02), [M04](#m04), [M06](#m06). Combines the objects into an exact task with meaningful approximate comparisons.

**Conceptual connection:** The project's $`F_e`$ is an overlap with a maximally entangled state. Watrous's fidelity uses the root convention, so that overlap is the square of the corresponding textbook fidelity.

**Sources:** [Watrous Sections 3.2.1 and 3.3.2](WATROUS_BRIDGE.md#source-map); [notation and fidelity conventions](WATROUS_BRIDGE.md#notation); [exact recovery](../science/model.md#5-exact-successful-recovery).

<a id="m12"></a>

### M12. Heralding and distinct resource counts

**Phase:** Background/preliminaries.

**Central question:** What must be counted when a perfect accepted output occurs only sometimes?

**Understanding target:** Distinguish calls per attempt $`m`$, success probability $`p`$, occupied order labels, their register encoding, and physical routing cost. Explain $`m/p`$ as mean noisy-channel calls under independent repetition on fresh entangled pairs. A failed unknown input cannot generally be reused, and a successful identity branch alone is not an unassisted-capacity theorem.

**Dependencies and story role:** [M06](#m06), [M08](#m08), [M11](#m11). Establishes the operational quantities against which the algebraic result will be assessed.

**Conceptual connection:** Encoding many order labels in a small number of qubits does not provide a small routing circuit. Minimum positive use count excludes a zero-call identity wire.

**Sources:** [Watrous instrument and capacity context](WATROUS_BRIDGE.md#source-map); [resource bridge](WATROUS_BRIDGE.md#resource-scaling); [resource accounting](../science/probability.md#resource-accounting); [fresh-pair interpretation](../science/interpretation.md#heralded-entanglement-generation).

<a id="m13"></a>

### M13. The complete recovery specification

**Phase:** Background/preliminaries. **Phase synthesis.**

**Central question:** What complete physical and mathematical problem do the preceding ideas define?

**Understanding target:** Assemble an arbitrary reference-target input, an independent controller and ancilla, independent ordinary CPTP slots, compulsory once-per-branch traversals, inaccessible environments, and a final control-only instrument. State the goal $`p\,\mathrm{id}_d`$ with $`p>0`$. Distinguish this architecture from protected input encodings, intermediate target gates, target-carrying bypass memory, and coherent environment access.

**Dependencies and story role:** [M01](#m01) through [M12](#m12), especially the joint description in M07, the exact criterion in M11, and the costs in M12. Integrates the background into one resource-conversion question.

**Conceptual connection:** A dimension-dependent universal construction and a lower bound for one completely known fixed channel have different quantifiers. The supplied generalized switch is itself a resource, not a consequence of the textbook channel axioms.

**Sources:** [Watrous-to-project conceptual gaps](WATROUS_BRIDGE.md#conceptual-gaps); [complete model](../science/model.md); [scientific question and proof dependencies](../science/README.md).

## Phase II: Narrative/innovation

<a id="m14"></a>

### M14. A valid switch independent of Kraus representation

**Phase:** Narrative/innovation.

**Central question:** Why does the order-controlled construction define a physical operation on channels themselves?

**Understanding target:** Explain how orthogonal order projectors and trace preservation of each complete composition give a trace-preserving switch. Follow why independent isometric changes of slot Kraus families combine into one tensor-product isometry on the history label, leaving the switch and its measurement branches unchanged.

**Dependencies and story role:** [M03](#m03), [M07](#m07), [M08](#m08), [M13](#m13). Establishes the physical validity of the resource before using interference to recover a state.

**Conceptual connection:** Each slot occurs once in every word. That common structure makes the construction insensitive to mathematical Kraus choices without granting access to the environments.

**Sources:** [Watrous Section 2.2.2, Corollary 2.23](WATROUS_BRIDGE.md#source-map); [normalized switch](../science/polynomial_realization.md#1-the-normalized-switch); [independence and Kraus representations](../science/polynomial_realization.md#independence-and-kraus-representations).

<a id="m15"></a>

### M15. Turning coefficients into control amplitudes

**Phase:** Narrative/innovation.

**Central question:** How can a formal signed or complex polynomial become an accepted quantum operation with correct normalization?

**Understanding target:** Connect control preparation, coefficient phases, and rank-one postselection to

```math
B_{\boldsymbol a}
=\frac{P(K^{(1)}_{a_1},\ldots,K^{(m)}_{a_m})}{L},
\qquad L=\sum_\pi|c_\pi|.
```

Explain why $`1/L`$ is the maximal scale for that fixed formal coefficient ray within the pure-preparation, single-rank-one-effect class. Distinguish this normalization statement from global optimization over polynomials, equivalent representations, or several accepted outcomes.

**Dependencies and story role:** [M06](#m06), [M09](#m09), [M14](#m14). Converts the algebraic object into a definite controller with a charged success probability.

**Conceptual connection:** Interference signs and normalization are physical parts of the protocol, not arbitrary rescalings of a formal identity.

**Sources:** [Watrous measurement foundations](WATROUS_BRIDGE.md#source-map); [polynomial realization](../science/polynomial_realization.md#2-every-multilinear-order-polynomial-is-realizable); [fixed-ray normalization](../science/probability.md#optimality-for-a-fixed-coefficient-ray).

<a id="m16"></a>

### M16. Centrality gives exact positive success

**Phase:** Narrative/innovation.

**Central question:** Why do inaccessible noise histories all give the same accepted target channel?

**Understanding target:** Explain the two steps of Theorem 1. Centrality makes every successful history operator a scalar multiple of identity, so their incoherent sum is $`p\,\mathrm{id}_d`$. If each slot has full Choi rank, Kraus spanning and multilinearity prevent every scalar from vanishing. The same controller works for heterogeneous full-rank slots and needs no channel description.

**Dependencies and story role:** [M05](#m05), [M10](#m10), [M11](#m11), [M15](#m15). Establishes achievability with both exactness and nonzero success.

**Conceptual connection:** The environment may retain the noise history. Recovery follows because all accepted histories act identically on the target. For arbitrary channels scalarity still holds, but positive success is not guaranteed outside the stated full-rank condition.

**Sources:** [Watrous Choi/Kraus background](WATROUS_BRIDGE.md#source-map); [Theorem 1 and its proof](../science/polynomial_realization.md#theorem-1-universal-exact-branch); [worked channel example](WATROUS_BRIDGE.md#worked-example).

<a id="m17"></a>

### M17. Why exact full-rank recovery forces centrality

**Phase:** Narrative/innovation.

**Central question:** Could a controller tailored to one known full-rank channel avoid the central-polynomial condition?

**Understanding target:** Follow the converse from a rank-one successful Choi operator to scalarity of every successful Kraus operator, then from a spanning Kraus family to centrality on arbitrary matrices. Understand why mixed preparations, independent control ancillas, general accepted effects, and final herald-dependent unitary corrections do not remove the lower bound.

**Dependencies and story role:** [M05](#m05), [M06](#m06), [M09](#m09), [M16](#m16). Completes the reverse implication needed for an exact characterization.

**Conceptual connection:** Positivity forbids cancellation between different Kraus histories even though order amplitudes can cancel within a history. The fixed-channel lower bound allows complete knowledge of the channel; it is not merely a limitation of one universal controller.

**Sources:** [Watrous positive Choi representation](WATROUS_BRIDGE.md#source-map); [rank-one converse](../science/query_threshold.md#2-the-rank-one-choi-argument); [general controls](../science/query_threshold.md#3-general-control-states-and-successful-events); [final corrections](../science/query_threshold.md#4-final-unitary-corrections-do-not-lower-the-count).

<a id="m18"></a>

### M18. The minimum degree becomes a call threshold

**Phase:** Narrative/innovation.

**Central question:** What exactly does $`m_\star(\mathcal N)=c(d)`$ determine?

**Understanding target:** Combine both theorem directions to identify the minimum positive number of compulsory calls for every fixed full-rank channel. Explain $`c(2)=4`$, $`c(3)=8`$, the upper bound $`c(4)\le13`$, and

```math
2d-1\le c(d)\le(d-1)^2+4,\qquad d>2.
```

Understand the lower-bound connection between central polynomials and polynomial identities through $`[P,X_{m+1}]`$. An exact minimum does not establish feasibility at every larger compulsory degree. For an entire $`k`$-qubit target, the lower bound is $`2^{k+1}-1`$ original calls.

**Dependencies and story role:** [M10](#m10), [M12](#m12), [M16](#m16), [M17](#m17). Turns the local construction and converse into the global query statement.

**Conceptual connection:** Algebraic degree becomes an operational invariant, while efficiency, success probability, and control complexity remain separate questions.

**Sources:** [Research ingredients beyond Watrous](WATROUS_BRIDGE.md#source-boundaries); [exact query characterization](../science/query_threshold.md#theorem-2-exact-query-characterization); [degree bounds](../science/query_threshold.md#7-a-dimension-dependent-lower-bound); [degree provenance](../references/RELATED_WORK.md#central-polynomials).

<a id="m19"></a>

### M19. Four calls for a qubit

**Phase:** Narrative/innovation.

**Central question:** How does the abstract mechanism become an explicit, minimal qubit protocol?

**Understanding target:** Connect the anticommutator of two commutators,

```math
P_2=\{[X_1,X_2],[X_3,X_4]\},
```

to eight signed order words and normalization $`L=8`$. Explain centrality using traceless $`2\times2`$ matrices and identify a nonzero Pauli evaluation. Distinguish the exact matrix-unit certificate excluding degrees one through three from sampled numerical checks.

**Dependencies and story role:** [M01](#m01), [M15](#m15), [M16](#m16), [M18](#m18). Supplies the smallest complete example of the characterization.

**Conceptual connection:** The protocol realizes an inherited central polynomial and a known qubit superswitch branch. Its explicit controller makes the general theorem tangible without turning this established special case into a separate novelty claim.

**Sources:** [Watrous operator foundations and worked example](WATROUS_BRIDGE.md#worked-example); [qubit construction and certificate](../science/qubit.md); [superswitch antecedent](../references/RELATED_WORK.md#kechrimparis2025).

<a id="m20"></a>

### M20. Eight calls for a qutrit

**Phase:** Narrative/innovation.

**Central question:** What establishes an explicit finite exact protocol beyond the qubit case?

**Understanding target:** Explain how the Drensky–Kasparian placement formula specifies 2,880 signed words of degree eight. Connect the exterior-algebra centrality argument, the scalar trace formula, and the exact nonzero matrix-unit witness to its channel realization. The polynomial is already multilinear, so no extra polarization factor or physical identity slot is introduced.

**Dependencies and story role:** [M10](#m10), [M15](#m15), [M16](#m16), [M18](#m18). Extends the explicit mechanism to a qutrit while keeping its proof dependencies visible.

**Conceptual connection:** The repository supplies the degree-eight construction and its centrality proof. Excluding all lower qutrit degrees is an imported algebraic theorem, not a conclusion of the finite coefficient or sampled-matrix checks. Controller support 2,880 is sufficient, not proven minimal.

**Sources:** [Boundary between textbook and research algebra](WATROUS_BRIDGE.md#source-boundaries); [qutrit construction](../science/qutrit.md); [imported qutrit minimum](../science/qutrit.md#5-why-eight-calls-are-optimal); [primary polynomial source](../references/RELATED_WORK.md#central-polynomials).

<a id="m21"></a>

### M21. Evaluating specified controllers

**Phase:** Narrative/innovation.

**Central question:** How are exact success probabilities obtained for the actual controller coefficients?

**Understanding target:** Relate $`p=L^{-2}\sum_{\boldsymbol a}|\lambda_{\boldsymbol a}|^2`$ to two useful calculation routes: qubit Pauli-vector covariances, and matrix-unit contraction counts for depolarizing slots. Explain why the general nonunital qubit formula contains a transpose rather than an adjoint. Understand the independent all-word and scalar-trace calculations of the qutrit probability.

**Dependencies and story role:** [M04](#m04), [M15](#m15), [M19](#m19), [M20](#m20). Quantifies the two explicit protocols after establishing their exactness.

**Conceptual connection:** The displayed controllers give $`p(\Delta_2)=3/64`$ and $`p(\Delta_3)=1/43740`$. These are exact values for specified preparations and effects, not global probability optima. Agreement of independent contractions checks normalization and transcription.

**Sources:** [Watrous Kraus and measurement foundations](WATROUS_BRIDGE.md#source-map); [general qubit probability](../science/probability.md#qubit-success-probability); [depolarizing contractions](../science/probability.md#exact-depolarizing-contractions); [qutrit probability](../science/probability.md#eight-call-qutrit-probability).

<a id="m22"></a>

### M22. An all-dimensional construction and its cost

**Phase:** Narrative/innovation.

**Central question:** What does an explicit construction in every dimension reveal beyond existence?

**Understanding target:** Connect the Regev polynomial's double alternation and determinant evaluation to an explicitly normalized controller of degree $`2d^2`$. Explain how Cauchy–Binet turns its Kraus-tuple probability into a constant times $`\det(J_{\mathcal N})^2`$. Full Choi rank is necessary and sufficient for this particular family's positive success, without classifying other controllers.

**Dependencies and story role:** [M04](#m04), [M10](#m10), [M15](#m15), [M18](#m18), [M21](#m21). Extends quantitative understanding from two small dimensions to a general family.

**Conceptual connection:** An all-dimensional exact formula can exhibit severe probability and control-support growth. This family is generally not query minimal and gives no practical-rate or uniform efficiency guarantee.

**Sources:** [Watrous matrix and channel foundations](WATROUS_BRIDGE.md#source-map); [Regev construction and probability](../science/probability.md#regev-construction); [direct Procesi determinant source and version qualification](../references/RELATED_WORK.md#regev-procesi).

<a id="m23"></a>

### M23. Two consequences of a CP noise floor

**Phase:** Narrative/innovation.

**Central question:** How can the same positive noise component guarantee quantum success and obstruct exact classical-order recovery?

**Understanding target:** Use $`\eta_j=d\,\lambda_{\min}(J_{\mathcal N_j})>0`$ to decompose each full-rank channel into a depolarizer and a residual channel. For a fixed exact controller, retaining the all-depolarizing contribution gives a positive lower bound on success. In a classical composition, every term containing a depolarizer is a replacer, giving a strictly positive gap from identity.

**Dependencies and story role:** [M04](#m04), [M08](#m08), [M11](#m11), [M16](#m16). Moves from the exact examples in M21 to structural guarantees and a comparison within the same access model.

**Conceptual connection:** Positivity controls both arguments, but coherent order and a classical mixture retain different conditional operations. The gap concerns accepted-channel exactness; it does not establish superior throughput or a no-go theorem for circuits with protected target memory.

**Sources:** [Watrous positive-map and distance foundations](WATROUS_BRIDGE.md#source-map); [CP noise-floor bound](../science/probability.md#cp-noise-floor); [quantitative classical-order gap](../science/interpretation.md#a-quantitative-gap-from-the-identity).

<a id="m24"></a>

### M24. Exact recovery near the identity boundary

**Phase:** Narrative/innovation.

**Central question:** Why can a minimum-call exact branch become rarer as the supplied channel approaches identity?

**Understanding target:** For $`\mathcal N_\epsilon=(1-\epsilon)\mathrm{id}+\epsilon\Delta`$, explain why a minimum-degree central polynomial vanishes when any argument is identity. Therefore, for $`0<\epsilon\le1`$,

```math
p_{\mathcal C}(\mathcal N_\epsilon)
=\epsilon^{c(d)}p_{\mathcal C}(\Delta),
\qquad m=c(d).
```

Relate the same law for the optimized probability to the common full-rank controller class. Understand the binomial-tail upper bound for arbitrary fixed $`m`$, without inferring achievability at every such degree.

**Dependencies and story role:** [M09](#m09), [M17](#m17), [M18](#m18), [M23](#m23). Reveals a singular resource consequence of requiring exactness.

**Conceptual connection:** At $`\epsilon=0`$, the supplied channel is identity and deterministic success is possible. The feasible controller class changes, so the optimized full-rank law does not extend to that endpoint.

**Sources:** [Watrous mixture and channel foundations](WATROUS_BRIDGE.md#source-map); [noise approaching identity](../science/probability.md#noise-approaching-the-identity); [resource bridge](WATROUS_BRIDGE.md#resource-scaling).

<a id="m25"></a>

### M25. Optimizing expected original-channel calls

**Phase:** Narrative/innovation.

**Central question:** Can spending more calls in one attempt reduce the total calls needed for a successful fresh pair?

**Understanding target:** Explain the matching fixed-dimension scaling

```math
\mathcal Q_d(\epsilon)
=\inf_{\text{exact fixed-batch protocols}}\frac{m}{p}
=\Theta_d(1/\epsilon),
\qquad 0<\epsilon\le1.
```

The binomial bound gives $`m/p\ge c(d)/\epsilon`$. Contiguous blocks of $`k=\lceil1/\epsilon\rceil`$ original calls give effective noise bounded away from zero, after which a central controller supplies the upper bound. Count all original calls, including failed attempts, rather than counting a block as one original use.

**Dependencies and story role:** [M12](#m12), [M18](#m18), [M23](#m23), [M24](#m24). Separates minimum calls per attempt from optimal mean calls.

**Conceptual connection:** The result concerns repeated fixed batches on fresh entangled pairs, known or comparably bounded $`\epsilon`$, and unmodified channel slots. Any explicit central degree suffices for a constructive upper bound. Dimension-dependent constants can be large; routing remains charged separately.

**Sources:** [Watrous bridge to resource scaling](WATROUS_BRIDGE.md#resource-scaling); [expected charged calls and proof](../science/probability.md#expected-charged-calls); [resource accounting](../science/probability.md#resource-accounting).

<a id="m26"></a>

### M26. Correlated and rank-deficient boundaries

**Phase:** Narrative/innovation.

**Central question:** Which conclusions fail when full rank or independence is removed?

**Understanding target:** Follow the correlated qubit example in which each marginal channel has full Choi rank, yet four slots share one random unitary, every order has the same product, and the overall target channel is depolarizing. Contrast rank-deficient examples: finite-order unitary channels can recover by ordinary composition, while repeated complete dephasing cannot be repaired by order control.

**Dependencies and story role:** [M05](#m05), [M07](#m07), [M16](#m16), [M17](#m17), [M23](#m23). Uses concrete alternatives to locate the assumptions of the general theorem.

**Conceptual connection:** Full-rank marginals do not replace joint independence. Conversely, full Choi rank is a sufficient hypothesis for the universal result and the converse, not a necessary condition for every recoverable individual channel.

**Sources:** [Watrous channel foundations and access boundaries](WATROUS_BRIDGE.md#conceptual-gaps); [exact correlated counterexample](../science/interpretation.md#correlated-instances); [rank-deficient boundaries](../science/interpretation.md#rank-deficient-channels).

<a id="m27"></a>

### M27. Randomizing wrappers change the resource

**Phase:** Narrative/innovation.

**Central question:** What does the all-channel wrapper corollary add to the full-rank result?

**Understanding target:** Explain how a finite unitary 1-design supplies an output depolarizer, giving $`\Delta\circ\mathcal N=\Delta`$ for every original CPTP channel. A central controller can then succeed on the composite slots. Every original slot is still traversed, but the attached randomizers are additional target operations and must be charged as resources.

**Dependencies and story role:** [M03](#m03), [M13](#m13), [M16](#m16), [M25](#m25), [M26](#m26). Separates an augmented corollary from the original conversion problem.

**Conceptual connection:** The original channel becomes informationally irrelevant in this construction. The unaugmented small-noise lower bound therefore does not survive the resource change. Randomizer environments must remain attached consistently to their slots across coherent orders.

**Sources:** [Watrous Sections 4.1.1–4.1.2 on mixed unitaries and Weyl operators](WATROUS_BRIDGE.md#source-map); [depolarization notation](WATROUS_BRIDGE.md#notation); [randomizing wrappers](../science/interpretation.md#randomizing-wrappers).

<a id="m28"></a>

### M28. Exact algebra with imperfect control

**Phase:** Narrative/innovation.

**Central question:** Why does universality over noisy channels not automatically imply robustness to control errors?

**Understanding target:** Distinguish changing an independent channel inside the model from changing routing, controller amplitudes, or target interactions. If the actual accepted operation differs from $`p\,\mathrm{id}`$ by diamond norm at most $`\delta<p`$, explain the bound $`2\delta/(p-\delta)`$ on the conditional trace norm error for every reference-inclusive input.

**Dependencies and story role:** [M11](#m11), [M12](#m12), [M16](#m16), [M21](#m21). Connects ideal exactness and small success probabilities to implementation accuracy.

**Conceptual connection:** Normalized success can depend nonlinearly on the input through its probability. The conditional estimate is pointwise, not a diamond distance between normalized channels. It motivates a separate control-accuracy analysis without establishing a fault-tolerance threshold or efficient physical router.

**Sources:** [Watrous Sections 3.2.1 and 3.3.2](WATROUS_BRIDGE.md#source-map); [conditional-error derivation](../science/interpretation.md#ideal-control-and-conditional-error); [routing and control resources](../science/probability.md#resource-accounting).

<a id="m29"></a>

### M29. Finite exact recovery in the communication landscape

**Phase:** Narrative/innovation.

**Central question:** Which operational capability survives once earlier mechanisms and their access assumptions are accounted for?

**Understanding target:** Locate finite exact identity heralding and its minimum-use characterization for each fixed full-rank channel relative to central-polynomial resetting with returned probes, the known qubit superswitch branch, finite-error cyclic transmission approaching identity, and compilation with coherent environment access. Explain why these statements can share algebra while supplying different channel resources or guarantees.

**Dependencies and story role:** [M13](#m13), [M18](#m18), [M20](#m20), [M23](#m23), [M25](#m25), [M26](#m26), [M28](#m28). Places the combined capability and its cost in a physically matched comparison.

**Conceptual connection:** The surviving consequence is a finite exact successful channel with a characterized minimum within the declared architecture. Its proof combines inherited identities with channel realization and a converse. Exact coefficients and numerical checks support specified implementations; they do not replace the general proofs or the imported qutrit lower bound.

**Sources:** [Watrous bridge and its research boundaries](WATROUS_BRIDGE.md#source-boundaries); [physical comparisons and primary sources](../references/RELATED_WORK.md#physical-comparisons); [proof and verification roles](../science/README.md#verification-route); [operational interpretation](../science/interpretation.md).

<a id="m30"></a>

### M30. The integrated operational theorem

**Phase:** Narrative/innovation. **Phase synthesis.**

**Central question:** How do the mechanism, converse, examples, costs, and boundaries form one complete scientific result?

**Understanding target:** Give a connected account of independent ordinary full-rank slots, compulsory traversals, an initially independent controller, and inaccessible environments. For repeated independent uses of any fixed full-rank channel $`\mathcal N`$, connect the exact successful branch and the minimum-use statement:

```math
\mathcal E_{\mathrm{succ}}=p\,\mathrm{id}_d,\qquad p>0,
\qquad
m_\star(\mathcal N)=c(d).
```

Connect this theorem to the explicit four- and eight-call protocols, the classical-order exactness gap, and the distinct expected-call optimization. Preserve the positive-minimum qualifier, the imported qutrit minimality, the unproved equality at $`d=4`$, and the $`\epsilon>0`$ domain of the small-noise laws.

**Dependencies and story role:** [M14](#m14) through [M29](#m29), integrated with the complete specification in [M13](#m13). Closes the narrative by relating every principal conclusion to the assumptions and proof ingredients that support it.

**Conceptual connection:** Perfect accepted entanglement from fresh trials is an operational consequence under the supplied switch access. Correlated uses, added target resources, control errors, routing overhead, and unknown-input failure prevent broader conclusions without further analysis. Minimum calls, positive probability, efficient implementation, and unassisted capacity remain different claims.

**Sources:** [Watrous-to-project proof route](WATROUS_BRIDGE.md#proof-route); [scientific overview](../science/README.md); [model and result map](../science/model.md#7-result-map); [query theorem](../science/query_threshold.md#theorem-2-exact-query-characterization); [probability and resources](../science/probability.md); [interpretation and boundaries](../science/interpretation.md).
