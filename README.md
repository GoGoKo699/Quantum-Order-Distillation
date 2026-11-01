# Quantum Order Distillation

**How many uses of a noisy quantum channel are needed to herald an exact identity channel by controlling only their order?**

**Sister project:** [Quantum Order Entanglement](https://github.com/GoGoKo699/Quantum-Order-Entanglement) studies retained controller–target entanglement from coherent order. This repository studies heralded channel recovery and its costs. The [shared-noise comparison](https://github.com/GoGoKo699/Quantum-Order-Entanglement/blob/main/science/noisy_separation.md) connects the two tasks while keeping their outputs, resources and conditioning explicit.

For a finite-dimensional channel with full Choi rank, the answer is an invariant from polynomial identity theory: the minimum degree of a nonidentity multilinear central polynomial. This repository develops the construction, proves its converse in a precise access model, and separates the minimum number of calls from the probability and cost of success.

The target traverses every independent channel slot once in each coherent order. The order control starts independently of the unknown input, and success is read from that control. The noise environments are inaccessible. The [model](science/model.md) specifies these assumptions and the scope of the supplied generalized switch.

## Principal result

Let $`c(d)`$ be the smallest positive degree of a multilinear polynomial that evaluates to a scalar matrix on every tuple in $`M_d(\mathbb C)`$ and is not identically zero there. For every full-Choi-rank channel $`\mathcal N`$ on a $`d`$-dimensional target,

```math
\boxed{
\min\{\text{positive channel uses for exact heralded identity}\}
=c(d).
}
```

The upper bound uses a controller determined only by the dimension and a chosen central polynomial. It also works for different full-rank channels in the individual slots. On success, the map is

```math
\mathcal S_{\mathrm{succ}}=p\,\mathrm{id}_d,\qquad p>0,
```

so it preserves an arbitrary input and its entanglement with any reference. The probability is independent of that input. The converse allows the controller to be optimized for the entire description of one fixed full-rank channel, including mixed controls, auxiliary control registers, general successful effects, and final herald-dependent unitary corrections. See the [construction](science/polynomial_realization.md) and [minimum-use proof](science/query_threshold.md).

### Explicit qubit and qutrit protocols

The algebraic minimum degrees are known for dimensions two and three. The associated protocols are explicit:

| Target | Calls per attempt | Occupied order labels | Success on complete depolarization | Mean channel calls per successful fresh entangled pair |
|---|---:|---:|---:|---:|
| [Qubit](science/qubit.md) | 4, minimum | 8 | $`3/64`$ | $`256/3\approx85.3`$ |
| [Qutrit](science/qutrit.md) | 8, minimum | 2,880 | $`1/43{,}740`$ | $`349{,}920`$ |

These probabilities belong to the stated preparation and success effect. They are not global probability optima. Occupied order labels measure controller support; encoding those labels in qubits does not establish an efficient routing circuit.

For $`d>2`$, the general bounds are

```math
2d-1\le c(d)\le(d-1)^2+4.
```

In particular, $`c(4)\le13`$; equality is not established here. For an entire $`k`$-qubit channel, the lower bound is $`2^{k+1}-1`$ calls. The finite-dimensional existence theorem therefore carries no efficiency guarantee in the number of target qubits.

## Probability and resource consequences

A positive success probability can still be small. The [probability account](science/probability.md) gives a lower bound from the smallest Choi eigenvalue, formulas for general qubit channels, two exact qutrit probability calculations, and a determinant formula for an all-dimensional Regev construction.

There is also a sharp distinction between fewest calls and expected calls. For

```math
\mathcal N_\epsilon=(1-\epsilon)\,\mathrm{id}_d+\epsilon\Delta_d,
\qquad 0\lt \epsilon\le1,
```

every minimum-call exact controller satisfies

```math
p_{\mathcal C}(\mathcal N_\epsilon)
=\epsilon^{c(d)}p_{\mathcal C}(\Delta_d).
```

Among exact protocols that repeat a fixed batch on **fresh entangled-pair inputs**, the batch size may instead be chosen as a function of the known noise strength. The optimal expected number of original noisy-slot calls per successfully heralded pair then scales as $`\Theta_d(1/\epsilon)`$. A binomial-tail bound gives the lower bound, and contiguous blocks of channel uses give the upper bound. Constants depend on dimension and can be large. These statements concern the unmodified channel slots; optional randomizing wrappers define a different resource model.

## Recovery with a noisy controller

The exact central-polynomial result assumes ideal order coherence. With unobserved controller noise, the useful question becomes the best corrected Bell fidelity at a specified accepted probability. The [shared noisy-routing model](science/noisy_routing_model.md) fixes four independent raw qubit slots, eight nested-switch orders, and their physical three-bit encoding. The channel measures along $`Z`$, prepares along a tilted axis, and includes a depolarizing admixture. Every slot is entanglement breaking and has full Choi rank when that admixture is positive.

For arbitrary controller preparations, controller ancillas, accepted instruments and final outcome-dependent target unitaries, we establish:

| Result | Exact scope and consequence |
|---|---|
| [Physical controller characterization](science/controller_optimization.md#controller-feasibility) | Accepted coefficient matrices are precisely $`H_x\succeq0`$, $`\sum_xH_x\preceq\mathrm{diag}(t)`$, $`t\ge0`$, $`\sum_it_i=1`$; ancillary controllers do not enlarge this set. |
| [Complete orthogonal-axis frontier](science/recovery_frontier.md#frontier-theorem) | For every channel admixture and every stated global visibility or equal local phase factor, three-$`X`$ readout with outcome-dependent corrections and acceptance attains the full probability–fidelity frontier. |
| [Robust tilted optimality](science/tilted_readout.md#tilted-optimality) | At admixture $`\epsilon=0.01`$ and physical local phase factor $`\gamma=0.8`$, the same readout is optimal for all $`\lvert\tan(\delta/2)\rvert\le1/5`$, despite noncommuting overlap matrices. |
| [Exact-recovery obstruction](science/noisy_recovery_obstruction.md#strict-schur-obstruction) | Independent full-Choi-rank slots and a strictly positive Schur-noise lower bound exclude any nonzero exact accepted unitary branch at finite call count in the stated access model. |

At the orthogonal point with $`\epsilon=0.01`$ and $`\gamma=0.8`$ (a 10% independent phase-flip probability per physical bit), the highest conditional Bell fidelity is approximately **0.758264**, available up to accepted probability **0.0619638**. Accepting every outcome instead gives optimized Bell fidelity **0.473296**. These optimize the stated Bell objective; a fidelity at or below one half does not classify every possible output's entanglement.

An arbitrary complex controller can have input-dependent acceptance. The [reduction proof](science/controller_optimization.md#real-controller-reduction) establishes that a real, input-independent representative attains the Bell-input optimum. The physical encoding is part of the local-noise theorem; merely relabeling the words in lexical binary order changes the noise model.

The [robust retained-$`C:T`$/accepted-$`R:T`$ comparison](https://github.com/GoGoKo699/Quantum-Order-Entanglement/blob/main/science/noisy_separation.md#uniform-separation) uses this raw channel family and per-bit noise with specified three-call retained and four-call accepted protocols. Their different bipartitions, resources and conditioning remain explicit.

## A complete three-call benchmark and randomized cost

For the fixed rank-deficient channel $`K_0=|+\rangle\langle0|`$, $`K_1=|-\rangle\langle1|`$, the [native three-call theorem](science/native_three_call_frontier.md#native-three-call-frontier) solves the complete success–half-diamond-error frontier using all six order labels and global logical controller dephasing. Four projector flags, final corrections and fractional acceptance attain the converse over arbitrary mixed, ancillary and complex controllers.

The [randomized-call theorem](science/randomized_call_cost.md#convex-call-cost) gives the exact optimum across fresh one-, two- and three-call trial types. It determines when an extra call lowers expected calls per accepted output at a fixed half-diamond error, including unequal logical visibilities. The table below compares its equal-visibility benchmark with physical local noise; the full account gives the attaining trial probabilities and success costs.

This benchmark averages corrected fine labels into the accepted target-only channel and requires input-independent aggregate acceptance. Fine outcomes may individually depend on the input. Its rank-two channel, six orders and global noise differ from the full-rank four-call physical-local-noise family above. The [unequal-noise criterion](science/randomized_call_cost.md#unequal-noise) and [conditional duration model](science/randomized_call_cost.md#duration-model) state exactly when the extra call remains useful; equal visibility across architectures is a model choice, rather than a hardware consequence.

## Physical encoding and the usefulness of an extra call

The [encoded local-noise analysis](science/local_controller_encoding.md) replaces one global visibility with independent phase damping on an explicit physical controller. Coherences decay with the Hamming distance between assigned codewords. For a pairwise-comparison encoding, it gives a matching high-coherence recovery frontier and a duration-dependent cost threshold.

The [all-assignment cost theorem](science/all_encoding_cost_bound.md) then answers the calibrated question over every injective assignment of the six orders to three-bit computational words, allowing nonuniform preparation, ancillary controllers and continuous final unitary corrections.

| Noise and resource model; half-diamond error at most $`3/8`$ | Supported conclusion |
|---|---|
| Global logical visibilities $`\eta_2=\eta_3=0.9`$ | Exact randomized cost $`183/62`$ versus $`101/32`$: approximately 6.48% saving. |
| Independent local phase damping, $`x=0.9`$, $`q=x^{3/2}`$, charged windows $`T_m=m\tau`$ | All 20,160 injective six-label assignments obey cost at least $`101/32`$. The one/two-call mixture attains it; any positive three-call trial share is strictly suboptimal. |
| Same local noise law and schedule, $`x=0.9604`$, specified alternative encoding | A certified construction costs approximately 2.823333 versus the optimal one/two-call baseline 2.889770: approximately 2.2990% saving. Global encoding optimality at this different point is not claimed. |

The exclusion uses 97 exact symmetry classes and 388 rational positive-matrix certificates. Its [standard-library verifier](verification/encoding_bound/README.md) covers every allowed controller strategy in that assignment model, rather than testing only a finite correction menu. The [local-controller checks](verification/local_controller/README.md) reproduce the encoding-specific frontier and positive construction.

These are rank-two-channel recovery and fresh-trial cost results. They retain their own scope alongside the full-rank central-polynomial theorem and four-call physical-noise frontier. The charged routing windows omit gate and readout overhead; neither cost model is a hardware validation. Retained controller–target entanglement and accepted reference–target recovery remain distinct tasks.

## Read the science

Start with the [scientific overview](science/README.md) for the connected argument, or choose a specific question:

| Question | Authoritative account |
|---|---|
| What operations and resources are supplied? | [Model and access](science/model.md) |
| How can inaccessible noise histories all produce the same successful channel? | [Polynomial realization](science/polynomial_realization.md) |
| Why is the minimum number of calls exactly a central-polynomial degree? | [Query threshold](science/query_threshold.md) |
| What are the actual four- and eight-call protocols? | [Qubit construction](science/qubit.md) and [qutrit construction](science/qutrit.md) |
| How likely is success, and what does it cost? | [Probability and resources](science/probability.md) |
| What is the optimal noisy-controller recovery tradeoff? | [Controller characterization](science/controller_optimization.md), [frontier](science/recovery_frontier.md), and [tilted theorem](science/tilted_readout.md) |
| When does controller noise exclude exact recovery? | [Strict Schur obstruction](science/noisy_recovery_obstruction.md) |
| How do physical codewords and duration change recovery cost? | [Encoded local noise](science/local_controller_encoding.md) and [all-assignment cost bound](science/all_encoding_cost_bound.md) |
| What is the physical meaning, and which assumptions are essential? | [Interpretation and boundaries](science/interpretation.md) |
| Which ingredients and capabilities were already known? | [Related work](references/RELATED_WORK.md) and [bibliography](references/references.bib) |

The [LLM reading guide](llms.txt) maps questions to the same authoritative files.

## Learning route

The single educational anchor is John Watrous, [*The Theory of Quantum Information* (2018)](https://jhwatrous.github.io/books-and-courses.html). The [Watrous bridge](learning/WATROUS_BRIDGE.md) selects the relevant sections, translates channel and Choi conventions, and explains the steps from ordinary quantum channels to interference between their orders. Central-polynomial identities and the exact query characterization are supplied by the project and its cited research sources.

The [30-milestone map](learning/MILESTONES.md) follows the scientific story through 13 background/preliminaries milestones and 17 narrative/innovation milestones. Each milestone states a central question, an understanding target, its earlier dependencies, and pointers to the relevant source and canonical scientific account. The map connects the exact central-polynomial mechanism and its cost to controller feasibility, noisy recovery frontiers, tilted optimality and the distinct retained-entanglement comparison. Its targets are defined by the project and apply to any reader.

Use the bridge and map together, following links to the full proofs when needed. The [scientific overview](science/README.md) remains the direct route to the complete argument.

## Reproduce the checks

From the repository root, with Python 3.11 or later:

```bash
python3 -m pip install -r verification/requirements.txt
python3 verification/run_checks.py
python3 verification/audit_verify.py --root verification/certificates --rebuild
```

The [verification guide](verification/README.md) explains the original exact algebra and channel checks, the complete noisy-frontier polynomial rebuild, and all 249 retained-state plus 83 recovery certificates. It also documents an optional exhaustive C++ probability calculation.

The analytical proofs establish the general claims. Finite checks verify transcription, normalization, concrete witnesses, and implementations; sampled matrix evaluations do not replace a centrality proof.

## Scope and attribution

With ideal order control, the immediate operational consequence is perfect heralded entanglement distribution under the supplied switch access. Failed attempts can disturb the target, so retries require fresh inputs. This does not establish an unassisted quantum-capacity theorem or a fault-tolerant implementation.

The minimum-use characterization is specific to compulsory traversals, independent slots, and control-only processing. Target encodings, intermediate target gates, target-carrying bypass memories, and accessible noise environments define larger models. The [interpretation](science/interpretation.md) proves a classical-order separation within the stated architecture and gives a correlated-noise counterexample.

Central polynomials, their known minimum degrees, the qubit superswitch branch, and central-polynomial quantum resetting are established antecedents. The [source comparison](references/RELATED_WORK.md) credits these inputs and distinguishes finite exact channel distillation from earlier communication activation, asymptotic heralded transmission, and resetting with returned probes.

## Purpose and contact

This repository serves as a record of the work and a guide for the author’s self-directed learning. For discussion or potential collaboration, please contact Ruge Lin at [gogoko699@gmail.com](mailto:gogoko699@gmail.com).
