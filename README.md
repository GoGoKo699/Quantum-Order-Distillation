# Quantum Order Distillation

**How many uses of a noisy quantum channel are needed to herald an exact identity channel by controlling only their order?**

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
\qquad 0<\epsilon\le1,
```

every minimum-call exact controller satisfies

```math
p_{\mathcal C}(\mathcal N_\epsilon)
=\epsilon^{c(d)}p_{\mathcal C}(\Delta_d).
```

Among exact protocols that repeat a fixed batch on **fresh entangled-pair inputs**, the batch size may instead be chosen as a function of the known noise strength. The optimal expected number of original noisy-slot calls per successfully heralded pair then scales as $`\Theta_d(1/\epsilon)`$. A binomial-tail bound gives the lower bound, and contiguous blocks of channel uses give the upper bound. Constants depend on dimension and can be large. These statements concern the unmodified channel slots; optional randomizing wrappers define a different resource model.

## Read the science

Start with the [scientific overview](science/README.md) for the connected argument, or choose a specific question:

| Question | Authoritative account |
|---|---|
| What operations and resources are supplied? | [Model and access](science/model.md) |
| How can inaccessible noise histories all produce the same successful channel? | [Polynomial realization](science/polynomial_realization.md) |
| Why is the minimum number of calls exactly a central-polynomial degree? | [Query threshold](science/query_threshold.md) |
| What are the actual four- and eight-call protocols? | [Qubit construction](science/qubit.md) and [qutrit construction](science/qutrit.md) |
| How likely is success, and what does it cost? | [Probability and resources](science/probability.md) |
| What is the physical meaning, and which assumptions are essential? | [Interpretation and boundaries](science/interpretation.md) |
| Which ingredients and capabilities were already known? | [Related work](references/RELATED_WORK.md) and [bibliography](references/references.bib) |

The [LLM reading guide](llms.txt) maps questions to the same authoritative files.

## Learning route

The single educational anchor is John Watrous, [*The Theory of Quantum Information* (2018)](https://jhwatrous.github.io/books-and-courses.html). The [Watrous bridge](learning/WATROUS_BRIDGE.md) selects the relevant sections, translates channel and Choi conventions, and explains the steps from ordinary quantum channels to interference between their orders. Central-polynomial identities and the exact query characterization are supplied by the project and its cited research sources.

The [30-milestone map](learning/MILESTONES.md) follows the scientific story through 13 background/preliminaries milestones and 17 narrative/innovation milestones. Each milestone states a central question, an understanding target, its earlier dependencies, and pointers to the relevant source and canonical scientific account. The map develops the mechanism, the converse, the explicit qubit and qutrit constructions, and the distinction between exactness and useful success rates. Its targets are defined by the project and apply to any reader.

Use the bridge and map together, following links to the full proofs when needed. The [scientific overview](science/README.md) remains the direct route to the complete argument.

## Reproduce the checks

From the repository root, with Python 3.11 or later:

```bash
python3 -m pip install -r verification/requirements.txt
python3 verification/run_checks.py
```

The [verification guide](verification/README.md) explains the exact coefficient checks, finite degree certificate, independent probability contractions, and numerical channel checks, including an optional exhaustive C++ calculation.

The analytical proofs establish the general claims. Finite checks verify transcription, normalization, concrete witnesses, and implementations; sampled matrix evaluations do not replace a centrality proof.

## Scope and attribution

The immediate operational consequence is perfect heralded entanglement distribution under the supplied switch access. Failed attempts can disturb the target, so retries require fresh inputs. This does not establish an unassisted quantum-capacity theorem or a fault-tolerant implementation.

The minimum-use characterization is specific to compulsory traversals, independent slots, and control-only processing. Target encodings, intermediate target gates, target-carrying bypass memories, and accessible noise environments define larger models. The [interpretation](science/interpretation.md) proves a classical-order separation within the stated architecture and gives a correlated-noise counterexample.

Central polynomials, their known minimum degrees, the qubit superswitch branch, and central-polynomial quantum resetting are established antecedents. The [source comparison](references/RELATED_WORK.md) credits these inputs and distinguishes finite exact channel distillation from earlier communication activation, asymptotic heralded transmission, and resetting with returned probes.

## Purpose and contact

This repository serves as a record of the work and a guide for the author’s self-directed learning. For discussion or potential collaboration, please contact Ruge Lin at [gogoko699@gmail.com](mailto:gogoko699@gmail.com).
