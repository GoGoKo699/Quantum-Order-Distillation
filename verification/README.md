# Verification

This suite checks the exact coefficient data, probability normalization, and finite-dimensional switch implementation used in the [scientific account](../science/probability.md). The construction and converse proofs remain the basis for the general claims. Random matrix checks test their transcription and implementation.

## Run the suite

Use **Python 3.11 or later**. The default suite requires only the Python standard library and NumPy. From the repository root, run:

```bash
python3 -m pip install -r verification/requirements.txt
python3 verification/run_checks.py
```

The requirements file pins NumPy to 2.3.5, the version used for the recorded results. The runner resolves its inputs relative to its own location, so it also works when invoked from another directory. It reads the coefficient export and canonical result files without changing them. Any failed equality, residual bound, or canonical-result comparison produces a nonzero exit status.

The default run includes the independent 96-term trace calculation of the qutrit probability. It does not require a C++ compiler. To additionally recompute the full histogram of all pairs of the 2880 qutrit words, use:

```bash
python3 verification/run_checks.py --exhaustive
```

That option requires `g++` with C++17 support. The helper is compiled in a temporary directory, which is removed after use. It computes the symmetric half of the pair table and doubles off-diagonal contributions, accounting exactly for all $`2880^2=8294400`$ ordered pairs. Its signed 64-bit counters are sufficient for this supplied coefficient list; arbitrary large coefficient lists would need a separate overflow check.

To save a machine-readable report of a run, choose a new output path:

```bash
python3 verification/run_checks.py --report latest_results.json
```

The committed [check_results.json](check_results.json) records the complete run, including the optional exhaustive histogram. The individual canonical JSON files below contain the detailed expected results. No source PDF or network access is required.

## What is checked

| Group | Implementation | Canonical output and checks |
|---|---|---|
| Qubit lower-degree exclusion | [qubit_minimum.py](qubit_minimum.py) | [qubit_minimum_results.json](qubit_minimum_results.json): exact matrix-unit evaluation certificates of ranks 1, 2, and 6, with determinants 1, -1, and -1, for degrees 1, 2, and 3. |
| Qutrit construction | [qutrit_degree8.py](qutrit_degree8.py) | [qutrit_results.json](qutrit_results.json) and [degree8_coefficients.json](degree8_coefficients.json): two independent coefficient generators agree on 2880 distinct signed words; the source tuple evaluates to $`-2I_3`$; 12 integer matrix tuples agree exactly with the scalar trace formula; one $`2\times2`$ tuple gives zero. |
| Prepared switch | [qubit_switch_checks.py](qubit_switch_checks.py) | [qubit_switch_results.json](qubit_switch_results.json): actual block-controlled switch Kraus matrices, trace preservation, identity Choi output, nonunital probability formula, untouched reference system, heralding from the full joint output, unequal complex controller phases, and independent rectangular Kraus isometries in every slot. |
| General control instruments and access assumptions | [check_access_model.py](check_access_model.py) | [access_model_checks.json](access_model_checks.json): a mixed state on order control and ancilla, a general success effect, the positive coefficient Gram representation, a classical-order infidelity bound, identity-mixture scaling, and a correlated-noise counterexample. |
| Independent qutrit probability | [qutrit_trace_crosscheck.py](qutrit_trace_crosscheck.py) | [qutrit_trace_results.json](qutrit_trace_results.json): 96 trace terms yield exactly $`1/43740`$ through 9216 ordinary-permutation contractions, independently of the 2880-word pair calculation. |
| Depolarizing word contractions | [depolarizing_probability.py](depolarizing_probability.py), optionally [cycle_histogram.cpp](cycle_histogram.cpp) | [depolarizing_results.json](depolarizing_results.json): the qubit histogram gives exactly $`3/64`$; the stored qutrit histogram obeys the normalization identities and gives $`1/43740`$ at dimension 3 and zero at dimensions 1 and 2. With `--exhaustive`, the entire qutrit histogram is regenerated from the committed coefficient list and compared exactly. |
| Regev probabilities and resources | [regev_probabilities.py](regev_probabilities.py) | [regev_results.json](regev_results.json): exact constants, rational probabilities, support sizes, and call counts for dimensions 2 through 5, plus finite-dimensional evaluations of the logarithmic asymptotic expansion. |
| Regev determinant identity | [regev_identity_checks.py](regev_identity_checks.py) | Recorded in [check_results.json](check_results.json): exact sparse summation on the matrix-unit bases in dimensions 2 and 3, six random complex qubit tuples, and the ordered-tuple Cauchy–Binet factor for a six-Kraus qubit channel. |

The scripts with `--output` or `--report` flags can also be run individually. The aggregate runner is the command that compares the results against the canonical files. It does not regenerate those files to make a failed comparison pass.

## Conventions and expected values

The exported qutrit word labels are one-based and specify **mathematical matrix multiplication from left to right**. Physical temporal order is reversed. Some Python routines convert labels to zero-based indices internally.

The main switch check uses unnormalized Choi matrices with output-then-reference row-major vectorization, so $`\operatorname{Tr}J_{\mathcal N}=d`$ for a channel and $`p=\operatorname{Tr}J_{\mathrm{succ}}/d`$. The independent access-model check uses column-major vectorization consistently within its own Choi calculations. This change of coordinate convention leaves its probabilities, eigenvalues, and matrix-norm comparisons unchanged.

The qubit prepared-switch check constructs the full eight-order channel before contracting the preparation and success vectors. Its inputs include a depolarizer, a full-rank Pauli channel, a generic nonunital channel, four distinct full-rank channels, amplitude damping, and a unitary channel. A three-dimensional reference is included in the state-preservation test. The Kraus-gauge test independently expands each of the four channels from four to six Kraus operators and compares the complete prepared output, including the unmeasured controller.

| Quantity | Exact expected value |
|---|---:|
| Four-call qubit depolarizing success | $`3/64`$ |
| Eight-call qutrit depolarizing success | $`1/43740`$ |
| Regev determinant constant, $`d=2`$ | $`-3`$ |
| Regev determinant constant, $`d=3`$ | $`720`$ |
| Eight-call Regev depolarizing success, $`d=2`$ | $`1/16384`$ |
| Eighteen-call Regev depolarizing success, $`d=3`$ | $`1/98411002933824`$ |

The nonzero entries of the exact qutrit word-pair histogram are

```math
(h_1,h_3,h_5,h_7,h_9)
=(124416,-259200,172224,-40320,2880).
```

The independent trace histogram has weights $`1536,-1920,384`$ at cycle counts $`4,6,8`$. These are different contractions of the same scalar probability and are compared through their exact rational answers.

## Reproducibility and limits

Integer coefficient dictionaries, certificate matrices, signed histograms, and rational probabilities are compared exactly. Canonical floating-point values use relative tolerance $`10^{-10}`$ and absolute tolerance $`10^{-12}`$. The individual calculations also enforce their scientific residual bounds, mostly $`10^{-11}`$ or $`10^{-12}`$, with stricter bounds for exact boundary examples. Small floating-point differences across numerical libraries are therefore allowed; discrepancies in exact data are not.

The qutrit integer examples and the main numerical probability checks use seed `20261028`. The independent access-model calculation uses seed `202610280314`. The report records both seeds and the Python and NumPy versions. Numerical examples are fixed reproducibility checks, not statistical evidence of a universal theorem.

The default run independently recomputes the qutrit probability from its short trace formula, validates the stored all-word histogram, and regenerates the coefficient list. Recomputing every coefficient pair is specifically the additional work performed by `--exhaustive`. The Regev identity is supported by the cited theorem in every dimension; the explicit polynomial evaluations in this suite cover dimensions 2 and 3. The suite does not simulate a physical routing apparatus or charge its gates, geometry, and error correction.

For source locations and the exact normalization derivations, see [probability and resources](../science/probability.md), the [qubit construction](../science/qubit.md), the [qutrit construction](../science/qutrit.md), and the [related-work account](../references/RELATED_WORK.md).
