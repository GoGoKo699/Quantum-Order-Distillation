# Related work and algebraic sources

The task studied here is heralded transmission through independent fresh channel slots in a generalized quantum switch. With ideal order control, the question is exact identity recovery. With noisy order control, the question is the attainable acceptance–corrected-Bell-fidelity frontier. Every order visits every prescribed slot once. The control is prepared independently of the target input and measured after the calls; the noise environments remain inaccessible. The [access model](../science/model.md) and [noisy routing model](../science/noisy_routing_model.md) specify these assumptions.

The [polynomial realization](../science/polynomial_realization.md) and [query threshold](../science/query_threshold.md) connect this task to an established algebraic invariant: the minimum positive degree $`c(d)`$ of a nonidentity multilinear central polynomial for $`M_d(\mathbb C)`$. For any fixed channel with full Choi rank, the minimum compulsory call count is $`c(d)`$. The algebraic identities, their degree bounds, and the use of central polynomials in quantum resetting are inherited inputs. The statements here concern their realization and converse for the specified channel-access model.

## Physical comparisons

| Work | Supplied resource and relevant result | Relation to this task |
|---|---|---|
| [Navascués (2018)](#navascues2018) | Controlled probes return after interacting with an uncontrolled target; probe interference realizes central polynomials. | Establishes the resetting mechanism in a model with accessible returned probes. |
| [Chiribella, Wilson and Chau (2021)](#cyclic2021) | Cyclic orders of completely depolarizing channels give heralded transmission with vanishing error as the number of calls grows. | Higher-dimensional communication activation is established; the displayed finite-call herald retains depolarization. |
| [Chiribella et al. (2021)](#two-slot-communication2021) | Two independent noisy-channel instances in a quantum switch enable perfect communication; an identical self-switch has a product-branch decomposition. | Communication activation and the distinction between retained-controller and reference entanglement are inherited. |
| [Das and Bandyopadhyay (2022)](#das-bandyopadhyay2022) | A quantum switch of two switches can improve the probability of heralded perfect qubit communication. | More switches improving communication is established; the fixed-channel three-call converse and randomized call-cost envelope are the narrower questions here. |
| [Kechrimparis, Moran and Kwon (2025)](#kechrimparis2025) | Recursive superswitches distill qubit channels; the paper also defines the general permutation switch. | Supplies the channel-level architecture and the inherited four-call qubit branch. |
| [Kechrimparis et al. (2024)](#kechrimparis2024) | Four-slot nested switching and a Pauli-specific branch reduction for state discrimination. | Supplies the eight-order architecture; its branch-exclusivity assumption must be checked before transfer to raw non-Pauli slots. |
| [Molitor et al. (2024)](#molitor2024) | Two-slot switching with dissipative control noise and outcome-dependent shielding. | Establishes the importance of control noise and rules out an unrestricted claim that every noisy controller destroys exactness. |
| [Regula and Takagi (2021)](#regula-takagi2021) | An accepted unwanted positive component obstructs perfect probabilistic distillation. | Supplies the obstruction principle; the routing-specific argument proves that nonzero acceptance cannot eliminate the diagonal component. |
| [Shen et al. (2026)](#shen2026) | A coherent Stinespring oracle with a specified environment basis permits Kraus-polynomial synthesis. | The returned environment is an additional resource; channel-only compilation is treated separately. |

<a id="two-slot-communication2021"></a>

## Communication baselines and the fixed-channel cost theorem

**Giulio Chiribella, Manik Banik, Some Sankar Bhattacharya, Tamal Guha, Mir Alimuddin, Arup Roy, Sutapa Saha, Sristy Agrawal and Guruprasad Kar**, *Indefinite causal order enables perfect quantum communication with zero capacity channels*, **New Journal of Physics 23**, 033039 (2021). [DOI](https://doi.org/10.1088/1367-2630/abe7a0); [arXiv:1810.10457v4](https://arxiv.org/abs/1810.10457v4), revised 5 September 2022.

Perfect communication through a coherent ordering of noisy independent channels is an established mechanism. Equation (3) gives the identical two-slot switch's commutator/anticommutator decomposition, with product controller–target branches. Its retained controller–target separability is compatible with preserving a message's entanglement with an untouched reference. These mechanisms and distinctions are inherited inputs to the present benchmark.

<a id="das-bandyopadhyay2022"></a>

**Debarshi Das and Somshubhro Bandyopadhyay**, *Quantum communication using a quantum switch of quantum switches*, **Proceedings of the Royal Society A 478**(2266), 20220231 (2022). [DOI](https://doi.org/10.1098/rspa.2022.0231); [arXiv:2111.08266v3](https://arxiv.org/abs/2111.08266v3), revised 6 October 2022 and identified by its authors as the published version; [institutional record](https://discovery.ucl.ac.uk/id/eprint/10156914/).

The paper coherently orders two quantum switches and gives cases with a higher probability of perfect heralded qubit transmission than the component switches, alongside cases without an advantage. Thus the broad phenomenon of improving heralded communication with additional switches is established prior work.

The [native three-call frontier](../science/native_three_call_frontier.md) addresses the fixed rank-deficient channel with Kraus operators $`K_0=|+\rangle\langle0|`$ and $`K_1=|-\rangle\langle1|`$, all six three-slot orders, and global unflagged logical order dephasing. It determines the full success–half-diamond-error frontier over the stated mixed, ancillary and complex controller model, with known final target unitary corrections. Four native controller flags attain a matching converse. Input independence is required for total acceptance; it is not assumed for every fine outcome.

The [randomized call-cost theorem](../science/randomized_call_cost.md) then gives the exact expected-call optimum across fresh one-, two- and three-call trials, with a matched protocol and convex-hull converse. Its unequal-noise and duration conclusions use the specified scalar logical-dephasing law. They neither price an unmeasured router implementation nor assert a physical multiqubit-noise advantage. The result is for this fixed channel and trial model, rather than a general channel-family or capacity theorem.

This three-call benchmark does not contradict the four-call minimum for full-Choi-rank qubit slots: its fixed channel has deficient Choi rank. It is also distinct from the four-call eight-order Bell-fidelity frontier with physical local controller noise. The call set, slot family, error criterion and noise model must remain attached to each statement. These comparisons establish the precise relation to the cited mechanisms without a claim of priority from a narrow search.

<a id="kechrimparis2025"></a>

## Exact qubit distillation with superswitches

**Spiros Kechrimparis, James Moran and Hyukjoon Kwon**, *Probabilistic channel distillation via indefinite causal order*, **Physical Review Research 7**, 043257 (2025). [DOI: 10.1103/mvrg-yj3l](https://doi.org/10.1103/mvrg-yj3l), [published full text](https://harvest.aps.org/v2/journals/articles/10.1103/mvrg-yj3l/fulltext).

References to sections and appendices use the journal version, published 5 December 2025. The earlier [arXiv:2501.13696v1](https://arxiv.org/abs/2501.13696v1) is dated 23 January 2025.

Appendix A, Eqs. (A1)–(A4), defines the generalized switch with all permutations and freely chosen control preparation and measurement. The successful second-order superswitch branch is the eight-word, degree-four qubit polynomial used in [the qubit construction](../science/qubit.md). Its Pauli-channel success probability is therefore an inherited special case.

Result 4 and Appendix I establish a higher-dimensional obstruction for the paper's recursive superswitch family. Arbitrary permutation coefficients extend beyond that family. The discussion in Section IX raises higher-dimensional switchlike distillation. The full-rank theorem addresses this capability within the stated order-only model; adding randomizing channel wrappers changes the supplied resources.

For the noisy-control question, Section IV gives the inherited four-call Pauli identity branch. Appendix A compares selected uniform-control permutation constructions, while the broader four-switch optimality question is posed as a conjecture. Section VII uses twirling to change general channels into the analyzed form; Section IX discusses noisy ancillary controls. These statements do not provide the arbitrary-controller corrected-fidelity frontier derived in [the controller reduction](../science/controller_optimization.md) and [the recovery frontier](../science/recovery_frontier.md).

<a id="kechrimparis2024"></a>

## Nested-switch algebra and raw-slot recovery

**Spiros Kechrimparis, James Moran, Athena Karsa, Changhyoup Lee and Hyukjoon Kwon**, *Enhancing quantum state discrimination with indefinite causal order*, **New Journal of Physics 26**, 123030 (2024). [DOI](https://doi.org/10.1088/1367-2630/ada05e); [arXiv:2406.19373v1](https://arxiv.org/html/2406.19373v1), 27 June 2024. The equation pointers here use this explicitly identified preprint version.

Section IV, Eqs. (28)–(30), defines the four-call nested process and expands its Kraus operators into commutators and anticommutators. The transition to Eqs. (31)–(32) assumes that each relevant Kraus pair cannot have both a nonzero commutator and a nonzero anticommutator. Pauli Kraus operators satisfy this condition. The raw measure–prepare family used here need not satisfy it. Appendix D studies a particular entangled-ancilla control construction; it does not optimize all allowed controller instruments.

Consequently, the eight-order router, its nested algebra, and the Pauli purification branch are inherited. The additional conclusions concern the exact feasible controller set, the complete corrected-fidelity frontier for the specified raw family and noise, and the certified [tilted readout optimum](../science/tilted_readout.md). The tilted proof uses a diagonal operator majorant even when the correction matrices do not commute. It is not a transplantation of Pauli branch exclusivity.

An output rotation after every supplied slot converts this particular raw family into a Pauli channel. Such wrappers alter the supplied slot resource and are excluded here. The allowed final herald-dependent target unitary acts after all routed calls and cannot generally perform those intermediate rotations. The [routing contract](../science/noisy_routing_model.md) makes this difference explicit. The optimization is over the specified eight orders, rather than all 24 four-slot permutations or arbitrary target circuits.

<a id="molitor2024"></a>

## Controller noise and its operational placement

**Otavio A. D. Molitor, André H. A. Malavazi, Roberto Dobal Baldijão, Alexandre C. Orthey Jr., Ismael L. Paiva and Pedro R. Dieguez**, *Quantum switch instabilities with an open control*, **Communications Physics 7**, 373 (2024). [DOI](https://doi.org/10.1038/s42005-024-01843-y); [published full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC11576514/); [arXiv:2405.09631v1](https://arxiv.org/abs/2405.09631v1).

This paper treats a two-slot switch whose control interacts with a thermal environment through a collisional model. Its conditional states and success probabilities respond differently to the bath: in the zero-temperature limit, a selected conditional state can be shielded while its probability decreases. Control noise therefore does not justify a universal impossibility claim by itself.

Here the noise is the stated unobserved Schur channel on the actual routing register, immediately before or after routing, with no accessible flag or replacement coherence. The [exact-recovery obstruction](../science/noisy_recovery_obstruction.md) requires a strictly positive Schur lower bound and full-Choi-rank independent slots. Protected sectors, accessible error flags, and other noise placements require a separate argument.

<a id="regula-takagi2021"></a>

## Accepted positive components and exact-recovery obstructions

**Bartosz Regula and Ryuji Takagi**, *Fundamental limitations on distillation of quantum channel resources*, **Nature Communications 12**, 4411 (2021). [DOI](https://doi.org/10.1038/s41467-021-24699-0); [published article](https://www.nature.com/articles/s41467-021-24699-0); [arXiv:2010.11942v3](https://arxiv.org/abs/2010.11942v3).

The probabilistic extension around main-text Eq. (34), with Supplementary Note 5, Theorem 21, identifies the obstruction caused by an unwanted component that survives the accepted transformation. The discussion allows exact probabilistic purification when postselection removes that component. Full Choi rank alone is therefore not a blanket prohibition on exact heralded recovery.

The model-specific step here is proving that every nonzero feasible controller effect retains positive diagonal routing weight when the Schur matrix has a strictly positive lower bound. That diagonal contribution is a positive mixture of definite-order compositions with full Choi rank. The accepted map consequently has full Choi rank, which precludes an exact unitary channel. The [proof](../science/noisy_recovery_obstruction.md) credits the positive-component principle while establishing its unavoidable-acceptance premise for this routing model.

<a id="matrix-factorization"></a>

## The matrix-factorization ingredient

**Moody T. Chu**, *On the first degree Fejér–Riesz factorization and its applications to X + A∗X⁻¹A = Q*, **Linear Algebra and its Applications 489**, 123–143 (2016). [DOI](https://doi.org/10.1016/j.laa.2015.09.051).

Matrix Fejér–Riesz factorization is an inherited tool for representing a positive matrix Laurent polynomial on the unit circle. The [controller dual](../science/controller_optimization.md) uses the degree-one case to express the correction-circle inequality through a finite positive-semidefinite block constraint. The feasibility theorem and the channel-specific attaining frontier require their separate arguments.

<a id="retained-and-accepted-entanglement"></a>

## Retained and accepted entanglement

The sister project [Quantum Order Entanglement](https://github.com/GoGoKo699/Quantum-Order-Entanglement) owns the complete [noisy protocol-separation proof](https://github.com/GoGoKo699/Quantum-Order-Entanglement/blob/main/science/noisy_separation.md). It compares a specified unconditional three-order controller–target output with the accepted reference–target output of the fixed four-call recovery branch. The channel family and per-bit noise strength agree; the bipartitions, call counts, controller sizes and conditioning differ. This is not an optimized equal-resource separation or a classification of all retained-controller protocols.

The proof uses **Asher Peres**, *Separability Criterion for Density Matrices*, **Physical Review Letters 77**, 1413–1415 (1996), [DOI](https://doi.org/10.1103/PhysRevLett.77.1413), and **Michał Horodecki, Paweł Horodecki and Ryszard Horodecki**, *Separability of mixed states: necessary and sufficient conditions*, **Physics Letters A 223**, 1–8 (1996), [DOI](https://doi.org/10.1016/S0375-9601(96)00706-2). PPT is sufficient on the retained occupied 3 × 2 support; the accepted two-qubit output is witnessed by Bell fidelity above one half. The latter is used as a sufficient witness, not as an entanglement classification for arbitrary complex controllers.

These comparisons identify inherited mechanisms and the additional conclusions at their stated scope. They do not establish publication priority from the absence of a matching result in a targeted search.

<a id="navascues2018"></a>

## Central-polynomial quantum resetting

**Miguel Navascués**, *Resetting Uncontrolled Quantum Systems*, **Physical Review X 8**, 031008 (2018). [DOI: 10.1103/PhysRevX.8.031008](https://doi.org/10.1103/PhysRevX.8.031008), [arXiv:1710.02470v4](https://arxiv.org/abs/1710.02470v4), revised 16 July 2018.

The target interacts with probes that return to a controlled laboratory. Joint probe processing realizes central polynomials in blocks of the interaction unitary. The main text and Appendix A establish exact resetting in arbitrary dimension; Appendix C relates resetting protocols to central polynomials and discusses success limitations. The qubit anticommutator of commutators used here also appears in footnote 16.

This supplies the quantum central-polynomial mechanism and a relation between polynomial degree and interaction count. In the present [realization](../science/polynomial_realization.md), interference combines orders within each fixed Kraus tuple, while distinct inaccessible noise histories are summed incoherently. Full Choi rank supplies the spanning property used in both directions of the channel theorem.

<a id="cyclic2021"></a>

## Cyclic orders and finite heralded noise

**Giulio Chiribella, Matt Wilson and H. F. Chau**, *Quantum and Classical Data Transmission through Completely Depolarizing Channels in a Superposition of Cyclic Orders*, **Physical Review Letters 127**, 190502 (2021). [DOI: 10.1103/PhysRevLett.127.190502](https://doi.org/10.1103/PhysRevLett.127.190502), [arXiv:2005.00618v5](https://arxiv.org/abs/2005.00618v5), revised 7 October 2021.

For $`N`$ cyclic orders of $`N`$ completely depolarizing channels, the normalized heralded channel in Eq. (8) is

```math
\mathcal E_0(\rho)
=\frac{N-1}{N-1+d^2}\rho
+\frac{d^2}{N-1+d^2}\frac{I_d}{d}.
```

The depolarizing weight is positive at every finite $`N`$ and tends to zero as $`N`$ increases. Thus this construction already gives arbitrarily high heralded fidelity in arbitrary dimension. Appendices E and F describe two-way-assisted communication consequences. The question considered here is finite-call exact identity heralding and its minimum call count. The [operational interpretation](../science/interpretation.md) distinguishes a successful identity branch from an unassisted capacity claim.

<a id="shen2026"></a>

## Coherent Kraus-polynomial compilation

**Zheyu Shen, Yusen Wu, Xiao Yuan, Xiao-Ming Zhang and Yukun Zhang**, *Quantum Algorithms for Multivariable Polynomial Transformations: From Efficient Synthesis to Quantum Channel Transformations*, [arXiv:2610.08714v1](https://arxiv.org/abs/2610.08714v1), 6 October 2026. [Full text](https://arxiv.org/html/2610.08714v1).

Section 5.2, Eq. (5.8), supplies a complete coherent Stinespring unitary, a fixed returned environment basis, and an identity sector. Corollary 5.5 compiles polynomial families in its Kraus amplitudes. Example 5.6 and the discussion following Remark 5.7 show that a change of Kraus basis can change the selected operation, including removing a state-preserving branch. Section 5.3 instead compiles a causal controller from supplied Choi data, allowing processing gates on signal and memory.

The order-only map in [the realization theorem](../science/polynomial_realization.md) is invariant under independent Kraus-representation changes in every channel slot. Its control acts on order labels, and the noise environments are unavailable. These distinct access assumptions determine which polynomial operations can be implemented.

<a id="central-polynomials"></a>

## Central polynomials and degree bounds

A central polynomial here has zero constant term, takes scalar-matrix values on every matrix substitution, and is not a polynomial identity. Over characteristic zero, taking a nonidentity multihomogeneous component and then fully polarizing preserves centrality and nonidentity without increasing the degree. Consequently the minimum degree can be expressed using multilinear polynomials, as required by the compulsory-slot architecture.

### Current degree statements

**Vesselin Drensky and Boyan Kostadinov**, *Central polynomials of minimal degree for matrices*, [arXiv:2601.07750v2](https://arxiv.org/abs/2601.07750v2), revised 4 July 2026. The introduction records the minimum degrees and original constructions; Section 2 gives the characteristic-zero linearization facts.

| Dimension | Algebraic input | Operational consequence in the stated model |
|---|---|---|
| $`d=2`$ | $`c(2)=4`$ | A four-call minimum; [explicit qubit branch](../science/qubit.md). |
| $`d=3`$ | $`c(3)=8`$ | An eight-call minimum; [explicit qutrit branch](../science/qutrit.md). |
| $`d=4`$ | $`c(4)\leq13`$ | A thirteen-call upper bound. Equality is not established by the cited results. |
| $`d>2`$ | $`c(d)\leq(d-1)^2+4`$ | The same general upper bound on the minimum call count. |

The new degree-twelve exclusion in Drensky and Kostadinov is restricted to **two-variable** polynomials for $`4\times4`$ matrices. It does not establish the unrestricted equality $`c(4)=13`$. The [query theorem](../science/query_threshold.md) characterizes a minimum; it does not assert that every larger number of compulsory slots is feasible.

### Explicit qutrit polynomial and historical minimality

**Vesselin Drensky and Azniv Kasparian**, *A new central polynomial for $`3\times3`$ matrices*, **Communications in Algebra 13**(3), 745–752 (1985). [DOI: 10.1080/00927878508823188](https://doi.org/10.1080/00927878508823188), [author-hosted copy](https://www.researchgate.net/profile/Azniv-Kasparian/publication/243045198_A_new_central_polynomial_for_3_3_matrices/links/0deec52050d137da71000000/A-new-central-polynomial-for-3-3-matrices.pdf).

Pages 746–747 give the degree-eight multilinear polynomial, pages 748–751 prove centrality and nonidentity, and page 752 gives its scalar trace expression. The [qutrit construction](../science/qutrit.md) uses this polynomial, reproduces its algebraic proof, and supplies the signed-word coefficients and exact evaluations needed for the channel realization. The central polynomial and its algebraic proof are inherited.

The original lower-bound attribution is **V. Drensky and A. Kasparian**, *Polynomial identities of eighth degree for $`3\times3`$ matrices*, **Annuaire de l'Université de Sofia, Faculté de Mathématiques et Mécanique, Livre 1, Mathématiques 77** (1983), 175–195. This metadata follows reference [19] of Drensky and Kostadinov v2. The historical article's full proof is not available in this account. The exclusion of degrees below eight is used as an established result reported in the current account by one of the original authors; its minimality proof is not independently reproduced here.

### General existence and the improved upper bound

**Edward Formanek**, *Central polynomials for matrix rings*, **Journal of Algebra 23**, 129–132 (1972), [DOI: 10.1016/0021-8693(72)90050-6](https://doi.org/10.1016/0021-8693%2872%2990050-6), supplies a degree-$`d^2`$ central polynomial. Polarization yields the multilinear form used for general existence in [the realization theorem](../science/polynomial_realization.md).

**Vesselin Drensky**, *New central polynomials for the matrix algebra*, **Israel Journal of Mathematics 92**, 235–248 (1995), [DOI: 10.1007/BF02762079](https://doi.org/10.1007/BF02762079), gives the improved degree $`(d-1)^2+4`$ for $`d>2`$. The interpretation of this construction as an upper bound, particularly at $`d=4`$, follows the current Drensky–Kostadinov account above.

### Polynomial identities and the elementary lower bound

**Vesselin Drensky**, *Weak polynomial identities and their applications*, **Communications in Mathematics 29**(2), 291–324 (2021), [DOI: 10.2478/cm-2021-0022](https://doi.org/10.2478/cm-2021-0022), [full text](https://dml.cz/bitstream/handle/10338.dmlcz/149196/ActaOstrav_29-2021-2_11.pdf), surveys the central-polynomial constructions and their degree bounds. Theorem 3 states the minimum polynomial-identity degree for matrix algebras.

The original source is **A. S. Amitsur and J. Levitzki**, *Minimal identities for algebras*, **Proceedings of the American Mathematical Society 1**(4), 449–463 (1950), [DOI: 10.1090/S0002-9939-1950-0036751-9](https://doi.org/10.1090/S0002-9939-1950-0036751-9). There is no nonzero polynomial identity of degree below $`2d`$ for $`M_d(\mathbb C)`$. For a positive-degree central polynomial $`P`$, the free polynomial $`[P,X_{m+1}]`$ is nonzero and is an identity on $`M_d`$. This gives $`c(d)\geq2d-1`$. The [query-threshold proof](../science/query_threshold.md) also gives a direct matrix-unit argument for the required identity lower bound.

<a id="regev-procesi"></a>

## Regev's polynomial and its determinant coefficient

**Claudio Procesi**, *A construction of swap or switch polynomials*, **Advances in Mathematics 400**, 108308 (2022), [DOI: 10.1016/j.aim.2022.108308](https://doi.org/10.1016/j.aim.2022.108308), gives the Regev determinant evaluation in Theorem 4.13. Equation references here use [arXiv:2102.10657v2](https://arxiv.org/abs/2102.10657v2), revised 17 September 2022, whose version note identifies a correction to the previous version. “Swap” and “switch” in this title refer to a tensor transposition operator.

An independent published statement of the same evaluation is **Claudio Procesi**, *A note on the Formanek Weingarten function*, **Note di Matematica 41**(1), 69–109 (2021), [DOI: 10.1285/i15900932v41n1p69](https://doi.org/10.1285/i15900932v41n1p69), [full text](https://ese-journals.unisalento.it/index.php/notemat/article/download/24154/20198). Theorem 7, Eq. (2.85), on page 99 states

```math
R_d(X;Y)=a_d\det[X_1,\ldots,X_{d^2}]
                    \det[Y_1,\ldots,Y_{d^2}]I_d,
\qquad
a_d=\frac{(-1)^{d-1}C_d^2}{(d!)^2(2d-1)},
\qquad
C_d=\frac{\prod_{j=1}^{d}(2j-1)!}{\prod_{j=1}^{d-1}j!}.
```

Here each determinant uses coordinates in the same matrix-unit basis. The polynomial is defined by unnormalized double alternation in Eq. (2.84). Equation (1.65), page 94, fixes the magnitude of $`C_d`$; its possible sign disappears after squaring. The proof on page 100 uses Lemma 2 and Proposition 10 on page 98, together with the full-cycle Weingarten value in Eq. (1.24), page 76.

A later [corrigendum](https://doi.org/10.1016/j.aim.2023.108970) to the 2022 paper appeared in **Advances in Mathematics 420**, 108970 (2023). Its full text is unavailable in this account, so no assertion is made that arXiv v2 incorporates that later correction. The determinant formula used here is independently supported by the 2021 theorem above.

The [probability account](../science/probability.md) separately supplies controller normalization, the sum over Kraus tuples, and the Choi-determinant calculation. These give success probabilities for specified constructions. Neither the Regev construction nor the explicit qubit and qutrit branches is asserted to maximize success probability over all allowed protocols.

## Bibliography

The accompanying [BibTeX bibliography](references.bib) records the journal metadata, persistent identifiers, and the versioned preprints used for the statements above.
