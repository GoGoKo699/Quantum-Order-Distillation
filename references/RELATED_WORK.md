# Related work and algebraic sources

The task studied here is heralded exact identity transmission through independent fresh channel slots in a generalized quantum switch. Every order visits every prescribed slot once. The control is prepared independently of the target input and measured after the calls; the noise environments remain inaccessible. The [access model](../science/model.md) specifies these assumptions.

The [polynomial realization](../science/polynomial_realization.md) and [query threshold](../science/query_threshold.md) connect this task to an established algebraic invariant: the minimum positive degree $`c(d)`$ of a nonidentity multilinear central polynomial for $`M_d(\mathbb C)`$. For any fixed channel with full Choi rank, the minimum compulsory call count is $`c(d)`$. The algebraic identities, their degree bounds, and the use of central polynomials in quantum resetting are inherited inputs. The statements here concern their realization and converse for the specified channel-access model.

## Physical comparisons

| Work | Supplied resource and relevant result | Relation to this task |
|---|---|---|
| [Navascués (2018)](#central-polynomial-quantum-resetting) | Controlled probes return after interacting with an uncontrolled target; probe interference realizes central polynomials. | Establishes the resetting mechanism in a model with accessible returned probes. |
| [Chiribella, Wilson and Chau (2021)](#cyclic-orders-and-finite-heralded-noise) | Cyclic orders of completely depolarizing channels give heralded transmission with vanishing error as the number of calls grows. | Higher-dimensional communication activation is established; the displayed finite-call herald retains depolarization. |
| [Kechrimparis, Moran and Kwon (2025)](#exact-qubit-distillation-with-superswitches) | Recursive superswitches distill qubit channels; the paper also defines the general permutation switch. | Supplies the channel-level architecture and the inherited four-call qubit branch. |
| [Shen et al. (2026)](#coherent-kraus-polynomial-compilation) | A coherent Stinespring oracle with a specified environment basis permits Kraus-polynomial synthesis. | The returned environment is an additional resource; channel-only compilation is treated separately. |

## Exact qubit distillation with superswitches

**Spiros Kechrimparis, James Moran and Hyukjoon Kwon**, *Probabilistic channel distillation via indefinite causal order*, **Physical Review Research 7**, 043257 (2025). [DOI: 10.1103/mvrg-yj3l](https://doi.org/10.1103/mvrg-yj3l), [published full text](https://harvest.aps.org/v2/journals/articles/10.1103/mvrg-yj3l/fulltext).

References to sections and appendices use the journal version, published 5 December 2025. The earlier [arXiv:2501.13696v1](https://arxiv.org/abs/2501.13696v1) is dated 23 January 2025.

Appendix A, Eqs. (A1)–(A4), defines the generalized switch with all permutations and freely chosen control preparation and measurement. The successful second-order superswitch branch is the eight-word, degree-four qubit polynomial used in [the qubit construction](../science/qubit.md). Its Pauli-channel success probability is therefore an inherited special case.

Result 4 and Appendix I establish a higher-dimensional obstruction for the paper's recursive superswitch family. Arbitrary permutation coefficients extend beyond that family. The discussion in Section IX raises higher-dimensional switchlike distillation. The full-rank theorem addresses this capability within the stated order-only model; adding randomizing channel wrappers changes the supplied resources.

## Central-polynomial quantum resetting

**Miguel Navascués**, *Resetting Uncontrolled Quantum Systems*, **Physical Review X 8**, 031008 (2018). [DOI: 10.1103/PhysRevX.8.031008](https://doi.org/10.1103/PhysRevX.8.031008), [arXiv:1710.02470v4](https://arxiv.org/abs/1710.02470v4), revised 16 July 2018.

The target interacts with probes that return to a controlled laboratory. Joint probe processing realizes central polynomials in blocks of the interaction unitary. The main text and Appendix A establish exact resetting in arbitrary dimension; Appendix C relates resetting protocols to central polynomials and discusses success limitations. The qubit anticommutator of commutators used here also appears in footnote 16.

This supplies the quantum central-polynomial mechanism and a relation between polynomial degree and interaction count. In the present [realization](../science/polynomial_realization.md), interference combines orders within each fixed Kraus tuple, while distinct inaccessible noise histories are summed incoherently. Full Choi rank supplies the spanning property used in both directions of the channel theorem.

## Cyclic orders and finite heralded noise

**Giulio Chiribella, Matt Wilson and H. F. Chau**, *Quantum and Classical Data Transmission through Completely Depolarizing Channels in a Superposition of Cyclic Orders*, **Physical Review Letters 127**, 190502 (2021). [DOI: 10.1103/PhysRevLett.127.190502](https://doi.org/10.1103/PhysRevLett.127.190502), [arXiv:2005.00618v5](https://arxiv.org/abs/2005.00618v5), revised 7 October 2021.

For $`N`$ cyclic orders of $`N`$ completely depolarizing channels, the normalized heralded channel in Eq. (8) is

```math
\mathcal E_0(\rho)
=\frac{N-1}{N-1+d^2}\rho
+\frac{d^2}{N-1+d^2}\frac{I_d}{d}.
```

The depolarizing weight is positive at every finite $`N`$ and tends to zero as $`N`$ increases. Thus this construction already gives arbitrarily high heralded fidelity in arbitrary dimension. Appendices E and F describe two-way-assisted communication consequences. The question considered here is finite-call exact identity heralding and its minimum call count. The [operational interpretation](../science/interpretation.md) distinguishes a successful identity branch from an unassisted capacity claim.

## Coherent Kraus-polynomial compilation

**Zheyu Shen, Yusen Wu, Xiao Yuan, Xiao-Ming Zhang and Yukun Zhang**, *Quantum Algorithms for Multivariable Polynomial Transformations: From Efficient Synthesis to Quantum Channel Transformations*, [arXiv:2610.08714v1](https://arxiv.org/abs/2610.08714v1), 6 October 2026. [Full text](https://arxiv.org/html/2610.08714v1).

Section 5.2, Eq. (5.8), supplies a complete coherent Stinespring unitary, a fixed returned environment basis, and an identity sector. Corollary 5.5 compiles polynomial families in its Kraus amplitudes. Example 5.6 and the discussion following Remark 5.7 show that a change of Kraus basis can change the selected operation, including removing a state-preserving branch. Section 5.3 instead compiles a causal controller from supplied Choi data, allowing processing gates on signal and memory.

The order-only map in [the realization theorem](../science/polynomial_realization.md) is invariant under independent Kraus-representation changes in every channel slot. Its control acts on order labels, and the noise environments are unavailable. These distinct access assumptions determine which polynomial operations can be implemented.

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
