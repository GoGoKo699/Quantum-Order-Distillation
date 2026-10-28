# Eight calls for a qutrit

## Construction and attribution

The construction in this document is explicit, multilinear, and exactly verified at the coefficient level. It uses a classical central polynomial due to **Vesselin Drensky and Azniv Kasparian**, published in 1985 [DK85]. The polynomial and its centrality are inherited mathematics. Their realization in the [order-only channel model](model.md) gives the protocol studied here.

The protocol uses eight independent channel slots, each occupied by the same full-Choi-rank qutrit channel. Its order control has support on **2,880 permutations**, with coefficients of equal magnitude and signs given below. Its success branch is exactly a scalar multiple of the identity channel. Eight channel calls are optimal within the stated order-only architecture because the minimum central-polynomial degree for 3×3 matrices is eight [DK85, DK26]. The claim is not an optimal success probability or a control-dimension lower bound.

No search over all 8! permutations is required to generate the protocol. The [Python module](../verification/qutrit_degree8.py) generates the coefficients from a compact placement formula and independently regenerates them from the original exterior-algebra construction. Their exact equality is checked. The [coefficient file](../verification/degree8_coefficients.json) lists all nonzero coefficients and states its ordering convention.

## 1. Compact polynomial specification

Let X₁,…,X₈ be noncommuting variables. Put

```math
A=[X_5,X_6],\qquad B=[X_7,X_8].
```

For five possibly noncommuting gap factors U₀,…,U₄, define the alternating placement operation

```math
\mathcal C_4(X_1,X_2,X_3,X_4;U_0,\ldots,U_4)
=\sum_{\sigma\in S_4}\mathrm{sgn}(\sigma)
U_0X_{\sigma(1)}U_1X_{\sigma(2)}U_2X_{\sigma(3)}U_3X_{\sigma(4)}U_4.
\qquad\text{(1)}
```

The slots are numbered **0 through 4**. For i≠j, define U⁽ⁱʲ⁾ by placing A in gap i, B in gap j, and identity factors in the other gaps. Define V⁽ᵏ⁾ by placing AB+BA in gap k and identities elsewhere. Let

```math
m_{ij}=\begin{cases}-1,&(i,j)=(1,3)\text{ or }(3,1),\\+1,&\text{otherwise},\end{cases}
\qquad
n_k=\begin{cases}-1,&k=1\text{ or }3,\\+1,&\text{otherwise}.\end{cases}
```

Then the Drensky–Kasparian polynomial is

```math
p_8(X_1,\ldots,X_8)
=\sum_{\substack{0\le i,j\le4\\i\ne j}}m_{ij}\,
\mathcal C_4(X_1,\ldots,X_4;U^{(ij)})
+\sum_{k=0}^{4}n_k\,
\mathcal C_4(X_1,\ldots,X_4;V^{(k)}).
\qquad\text{(2)}
```

This is the explicit formula on printed page 747 of [DK85], with notation changed for readability. Every monomial contains each of X₁,…,X₈ exactly once. The identity gap factors in (1) are simply omitted when writing a monomial; they require no extra physical slot or target operation.

### Support and normalization

The first sum has 20 ordered placements, 24 permutations of the first four variables, and four terms from expanding the two commutators. It contributes 1,920 distinct words. The second sum has five gaps, two orders AB and BA, 24 permutations, and four commutator expansions, contributing 960 distinct words. The two types cannot overlap: in the first, at least one of X₁,…,X₄ separates the A and B blocks; in the second, their blocks are adjacent. Within either type the placement, variable permutation, and internal commutator order are uniquely recoverable from the word.

Consequently the support size is exactly

```math
20\cdot24\cdot4+5\cdot2\cdot24\cdot4=2880.
```

Every nonzero coefficient is ±1. Alternation in the first four variables pairs positive and negative coefficients, so there are 1,440 of each. The coefficient norm used by the switch compiler is therefore

```math
L=\sum_\pi|c_\pi|=2880.
\qquad\text{(3)}
```

The control dimension 2,880 is a **sufficient explicit implementation**, not a minimality statement. It fits within a 12-qubit Hilbert space, but that observation does not supply an efficient routing circuit or account for the gates needed to prepare and measure the prescribed states.

## 2. Centrality and its scalar value

For arbitrary complex 3×3 matrices X₁,…,X₈, the polynomial obeys

```math
p_8(X_1,\ldots,X_8)=\lambda(X_1,\ldots,X_8)I_3,
\qquad\text{(4)}
```

where

```math
\lambda(X_1,\ldots,X_8)
=2\sum_{\sigma\in S_4}\mathrm{sgn}(\sigma)
\mathrm{tr}\!\left(X_{\sigma(1)}X_{\sigma(2)}[X_5,X_6]\right)
\mathrm{tr}\!\left(X_{\sigma(3)}X_{\sigma(4)}[X_7,X_8]\right).
\qquad\text{(5)}
```

Equation (5), equivalent to the scalar formula on printed page 752 of [DK85], is convenient for exact evaluation and probability calculation. It is not a replacement for the actual order polynomial in a physical implementation: the switch implements (2), while (5) evaluates its resulting scalar.

### 2.1 Exterior-algebra definition

Let G be the exterior algebra on anticommuting generators g₁,…,g₈, and define matrices with odd exterior entries

```math
X=\sum_{i=1}^{4}X_i g_i,\qquad
Y=X_5g_5+X_6g_6,\qquad
Z=X_7g_7+X_8g_8.
```

Set

```math
P=X^2+Y^2+Z^2,\qquad Q=X(Y^2+Z^2)X,
```

and

```math
F=P^4-XP^3X+P^2Q+PQP+QP^2.
\qquad\text{(6)}
```

The coefficient of g₁g₂⋯g₈ in F is p₈. In particular, it is already multilinear; **no further polarization or factorial normalization is needed**.

The component contributing to the top exterior coefficient has four X factors, two Y factors, and two Z factors. It consists of 12 patterns from P⁴, six patterns with negative sign from XP³X, and 12 patterns from the three Q terms. Assigning the four X labels in all orders and the two Y and Z labels in both orders, with exterior signs, gives precisely (2). The independent generators in the accompanying code check this equality as an equality in the free algebra, not merely on sampled matrices.

### 2.2 Self-contained centrality proof

The even subalgebra of G is commutative. Hence P and Q are 3×3 matrices over a commutative algebra, to which Cayley–Hamilton applies. For an even-length word of odd matrices, moving its first odd factor around the trace changes the sign. In particular,

```math
\mathrm{tr}(X^{2k})=
\mathrm{tr}(Y^{2k})=
\mathrm{tr}(Z^{2k})=0.
```

A word with more than four X factors, more than two Y factors, or more than two Z factors vanishes, because it repeats an exterior generator from the corresponding group. Define the even scalar elements

```math
\begin{aligned}
a&=\mathrm{tr}(X^2Y^2),& b&=\mathrm{tr}(X^2Z^2),&
c&=\mathrm{tr}(Y^2Z^2),\\
u&=\mathrm{tr}(X^4Y^2),& v&=\mathrm{tr}(X^4Z^2),&
w&=\mathrm{tr}(X^2Y^2Z^2)+\mathrm{tr}(X^2Z^2Y^2).
\end{aligned}
```

Graded cyclicity and the degree restrictions yield

```math
\begin{gathered}
\mathrm{tr}P=0,\qquad \mathrm{tr}Q=-a-b,\qquad
\tfrac12\mathrm{tr}P^2=a+b+c,\\
\mathrm{tr}(PQ)=-u-v,\qquad
\tfrac13\mathrm{tr}P^3=u+v+w,\qquad
\mathrm{tr}(P^2Q)=0.
\qquad\text{(7)}
\end{gathered}
```

For clarity, the last cancellation is especially easy to track after discarding forbidden degrees. The remaining four traces are

```math
\mathrm{tr}(X^2Y^2XZ^2X),\quad
\mathrm{tr}(Y^2X^3Z^2X),\quad
\mathrm{tr}(X^2Z^2XY^2X),\quad
\mathrm{tr}(Z^2X^3Y^2X).
```

Moving the last X of the first and third terms to the front produces a minus sign; moving the initial Y² or Z² of the second and fourth terms to the end produces no sign. They therefore cancel in pairs. Similarly, the mixed terms in tr(PQ) cancel because Y²X and Z²X both have odd degree.

Cayley–Hamilton for traceless P gives

```math
P^3=\tfrac12\mathrm{tr}(P^2)P+\tfrac13\mathrm{tr}(P^3)I.
\qquad\text{(8)}
```

The component linear in Q of Cayley–Hamilton for P+tQ gives

```math
\begin{aligned}
P^2Q+PQP+QP^2={}&\mathrm{tr}(Q)P^2+
\mathrm{tr}(PQ)P+\tfrac12\mathrm{tr}(P^2)Q\\
&+\left[\mathrm{tr}(P^2Q)-
\tfrac12\mathrm{tr}(P^2)\mathrm{tr}(Q)\right]I.
\end{aligned}
\qquad\text{(9)}
```

Using (7)–(8) and retaining only allowed exterior degrees,

```math
\begin{aligned}
P^4-XP^3X={}&
a(X^2Z^2+Z^2X^2-XZ^2X)\\
&+b(X^2Y^2+Y^2X^2-XY^2X)+uZ^2+vY^2.
\end{aligned}
\qquad\text{(10)}
```

The nonscalar terms on the right of (9) are the negatives of (10). Its scalar term is

```math
(a+b+c)(a+b)I=2abI,
```

since a²=b²=ac=bc=0 by the degree restrictions. Thus

```math
F=2\mathrm{tr}(X^2Y^2)\mathrm{tr}(X^2Z^2)I_3.
\qquad\text{(11)}
```

Taking the top exterior coefficient proves (4)–(5). This proof follows the mechanism in [DK85], with the cancellations organized explicitly for this application.

## 3. Exact nonidentity witness

Take the one-based matrix-unit tuple from printed page 751 of [DK85]:

```math
(X_1,\ldots,X_8)=
(E_{11},E_{12},E_{23},E_{31},E_{12},E_{23},E_{32},E_{21}).
\qquad\text{(12)}
```

Then A=E₁₃ and B=E₃₁. In (5), the first trace is nonzero only for the ordered pair (X₄,X₁), and the second trace for the remaining ordered pair (X₂,X₃). Both traces equal one. The contributing permutation is (4,1,2,3), whose sign is −1. Consequently

```math
p_8(E_{11},E_{12},E_{23},E_{31},E_{12},E_{23},E_{32},E_{21})=-2I_3.
\qquad\text{(13)}
```

This exact evaluation proves that p₈ is not a polynomial identity. Repetition of some matrix values across different slots is permitted; the polynomial remains multilinear in its eight separate variables.

## 4. Compilation into a channel protocol

Let $`\mathcal N(\rho)=\sum_aK_a\rho K_a^\dagger`$ be a qutrit channel of Choi rank nine. Use eight independent fresh copies in the standard eight-switch. Restrict its order control to the 2,880 words in the support of p₈. All word labels here refer to matrix-product order, left to right; the corresponding physical chronological order is reversed.

Prepare and postselect on

```math
|\alpha\rangle=\frac1{\sqrt{2880}}\sum_{\pi\in\mathrm{supp}(p_8)}|\pi\rangle,
\qquad
|\beta\rangle=\frac1{\sqrt{2880}}\sum_{\pi\in\mathrm{supp}(p_8)}c_\pi|\pi\rangle.
\qquad\text{(14)}
```

Both are normalized. The success Kraus operator for tuple $`\mathbf a=(a_1,\ldots,a_8)`$ is

```math
A_{\mathbf a}=\frac{p_8(K_{a_1},\ldots,K_{a_8})}{2880}
=\frac{\lambda(K_{a_1},\ldots,K_{a_8})}{2880}I_3.
\qquad\text{(15)}
```

Hence the unnormalized success map is exactly $`q(\mathcal N)\,\mathrm{id}_3`$, with

```math
q(\mathcal N)=\frac1{2880^2}
\sum_{a_1,\ldots,a_8}\left|
\lambda(K_{a_1},\ldots,K_{a_8})\right|^2.
\qquad\text{(16)}
```

Equation (16) is input independent, including for an input entangled with a reference. Because the Kraus family spans all 3×3 matrices, multilinearity and (13) imply that at least one Kraus tuple has nonzero scalar. Thus $`q(\mathcal N)>0`$ for every full-rank qutrit channel. The scalar-tuple expression is Kraus-gauge invariant, as also follows directly from the channel-defined switch.

The control preparation and postselection are channel independent; no channel tomography, environment measurement, or Kraus-label access is part of the protocol.

The [probability analysis](probability.md) obtains the exact depolarizing-channel value

```math
q(\Delta_3)=\frac1{43740}\approx2.28624\times10^{-5}.
\qquad\text{(17)}
```

Its two independent certificates are the complete 2,880-word contraction and the 96-term expansion of (5). Their derivations and exact cycle histograms are given in the probability account. The [trace-contraction script](../verification/qutrit_trace_crosscheck.py) and its [exact output](../verification/qutrit_trace_results.json) reproduce the second calculation. This is the probability for the explicitly stated control preparation and postselection; it is not an optimized probability claim.

## 5. Why eight calls are optimal

The [query-threshold theorem](query_threshold.md) identifies the minimum number of calls with the minimum positive central-polynomial degree. Drensky and Kasparian's degree-eight construction and the absence of lower-degree central polynomials for three-by-three matrices are established algebraic results. The current primary-author account [DK26, Introduction] confirms both statements for arbitrary numbers of variables.

Thus the explicit polynomial above attains the exact eight-call threshold. Its finite verification does not reprove the imported representation-theoretic lower bound. General degree bounds and the unresolved four-dimensional minimum are stated in [the query account](query_threshold.md#5-known-degrees-and-their-quantum-meaning).

## 6. Verification record and reproducibility

The accompanying script uses only Python's standard library and exact integers. It verifies:

1. Equality of the Capelli-placement coefficient dictionary and an independently generated exterior-algebra coefficient dictionary.
2. Degree eight and exactly one occurrence of every variable in every word.
3. Exact support 2,880, coefficient histogram ±1 with 1,440 of each, and coefficient norm L=2,880.
4. The original matrix-unit nonidentity witness (13), both from the full polynomial and from the trace formula.
5. Twelve deterministic random integer 3×3 tuples for which full polynomial evaluation agrees exactly with (5)I₃.
6. A lower-dimensional vanishing check on one integer 2×2 tuple, consistent with embedding M₂ as a proper corner of M₃.

The sampled matrix checks catch implementation and transcription mistakes; they are not presented as a proof of centrality. Centrality is proved in Section 2, and nonidentity by Section 3. The two coefficient-generation routes agree before any matrix substitution, providing an exact algebraic transcription check.

Reproduce the exports and report with:

```bash
python3 verification/qutrit_degree8.py --export /tmp/qutrit_coefficients.json --report /tmp/qutrit_verification.json
```

The JSON words use labels 1,…,8 in mathematical matrix-product order. The Python module uses zero-based labels internally. Keep this conversion explicit when compiling routing instructions.

## 7. Scientific limits

The construction establishes exact finite-query capability and query optimality in the declared architecture. It does not establish an optimal coefficient norm, success probability, number of coherent orders, routing complexity, or robustness to control noise. The success probability can be small, and a heralded identity branch is not by itself a positive unassisted quantum-capacity theorem. Its immediate operational use is heralded perfect entanglement distribution from fresh attempts under the supplied switch-access assumptions.

The eight-call result concerns fresh channel copies with inaccessible environments. It does not extend automatically to correlated repeated uses of one device, whose Kraus tuples need not factor independently. It also does not constrain a more general architecture with target encodings, interleaved operations, or a target-carrying bypass memory.

## Sources

**[DK85]** V. Drensky and A. Kasparian, *A new central polynomial for 3×3 matrices*, Communications in Algebra **13**(3), 745–752 (1985), [DOI:10.1080/00927878508823188](https://doi.org/10.1080/00927878508823188). The placement formula is on p.747, the centrality argument on pp.748–751, the nonzero evaluation on p.751, and the scalar trace expression on p.752.

**[DK26]** V. Drensky and B. Kostadinov, *Central polynomials of minimal degree for matrices*, [arXiv:2601.07750v2](https://arxiv.org/abs/2601.07750v2), revised 4 July 2026. The introduction records the degree minima and their original sources. Section 2 explains characteristic-zero linearization.

The [related-work account](../references/RELATED_WORK.md#central-polynomials-and-degree-bounds) and [bibliography](../references/references.bib) give the broader attribution and distinguish the inherited algebra from the channel realization.
