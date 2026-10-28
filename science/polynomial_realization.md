# Turning order interference into a channel

The controller in the [channel model](model.md) changes the order of complete channel uses. Its measurement can make different products of Kraus operators interfere while each inaccessible environment history remains a separate Kraus label. This is the link between a generalized switch and a noncommutative polynomial.

## 1. The normalized switch

Let $`\mathcal N_j`$ have Kraus operators $`K^{(j)}_a`$, for $`j=1,\ldots,m`$. A word $`\pi=(\pi_1,\ldots,\pi_m)`$ lists the operators in matrix-product order, so the rightmost slot acts first. Define

```math
A_{\boldsymbol a,\pi}
=K^{(\pi_1)}_{a_{\pi_1}}\cdots K^{(\pi_m)}_{a_{\pi_m}},
\qquad
W_{\boldsymbol a}
=\sum_{\pi\in S_m}A_{\boldsymbol a,\pi}
 \otimes\lvert\pi\rangle\langle\pi\rvert.
```

These are Kraus operators of a channel on the target and its order control. Orthogonality of the control projectors gives

```math
\sum_{\boldsymbol a}W_{\boldsymbol a}^{\dagger}W_{\boldsymbol a}
=\sum_{\pi\in S_m}
 \left(\sum_{\boldsymbol a}
 A_{\boldsymbol a,\pi}^{\dagger}A_{\boldsymbol a,\pi}\right)
 \otimes\lvert\pi\rangle\langle\pi\rvert
=I\otimes I.
```

The inner sum is the trace-preservation identity for a composition of the input channels. This is the standard generalized-switch construction used in [Kechrimparis, Moran and Kwon, Appendix A](../references/RELATED_WORK.md#kechrimparis2025). It is an assumed higher-order resource. Counting its channel slots does not give a gate-level implementation cost.

### Independence and Kraus representations

Suppose a Kraus family is changed by an isometry,
$`\widetilde K^{(j)}_b=\sum_a u^{(j)}_{ba}K^{(j)}_a`$.
Every order uses the operator from slot $`j`$ exactly once. Scalar coefficients therefore commute through the word, and

```math
\widetilde W_{\boldsymbol b}
=\sum_{\boldsymbol a}
\left(\prod_{j=1}^m u^{(j)}_{b_j a_j}\right)
W_{\boldsymbol a}.
```

The tensor product of these isometries is again an isometry. The switch channel and every control-measurement branch are consequently independent of the chosen Kraus representations.

The tuple $`\boldsymbol a`$ labels independent environments. It is never measured or coherently selected by the protocol. Replacing independent environments by a shared noise process changes the model; a [full-rank marginal counterexample](interpretation.md#correlated-instances) shows why the distinction matters.

## 2. Every multilinear order polynomial is realizable

Let

```math
P(X_1,\ldots,X_m)
=\sum_{\pi\in S_m}c_\pi X_{\pi_1}\cdots X_{\pi_m},
\qquad
L=\sum_\pi |c_\pi|>0.
```

Coefficients may vanish, in which case their order states need not be populated. Prepare a normalized control state $`\lvert\alpha\rangle`$ and postselect the rank-one effect $`\lvert\beta\rangle\langle\beta\rvert`$, with

```math
\alpha_\pi=\sqrt{\frac{|c_\pi|}{L}},
\qquad
\beta_\pi=e^{-i\arg(c_\pi)}
           \sqrt{\frac{|c_\pi|}{L}}.
```

For zero coefficients both amplitudes are zero. Each vector has squared norm one, and its successful target Kraus operator is

```math
B_{\boldsymbol a}
=\langle\beta|W_{\boldsymbol a}|\alpha\rangle
=\frac{1}{L}
P\!\left(K^{(1)}_{a_1},\ldots,K^{(m)}_{a_m}\right).
\qquad\text{(1)}
```

Thus any multilinear polynomial using every variable once in each monomial has a physical realization as a trace-nonincreasing branch. The coefficient normalization is part of the construction.

For this fixed coefficient representation and one pure preparation/one rank-one postselection, the scale $`1/L`$ is optimal. The [probability account](probability.md#coefficient-normalization) proves that normalization bound and distinguishes it from optimization across polynomial representatives or several successful outcomes.

## 3. Centrality gives exact recovery

A positive-degree multilinear polynomial is **central for $`M_d(\mathbb C)`$** when its value on every tuple of $`d\times d`$ matrices is a scalar matrix. It is **nonidentity** when at least one such value is nonzero.

### Theorem 1. Universal exact branch

Let $`P`$ be a nonidentity multilinear central polynomial for $`M_d(\mathbb C)`$ of degree $`m`$. The controller above gives

```math
\mathcal S_{\mathrm{succ}}
(\mathcal N_1,\ldots,\mathcal N_m)
=p(\mathcal N_1,\ldots,\mathcal N_m)\,\mathrm{id}_d
```

for arbitrary $`d`$-dimensional channels. If every input channel has full Choi rank, then $`p>0`$. The same control preparation and measurement work for all those channels.

**Proof.** Centrality supplies scalars $`\lambda_{\boldsymbol a}`$ such that
$`P(K^{(1)}_{a_1},\ldots,K^{(m)}_{a_m})=\lambda_{\boldsymbol a}I`$.
Equation (1) therefore gives

```math
\mathcal S_{\mathrm{succ}}(\rho)
=\sum_{\boldsymbol a}B_{\boldsymbol a}\rho B_{\boldsymbol a}^{\dagger}
=\left(\frac{1}{L^2}
\sum_{\boldsymbol a}|\lambda_{\boldsymbol a}|^2\right)\rho,
\qquad
p=\frac{\sum_{\boldsymbol a}|\lambda_{\boldsymbol a}|^2}{L^2}.
\qquad\text{(2)}
```

This is an identity between linear maps, so for every reference system $`R`$ and every joint state $`\rho_{RT}`$,

```math
(\mathrm{id}_R\otimes\mathcal S_{\mathrm{succ}})(\rho_{RT})
=p\,\rho_{RT}.
```

The success probability is independent of the input. Complete positivity and trace nonincrease give $`0\le p\le1`$.

Use the unnormalized Choi convention

```math
J(\mathcal N_j)
=\sum_a\lvert K^{(j)}_a\rangle\!\rangle
       \langle\!\langle K^{(j)}_a\rvert,
\qquad \mathrm{Tr}J(\mathcal N_j)=d.
```

The range of this positive matrix is the span of its Kraus vectors. Full Choi rank is therefore equivalent to
$`\mathrm{span}\{K^{(j)}_a\}=M_d(\mathbb C)`$.

If every scalar in (2) vanished, write arbitrary $`X_j`$ as a linear combination of the Kraus operators of slot $`j`$. Multilinearity would then force $`P(X_1,\ldots,X_m)=0`$ for every matrix tuple, contradicting nonidentity. At least one scalar is nonzero, and hence $`p>0`$. No channel description enters the controller coefficients. $`\square`$

Full rank is a sufficient condition for this universal positivity statement, not a necessary condition for every possible channel to be distillable. The [boundary examples](interpretation.md#rank-deficient-channels) distinguish these questions.

## 4. Why finite-dimensional constructions exist

Existence of central polynomials is inherited algebra. Formanek's construction supplies a positive-degree homogeneous central polynomial of degree $`d^2`$. Later constructions give degree $`(d-1)^2+4`$ for $`d>2`$. The exact source and degree distinctions are collected in [the algebraic references](../references/RELATED_WORK.md#central-polynomials).

For completeness, the passage from those polynomials to order control needs no repetition of a channel variable inside one slot. It uses characteristic-zero linearization.

Suppose $`Q(Y_1,\ldots,Y_r)`$ is multihomogeneous of degrees $`m_1,\ldots,m_r`$, with total degree $`m`$. Replace

```math
Y_j\longmapsto \sum_{k=1}^{m_j}t_{jk}X_{jk}
```

and take the coefficient of $`\prod_{j,k}t_{jk}`$. Call the resulting polynomial $`\widetilde Q`$. Its $`m`$ variables occur once in every word. Because $`Q`$ is scalar on matrices for all scalar parameters, every coefficient is scalar, so $`\widetilde Q`$ is central. Setting all $`X_{jk}=Y_j`$ yields

```math
\widetilde Q(\{Y_j\})
=\left(\prod_j m_j!\right)Q(Y_1,\ldots,Y_r).
```

The multiplier is nonzero over $`\mathbb C`$, so nonidentity is preserved. If an initial homogeneous polynomial is not multihomogeneous, scale its variables independently and select a nonidentity multihomogeneous component first.

The [qutrit construction](qutrit.md) is already multilinear. Its eight slots require no additional linearization. The [Regev construction](probability.md#regev-construction) offers a different explicit family with a closed success formula, at a larger call count.

## 5. Meaning of the result

The environment can retain which Kraus tuple occurred. Exact recovery is possible because **every** successful tuple operator is scalar, so all tuple branches induce the same target channel. The protocol does not erase an environment record by accessing it.

The theorem supplies a finite exact identity branch, with its probability and controller resources charged separately. [The converse](query_threshold.md) identifies when such a branch can first appear. [The resource account](probability.md) explains why a minimum query count alone does not establish an efficient communication rate.

