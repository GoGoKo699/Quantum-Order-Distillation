# The minimum number of channel uses

The construction in [polynomial realization](polynomial_realization.md) turns a central polynomial into a successful switch branch. Full Choi rank also forces the converse: a successful branch must come from a central polynomial.

This section concerns the precise [order-only architecture](model.md). The target visits every supplied slot once in each order, its control starts independently of the target, and success is determined by a final control measurement. There are no interleaved target operations or target-carrying bypass registers.

## 1. The algebraic quantity

For $`d\ge2`$, define

```math
c(d)=\min\left\{
m\ge1:
\begin{array}{l}
\text{a nonidentity multilinear central polynomial}\\
\text{for }M_d(\mathbb C)\text{ has degree }m
\end{array}
\right\}.
```

The minimum is finite by the central-polynomial existence results. Characteristic-zero linearization preserves the minimum positive degree, so this agrees with the minimum degree in the usual multihomogeneous formulation.

For one fixed channel $`\mathcal N`$, let $`m_\star(\mathcal N)`$ be the minimum **positive** number of independent copies that can herald $`p\,\mathrm{id}_d`$ with $`p>0`$ in this architecture. The positive-use condition specifies the conversion task: a zero-slot identity wire is not a use of the noisy channel.

### Theorem 2. Exact query characterization

For every full-Choi-rank channel $`\mathcal N:M_d(\mathbb C)\to M_d(\mathbb C)`$,

```math
\boxed{m_\star(\mathcal N)=c(d).}
```

The controller in the lower bound may depend on the entire description of $`\mathcal N`$. The upper bound has one controller determined only by $`d`$.

More precisely, for a fixed positive number $`m`$, a successful branch exists if and only if a nonidentity multilinear central polynomial of degree $`m`$ exists. This does not assert feasibility for every integer $`m\ge c(d)`$: extra compulsory slots cannot simply be omitted.

## 2. The rank-one Choi argument

First consider a pure initial control $`\lvert\alpha\rangle`$ and a rank-one successful effect $`\lvert\beta\rangle\langle\beta\rvert`$. Their target Kraus operators have the form

```math
B_{\boldsymbol a}
=P(K_{a_1},\ldots,K_{a_m}),
\qquad
P(X_1,\ldots,X_m)
=\sum_\pi \alpha_\pi\overline{\beta_\pi}
 X_{\pi_1}\cdots X_{\pi_m}.
```

Suppose the successful map is $`p\,\mathrm{id}_d`$, with $`p>0`$. Its unnormalized Choi operator is

```math
\sum_{\boldsymbol a}
 \lvert B_{\boldsymbol a}\rangle\!\rangle
 \langle\!\langle B_{\boldsymbol a}\rvert
=p\,\lvert I\rangle\!\rangle\langle\!\langle I\rvert.
\tag{3}
```

For every vector $`v`$ orthogonal to $`\lvert I\rangle\!\rangle`$, taking the quadratic form of (3) gives a sum of nonnegative terms,

```math
\sum_{\boldsymbol a}
|\langle v\mid B_{\boldsymbol a}\rangle\!\rangle|^2=0.
```

Every term is zero. Consequently each $`B_{\boldsymbol a}`$ is a scalar multiple of $`I`$.

Full Choi rank means the Kraus operators span all $`d\times d`$ matrices. For arbitrary $`X_1,\ldots,X_m`$, expand each $`X_j`$ in that spanning family. Multilinearity expresses $`P(X_1,\ldots,X_m)`$ as a linear combination of the scalar values $`B_{\boldsymbol a}`$. It is therefore scalar. Since $`p>0`$, at least one $`B_{\boldsymbol a}`$ is nonzero. The polynomial is central and nonidentity, so $`m\ge c(d)`$.

At $`m=c(d)`$, Theorem 1 supplies the opposite inequality. This proves the characterization for pure controls and rank-one effects.

## 3. General control states and successful events

Let the initial control density operator and successful effect have spectral decompositions

```math
\omega=\sum_r w_r\lvert\alpha_r\rangle\langle\alpha_r\rvert,
\qquad
E=\sum_s e_s\lvert\beta_s\rangle\langle\beta_s\rvert,
\qquad 0\le E\le I.
```

The successful target map has a Kraus decomposition indexed by $`(\boldsymbol a,r,s)`$,

```math
B_{\boldsymbol a,r,s}
=\sqrt{w_r e_s}
\sum_\pi
(\alpha_r)_\pi\overline{(\beta_s)_\pi}
K_{a_{\pi_1}}\cdots K_{a_{\pi_m}}.
```

If the total Choi operator has rank one, the same nonnegative-sum argument forces every contributing Kraus vector into its one-dimensional support. Positive success ensures that at least one pair $`(r,s)`$ supplies a nonidentity central polynomial.

An additional control ancilla initialized independently of the target does not alter the proof. For pure components on control and ancilla, the coefficient of order $`\pi`$ is the inner product of the corresponding ancilla amplitude vectors. Each induced target Kraus operator remains a scalar-coefficient multilinear order polynomial. Arbitrary processing and discarding of the controller **after** the switch are absorbed into its final instrument.

This reduction concerns the controller of the supplied generalized switch. It does not replace a circuit with arbitrary target-control interactions at intermediate times by the same model.

## 4. Final unitary corrections do not lower the count

One can enlarge the successful event to include a known target unitary correction depending on the herald. The threshold for full-rank channels remains $`c(d)`$.

Let $`h`$ label a successful herald and $`U_h`$ its final correction. Decompose its uncorrected branch into Kraus polynomials $`B_{\boldsymbol a,r,s,h}`$. The sum of the corrected branches has Choi matrix $`p\lvert I\rangle\!\rangle\langle\!\langle I\rvert`$. Positivity therefore forces every $`U_hB_{\boldsymbol a,r,s,h}`$ to be scalar. For each fixed $`r,s,h`$, spanning and multilinearity put the polynomial's values in $`\mathrm{span}(U_h^\dagger)`$ on arbitrary matrix tuples. Positive total success guarantees a nonzero such polynomial.

To see why this line must be scalar, consider a nonzero polynomial whose values all belong to the line spanned by a unitary $`U`$. Simultaneous conjugation of its arguments obeys

```math
P(VX_1V^\dagger,\ldots,VX_mV^\dagger)
=V P(X_1,\ldots,X_m)V^\dagger.
```

Choose an evaluation with nonzero value. The line spanned by $`U`$ must then be invariant under conjugation by every unitary $`V`$. Differentiating at $`V=e^{itH}`$ gives $`[H,U]=a_HU`$ for every Hermitian $`H`$. However,

```math
\mathrm{Tr}\!\left(U^\dagger[H,U]\right)
=\mathrm{Tr}(U^\dagger H U)-\mathrm{Tr}H=0.
```

Thus $`a_H=0`$. The unitary $`U`$ commutes with every Hermitian matrix and is scalar. Applying this argument separately to every nonzero corrected outcome recovers the central-polynomial lower bound.

The statement is specific to full-rank slots and order polynomials. Rank-deficient channels can have useful nonscalar unitary branches.

## 5. Known degrees and their quantum meaning

The algebraic results and their provenance are recorded in [Related work](../references/RELATED_WORK.md#central-polynomials).

| Target dimension | Minimum positive uses in the stated architecture | What is established |
|---|---:|---|
| $`d=2`$ | $`4`$ | Exact; an eight-order commutator construction is explicit |
| $`d=3`$ | $`8`$ | Exact; an explicit degree-eight polynomial is supplied |
| $`d=4`$ | At most $`13`$ | A construction, with no assertion that 13 is minimal |
| $`d>2`$ | At most $`(d-1)^2+4`$ | General inherited central-polynomial construction |

The qubit and qutrit minima are inherited algebraic degrees transported through Theorem 2. The [qubit certificate](qubit.md) also gives a direct finite linear-algebra verification that degrees one, two and three cannot work. The [qutrit account](qutrit.md) states the imported lower-bound theorem separately from its explicit construction.

The degree-four-matrix question remains open in the general formulation used here. In particular, a bound restricted to polynomials in **two variables** cannot establish a lower bound for all multilinear order polynomials.

## 6. What is optimized

The theorem optimizes the number of independent noisy slots over every allowed control preparation and successful control effect. It does not minimize the controller dimension, the number of coherent paths, the probability-normalized number of uses, or gate complexity.

For example, the explicit qutrit construction uses eight slots and 2880 occupied order states. Its [success probability](probability.md) is positive but much smaller than one. These are separate resources, and both remain visible in the [interpretation](interpretation.md).


## 7. A dimension-dependent lower bound

There is a useful general bound independent of the exact unknown degrees:

```math
2d-1\le c(d)\le(d-1)^2+4,\qquad d>2.
```

For the lower bound, first consider a nonzero multilinear polynomial of degree $`r\le2d-1`$. Choose one word with a nonzero coefficient. Assign its consecutive variables the first $`r`$ matrix units in

```math
E_{11},E_{12},E_{22},E_{23},\ldots,
E_{d-1,d-1},E_{d-1,d},E_{dd}.
```

A nonzero product must follow the directed path from the first vertex through the available edges. Every loop must be used before the edge leaving its vertex. Exactly the chosen word has a nonzero product, so the polynomial cannot be an identity. This is the elementary matrix-unit lower-bound argument associated with the minimal-degree part of the Amitsur–Levitzki theorem.

If $`P`$ is a nonidentity multilinear central polynomial of degree $`m`$, the free polynomial $`[P,X_{m+1}]`$ is a matrix identity. It is a nonzero free polynomial: in its two groups of words, the new variable occurs at opposite ends, so the supports cannot cancel. Its degree $`m+1`$ must be at least $`2d`$. Thus $`m\ge2d-1`$.

For a target register of $`k`$ qubits, $`d=2^k`$, this gives

```math
m_\star(\mathcal N)\ge2^{k+1}-1.
```

The general finite-dimensional capability is therefore not an efficient-in-qubit-number guarantee in this architecture. The query model charges each call to the entire $`d`$-dimensional channel. Quantum error-correcting encodings and protected target memories are outside that model.

