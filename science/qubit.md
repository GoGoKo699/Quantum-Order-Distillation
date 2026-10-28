# Four calls for a qubit

For arbitrary $`2\times2`$ matrices, define

```math
P_2(X_1,X_2,X_3,X_4)
=\{[X_1,X_2],[X_3,X_4]\}.
```

This central polynomial appears in [Navascués's resetting construction](../references/RELATED_WORK.md#navascues2018). Its eight-order channel realization is the known successful second-order superswitch branch of [Kechrimparis, Moran and Kwon](../references/RELATED_WORK.md#kechrimparis2025). It gives an explicit example of the general [construction and converse](query_threshold.md).

## 1. The eight orders

Words list matrices from left to right. Chronological routing follows the reverse word.

| Word | Coefficient | Word | Coefficient |
|---|---:|---|---:|
| 1234 | +1 | 3412 | +1 |
| 1243 | -1 | 3421 | -1 |
| 2134 | -1 | 4312 | -1 |
| 2143 | +1 | 4321 | +1 |

There are eight coefficients of unit magnitude, so $`L=8`$. Prepare the uniform superposition of these orders and postselect on the uniform superposition with the listed signs. Each successful noise-history Kraus operator is $`P_2/8`$.

The controller occupies eight orthogonal order states, which fit in three qubits. A physical router must implement the corresponding coherent paths.

## 2. Centrality and a nonzero evaluation

Every commutator is traceless. For a traceless $`2\times2`$ matrix $`C`$, Cayley–Hamilton gives

```math
C^2=\frac12\mathrm{Tr}(C^2)I.
```

Apply this identity to $`C+D`$ and subtract its versions for $`C`$ and $`D`$. For traceless $`C,D`$,

```math
CD+DC=\mathrm{Tr}(CD)I,
\qquad
P_2=\mathrm{Tr}\!\left([X_1,X_2][X_3,X_4]\right)I.
```

The polynomial is nonidentity:
$`P_2(\sigma_x,\sigma_y,\sigma_x,\sigma_y)=-8I`$.
The spanning argument therefore makes its probability positive for every full-Choi-rank qubit channel, including nonunital channels.

For a Pauli channel
$`\mathcal N(\rho)=\sum_{\mu=0}^3p_\mu\sigma_\mu\rho\sigma_\mu`$,
the probability is the established superswitch expression

```math
p_{\mathrm{succ}}
=4(p_1^2p_2^2+p_1^2p_3^2+p_2^2p_3^2).
```

For complete depolarization, $`p_\mu=1/4`$, it is $`3/64`$. The [probability account](probability.md) derives this normalization and gives a formula for arbitrary, possibly different and nonunital, qubit slots.

## 3. An exact certificate excluding fewer calls

The converse reduces a full-rank qubit identity branch with $`m`$ calls to a multilinear central polynomial
$`P=\sum_{\pi\in S_m}c_\pi X_{\pi_1}\cdots X_{\pi_m}`$.

A $`2\times2`$ matrix is scalar exactly when its two off-diagonal entries and its diagonal difference vanish. By multilinearity, these conditions need only be imposed on tuples from $`E_{00},E_{01},E_{10},E_{11}`$. They are linear equations in the $`m!`$ coefficients. For degrees one, two and three, a square subsystem already has full rank:

| Degree | Coefficients | Certificate rank | Determinant |
|---|---:|---:|---:|
| 1 | 1 | 1 | 1 |
| 2 | 2 | 2 | -1 |
| 3 | 6 | 6 | -1 |

For degree three, order the coefficient columns as $`123,132,213,231,312,321`$. One certificate matrix is

```math
\begin{pmatrix}
1&1&1&1&1&1\\
1&0&1&0&0&0\\
0&0&0&0&1&1\\
0&1&0&0&1&0\\
1&0&0&1&-1&0\\
0&0&1&1&0&0
\end{pmatrix}.
```

Its rows arise from these scalarity conditions:

| Matrix tuple | Coordinate |
|---|---|
| $`(E_{00},E_{00},E_{00})`$ | $`M_{00}-M_{11}`$ |
| $`(E_{00},E_{00},E_{01})`$ | $`M_{01}`$ |
| $`(E_{00},E_{00},E_{10})`$ | $`M_{10}`$ |
| $`(E_{00},E_{01},E_{00})`$ | $`M_{01}`$ |
| $`(E_{00},E_{01},E_{10})`$ | $`M_{00}-M_{11}`$ |
| $`(E_{00},E_{10},E_{00})`$ | $`M_{10}`$ |

Its nonzero determinant forces every complex coefficient to vanish. The smaller certificates and their generating tuples are included in the [exact output](../verification/qubit_minimum_results.json). Together with the degree-four construction, this proves $`c(2)=4`$ and the corresponding four-call threshold.

## 4. Reproduction

From the repository root:

```bash
python3 verification/qubit_minimum.py
```

The script regenerates the scalarity rows and exact rational determinants, then compares them with the committed certificate. This is a finite algebraic proof certificate. The separate [channel checks](../verification/README.md) verify the physical compiler, reference-system behavior, and probability formulas.

