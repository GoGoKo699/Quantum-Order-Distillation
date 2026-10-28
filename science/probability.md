# Success probability and resources

The [order-only model](model.md) supplies coherent control of independent channel slots while their noise environments remain inaccessible. The [realization theorem](polynomial_realization.md) establishes exact successful identity branches; the [query threshold](query_threshold.md) gives their minimum positive call count. This page calculates success probabilities, normalizes the controller, and charges channel uses and order labels. Exact recovery with positive probability does not establish an efficient implementation or a positive unassisted communication rate.

## Coefficient normalization

The target has dimension $`d`$, and each attempt uses a positive number $`m\ge1`$ of channel slots. Slot $`j`$ contains an independent CPTP channel $`\mathcal N_j`$ with Kraus operators $`K^{(j)}_{a_j}`$. A word $`\pi=(\pi_1,\ldots,\pi_m)`$ denotes **mathematical matrix-product order**,

```math
K_{\boldsymbol a,\pi}
=K^{(\pi_1)}_{a_{\pi_1}}\cdots K^{(\pi_m)}_{a_{\pi_m}}.
```

The physical temporal order is reversed. The switch Kraus operator is

```math
W_{\boldsymbol a}=\sum_{\pi}K_{\boldsymbol a,\pi}\otimes|\pi\rangle\langle\pi|.
\tag{1}
```

This is the arbitrary-order construction in Kechrimparis–Moran–Kwon, published Appendix A, Eqs. (A1)–(A4). Their controller preparation and measurement basis are freely chosen. [Kechrimparis, Moran and Kwon](../references/RELATED_WORK.md#kechrimparis2025)

Let

```math
P(X_1,\ldots,X_m)=\sum_\pi c_\pi X_{\pi_1}\cdots X_{\pi_m},
\qquad L=\sum_\pi|c_\pi|>0.
```

Identical formal words must first be combined. Initialize and postselect the controller in

```math
\alpha_\pi=\sqrt{|c_\pi|/L},\qquad
\beta_\pi=e^{-i\arg c_\pi}\sqrt{|c_\pi|/L}.
\tag{2}
```

Both vectors have unit norm. The successful target Kraus operators are exactly

```math
B_{\boldsymbol a}=\langle\beta|W_{\boldsymbol a}|\alpha\rangle
=\frac{P(K^{(1)}_{a_1},\ldots,K^{(m)}_{a_m})}{L}.
\tag{3}
```

If $`P`$ is central, write its value as $`\lambda_{\boldsymbol a}I_d`$. The successful operation and its probability are

```math
\mathcal S(\rho)=p\rho,
\qquad
p=\frac{1}{L^2}\sum_{\boldsymbol a}|\lambda_{\boldsymbol a}|^2.
\tag{4}
```

The probability is independent of the input and $`0\le p\le1`$, because the construction is a branch of a trace-preserving switch followed by a valid measurement. For every reference system $`R`$,
$`(\mathcal S\otimes\mathrm{id}_R)(\rho_{TR})=p\rho_{TR}`$.

### Optimality for a fixed coefficient ray

Suppose a pure preparation and one rank-one success effect implement coefficients $`t c_\pi`$. The effect vector can have norm at most one. Then

```math
|t|L=\sum_\pi|\alpha_\pi\beta_\pi^*|
\le\|\alpha\|_2\|\beta\|_2\le1.
\tag{5}
```

Equation (2) saturates the bound. Thus $`|t|=1/L`$ is optimal for this **fixed formal coefficient ray**, and Eq. (4) is its maximal single-outcome probability within this preparation/postselection class. This does not optimize over other central polynomials, alternative representatives differing by matrix identities, mixed controllers, or several accepted outcomes. In particular, the explicit qubit and qutrit probabilities below are not asserted to be globally optimal over all allowed controllers.

### Choi convention and channel independence

The [realization proof](polynomial_realization.md#independence-and-kraus-representations) establishes Kraus-gauge invariance for the complete switch. Independent isometric changes of the slot Kraus families induce their tensor-product isometry on the complete noise-history label. The controller therefore acts on channels, with no Kraus representation supplied as an additional resource.

Use the **unnormalized** Choi matrix

```math
J_{\mathcal N}=(\mathcal N\otimes\mathrm{id})(|\Omega_d\rangle\langle\Omega_d|)
=\sum_a|K_a\rangle\!\rangle\langle\!\langle K_a|,
\quad |\Omega_d\rangle=\sum_{i=1}^d|i\rangle|i\rangle.
\tag{7}
```

Thus $`\mathrm{Tr}J_{\mathcal N}=d`$, $`J_\Delta=I_{d^2}/d`$ for the completely depolarizing channel, and
$`J_{\mathcal S}=p|\Omega_d\rangle\langle\Omega_d|`$.
The probability is $`\mathrm{Tr}J_{\mathcal S}/d`$. This factor of $`d`$ is essential.

## Full-rank noise floor and small-noise limits

### CP noise floor

For a full-Choi-rank channel define

```math
\eta_j=d\,\lambda_{\min}(J_{\mathcal N_j})\in(0,1].
\tag{8}
```

Then $`J_{\mathcal N_j}-\eta_jJ_\Delta\ge0`$, so

```math
\mathcal N_j=\eta_j\Delta+(1-\eta_j)\mathcal R_j
\tag{9}
```

with $`\mathcal R_j`$ CPTP when $`\eta_j<1`$. The case $`\eta_j=1`$ is $`\mathcal N_j=\Delta`$. With normalized Choi states $`\rho_{\mathcal N}=J_{\mathcal N}/d`$, the same parameter is $`d^2\lambda_{\min}(\rho_{\mathcal N})`$.

The switch branch is multilinear in independent slots. Expanding Eq. (9) produces a sum of completely positive operations with nonnegative weights. Retaining the all-depolarizing term proves

```math
\boxed{\quad
p(\mathcal N_1,\ldots,\mathcal N_m)
\ge p(\Delta,\ldots,\Delta)\prod_{j=1}^m\eta_j.
\quad}
\tag{10}
```

The bound applies to a fixed valid controller, including a general controller instrument whose successful map is an identity branch. For a nonidentity central-polynomial controller, its depolarizing probability is strictly positive because the depolarizer's matrix-unit Kraus operators span the full matrix algebra. The inequality is a uniform guarantee on classes with a Choi eigenvalue floor, not a claim of an appreciable rate near the Choi boundary.

### Noise approaching the identity

Let $`c(d)`$ be the least positive degree of a nonidentity multilinear central polynomial for $`M_d(\mathbb C)`$, with $`d\ge2`$. Assume the order-only model and the full-rank converse established in the model/theorem account. Consider

```math
\mathcal N_\epsilon=(1-\epsilon)\mathrm{id}+\epsilon\Delta,
\qquad 0<\epsilon\le1.
\tag{11}
```

For $`m=c(d)`$, every exact successful controller satisfies

```math
\boxed{\quad
p_{\mathcal C}(\mathcal N_\epsilon,\ldots,\mathcal N_\epsilon)
=\epsilon^{c(d)}p_{\mathcal C}(\Delta,\ldots,\Delta).
\quad}
\tag{12}
```

**Proof.** Resolve mixed preparation and a general success effect into pure components. Choi positivity and the full-rank converse imply that each component is a central polynomial, or a matrix identity, of degree $`m`$. Substituting $`I`$ for any one variable leaves a multilinear central polynomial of degree $`m-1`$. By minimality it must vanish identically. In the Kraus representation consisting of $`\sqrt{1-\epsilon}I`$ and the depolarizing Kraus operators multiplied by $`\sqrt\epsilon`$, every history containing an identity Kraus operator therefore contributes zero. The all-depolarizing histories have amplitudes multiplied by $`\epsilon^{m/2}`$, proving Eq. (12) componentwise and after summation.

The same law holds for the optimized probability, with optimization interpreted as a supremum:

```math
p_{\mathrm{opt},c(d)}(\mathcal N_\epsilon)
=\epsilon^{c(d)}p_{\mathrm{opt},c(d)}(\Delta),
\qquad \epsilon>0.
\tag{13}
```

Indeed, every $`\mathcal N_\epsilon`$ with $`\epsilon>0`$ and $`\Delta`$ has full Choi rank. The converse identifies the same feasible exact-identity controller class for all these channels. Equation (12) holds for every member of that common class, so it survives taking the supremum; existence of a maximizing controller is unnecessary.

This is specifically an identity-channel admixture statement. Replacing the identity by an arbitrary unitary does not follow from the degree argument, because substitution then introduces a non-scalar matrix constant. At $`\epsilon=0`$ the feasible controller class changes: the identity channel itself permits deterministic success, whereas any fixed central controller appearing in Eq. (12) has zero success. Thus Eq. (13) is not an equality at that endpoint.

For completeness, a related upper bound holds for any fixed number $`m`$ of calls:

```math
p_{\mathrm{opt},m}(\mathcal N_\epsilon)
\le\Pr\{\mathrm{Binomial}(m,\epsilon)\ge c(d)\}.
\tag{14}
```

Here $`p_{\mathrm{opt},m}=0`$ when no positive-probability exact identity controller exists, including $`m<c(d)`$.

After choosing which slots contribute $`\Delta`$, identity substitution leaves a central polynomial in the remaining $`k`$ variables. It vanishes for $`1\le k<c(d)`$. It also vanishes for $`k=0`$: centrality with one arbitrary variable and all others equal to $`I`$ forces $`\sum_\pi c_\pi=0`$. Every conditional branch with $`k\ge c(d)`$ has probability at most one. Summing its mixture weights gives Eq. (14). This is a bound, not an achievability result for larger $`m`$.

### Expected charged calls

The minimum-call probability law does not imply that $`\epsilon^{-c(d)}`$ attempts are the best use of the supplied channels. Allow more calls per attempt, charge every call including failed attempts, and repeat independent attempts on freshly prepared entangled pairs. In the **unaugmented order-only model**, define

```math
\mathcal Q_d(\epsilon)
=\inf_{\substack{m\ge1,\ \mathcal C,\ p_{\mathcal C}>0\\
\mathcal S_{\mathcal C}(\mathcal N_\epsilon,\ldots,\mathcal N_\epsilon)
=p_{\mathcal C}\,\mathrm{id}_d}}
\frac{m}{p_{\mathcal C}}.
\tag{14a}
```

Here each attempt has a fixed batch size $`m`$; the infimum may choose that size and controller as functions of the known parameter $`\epsilon`$. Additional depolarizing gadgets or other target processing are outside this model and outside the following lower bound.

Let $`c=c(d)`$, and choose any degree-$`c`$ central controller with depolarizing success $`p_D>0`$. Then, for $`0<\epsilon\le1`$,

```math
\boxed{\quad
\frac{c}{\epsilon}\le\mathcal Q_d(\epsilon)
\le\frac{c(1+\epsilon)}{p_D(1-e^{-1})^c\,\epsilon}
\le\frac{2c}{p_D(1-e^{-1})^c\,\epsilon}.
\quad}
\tag{14b}
```

In particular, $`\mathcal Q_d(\epsilon)=\Theta_d(\epsilon^{-1})`$ for each fixed $`d`$. The constants are dimension dependent and can be large; no optimal prefactor or uniform-in-d efficiency is asserted.

**Lower bound.** Markov's inequality applied to Eq. (14) gives
$`p\le\Pr\{\mathrm{Binomial}(m,\epsilon)\ge c\}\le m\epsilon/c`$.
Every positive-probability fixed-batch protocol therefore has $`m/p\ge c/\epsilon`$.

**Construction.** Replace each of the base controller's $`c`$ slots by a contiguous block of $`k=\lceil1/\epsilon\rceil`$ supplied channel uses, always in the same internal order. Superpose only the original block permutations. The target still visits every one of the $`ck`$ supplied boxes once, and there is no intervening target operation. Algebraically the new polynomial is

```math
P(X_1\cdots X_k,\ X_{k+1}\cdots X_{2k},\ldots,
  X_{(c-1)k+1}\cdots X_{ck}),
\tag{14c}
```

a degree-$`ck`$ multilinear central polynomial implemented by a subset of the $`(ck)`$-switch orders. Its active coefficient support and controller vectors are the same as those of the base polynomial; the longer paths still have to be physically implemented and charged.

Because $`\Delta^2=\Delta`$, each independent composite block is
$`\mathcal N_\epsilon^k=\mathcal N_\delta`$, where
$`\delta=1-(1-\epsilon)^k\ge1-e^{-1}`$.
Equation (12) applied to the $`c`$ block channels gives success $`\delta^c p_D`$. Using $`k\le(1+\epsilon)/\epsilon`$ proves the upper bound in Eq. (14b). Choosing the batch size assumes $`\epsilon`$ is known, or bounded to within a fixed multiplicative factor. It is not a channel-estimation protocol.

The constructive upper bound does not require solving the minimal-degree problem in an arbitrary dimension. Take any explicit nonidentity central polynomial of known degree $`r`$ and depolarizing success $`p_D>0`$, such as the degree-$`2d^2`$ Regev construction below. The same blocks and the CP-floor bound give success at least $`\delta^r p_D`$, hence cost at most
$`2r/[p_D(1-e^{-1})^r\epsilon]`$.
Together with the unchanged lower bound $`c(d)/\epsilon`$, this proves the same $`\Theta_d(\epsilon^{-1})`$ scaling by a fully explicit construction even when $`c(d)`$ is unknown.

For qubits, [Kechrimparis, Moran and Kwon, Secs. V–VI and Appendices F, J, K](../references/RELATED_WORK.md#kechrimparis2025) already studies improved success under higher-order superswitches and a positive large-depth success limit at each fixed nonzero depolarizing noise. That qubit benefit from additional uses is inherited. Equation (14b) supplies a matching charged-query scaling and an architecture-wide lower bound in every fixed dimension; it is not a claim that increasing the qubit superswitch order was previously unknown.

## Qubit success probability

The [four-call qubit construction](qubit.md) uses

```math
P_2=\{[X_1,X_2],[X_3,X_4]\},\qquad L=8.
\tag{15}
```

Its eight order words and centrality proof are given on that page. It is the known qubit superswitch branch of [Kechrimparis, Moran and Kwon](../references/RELATED_WORK.md#kechrimparis2025). The calculation here extends its probability formula to general, potentially nonunital qubit channels and to four different input channels.

Expand each Kraus operator in the usual Pauli basis:

```math
K^{(j)}_a=k^{(j)}_{a0}I+\boldsymbol k^{(j)}_a\cdot\boldsymbol\sigma,
\qquad k^{(j)}_{ai}=\tfrac12\mathrm{Tr}(\sigma_iK^{(j)}_a).
```

The vector components are generally complex. Define the positive semidefinite covariance

```math
Q_j=\sum_a\boldsymbol k^{(j)}_a\boldsymbol k^{(j)\dagger}_a.
\tag{16}
```

The polynomial value and normalized Kraus scalar are

```math
P_2(K_1,K_2,K_3,K_4)
=-8[(\boldsymbol k_1\times\boldsymbol k_2)^T
      (\boldsymbol k_3\times\boldsymbol k_4)]I,
\qquad
B=-[(\boldsymbol k_1\times\boldsymbol k_2)^T
      (\boldsymbol k_3\times\boldsymbol k_4)]I.
\tag{17}
```

Define

```math
S(Q_1,Q_2)_{ij}
=\sum_{k,l,m,n}\varepsilon_{ikl}\varepsilon_{jmn}
 (Q_1)_{km}(Q_2)_{ln}.
\tag{18}
```

It is the covariance of the cross products, so it is positive semidefinite. Summing Eq. (17) gives the exact formula for four possibly different channels:

```math
p=\mathrm{Tr}\big[S(Q_1,Q_2)^T S(Q_3,Q_4)\big].
\tag{19}
```

For identical channels, $`S(Q,Q)=2\mathrm{adj}(Q)^T`$, yielding

```math
\boxed{\quad p=4\mathrm{Tr}\big[\mathrm{adj}(Q)
                                  \mathrm{adj}(Q)^T\big].\quad}
\tag{20}
```

The transpose in these formulas is **not** an adjoint. It arises from the bilinear Pauli-vector dot product. Replacing it with an adjoint is incorrect for general nonunital channels. Positivity follows because $`\mathrm{adj}(Q)`$ and its transpose are positive semidefinite. Full Choi rank implies $`Q>0`$, making Eq. (20) strictly positive.

For a Pauli channel $`\mathcal N(\rho)=\sum_{i=0}^3p_i\sigma_i\rho\sigma_i`$, $`Q=\mathrm{diag}(p_1,p_2,p_3)`$, and

```math
p=4(p_1^2p_2^2+p_1^2p_3^2+p_2^2p_3^2).
\tag{21}
```

This recovers [Kechrimparis, Moran and Kwon, Eq. (5)](../references/RELATED_WORK.md#kechrimparis2025). In particular,

```math
p(\Delta_2)=\frac3{64},\qquad
p((1-\epsilon)\mathrm{id}+\epsilon\Delta_2)=\frac{3\epsilon^4}{64}.
\tag{22}
```

Full Choi rank is sufficient, not necessary: some rank-deficient Pauli channels have positive Eq. (21). Conversely, a unitary channel and ordinary amplitude damping make this particular central polynomial vanish. The numerical suite includes those boundary examples to test the complex-transpose issue and avoid overgeneralizing the full-rank theorem.

## Exact depolarizing contractions

An independent way to calculate success probabilities avoids enumerating all Kraus histories. For $`m`$ identical completely depolarizing slots, use Kraus operators
$`K_{ab}=|a\rangle\langle b|/\sqrt d`$.
For a word $`\pi`$, define the augmented cycle on $`\{0,1,\ldots,m\}`$,

```math
C_\pi=(0,\pi_1,\ldots,\pi_m).
```

Let $`\#\mathrm{cyc}(g)`$ count all cycles of a permutation, including fixed points. The normalized Hilbert–Schmidt Gram entry is

```math
g_{\pi\tau}:=\frac1d\sum_{\boldsymbol a}
 \mathrm{Tr}(K_{\boldsymbol a,\pi}^\dagger K_{\boldsymbol a,\tau})
=d^{\#\mathrm{cyc}(C_\pi C_\tau^{-1})-m-1}.
\tag{23}
```

**Derivation.** The $`\pi`$-product of matrix units supplies $`m-1`$ internal Kronecker constraints. After imposing them, label the $`m+1`$ remaining wire indices by $`x_0,\ldots,x_m`$, where $`x_j`$ for $`j\ge1`$ is the outgoing index of slot $`j`$, and $`x_0`$ is the first incoming index. The $`\tau`$-product and the trace identify
$`x_{C_\pi^{-1}(j)}=x_{C_\tau^{-1}(j)}`$ for every $`j=0,\ldots,m`$.
There is one freely summed index for each cycle of $`C_\tau^{-1}C_\pi`$, which has the same cycle count as $`C_\pi C_\tau^{-1}`$. Each contributes a factor $`d`$; the Kraus normalization supplies $`d^{-m}`$, and the normalized trace supplies $`d^{-1}`$. This proves Eq. (23).

For any coefficient vector,

```math
p_{I/d}=\frac1{L^2}\sum_{\pi,\tau}\overline{c_\pi}c_\tau
d^{\#\mathrm{cyc}(C_\pi C_\tau^{-1})-m-1}.
\tag{24}
```

For a central polynomial this is its input-independent identity success probability. For a noncentral polynomial it is only the success probability on the maximally mixed input.

If coefficients are real integers, collect the exact signed histogram

```math
h_k=\sum_{\pi,\tau:\ \#\mathrm{cyc}(C_\pi C_\tau^{-1})=k}c_\pi c_\tau.
\tag{25}
```

Then $`p=L^{-2}\sum_k h_kd^{k-m-1}`$. Checks include
$`h_{m+1}=\sum_\pi c_\pi^2`$ and $`\sum_kh_k=(\sum_\pi c_\pi)^2`$.
The qubit eight-word polynomial gives $`h_1=16,h_3=-24,h_5=8`$, recovering $`3/64`$ at $`d=2`$.

### Eight-call qutrit probability

The [explicit qutrit construction](qutrit.md) uses the Drensky–Kasparian degree-eight central polynomial from [Drensky and Kasparian, pp. 746–747](../references/RELATED_WORK.md#central-polynomials). Its [coefficient list](../verification/degree8_coefficients.json) contains 2880 distinct words with coefficients $`\pm1`$, hence $`L=2880`$. The construction and centrality proof are given on that page. Applying Eq. (25) to the coefficient list gives

| Cycle count $`k`$ | $`h_k`$ |
|---:|---:|
| 1 | 124416 |
| 3 | -259200 |
| 5 | 172224 |
| 7 | -40320 |
| 9 | 2880 |

Thus

```math
\boxed{p(\Delta_3)=\frac1{43740}
\simeq2.286236854\times10^{-5}.}
\tag{26}
```

The same exact contraction expression vanishes at dimensions one and two, as required when the qutrit central polynomial is restricted to smaller matrix algebras.

There is also a structurally independent check. On $`3\times3`$ matrices the scalar value is, by the printed trace formula in [Drensky and Kasparian, p. 752](../references/RELATED_WORK.md#central-polynomials),

```math
\lambda=2\sum_{\sigma\in S_4}\mathrm{sgn}(\sigma)
\mathrm{Tr}(X_{\sigma_1}X_{\sigma_2}[X_5,X_6])
\mathrm{Tr}(X_{\sigma_3}X_{\sigma_4}[X_7,X_8]).
\tag{27}
```

Expanding this yields 96 trace terms. Their pairwise matrix-unit contractions use ordinary permutations on eight symbols, without the augmented boundary index. The exact weights at cycle counts $`4,6,8`$ are $`1536,-1920,384`$. Therefore

```math
\sum_{\boldsymbol a}|\lambda_{\boldsymbol a}|^2
=384\left(1-\frac5{3^2}+\frac4{3^4}\right)
=\frac{5120}{27}.
```

Dividing by $`2880^2`$ independently gives Eq. (26). At the minimum degree $`c(3)=8`$, Eq. (12) also gives
$`p(\mathcal N_\epsilon)=\epsilon^8/43740`$ for this specific controller.

## Regev construction

The Regev polynomial supplies a simple all-dimensional quantitative formula, though its $`2d^2`$ calls are generally not query minimal. Let $`r=d^2`$, let $`X`$ and $`Y`$ each contain $`r`$ distinct variables, and define

```math
m_k(X)=X_{(k-1)^2+1}\cdots X_{k^2},
\qquad
R_d=\mathrm{Alt}_X\mathrm{Alt}_Y
  [m_1(X)m_1(Y)\cdots m_d(X)m_d(Y)].
\tag{28}
```

Alternation here is the unnormalized signed sum. There are exactly
$`s=L=(r!)^2`$ distinct signed words. [Procesi's published 2021 account](../references/RELATED_WORK.md#central-polynomials), Theorem 7 and Eq. (2.85) on page 99, gives the determinant evaluation, with the constant from Eq. (1.65) on page 94. Its proof on page 100 follows Formanek. The same formula appears in the revised 2022 account, Theorem 4.13 and Eq. (37):

```math
R_d(X;Y)=a_d\det[\mathrm{vec}X_1,\ldots,\mathrm{vec}X_r]
                \det[\mathrm{vec}Y_1,\ldots,\mathrm{vec}Y_r]I_d,
```

```math
a_d=\frac{(-1)^{d-1}C_d^2}{(d!)^2(2d-1)},
\qquad
C_d=\frac{\prod_{k=1}^d(2k-1)!}{\prod_{k=1}^{d-1}k!}.
\tag{29}
```

The sign convention of $`\mathrm{vec}`$ does not affect probabilities. The same coordinate basis must be used in both determinants.

For identical slots $`\mathcal N`$, form the matrix $`V`$ with columns $`\mathrm{vec}K_a`$, so $`VV^\dagger=J_{\mathcal N}`$. Cauchy–Binet, including all orderings of each set of distinct columns, gives

```math
\sum_{a_1,\ldots,a_r}
 |\det[\mathrm{vec}K_{a_1},\ldots,\mathrm{vec}K_{a_r}]|^2
=r!\det J_{\mathcal N}.
\tag{30}
```

Equations (4), (29), and (30) yield

```math
\boxed{p_R(\mathcal N)
=\frac{|a_d|^2\det(J_{\mathcal N})^2}{(r!)^2}.}
\tag{31}
```

This particular protocol succeeds precisely when $`J_{\mathcal N}`$ has full rank. That equivalence does not characterize distillability by other polynomials or controllers.

Since $`J_{\mathcal N}\ge0`$ and $`\mathrm{Tr}J_{\mathcal N}=d`$, the arithmetic–geometric mean inequality gives $`\det J_{\mathcal N}\le d^{-d^2}`$, with equality at the depolarizer. Consequently

```math
p_R(\mathcal N)\le p_R(\Delta_d)
=\frac{|a_d|^2}{(d^2!)^2d^{2d^2}}.
\tag{32}
```

For example, $`a_2=-3`$, $`a_3=720`$, and
$`p_R(\Delta_2)=1/16384`$,
$`p_R(\Delta_3)=1/98411002933824`$.

Stirling summation gives, with natural logarithms,

```math
\log p_R(\Delta_d)
=-4d^2\log d+(4\log2-1)d^2+O(d\log d).
\tag{33}
```

One way to check the expansion is the identity
$`C_d=2^{d^2}\pi^{-d/2}\prod_{k=1}^d\Gamma(k+1/2)`$.
Stirling's expansion and Euler–Maclaurin summation give

```math
\log C_d=\tfrac12d^2\log d+(\log2-\tfrac34)d^2
+\tfrac12d\log d+(\tfrac12\log2-\tfrac12)d
+\tfrac1{24}\log d+O(1).
```

Substituting this into
$`\log p_R(\Delta_d)=4\log C_d-4\log(d!)-2\log(2d-1)-2\log((d^2)!)-2d^2\log d`$
gives the sharper expression

```math
\log p_R(\Delta_d)
=-4d^2\log d+(4\log2-1)d^2
-2d\log d+(2\log2+2)d-\frac{35}{6}\log d+O(1).
\tag{34}
```

The exact formula, rather than the asymptotic approximation, is used for all quoted finite-dimensional numbers. Equation (33) exhibits severe probability suppression; this construction is an existence proof with explicit constants, not a practical-rate proposal.

## Resource accounting

The support size $`s`$ counts coherently addressed order labels in these explicit constructions. A register can encode them in $`\lceil\log_2s\rceil`$ qubits, but this says nothing by itself about the gate count, spatial routing, depth, or fault tolerance needed to implement the corresponding generalized switch.

| Protocol | Target dimension | Calls per attempt | Active order labels | Controller qubits to encode labels | Success on $`\Delta`$ |
|---|---:|---:|---:|---:|---:|
| Qubit commutator/anticommutator | 2 | 4 | 8 | 3 | $`3/64`$ |
| Drensky–Kasparian | 3 | 8 | 2880 | 12 | $`1/43740`$ |
| Regev | 2 | 8 | 576 | 10 | $`1/16384`$ |
| Regev | 3 | 18 | 131681894400 | 37 | $`1/98411002933824`$ |

For the Regev family,
$`\log_2s=4d^2\log_2d-2(\log_2e)d^2+O(\log d)`$.
Each attempt uses every independent channel slot once in every order branch. The noise environments are neither controlled nor measured. These facts do not make the coherent routing or independent controller resources free.

For independent attempts on freshly prepared Bell pairs, success probability $`p`$ gives mean $`1/p`$ attempts and $`m/p`$ channel uses per heralded pair, conditional on the supplied ideal switch resource. Thus the explicit depolarizing qutrit example has mean 43740 attempts and 349920 channel calls per successful pair before charging routing and control costs. This accounting is distinct from an ordinary one-way channel capacity. Failure generally disturbs the target, so an unknown quantum state cannot simply be reused on the next attempt. The result supports heralded entanglement generation with the appropriate supplied routing and classical communication, and does not establish positive unassisted quantum capacity.

When the noise has the identity-mixture form of Eq. (11), Eq. (14b) describes how optimizing the fixed batch size changes this charged-call cost. It applies to the unmodified supplied channels. Allowing an additional depolarizing wrapper would change the access model and remove the small-noise lower bound; that augmented model must be accounted for separately.

The number of calls is query minimal for the qubit and qutrit examples because the inherited central-degree values are $`c(2)=4`$ and $`c(3)=8`$; their support sizes and success probabilities have not been optimized globally. [Drensky and Kostadinov](../references/RELATED_WORK.md#central-polynomials)

## Verification

Run the [documented suite](../verification/README.md) from the repository root:

```bash
python3 -m pip install -r verification/requirements.txt
python3 verification/run_checks.py
```

The default suite needs Python 3.11 or later and NumPy. It includes the exact qubit lower-degree certificates, both exact qutrit coefficient generators, the source's nonzero witness, and the independent 96-term qutrit trace calculation. Actual full switch Kraus matrices test trace preservation, generic nonunital probabilities, reference preservation, controller phases, and rectangular Kraus changes in every slot. An independent access-model calculation tests mixed ancillary controls, general success effects, the classical-order bound, and a correlated-noise counterexample. Regev matrix-unit evaluations and random tuples check the determinant constant and Cauchy–Binet normalization.

The optional `--exhaustive` flag also compiles the supplied C++17 helper and recomputes every pair of the 2880 qutrit words. The all-word histogram and trace contraction both give exactly $`1/43740`$. Integer histograms and rational probabilities are compared exactly; floating results are compared with explicit tolerances and must satisfy the residual bounds in the scripts. The [recorded report](../verification/check_results.json) states which checks ran.

These calculations test transcription, normalization, and numerical implementation. They do not replace the construction and converse proofs or the cited algebraic theorems. The [verification README](../verification/README.md) identifies each script and its canonical output.

## References

| Source | Locations used here |
|---|---|
| [Kechrimparis, Moran and Kwon (2025)](../references/RELATED_WORK.md#kechrimparis2025) | Published Eq. (5) gives the Pauli probability; Appendix A, Eqs. (A1)–(A6), gives the generalized switch and controller choices; Secs. V–VI and Appendices F, J, K give the qubit higher-order analysis. |
| [Drensky and Kasparian (1985)](../references/RELATED_WORK.md#central-polynomials) | Pages 746–747 define the qutrit polynomial; printed page 752 gives the scalar trace formula used in the independent probability calculation. |
| [Procesi (2021)](../references/RELATED_WORK.md#central-polynomials) | *A note on the Formanek Weingarten function*: Eq. (1.65), printed page 94, gives $`C_d`$; Eq. (2.84) defines the polynomial; Theorem 7 and Eq. (2.85), page 99, give the exact determinant evaluation, proved on page 100. This is the direct source for the constant used here. |
| [Procesi (2022), revised arXiv version 2](../references/RELATED_WORK.md#central-polynomials) | A parallel formulation appears in Eq. (26), Eqs. (35)–(36), and Theorem 4.13, Eq. (37). The 2021 theorem independently supplies the formula required here. |
| [Drensky and Kostadinov (2026)](../references/RELATED_WORK.md#central-polynomials) | The introduction records the inherited qubit and qutrit minimum degrees. The [threshold account](query_threshold.md) separates those results from the general channel converse. |

The [related-work account](../references/RELATED_WORK.md) and [bibliography](../references/references.bib) give full source records and distinguish the inherited central-polynomial constructions from their order-only channel interpretation.
