# A cost exclusion for every six-label computational encoding

At two-call visibility $`x=0.9`$, independent local phase damping and proportional routing windows, assigning the six three-slot orders to different three-bit words cannot recover the extra-call saving of the global logical-noise model. An exact finite certificate covers every injective six-label computational assignment and every allowed controller strategy.

This is the fixed rank-two channel $`K_0=|+\rangle\langle0|`$, $`K_1=|-\rangle\langle1|`$, with $`\mathcal M^2=\Delta`$. It is separate from the [full-rank four-call theorem](query_threshold.md), the [global logical-dephasing benchmark](native_three_call_frontier.md), and the [comparison-encoding local frontier](local_controller_encoding.md).

<a id="all-encoding-cost-theorem"></a>
## 1. Exact benchmark and scope

Use the [same access and accepted target-only output contract](local_controller_encoding.md#1-access-encoding-and-noise): independent compulsory slots; one traversal of every slot in each occupied order; initially target-independent mixed/ancillary controls; inaccessible environments; no target bypass, intermediate gates or coding across failed targets; arbitrary final controller instruments and known outcome-dependent target unitaries. Only the corrected accepted aggregate must have input-independent acceptance $`\mathcal E=p\Lambda`$. Fine and architecture-component probabilities may depend on input.

A three-call encoding injectively maps $`012,021,102,120,201,210`$ into six of the eight physical three-bit computational words. The other two words are unoccupied. Its noise kernel is $`L_{ij}=q^{d_H(c_i,c_j)}`$. Two orders use one physical control bit. The declared exposure and charged routing windows are $`T_m=m\tau`$, with

```math
x=e^{-2\kappa\tau}=\frac9{10},\qquad
q=e^{-3\kappa\tau}=x^{3/2}=\sqrt{\frac{729}{1000}}.
```

Write $`\epsilon=\frac12\lVert\Lambda-\mathrm{id}\rVert_\diamond`$ and $`C=\mathbb E[m]/p`$. These are fresh-target trial costs. Gate-level routing, preparation/readout overhead, loss and gate faults are uncharged and unverified.

**Theorem.** Under this contract, the exact optimum at error budget $`3/8`$ is

```math
C_{123,*}(3/8)=C_{12,*}(3/8)=101/32=3.15625.
```

More strongly, if $`t_3`$ is the raw probability of choosing a three-call trial, every feasible policy satisfies

```math
\mathbb E[m]\ge\frac{101}{32}p+\frac3{10}t_3.
```

Every positive three-call share is therefore strictly suboptimal at this benchmark. Charged routing time has optimum $`(101/32)\tau`$. The theorem includes classical randomization over encodings and architectures, arbitrary nonuniform preparations, complex effects and controller ancillas, and continuous target-unitary corrections. It does not cover redundant eight-word encodings, noncomputational/protected codes or a different controller channel.

## 2. A support functional for arbitrary controllers

The [coefficient characterization](controller_optimization.md#controller-feasibility) gives accepted matrices

```math
G_y\succeq0,\qquad \sum_yG_y\preceq D,
\qquad D\succeq0\text{ diagonal},\qquad\mathrm{Tr}D=1.
```

Noise acts by $`G_y\circ L`$. This reduction covers mixed/ancillary preparations and the rejected part of the instrument. It does not set $`D=I/6`$.

For the normalized-input Choi map, let $`H`$ be the trace kernel and $`M_Q`$ the Bell-overlap kernel for a unitary Bell vector $`Q`$. Inverting the final correction merely renames this unrestricted unitary. With

```math
B=69/4,\qquad c=1-B/2=-61/8,
\qquad W_Q=L\circ(BM_Q+cH),
```

the corrected support of a three-call strategy is

```math
S_3=cP_3+BF_3=\sum_y\mathrm{Tr}(G_yW_{Q_y}).
```

Here $`P_3`$ is the accepted Choi trace and $`F_3`$ its unnormalized corrected Bell overlap, without assuming that component acceptance is input independent.

**Certificate lemma.** For every permitted encoding, there is a positive semidefinite $`K`$ satisfying $`K\succeq W_Q`$ for every unitary $`Q`$ and $`K_{ii}\le27/10`$. Hence

```math
S_3\le\mathrm{Tr}\!\left(K\sum_yG_y\right)
\le\mathrm{Tr}(KD)\le27/10<3.
```

Positivity of $`K`$ is essential for the middle inequality, because the accepted effects need not consume the full diagonal budget. Sections 3–5 supply the complete finite certificate, independent of solver status.

## 3. Exact kernels and continuous target corrections

The [native kernel expansion](native_three_call_frontier.md#exact-bell-kernels) gives, in the $`I,U,Y,V`$ Bell basis,

```math
H=3T+\tfrac34(O+E),\qquad M_I=\tfrac32T,
\qquad M_U=\tfrac34O,
```

```math
M_Y=\tfrac12(T+E),\qquad M_V=T+\tfrac14E.
```

All target-Bell cross kernels vanish except $`U,V`$. For $`Q_\theta=\cos\theta\,U+\sin\theta\,V`$,

```math
M_\theta=\cos^2\theta M_U+\sin^2\theta M_V
+\cos\theta\sin\theta C_{\mathrm{rot}},
```

where the real symmetric cross kernel in the stated six-word order is

```math
C_{\mathrm{rot}}=\frac14
\begin{pmatrix}
0&1&-1&0&0&0\\
1&0&0&0&-1&0\\
-1&0&0&1&0&0\\
0&0&1&0&0&-1\\
0&-1&0&0&0&1\\
0&0&0&-1&1&0
\end{pmatrix}.
```

Every target unitary, up to phase, has form $`aI+i(bU+dY+eV)`$ with real coefficients whose squares sum to one. The $`U,V`$ coefficients share a phase, so only their real cross kernel enters the overlap. Its witness is a convex combination of the $`I`$, $`Y`$ and real-angle witnesses above. This statement holds for arbitrary complex controller matrices. It suffices to dominate $`W_I,W_Y`$ and the entire angle family.

Set $`\phi=2\theta`$ and

```math
W_0=L\circ\bigl(B(M_U+M_V)/2+cH\bigr),
```

```math
W_c=L\circ B(M_U-M_V)/2,\qquad
W_s=L\circ BC_{\mathrm{rot}}/2.
```

Then $`W_\theta=W_0+W_c\cos\phi+W_s\sin\phi`$. For Hermitian $`X`$ define

```math
A_1=-(W_c-iW_s)/2,\qquad
Z=\begin{pmatrix}X&A_1\\A_1^\dagger&K-W_0-X\end{pmatrix}.
```

If $`Z\succeq0`$, its congruence with the stacked matrix having blocks $`I`$ and $`e^{i\phi}I`$ is $`K-W_\theta`$. Thus the finite constraints

```math
K\succeq0,\qquad K-W_I\succeq0,\qquad
K-W_Y\succeq0,\qquad Z\succeq0,
\qquad K_{ii}\le27/10
```

prove domination for every angle. No correction grid or necessity theorem for matrix spectral factorization is used.

## 4. Every encoding is represented

There are $`8!/2!=20{,}160`$ injective assignments. The 48 cube isometries preserve Hamming distances. Physical relabeling of the three identical slots reindexes Kraus histories and preserves every full target-Choi coordinate kernel, including the nonsymmetric $`U,V`$ cross kernel. These transformations transport $`G,D,K,X`$ by the same label permutation and preserve arbitrary nonuniform preparations.

Exact canonicalization gives 420 cube classes and 97 combined classes. Their orbit sizes are

```math
3\cdot48+49\cdot144+45\cdot288=20{,}160.
```

An independent Burnside check gives fixed-point sum $`20{,}160+7{,}488+288=27{,}936`$ for the group of order $`288`$, hence $`97`$ orbits. The [coverage verifier](../verification/encoding_bound/README.md) reconstructs a symmetry witness for every assignment, checks the kernel invariances exactly, and verifies that the rational certificate contains each class once.

## 5. Exact positive-matrix certificates

The [delivered certificate data](../verification/encoding_bound/README.md) contain rational $`K,X`$ for all 97 classes and auxiliary rational Gram factors. Their entries are integer arrays divided by $`10^6`$. Numerical optimization located candidates; the proof uses their subsequently checked rational bytes.

Let $`\widehat q=170763/200000`$. Exact rational squaring establishes

```math
0<\widehat q-10^{-7}<\sqrt{729/1000}<\widehat q<1.
```

At $`\widehat q`$, each of $`K,K-W_I,K-W_Y`$ and the real 24-dimensional representation of $`Z`$ has a supplied rational decomposition

```math
\widehat M=FF^T+R,
\qquad r_i=R_{ii}-\sum_{j\ne i}|R_{ij}|>0.
```

The verifier multiplies the integer Gram factors, clears denominators, and checks every residual row. Strict diagonal dominance proves $`\widehat M\succeq(\min_i r_i)I`$. Across all $`388`$ positive-matrix checks, the smallest exact margin is

```math
\mu_{\min}=\frac{503251611580961}{512000000000000000}
\approx0.000982913303869.
```

Each variable matrix entry is a coefficient of magnitude at most $`61/16<4`$ times $`q^j`$, with $`0\le j\le3`$. Since $`q,\widehat q\in[0,1]`$,

```math
\lVert M(q)-M(\widehat q)\rVert_{\mathrm{op}}
<24\cdot4\cdot3\cdot10^{-7}
<3/100000<\mu_{\min}.
```

The first bound uses the maximum absolute row-sum norm; the realification of a Hermitian matrix has its eigenvalues duplicated. This proves positivity at the actual algebraic visibility. Hermiticity and $`K_{ii}\le27/10`$ are checked exactly. No optimizer, floating-point eigenvalue, Cholesky routine or numerical success flag enters the standard-library verifier. This completes the certificate lemma.

## 6. From support to the optimal mixed call cost

For one call the normalized Choi spectrum has maximum $`1/2`$. Hence $`f_1\le p_1/2`$ and $`cp_1+Bf_1\le1`$.

For two calls at $`x=9/10`$, direct Kraus multiplication gives

```math
H_2=\begin{pmatrix}1&x/2\\x/2&1\end{pmatrix},
\qquad M_I=M_Y=M_V=\frac14\begin{pmatrix}1&x\\x&1\end{pmatrix},
```

```math
M_U=\frac14\begin{pmatrix}1&-x\\-x&1\end{pmatrix}.
```

All target-Bell cross kernels vanish. Every unitary kernel is therefore bounded by $`(19/40)I`$. If $`P_\pm`$ are the two order-parity projectors,

```math
BM_Q+cH_2\preceq B(19/40)I+cH_2
=4P_- -(229/80)P_+\preceq K_2:=4P_-.
```

The positive $`K_2`$ has diagonal two, so the same arbitrary diagonal-budget argument proves $`cp_2+Bf_2\le2`$. This is not a restriction to uniform two-call preparation.

Average these bounds and the three-call support bound over raw architecture probabilities $`t_1,t_2,t_3`$. For aggregate normalized-input Choi trace $`p`$ and Bell numerator $`f`$,

```math
cp+Bf\le t_1+2t_2+(27/10)t_3
=\mathbb E[m]-(3/10)t_3.
```

Only now use aggregate scaled trace preservation. The Bell input and Bell-projector measurement give $`\epsilon\ge1-f/p`$. At error at most $`3/8`$, $`f\ge5p/8`$, and the theorem's strict lower bound follows.

### Matching attainment

Choose a Hadamard-corrected one-call trial with raw probability $`21/61`$ and the balanced two-call minus-herald trial, corrected by $`U`$, with raw probability $`40/61`$. The latter has $`p_-=11/40`$ and $`f_-=19/80`$. Thus

```math
p=32/61,\qquad f=20/61,\qquad
\mathbb E[m]=101/61.
```

Both accepted maps are positive mixtures of the identity and unitary error channels whose Choi vectors are orthogonal to $`\Phi`$. The Bell lower bound and unitary triangle upper bound coincide: the aggregate half-diamond error is exactly $`3/8`$, with cost $`101/32`$. Its accepted reference–target Choi state is entangled because its Bell overlap is $`5/8>1/2`$, the separable-state maximum. This does not identify retained controller–target entanglement.

### Common-error corollary

The same support inequality and the one/two-call attaining mixtures give, for every $`3/22\le\epsilon\le1/2`$,

```math
C_{123,*}(\epsilon)=1+\frac{69}{4}(1/2-\epsilon).
```

The certificate does not exclude an accuracy extension below $`3/22`$ for arbitrary encodings. The promised fixed benchmark and this common-error interval are the matched cost conclusions.

## 7. Interpretation and reproducibility

The [verification guide](../verification/encoding_bound/README.md) gives a one-command standard-library exact check, the independent baseline/reduction/orbit audits and the full certificate inventory.

The earlier [2.299% alternative-encoding gain](local_controller_encoding.md#6-an-exact-advantage-from-another-encoding) at $`x=0.9604`$ remains valid at a different noise point. This all-assignment exclusion fixes $`x=0.9`$ and proportional exposure/routing duration; it is not an all-noise optimal phase diagram, a general quantum-encoding theorem, a communication-capacity statement or a hardware implementation result.
