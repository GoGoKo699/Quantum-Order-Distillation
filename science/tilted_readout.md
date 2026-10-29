# Optimal simple readout in a noncommuting neighborhood

The [orthogonal-family proof](recovery_frontier.md) works because its acceptance and four endpoint overlap matrices share the physical Walsh basis. Tilting the channel can destroy that simultaneous diagonalization. A diagonal operator envelope nevertheless proves that individual controller-qubit $`X`$ measurements still attain the entire recovery frontier on an explicit interval.

Use the [physical eight-order raw-slot model](noisy_routing_model.md), fix $`\epsilon=1/100`$ and local visibility $`\gamma=4/5`$, and set

```math
\theta=\pi/2+\delta,
\qquad t=\tan(\delta/2).
```

The two pure measure–prepare output columns are

```math
\frac{(1-t,1+t)^T}{\sqrt{2(1+t^2)}},
\qquad
\frac{(-1-t,1-t)^T}{\sqrt{2(1+t^2)}}.
```

The raw-slot restriction remains essential: the allowed target action is one herald-dependent unitary after the complete noisy switch, not a wrapper after each slot.

<a id="tilted-optimality"></a>
## 1. Certified full-frontier theorem

**Theorem.** For

```math
|t|\le\frac15,
\qquad\text{equivalently}\qquad
|\delta|\le2\arctan(1/5)
=0.39479111969976\ldots,
\qquad\text{(1)}
```

uniform controller preparation and the three physical controller $`X`$ measurements, with each outcome's final target correction optimized, attain the complete input-independent acceptance–Bell-fidelity frontier on the eight occupied orders. The competing class includes arbitrary finite ancillary controllers, mixed preparations, and collective accepted effects under the same access contract.

Let $`b_j(t)`$ be the acceptance eigenvalues in the physical Walsh basis, and let $`m_j(t)=\max_U(A^U_\Gamma(t))_{jj}`$ be each individual mode's best corrected overlap. Then

```math
u_*(p)=pF_*(p)
=\max\left\{
\sum_jm_j(t)x_j:
0\le x_j\le\frac18,
\quad\sum_jb_j(t)x_j=p
\right\}.
\qquad\text{(2)}
```

The allocation accepts modes in decreasing order of $`m_j(t)/b_j(t)`$, omits zero-capacity modes, and partially accepts the last needed mode. Outcome-dependent unitary optimization and partial acceptance are part of the attaining simple-readout protocol.

<a id="envelope-lemma"></a>
## 2. A sufficient diagonal-envelope condition

**Lemma.** Suppose $`B_\Gamma`$ is diagonal in a controller basis and

```math
D=\mathrm{diag}(m_j),
\qquad m_j=\max_U(A^U_\Gamma)_{jj},
\qquad D\succeq A^U_\Gamma
\text{ for every final }U.
\qquad\text{(3)}
```

Under the uniform coefficient budget, the complete controller frontier equals the mode allocation (2).

**Proof.** The [controller symmetry theorem](controller_optimization.md#uniform-preparation) gives $`\sum_xH_x\preceq I/8`$. Apply $`\mathrm{Tr}[(D-A^{U_x}_\Gamma)H_x]\ge0`$ to every outcome. The resulting bound depends only on the diagonal entries of $`\sum_xH_x`$ and has precisely the constraints in (2). Measuring the stated controller basis, optimizing its individual corrections and thinning its outcomes attains every allowed diagonal allocation. Here that basis is the physical $`X`$ basis, so these effects are implementable. The [real-controller reduction](controller_optimization.md#real-controller-reduction) gives constant-success attaining maps. $`\square`$

<a id="exact-tilted-arrays"></a>
## 3. Reconstructing the tilted objective matrices exactly

The finite identities below are part of a computer-assisted proof, rather than numerical samples. They are recorded in [tilt_exact_polynomials.npz](../verification/certificates/tilt_exact_polynomials.npz) and independently reconstructed by [audit_verify.py](../verification/audit_verify.py).

For a tuple with $`k`$ pure measure–prepare labels, use the two polynomial matrices

```math
L_0(t)=\begin{pmatrix}1-t&0\\1+t&0\end{pmatrix},
\qquad
L_1(t)=\begin{pmatrix}0&-1-t\\0&1-t\end{pmatrix},
```

and the four constant matrix units. Form all eight word products $`P_{\boldsymbol a,\pi}(t)`$ for every one of the $`6^4`$ tuples. The tuple's common squared physical scale is

```math
\frac{99^k}{16\,100^4(1+t^2)^k}.
\qquad\text{(4)}
```

Let $`W_{\pi j}=(-1)^{b_\pi\cdot j}`$ be the unnormalized physical Walsh matrix and

```math
(G_n)_{\pi\sigma}
=4^{d_H(b_\pi,b_\sigma)}5^{3-d_H(b_\pi,b_\sigma)},
\qquad\Gamma=G_n/125.
```

Use the real matrix basis $`(M_0,M_1,M_2,M_3)=(I,iY,X,Z)`$. Define tuple matrices

```math
(T_{ab})_{\pi\sigma}
=\overline{\mathrm{Tr}(M_aP_{\boldsymbol a,\pi})}
\mathrm{Tr}(M_bP_{\boldsymbol a,\sigma}),
\qquad
S_{\pi\sigma}
=\mathrm{Tr}(P_{\boldsymbol a,\pi}^TP_{\boldsymbol a,\sigma}).
```

The exact polynomial numerators in the normalized Walsh basis are

```math
G_{ab}^{\mathrm{num}}(t)
=\sum_{\boldsymbol a}99^k(1+t^2)^{4-k}
W^T(G_n\circ T_{ab})W,
```

```math
B^{\mathrm{num}}(t)
=\sum_{\boldsymbol a}99^k(1+t^2)^{4-k}
W^T(G_n\circ S)W.
\qquad\text{(5)}
```

Both have degree at most eight. The common overlap denominator and the acceptance normalization are

```math
D(t)=64\,100^4\,125\,8\,(1+t^2)^4
=6400000000000(1+t^2)^4,
```

```math
G_{ab}(t)=G_{ab}^{\mathrm{num}}(t)/D(t),
\qquad B_\Gamma(t)=B^{\mathrm{num}}(t)/(D(t)/2).
\qquad\text{(6)}
```

The archive entries `Gx[p,a,b,i,j]` and `Bx[p,i,j]`, for powers `p=0,...,8`, are the coefficient arrays of these numerators. The factor 8 in $`D`$ is the Walsh basis normalization; the separate factor $`1/8`$ in an actual accepted outcome comes from uniform input preparation. At every $`t`$, $`\mathrm{Tr}B_\Gamma=8`$, consistent with total acceptance one.

On either correction circle $`A(\varphi)=C_0+C_c\cos(2\varphi)+C_s\sin(2\varphi)`$. Its matrices are recovered from the corresponding component pair $`(a,b)=(0,1)`$ or $`(2,3)`$ by

```math
C_0=(G_{aa}+G_{bb})/2,
\quad C_c=(G_{aa}-G_{bb})/2,
\quad C_s=(G_{ab}+G_{ba})/2.
\qquad\text{(7)}
```

Consequently each mode's reward in (2) is explicitly

```math
m_j(t)=\max_{\ell\in\{\parallel,\perp\}}
\left\{(C_{0,\ell})_{jj}
+\sqrt{(C_{c,\ell})_{jj}^2+(C_{s,\ell})_{jj}^2}\right\}.
\qquad\text{(8)}
```

This formula, (5)–(7), and the finite arrays specify the complete continuum frontier without a numerical SDP fit.

<a id="block-majorant-proof"></a>
## 4. The all-unitary block majorant

The [correction-sector reduction](controller_optimization.md#real-controller-reduction) applies throughout this family, so a common majorant for its two real circles covers all final qubit unitaries, including complex controller coefficients.

Exact coefficient reconstruction gives the following support identities in the physical Walsh basis:

- $`B_\Gamma`$ is diagonal.
- Every matrix on the $`I/iY`$ circle is diagonal.
- On the $`iX/iZ`$ circle, $`C_0`$ is diagonal and off-diagonal harmonics occur only in the mode pairs $`(1,5)`$ and $`(2,6)`$.

For one nontrivial block write $`A(\alpha)=C_0+C_c\cos\alpha+C_s\sin\alpha`$, where $`\alpha=2\varphi`$. Let $`h_1,h_2`$ be its diagonal harmonic vectors $`(C_{c,jj},C_{s,jj})`$, and let $`k`$ be its off-diagonal harmonic vector. The coefficient identities also give

```math
h_1\times h_2=0,
\qquad h_1\cdot k=h_2\cdot k=0.
\qquad\text{(9)}
```

Here the two-dimensional cross product means the scalar determinant. The two mode pairs have identical diagonal vectors and opposite off-diagonal vectors, so they have the same positivity margin.

If

```math
-h_1\cdot h_2-\lVert k\rVert^2>0,
\qquad\text{(10)}
```

then $`h_1,h_2`$ are nonzero and antiparallel. Rotate the harmonic angle to put the block in the form

```math
\begin{pmatrix}a_0&0\\0&a'_0\end{pmatrix}
+\begin{pmatrix}r_1&0\\0&-r_2\end{pmatrix}\cos\alpha
+\kappa\begin{pmatrix}0&1\\1&0\end{pmatrix}\sin\alpha,
\qquad r_1,r_2>0.
\qquad\text{(11)}
```

Its independent diagonal maxima are $`a_0+r_1`$ and $`a'_0+r_2`$. Their diagonal majorant minus (11) has nonnegative diagonal entries $`r_1(1-\cos\alpha)`$ and $`r_2(1+\cos\alpha)`$, and determinant

```math
(r_1r_2-\kappa^2)\sin^2\alpha\ge0.
\qquad\text{(12)}
```

It is positive semidefinite for every correction angle. Taking each diagonal maximum also over the diagonal $`I/iY`$ sector only increases the majorant. Thus (10) proves condition (3), even though the objective matrices need not commute.

<a id="bernstein-certificate"></a>
## 5. A positive rational certificate on the whole interval

The doubled harmonic numerators in (7) share denominator $`2D(t)`$. Clearing their common factors gives exactly

```math
-h_1\cdot h_2-\lVert k\rVert^2
=\frac{231234600960000\,q_{16}(t)}{4D(t)^2}.
\qquad\text{(13)}
```

Write $`q_{16}(t)=\sum_{j=0}^{16}q_jt^j`$. Its coefficients are palindromic, $`q_j=q_{16-j}`$, and hence specified completely by:

| $`j`$ | $`q_j`$ |
|---:|---:|
| 0 | 11150726384 |
| 1 | 69699520000 |
| 2 | 89985970672 |
| 3 | 60943517800 |
| 4 | 52161035996 |
| 5 | -671075691000 |
| 6 | -1364298119321 |
| 7 | -1830048022000 |
| 8 | -2675247822058 |

Map $`t=(2v-1)/5`$ with $`v\in[0,1]`$, and let $`c_j`$ be the resulting monomial coefficients. The coefficients in the degree-16 Bernstein basis are

```math
\beta_k=\sum_{j=0}^kc_j
\frac{\binom{k}{j}}{\binom{16}{j}}.
```

All 17 exact rational coefficients are positive, with

```math
\min_k\beta_k
=\frac{84073951440936790684}{152587890625}>0.
\qquad\text{(14)}
```

The Bernstein basis is nonnegative and sums to one, so (14) proves $`q_{16}(t)>0`$ on $`[-1/5,1/5]`$. Equations (9) and (13) then establish the all-unitary diagonal envelope. The lemma proves the full-frontier theorem (1)–(2).

The bound is uniform. Combining the minimum in (14) with the largest denominator on this interval gives

```math
-h_1\cdot h_2-\lVert k\rVert^2
>0.0005682082787.
```

The verifier checks the support identities, vector orthogonality, equal pair margins, all polynomial coefficients, and all Bernstein coefficients. With `--rebuild`, it regenerates every coefficient using Python integers. It also bounds the original fixed-width absolute intermediate sums by $`276540930048000`$ for overlap and $`305648793600000`$ for acceptance, both below $`2^{63}`$. An independent direct rational Kraus contraction at $`t=1/10`$ checks the physical normalization in (6). These finite reconstructions support the exact proof obligations; floating-point optimization is not used to prove the interval.

## 6. Evaluated optima and their interpretation

Representative numerical evaluations of the attaining formulas are:

| Tilt $`\delta`$ in radians | Maximum $`F`$ | Acceptance at maximum $`F`$ | $`F`$ at $`p=.1245025`$ |
|---:|---:|---:|---:|
| 0 | 0.758264363 | 0.061963764 | 0.688117779 |
| 0.05 | 0.768570252 | 0.067010657 | 0.694780568 |
| 0.20 | 0.791654670 | 0.083660370 | 0.718752556 |

At these points, the best-fidelity herald is still the $`P_2`$ controller outcome, with its final correction optimized. Holding the zero-tilt correction fixed would compare against a weaker protocol than the theorem permits.

Noncommutation of objective matrices alone therefore does not establish an advantage for a collective controller measurement. The stronger operator envelope can make the physical product readout globally optimal. This conclusion covers the stated eight-order architecture, noise point and tilt interval; it is not an assertion about every noisy switch or all 24 orders.

For the [sister entanglement separation](https://github.com/GoGoKo699/Quantum-Order-Entanglement/blob/main/science/noisy_separation.md), a smaller interval $`|\delta|\le.2`$ is certified using a fixed identity correction. That comparison is protocol-specific and concerns distinct controller–target and reference–target bipartitions. It is independent of a converse claim that fidelity at or below $`1/2`$ implies separability.
