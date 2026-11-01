# Recovery with an encoded local phase-noise controller

Independent dephasing of physical control bits assigns different decay factors to different order pairs. For one natural six-word encoding, the complete three-call recovery frontier remains a four-flag problem in a proved coherence regime. Another encoding gives a certified lower-noise cost advantage. The encoding and the exposure window are therefore part of the resource model.

This page concerns the fixed rank-two channel from the [native benchmark](native_three_call_frontier.md), with $`K_0=|+\rangle\langle0|`$, $`K_1=|-\rangle\langle1|`$ and $`\mathcal M^2=\Delta`$. Its three-call results are distinct from both global logical dephasing and the [four-call full-rank physical-noise family](recovery_frontier.md).

## 1. Access, encoding and noise

Use independent compulsory channel slots. Every occupied order traverses every slot once. The target has no intermediate gates, protected encoding or bypass, accessible environment, or coding across failed trials. The controller starts independently of target and reference; arbitrary finite mixed and ancillary controls and final controller instruments are allowed. Fine outcomes may select known target unitaries. Noise acts before all final controller processing, with no protected transfer before noise.

The evaluated corrected accepted target-only aggregate is $`\mathcal E=p\Lambda`$, with aggregate input-independent $`p>0`$ and $`\Lambda`$ trace preserving. Individual fine outcomes and trial components may have input-dependent acceptance. Their corrected labels are averaged. Write

```math
\epsilon=\tfrac12\lVert\Lambda-\mathrm{id}\rVert_\diamond,
\qquad f=\langle\Phi|J(\mathcal E)|\Phi\rangle,
\qquad |\Phi\rangle=(|00\rangle+|11\rangle)/\sqrt2.
```

All Choi matrices use this normalized Bell input. The matrix-product word's rightmost slot acts first.

| Word | Pairwise-comparison encoding | Alternative encoding used in Section 6 |
|---|---|---|
| 012 | 000 | 000 |
| 021 | 001 | 011 |
| 102 | 100 | 101 |
| 120 | 110 | 010 |
| 201 | 011 | 100 |
| 210 | 111 | 001 |

In the comparison encoding, each bit says whether the first member of the named pair $`(0,1),(0,2),(1,2)`$ occurs after the second. The occupied words are the six transitive tournaments; the excluded words 010 and 101 are opposite cube vertices. Identical independent phase damping gives

```math
L_{\pi\sigma}(q)=q^{d_H(c_\pi,c_\sigma)},\qquad
q=e^{-\gamma T_3},\qquad
\mathcal L_C=\tfrac\gamma2\sum_{k=1}^3(\mathrm{Ad}_{Z_k}-\mathrm{id}).
```

The two-order controller uses one physical bit, with $`x=e^{-\gamma T_2}`$. The routing primitive is diagonal in order labels, so the declared phase noise commutes with routing and can be placed immediately after it. In the comparison encoding, Hamming distance is Kendall distance. Slot relabeling permutes and complements comparison bits, preserving that metric; reversing a word complements all three bits. These symmetries are specific to this encoding.

The charged serial schedule is $`T_m=m\tau`$. Its routing windows are also the controller exposure windows. Preparation, measurement and correction are ideal and their overhead is uncharged. Gate faults, loss, leakage, a gate-level router and measured laboratory elapsed time remain unverified. Calls and charged routing time per accepted fresh target are $`\mathbb E[m]/p`$ and $`\mathbb E[T_m]/p`$ respectively.

<a id="local-four-flag-frontier"></a>
## 2. Exact four-flag theorem

Use the [native projectors](native_three_call_frontier.md#2-four-native-controller-flags) $`T,A,O,E`$, with $`A`$ the rank-one sign projector. The verification programs call this same sign projector $`S`$. Put $`U=(X-Z)/\sqrt2`$ and $`V=(X+Z)/\sqrt2`$.

Let $`q_0`$ be the unique root in $`(0,1)`$ of

```math
8q_0^4+2q_0^3+q_0^2+q_0-6=0,
\qquad q_0=0.811155656051003\ldots.
```

**Theorem.** For the comparison encoding and $`q\in[q_0,1]`$, uniform control, this PVM, the corrections below, and fractional flag acceptance attain the complete success–half-diamond-error frontier over the controller model in Section 1.

| Flag | Probability $`p_z`$ | Corrected Bell numerator $`f_z`$ | Correction |
|---|---|---|---|
| $`O`$ | $`(4+2q-q^2-2q^3)/12`$ | $`(2+q+q^2+2q^3)/24`$ | $`U`$ |
| $`E`$ | $`(4-2q-q^2+2q^3)/12`$ | $`(1+q^3)/12`$ | $`Y`$ |
| $`T`$ | $`(2+2q+q^2+q^3)/12`$ | $`(1+2q+2q^2+q^3)/24`$ | $`I`$ |
| $`A`$ | $`(2-2q+q^2-q^3)/12`$ | $`(1-q^3)/24`$ | $`Y`$ |

In particular, the sign flag uses $`Y`$, whereas the global-dephasing native table uses $`I`$. The exact optimum is

```math
f_{3,*}(p,q)=\max\left\{\sum_zf_za_z:
0\le a_z\le1,\quad\sum_zp_za_z=p\right\},
\qquad \epsilon_{3,*}(p,q)=1-f_{3,*}(p,q)/p.
```

The decreasing fidelity order is $`O,E,\{T,A\}`$, with

```math
\frac{f_T}{p_T}=\frac{f_A}{p_A}
=\frac{1+q+q^2}{2(2+q^2)}.
```

Omit the zero-probability sign flag at $`q=1`$. Below $`q_0`$ the same native maps remain feasible, but their global optimality is not asserted.

For the ratio order, the nontrivial positive-denominator cross-products reduce to $`q(-4+2q+5q^2-7q^3+2q^4+8q^5)`$ for $`O,E`$ and $`q(-2+q+5q^2-q^3)`$ for $`E,T`$. Both bracketed polynomials are positive and increasing on $`[4/5,1]`$, which contains the theorem's interval.

## 3. Controller-independent converse

The [Gram reduction](controller_optimization.md#controller-feasibility) supplies $`G_y\succeq0`$ with $`\sum_yG_y\preceq D`$, where $`D`$ is diagonal, positive and has trace one. Local noise acts by $`G_y\mapsto G_y\circ L(q)`$. Set

```math
(t_G,a_G,o_G,e_G)
=(\mathrm{Tr}(TG),\mathrm{Tr}(AG),\mathrm{Tr}(OG),\mathrm{Tr}(EG)).
```

Here $`e_G`$ is an even-standard-sector mass, not the error. Define

```math
\ell_T=1+2q+2q^2+q^3,\qquad
\ell_A=1-2q+2q^2-q^3,
```

```math
\ell_O=1+q-q^2-q^3,\qquad
\ell_E=1-q-q^2+q^3.
```

Slot covariance and reversal symmetry make the target kernels diagonal between the trivial, sign and standard sectors. The only target-Bell cross block is $`U,V`$. Direct contraction of the [ideal Bell kernels](native_three_call_frontier.md#exact-bell-kernels) gives the diagonal coefficients:

```math
J_I=\tfrac14(\ell_Tt_G+\ell_Aa_G+\ell_Oo_G+\ell_Ee_G),
```

```math
J_U=\tfrac{\ell_O}{4}t_G+\tfrac{\ell_E}{4}a_G+u_Oo_G+u_Ee_G,
```

```math
J_Y=\tfrac{\ell_T+2\ell_E}{12}t_G
+\tfrac{1-q^3}{4}(a_G+o_G)+\tfrac{1+q^3}{4}e_G,
```

```math
J_V=\tfrac{2\ell_T+\ell_E}{12}t_G
+\tfrac{2\ell_A+\ell_O}{12}a_G+v_Oo_G+v_Ee_G,
```

where

```math
u_O=(2+q+q^2+2q^3)/8,\qquad
v_O=(2+q-q^2-2q^3)/8,
```

```math
u_E=(2-q+q^2-2q^3)/8,\qquad
v_E=(2-q-q^2+2q^3)/8.
```

For example, $`T\circ L=L/6`$ has eigenvalues $`\ell_z/6`$ on these four sectors. Also $`A\circ L`$ is the sign conjugation of $`L/6`$, while $`O\circ L=(I-q^3R)/2-A\circ L`$ and $`E\circ L=(I+q^3R)/2-T\circ L`$. These identities derive the displayed coefficients from the ideal kernels. The exact verifier independently reconstructs them from all Kraus histories.

### The complex cross term

The real symmetric $`U,V`$ kernel is supported on Hamming-distance-one entries and its imaginary antisymmetric part on distance two. Thus, for the standard-sector Choi cross term,

```math
\chi_q=q\mathrm{Re}\chi_1+iq^2\mathrm{Im}\chi_1,
\qquad |\chi_q|^2\le\frac{3q^2}{16}o_Ge_G.
```

The second inequality follows from the ideal standard-sector positive $`U,V`$ block, whose diagonal entries are $`3o_G/4,e_G/4`$. This covers complex $`G`$.

The largest trivial and sign column coefficients are respectively $`c_T=\ell_T/4`$ and $`c_A=(1-q^3)/4`$. For $`q\ge q_0`$, the standard columns have maxima $`c_O=u_O`$ and $`c_E=(1+q^3)/4`$. To bound every eigenvector of the standard $`U,V`$ block, put $`b=c_Oo_G+c_Ee_G`$. The diagonal entries of $`bI-J_{UV}`$ are

```math
(c_E-u_E)e_G,
\qquad (c_O-v_O)o_G+(c_E-v_E)e_G.
```

Their factors are nonnegative, and

```math
(c_O-v_O)(c_E-u_E)-\frac{3q^2}{16}
=\frac{q^2}{32}(8q^4+2q^3+q^2+q-6)\ge0.
```

The determinant is therefore nonnegative by the cross-term bound. The isolated $`I,Y`$ entries obey the same standard-sector bound. Positivity and subadditivity across the three isotypic compressions now give

```math
\lambda_{\max}J(G\circ L)
\le c_Tt_G+c_Aa_G+c_Oo_G+c_Ee_G.
```

A unitary-corrected Bell overlap cannot exceed this largest eigenvalue. Constant projector diagonals give accepted masses at most $`1/6,1/6,1/3,1/3`$ for $`T,A,O,E`$. The trace kernel is diagonal in these same four sectors and has coefficients $`6p_T,6p_A,3p_O,3p_E`$. Rescaling masses by their capacities yields exactly the theorem's acceptance box and trace constraint.

### Attainment and the metric

The native coefficient matrices are $`G_z=a_zP_z/6`$. Their corrected maps are positive mixtures of the unitary channels $`I,U,Y,V`$, with identity Choi weight $`f/p`$. They are therefore scaled trace-preserving maps. All nonidentity Bell vectors are orthogonal to $`\Phi`$. The Bell-input lower bound and the unitary-channel triangle upper bound agree, so their half-diamond error is exactly $`1-f/p`$. This proves matching attainment.

## 4. Mixing call budgets and charging duration

Combine the local flag table with the one-call point $`(p,f)=(1,1/2)`$ and two-call flags $`((2\mp x)/4,(1+x)/8)`$. The [accepted-flow LP and convex-hull proof](randomized_call_cost.md#convex-call-cost) applies unchanged: use any specified per-trial prices, and average corrected target-only outcomes. Mixtures with one-call $`Z`$ dephasing still have error Choi support orthogonal to $`\Phi`$, so Bell deficit remains the exact half-diamond error.

At $`x=0.9`$, put $`B=69/4`$. A three-call policy can beat the affine one/two-call cost support only if

```math
P+B(F-P/2)>3.
```

At $`q=0.9`$, the exact four-flag support maximum is

```math
\sum_z[p_z+B(f_z-p_z/2)]_+
=39529/16000=2.4705625<3.
```

Here all four summands are positive. The three-call minimum error is the $`O`$ error $`237/883`$, larger than the two-call floor $`3/22`$. Additional dephasing can be applied before final measurement, so both exclusions extend to all $`q\le0.9`$ for this encoding, including below the regime where the native frontier is proved. In particular, every ratio $`T_3/T_2\ge1`$ excludes query activation under the calibrated exponential exposure law.

For a separate time-price sensitivity with $`T_1/T_2=1/2`$ and $`r=T_3/T_2`$, let

```math
F_{\mathrm{all}}(q)=\frac{6+3q+3q^2+4q^3}{24},
\qquad q=0.9^r.
```

A strict routing-time saving at an interior common error exists precisely for $`1\le r<r_*`$, where

```math
r_*=\frac12+\frac{69}{8}\bigl(F_{\mathrm{all}}(0.9^{r_*})-\tfrac12\bigr),
```

```math
1.14647049186791<r_*<1.14647049186792.
```

The support difference is strictly decreasing in $`r`$; stronger dephasing and a larger three-call price preserve the exclusion beyond the crossing. Ratios below $`3/2`$ here are conditional prices, not implementations of three identical serial slots of duration $`\tau`$ when $`T_2=2\tau`$. The serial schedule has $`r=3/2`$, so the fixed $`x=0.9`$ comparison encoding saves neither calls nor charged routing time. The [all-assignment converse](all_encoding_cost_bound.md) strengthens this fixed-benchmark conclusion to every permitted injective encoding.

## 5. Exact comparison-encoding activation at lower noise

For $`T_m=m\tau`$, write $`t=e^{-\gamma\tau}`$, $`x=t^2`$ and $`q=t^3`$. On $`x\in[0.9,1]`$, all four support coefficients against the optimal one/two-call affine bound are nonnegative. Their sum exceeds three exactly when

```math
P(t)=4t^{11}+24t^9+3t^8+18t^6+3t^5+18t^3-54t^2-12>0.
```

Indeed, the sum is $`1+B(x)(F_{\mathrm{all}}(q)-1/2)`$, where $`B(x)=2(6+x)/(2x-1)`$. For $`T,A`$ positivity follows from $`2(2+q^2)-B(x)(1-q)>0`$ using $`B(x)\le69/4`$ and $`q>4/5`$; the $`O,E`$ fidelity-surplus terms are positive there. The crossing is the displayed polynomial. Its derivative is positive on $`[0.9,1]`$: already $`216t^8+54t^2-108t>0`$ there, and the other derivative terms are positive. Exact rational sign brackets give

```math
t_*=0.988094367112997\ldots,
\qquad x_*=t_*^2=0.976330478320434\ldots.
```

This is also the criterion for any comparison-encoding activation on that interval: its $`O`$ error cannot improve the two-call floor, since

```math
\epsilon_O(q)-\frac{3(1-q)}{2(2-q)}
=\frac{3q(1-q^2)}{(2-q)(4+2q-q^2-2q^3)}\ge0,
```

and $`q=t^3\le t^2=x`$. Thus $`x>x_*`$ gives a cost benefit at common feasible errors; it is not an accuracy-floor extension. The same criterion holds for the proportional time price. It is not a universal encoding threshold.

## 6. An exact advantage from another encoding

Keep the same local noise and $`T_m=m\tau`$, but use the alternative code in Section 1 and choose $`t=49/50`$, so $`x=2401/2500=0.9604`$ and $`q=117649/125000=0.941192`$. Uniform control and the same PVM, rejecting $`A`$ and correcting $`T,O,E`$ by $`I,U,Y`$, give an exact constructive advantage when mixed with one-call Hadamard-corrected trials.

| Quantity | Exact value |
|---|---|
| Raw one-call probability | $`1505323412240983/13224073412240983`$ |
| Raw three-call probability | $`11718750000000000/13224073412240983`$ |
| Aggregate acceptance | $`12985210999106966/13224073412240983`$ |
| Attained calls per accepted target | $`36661573412240983/12985210999106966\approx2.823332899`$ |
| Optimal one/two-call baseline | $`26609/9208\approx2.889769765`$ |
| Relative saving | Approximately 2.2990% |

The [exact protocol verifier](../verification/local_controller/README.md) reconstructs every branch from Kraus words, checks all Bell offdiagonals and partial traces, and verifies the raw mixture and strict baseline inequality with rational arithmetic. The corrected aggregate's identity Bell weight is exactly $`5/8`$, so its half-diamond error is exactly $`3/8`$. Its accepted reference–target Choi state is entangled because $`5/8>1/2`$, the separable-state Bell bound. This is distinct from retaining controller–target entanglement.

This is an attaining protocol against a globally optimal one/two-call comparator. It is not a global optimum over alternative-code controls, encodings or layouts. At this same noise point the comparison encoding has no activation, demonstrating encoding dependence without changing the serial schedule. The all-assignment exclusion at $`x=0.9`$ concerns a different operating point and is compatible with this lower-noise example.

## Verification and boundaries

The [local-controller verification guide](../verification/local_controller/README.md) covers exact spectral/kernel identities, native branch polynomials, threshold brackets and the alternative construction. The [all-encoding certificate](all_encoding_cost_bound.md) covers the fixed calibrated benchmark over all injective six-label bit assignments.

The high-coherence frontier is specific to the comparison encoding. The all-encoding result fixes one noise point and does not supply an optimal phase diagram. Neither result certifies redundant eight-word encodings, protected/noncomputational control codes, a different noise map, target filters, extra calls, gate errors or hardware routing overhead.
