# Exact cost of randomized one-, two- and three-call trials

A trial may choose its call count classically, use fresh target and reference inputs, and accept selected corrected controller outcomes. Charging both successful and unsuccessful trials gives an exact finite linear program for the best expected noisy calls per accepted fresh output. Its feasible native flags come from the [complete three-call theorem](native_three_call_frontier.md#native-three-call-frontier) and the one- and two-call proofs below.

Write $`\mathrm{Ad}_W(\rho)=W\rho W^\dagger`$. The raw channel is fixed throughout: $`K_0=|+\rangle\langle0|`$, $`K_1=|-\rangle\langle1|`$. Slots are independent and compulsory; the controller begins independently of target and reference; environments, protected bypasses and intermediate target operations are unavailable. Global logical order dephasing acts after routing and before final controller processing. Let $`x=\eta_2`$ and $`z=\eta_3`$ be its two- and three-call visibilities. They are independent architecture parameters unless a duration law is explicitly imposed.

Fine controller outcomes, final correction labels and trial-type labels are averaged in the accepted target-only channel. The total accepted map must be $`p\Lambda`$ with input-independent $`p\gt0`$. Its half-diamond error is $`e=\frac12\lVert\Lambda-\mathrm{id}\rVert_\diamond`$. No batch coding or collective distillation across failed outputs is included. Call cost is distinct from routing, clock time and controller-hardware cost.

## 1. Complete one- and two-call baselines

One call followed by a Hadamard gives $`D_Z=(\mathrm{id}+\mathrm{Ad}_Z)/2`$. It has acceptance one, Bell fidelity $`1/2`$ and half-diamond error $`1/2`$. With one routing label, a controller independent of the input can only scale the accepted map. A final target unitary cannot exceed the original Choi largest eigenvalue $`1/2`$, so no one-call strategy improves that error. Partial acceptance only increases the call cost per accepted trial.

For two calls, take orders $`01,10`$ and prepare the uniform order qubit. The symmetric and antisymmetric controller projectors are $`P_\pm=(I\pm S)/2`$, where $`S`$ exchanges those words. Before corrections, direct multiplication of the four environmental histories gives

```math
\mathcal E_-^{\mathrm{ideal}}=\tfrac14\mathrm{Ad}_U,
\qquad
\mathcal E_+^{\mathrm{ideal}}
=\tfrac14(\mathrm{id}+\mathrm{Ad}_Y+\mathrm{Ad}_V),
\qquad
U=(X-Z)/\sqrt2,\quad V=(X+Z)/\sqrt2.
```

Each definite order implements $`\mathcal M^2=\Delta`$. Global dephasing therefore gives

```math
\mathcal E_\pm(x)=x\mathcal E_\pm^{\mathrm{ideal}}
+\frac{1-x}{2}\Delta.
```

Correct the minus outcome by $`U`$ and the plus outcome by $`I`$. Their exact data are:

| Flag | Acceptance | Corrected Bell overlap | Calls per attempt |
|---|---|---|---:|
| Minus | $`p_-=(2-x)/4`$ | $`f_-=(1+x)/8`$ | 2 |
| Plus | $`p_+=(2+x)/4`$ | $`f_+=(1+x)/8`$ | 2 |

### Why this is the full two-call optimum

This table is not a restriction to a chosen measurement. Write $`\mathcal S_{ij}`$ for the two-order target cross maps. Identical independent slots imply $`\mathcal S_{00}=\mathcal S_{11}=\Delta`$ and $`\mathcal S_{01}=\mathcal S_{10}`$: exchanging the two environmental indices turns either cross sum into the other. Thus, for any complex positive accepted coefficient matrix $`G`$, the full target map depends only on $`\mathrm{Tr}G`$ and the real part of $`G_{01}`$. In the $`P_+,P_-`$ basis it depends only on the nonnegative mode masses $`h_\pm=\mathrm{Tr}(P_\pm G)`$ and equals

```math
\mathcal E_x(G)=2h_-\mathcal E_-(x)+2h_+\mathcal E_+(x).
```

The [general controller budget](controller_optimization.md#controller-feasibility) gives $`\sum_yG_y\preceq D`$ for a diagonal trace-one $`D`$. Since both projector diagonals equal $`1/2`$,

```math
\sum_yh_{-,y}\le1/2,\qquad\sum_yh_{+,y}\le1/2.
```

Known final target unitaries preserve Choi spectra. The largest eigenvalue of either uncorrected noisy native Choi matrix is $`(1+x)/8`$; the stated corrections attain it. Largest-eigenvalue subadditivity bounds every accepted outcome by the same positive combination of native rewards. Rescaling the two mode masses therefore gives exactly fractional acceptance of the table, with no improvement from mixed, ancillary or complex controllers.

Every two-call accepted branch is a scalar multiple of a trace-preserving channel, including arbitrary fine branches before or after their final corrections. Native corrected mixtures have an identity Choi component and a positive nonidentity component orthogonal to it, so their half-diamond error equals their Bell deficit, as in the three-call attainment proof. Hence the table gives the full two-call success–diamond-error frontier, not just a Bell-input bound.

### The optimal randomized one-/two-call baseline

For $`x\gt1/2`$, the minimum feasible error and its cost are the fully accepted minus flag:

```math
e_2=\frac{3(1-x)}{2(2-x)},\qquad c_2=\frac8{2-x}.
```

The deterministic two-call point has $`e=(3-x)/4\ge1/2`$ and cost two, and is dominated by one call. Between cumulative native vertices, cost is affine in error. The exact one-/two-call lower cost hull therefore joins $`(e_2,c_2)`$ to $`(1/2,1)`$:

```math
C_{12}(e;x)=1+B(1/2-e),\qquad
B=\frac{2(6+x)}{2x-1},\qquad e_2\le e\le1/2.
```

It is infeasible below $`e_2`$ and equals one for $`e\ge1/2`$. For $`x\le1/2`$, one call dominates every two-call strategy, and no one-/two-call mixture reaches error below $`1/2`$.

<a id="convex-call-cost"></a>
## 2. The exact finite linear program

For each call count $`m`$, use the native flag pairs $`(p_{mi},f_{mi})`$. One call has the single pair $`(1,1/2)`$; two calls use the preceding table. Three calls use:

| Flag | Acceptance | Corrected Bell overlap |
|---|---|---|
| $`O`$ | $`(4-z)/12`$ | $`(1+2z)/12`$ |
| $`E`$ | $`(4-z)/12`$ | $`(1+z)/12`$ |
| $`T`$ | $`(1+2z)/6`$ | $`(1+5z)/24`$ |
| $`A`$ | $`(1-z)/6`$ | $`(1-z)/24`$ |

**Theorem.** The exact optimum among all randomized one-, two- and three-call fresh trials is

```math
\begin{aligned}
\mathrm{minimize}\quad&\sum_{m=1}^3m v_m,\\
\mathrm{subject\ to}\quad&v_m\ge0,\quad0\le u_{mi}\le v_m,\\
&\sum_{m,i}p_{mi}u_{mi}=1,\\
&\sum_{m,i}f_{mi}u_{mi}\ge1-e.
\end{aligned}
```

Here $`v_m`$ is the mass of attempts of type $`m`$ per accepted output, and $`u_{mi}`$ is that mass times the acceptance fraction of native flag $`i`$. Equivalently, the answer is the lower convex envelope in $`(e,C)`$ of the one-call point and the two- and three-call cumulative-prefix vertices

```math
P_k=\sum_{i=1}^kp_i,\qquad F_k=\sum_{i=1}^kf_i,
\qquad (e_k,C_k)=(1-F_k/P_k,m/P_k).
```

Flags are ordered by decreasing $`f_i/p_i`$. Omit zero-success prefixes and keep the smaller cost at duplicate errors.

**Converse.** For raw trial probabilities $`t_j`$ and trial acceptance $`p_j`$, the accepted weights are $`\alpha_j=t_jp_j/\sum_\ell t_\ell p_\ell`$. Expected calls per accepted output are $`\sum_j\alpha_jm_j/p_j`$. Bell overlap is linear in unnormalized accepted Choi matrices, and the final half-diamond error obeys $`e\ge1-F_\Phi`$.

The one- and two-call converse above applies to all their branches. Their aggregate contributions are scaled trace preserving. Subtracting them from the input-independent total accepted map shows that its aggregate three-call contribution is scaled trace preserving too. All random three-call choices can be grouped into one allowed mixed or ancillary three-call controller. The complete native three-call converse therefore applies, even if its fine outcomes have input-dependent success. Replacing each grouped contribution by a native strategy with the same acceptance and at least its Bell overlap gives the displayed LP after scaling total success to one. No conditional channel is assigned to an input-dependent fine piece.

**Attainment.** For any LP point, choose raw call-count probabilities $`t_m=v_m/\sum_nv_n`$ and flag acceptance fractions $`u_{mi}/v_m`$, omitting a type with $`v_m=0`$. Overall success is $`1/\sum_nv_n`$, so expected calls per acceptance are exactly the objective.

The native two- and three-call channels have traceless nonidentity unitary components, whose Bell vectors are orthogonal to the identity Bell vector. The one-call nonidentity unitary is $`Z`$, also traceless. It need not be orthogonal to each other error vector: their positive combined error Choi matrix still has support orthogonal to $`\Phi`$. Thus the mixture's Bell lower bound and channel triangle upper bound coincide, giving exact half-diamond error equal to its Bell deficit. This realizes every LP point. For a fixed call family, each cost segment is affine in error, so cumulative-prefix vertices suffice for its lower hull. $`\square`$

For accepted hull weights $`\alpha_i`$, the actual raw trial probabilities are

```math
t_i=\frac{\alpha_i/p_i}{\sum_j\alpha_j/p_j}.
```

Expected calls per accepted fresh trial count failed attempts. They do not describe replaying the same unknown input after a failed run.

## 3. Exact equal-visibility examples

At ideal coherence $`x=z=1`$, the lower hull has vertices $`(e,C)=(0,8),(1/3,3),(1/2,1)`$. Hence

```math
C_{123}(e;1,1)=
\begin{cases}
8-15e,&0\le e\le1/3,\\
7-12e,&1/3\le e\le1/2,\\
1,&e\ge1/2.
\end{cases}
```

At $`x=z=9/10`$, the exact retained hull is:

| Trial protocol | Error | Calls per accepted trial |
|---|---|---|
| Three calls, $`O`$ only | $`3/31`$ | $`360/31`$ |
| Two calls, minus only | $`3/22`$ | $`80/11`$ |
| Three calls, $`O,E,T`$ | $`87/236`$ | $`180/59`$ |
| One call | $`1/2`$ | 1 |

The $`OE`$ prefix and deterministic $`OETA`$ point are dominated. The full cost is

```math
C_{123}(e;9/10,9/10)=
\begin{cases}
\infty,&e\lt3/31,\\
200/9-(2960/27)e,&3/31\le e\le3/22,\\
1960/201-(10960/603)e,&3/22\le e\le87/236,\\
273/31-(484/31)e,&87/236\le e\le1/2,\\
1,&e\ge1/2.
\end{cases}
```

For the fixed external tolerance $`e=3/8`$, choose a one-call trial with raw probability $`1/21`$ and an $`OET`$ three-call trial with probability $`20/21`$. The latter's success is $`59/60`$, so total success is $`62/63`$, mean calls per attempt are $`61/21`$, and cost is $`183/62`$. Accepted weights are $`3/62`$ and $`59/62`$. The optimal one-/two-call cost at the same tolerance is $`101/32`$, giving relative saving

```math
1-\frac{183/62}{101/32}=\frac{203}{3131}
\simeq6.48355\%.
```

The largest relative saving on this visibility slice is $`203/3083\simeq6.58450\%`$, attained at error $`87/236`$.

<a id="unequal-noise"></a>
## 4. Independent logical visibilities

An accuracy extension and a cost saving at a common error are different questions. The following criteria retain that distinction.

### Strict saving at an interior common error

For $`x\gt1/2`$, use $`B=2(6+x)/(2x-1)`$ and define

```math
k_O=\frac{8-2z+B(5z-2)}{24},\qquad
k_E=\frac{8-2z+B(3z-2)}{24},\qquad
k_T=\frac{4+8z+B(z-1)}{24}.
```

There is a strict cost improvement at an interior common error $`e\in(e_2,1/2)`$ exactly when

```math
[k_O]_+ +[k_E]_+ +[k_T]_+\gt3,
\qquad [u]_+=\max\{u,0\}.
```

**Proof.** A three-call native accepted map with totals $`P,F`$ lies below the one-/two-call support line precisely when $`P+B(F-P/2)\gt3`$. The four independent flag contributions are the displayed $`k_i`$ and $`k_A=(1-z)(4-B)/24\le0`$. Maximizing fractional acceptance chooses each positive contribution, giving the criterion. The full native converse extends this support inequality to every allowed controller. If no flag choice violates it, no mixture improves the baseline. If one point does, mixing it with the baseline endpoints gives strict savings at every interior common tolerance. At $`e=1/2`$, the minimum possible call cost one cannot be improved. $`\square`$

For equal visibilities $`x=z=\eta`$, this common-error cost region is re-entrant:

```math
\eta\in
\left(\frac12,\frac{35-\sqrt{985}}6\right)
\ \cup\
\left(\frac{5+\sqrt{233}}{26},1\right].
```

The two finite thresholds are approximately $`0.6025483911749282`$ and $`0.779397597018221`$. The intermediate region still has an accuracy extension.

### Any three-call activation

Call three calls *active* when they either enable a lower error or reduce the cost at an already feasible error. The exact criterion is

```math
\begin{aligned}
x\le\tfrac12:\quad&z\gt\tfrac25,\\
x\gt\tfrac12:\quad&
z\gt\min\left\{\frac{2x}{3-x},\frac{57x+4}{13x+52}\right\}.
\end{aligned}
```

**Proof.** Below error $`1/2`$, every relevant three-call minimum is the $`O`$ flag, with error $`3(1-z)/(4-z)`$. It extends the one-call floor exactly for $`z\gt2/5`$ if $`x\le1/2`$, or the one-/two-call floor exactly for $`z\gt q(x)=2x/(3-x)`$ otherwise.

For $`x\gt1/2`$ suppose instead $`z\le q(x)`$. Every three-call cost segment is affine in error, so it suffices to compare its cumulative vertices with the baseline support line. The exact crossing thresholds are:

| Three-call cumulative flags | Threshold on $`z`$ to lie below the support line |
|---|---|
| $`O`$ | $`(66x-20)/(3x+31)`$ |
| $`O,T`$ | $`(21x-4)/(4x+11)`$ |
| $`O,E`$ | $`(30x-2)/(2x+25)`$ |
| $`O,E,T`$ | $`(57x+4)/(13x+52)`$ |

Only the appropriate second prefix occurs in each flag ordering. The first three thresholds exceed $`q(x)`$ throughout $`1/2\lt x\le1`$, as follows by clearing positive denominators and checking the resulting quadratic signs. A $`T`$-first prefix occurs only at error at least $`1/2`$ and cannot improve this range. The $`OET`$ vertex has

```math
e_{OET}=\frac{15-7z}{4(5+z)},\qquad
c_{OET}=\frac{18}{5+z}.
```

Its threshold $`r(x)=(57x+4)/(13x+52)`$ is at least $`5/9`$. Whenever it matters, the $`O,E,T`$ ordering and error below $`1/2`$ therefore hold. The fully accepted $`OETA`$ point adds no new activation. If its error is below $`1/2`$, the slope of its final sign-acceptance segment is $`-36/(5z)`$, larger than the slope $`-24/(5z-3)`$ from that full point to one call; mixing $`OET`$ with one call dominates it. If its error is at least $`1/2`$, one call directly dominates it. This exhausts the vertices. Taking the union of $`z\gt q(x)`$ and the remaining support crossing $`z\gt r(x)`$ gives the criterion. $`\square`$

The two thresholds cross at the positive root of $`83x^2-63x-12=0`$:

```math
x_*=(63+\sqrt{7953})/166\simeq0.9167445426413625.
```

For ideal two-call coherence $`x=1`$, any three-call activation requires $`z\gt61/65`$.

### A fixed external error on the priority slice

At $`x=9/10`$, accuracy extends exactly for $`z\gt6/7`$, while common-error cost improves exactly for $`z\gt79/91`$. At the fixed tolerance $`e=3/8`$,

```math
C_{123}(3/8;9/10,z)=
\begin{cases}
101/32,&z\le79/91,\\
(789-680z)/(135-83z),&79/91\lt z\lt15/17,\\
(17z+3)/(18z-10),&15/17\le z\le1.
\end{cases}
```

The middle branch mixes two-call minus with three-call $`OET`$; the last mixes one call with $`OET`$. The external tolerance remains $`3/8`$ throughout this comparison.

<a id="duration-model"></a>
## 5. Conditional duration and other attempt costs

Suppose logical order dephasing is generated by $`\kappa(\mathrm{Diag}-\mathrm{id})`$ for duration $`T_m`$, so $`\eta_m=e^{-\kappa T_m}`$. This is an additional scalar logical-noise model; $`\kappa`$ is a rate, not the physical local-qubit visibility $`\gamma`$ used in the four-call analysis.

At $`x=9/10`$, the fixed-error saving above persists exactly when

```math
\kappa(T_3-T_2)\lt\log(819/790)
=0.03605113839200222\ldots,
\qquad T_3/T_2\lt1.3421693427268673\ldots.
```

Call-proportional duration $`T_3/T_2=3/2`$ gives $`z=(9/10)^{3/2}\lt6/7`$, so three calls add no advantage on this particular visibility slice.

For a general $`x`$ with that duration ratio, put $`t=\sqrt x`$, so $`z=t^3`$. Accuracy never extends: for $`0\lt t\lt1`$, $`3t-t^3\lt2`$, equivalently $`(t-1)^2(t+2)\gt0`$, gives $`x^{3/2}\lt2x/(3-x)`$. For $`x\le1/2`$, also $`z\le(1/2)^{3/2}\lt2/5`$. Any activation is therefore precisely the remaining $`OET`$ cost condition

```math
13t^5+52t^3-57t^2-4\gt0.
```

On $`t\in[\sqrt{1/2},1]`$, its derivative $`t(65t^3+156t-114)`$ is positive, and the polynomial changes sign between the endpoints. There is one threshold $`x_c=t_c^2`$. Equivalently, $`x_c`$ is the root in this interval of

```math
169x^5+1352x^4+2704x^3-3249x^2-456x-16=0.
```

Exact rational sign substitution certifies

```math
0.9199910805350899\lt x_c\lt0.9199910805350900.
```

Thus duration ratio $`3/2`$ loses activation at $`x=.9`$ but retains a call-cost gain for $`x\gt x_c`$.

If elapsed time or another resource is charged per attempt, use its explicit weights $`a_m`$ and replace the LP objective by $`\sum_ma_mv_m`$ and vertex costs by $`a_m/P_k`$. Recompute the hull and thresholds for those weights. Local noise on an encoded multiqubit controller generally has pair-dependent coherence factors and needs a different analysis; the scalar duration law and call-cost thresholds do not certify a hardware advantage.

## 6. Verification and attribution boundaries

The [exact invariance audit](../verification/native_benchmark/verify_exact_invariance.py) checks the finite complex-controller kernels and the unequal-noise activation predicate using rational arithmetic. The [closure checker](../verification/native_benchmark/verify_complete_closure.py) reproduces the finite hulls, trial plans, phase calculations and figures. The exact constants are in [exact_cost_and_noise_constants.json](../verification/native_benchmark/exact_cost_and_noise_constants.json). These checks support the displayed finite algebra and examples; the all-controller converse and the convex closure use the proofs on this page and the [native three-call theorem](native_three_call_frontier.md).

Numerical finite-correction-menu SDP comparisons are corroboration only. In particular, an inaccurate ideal-coherence, zero-error numerical solve is not used as a feasibility certificate. The analytic native maps give feasibility at that boundary. Closest-work attribution and the limits of the literature comparison belong to the [related-work account](../references/RELATED_WORK.md); neither this hull calculation nor a search establishes publication priority, channel capacity, or device-independent causal nonseparability.
