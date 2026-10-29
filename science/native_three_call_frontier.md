# The native three-call success–diamond-error frontier

For the rank-deficient qubit channel $`\mathcal M(\rho)=K_0\rho K_0^\dagger+K_1\rho K_1^\dagger`$, with $`K_0=|+\rangle\langle0|`$ and $`K_1=|-\rangle\langle1|`$, four controller flags determine the complete three-call recovery frontier. Fractional acceptance of those flags attains the optimum over arbitrary mixed or ancillary controllers, complex accepted effects, and known final target-unitary corrections.

This is the zero-depolarizing-admixture point of the [common channel family](noisy_routing_model.md#channel-family). It uses all six three-slot orders and **global logical order dephasing**. It is distinct from the [four-call full-rank family with physical local noise](recovery_frontier.md). In particular, the [four-call full-rank minimum](query_threshold.md) does not apply to this rank-deficient channel.

## 1. The task and its output

Use three independent compulsory ordinary-channel slots and the matrix-product words $`012,021,102,120,201,210`$. The rightmost slot acts first. The controller and any finite ancilla begin independently of target and reference. Slot environments are inaccessible; no protected encoding, target-carrying bypass or intermediate target gate is supplied. After routing, unflagged global dephasing $`\mathcal D_\eta=\eta\,\mathrm{id}+(1-\eta)\mathrm{Diag}`$ acts on the six actual order labels, before any final controller processing. Each final instrument outcome may select a known target unitary.

The evaluated output is the corrected accepted **target-only** channel. Fine controller outcomes and correction labels are averaged, while the overall acceptance flag remains available. Require aggregate input-independent acceptance $`\mathcal E=p\Lambda`$, with $`p\gt0`$ and $`\Lambda`$ trace preserving. Individual fine outcomes of an arbitrary three-call strategy may have input-dependent probabilities; the converse does not normalize them separately.

The error metric is the half-diamond distance

```math
e=\frac12\lVert\Lambda-\mathrm{id}\rVert_\diamond.
```

Choi matrices use the normalized Bell input $`|\Phi\rangle=(|00\rangle+|11\rangle)/\sqrt2`$. Write $`f=\langle\Phi|J(\mathcal E)|\Phi\rangle`$ for unnormalized corrected Bell overlap and $`F_\Phi=f/p`$. Set $`U=(X-Z)/\sqrt2`$ and $`V=(X+Z)/\sqrt2`$.

## 2. Four native controller flags

Let $`R`$ reverse each word, $`\boldsymbol1=(1,1,1,1,1,1)^T`$, and $`s=(1,-1,-1,1,1,-1)^T`$. Define

```math
T=\frac{\boldsymbol1\boldsymbol1^T}{6},\qquad
A=\frac{ss^T}{6},\qquad
O=\frac{I-R}{2}-A,\qquad
E=\frac{I+R}{2}-T.
```

They are mutually orthogonal projectors of ranks $`1,1,2,2`$, sum to $`I_6`$, and have constant diagonals $`T_{ii}=A_{ii}=1/6`$ and $`O_{ii}=E_{ii}=1/3`$. Prepare the uniform six-order controller and measure this PVM. Correct by $`I`$ on $`T`$, $`U`$ on $`O`$, $`Y`$ on $`E`$, and $`I`$ on $`A`$.

| Flag | Acceptance $`p_z`$ | Corrected Bell overlap $`f_z`$ | Correction |
|---|---|---|---|
| $`O`$ | $`(4-\eta)/12`$ | $`(1+2\eta)/12`$ | $`U`$ |
| $`E`$ | $`(4-\eta)/12`$ | $`(1+\eta)/12`$ | $`Y`$ |
| $`T`$ | $`(1+2\eta)/6`$ | $`(1+5\eta)/24`$ | $`I`$ |
| $`A`$ | $`(1-\eta)/6`$ | $`(1-\eta)/24`$ | $`I`$ |

The probabilities sum to one. Zero-probability flags are omitted when forming fidelity ratios.

<a id="native-three-call-frontier"></a>
## 3. The complete frontier

**Theorem.** For every $`0\lt p\le1`$ and $`0\le\eta\le1`$, define

```math
f_{3,*}(p,\eta)=
\max\left\{\sum_zf_zx_z:
0\le x_z\le1,\quad\sum_zp_zx_z=p\right\}.
```

The exact optimum over the full controller model in Section 1 is

```math
e_{3,*}(p,\eta)=1-\frac{f_{3,*}(p,\eta)}p.
```

It is attained by the native PVM, accepting flag $`z`$ with probability $`x_z`$. Continuous rotated measurements, complex preparations and effects, controller ancillas, and arbitrary known final target-unitary corrections cannot improve this frontier.

### Universal coefficient budget

The [controller-feasibility proof](controller_optimization.md#controller-feasibility) does not require full Choi rank. Each accepted fine outcome $`y`$ has a pre-noise coefficient matrix $`G_y\succeq0`$ and target operation

```math
\mathcal E_{1,y}(\rho)=
\sum_{\boldsymbol a,\pi,\sigma}(G_y)_{\pi\sigma}
A_{\boldsymbol a,\pi}\rho A_{\boldsymbol a,\sigma}^\dagger.
```

Including the rejected effect gives

```math
\sum_yG_y\preceq D,\qquad
D=\mathrm{diag}(d_\pi),\qquad d_\pi\ge0,\qquad
\mathrm{Tr}D=1.
```

This includes arbitrary mixed and ancillary controllers without assuming uniform preparation. For a positive $`G`$, put $`t=\mathrm{Tr}(TG)`$, $`a=\mathrm{Tr}(AG)`$, $`o=\mathrm{Tr}(OG)`$, and $`v=\mathrm{Tr}(EG)`$. They are nonnegative and sum to $`q=\mathrm{Tr}G`$. The symbol $`v`$ here is the even-standard-sector mass, distinct from the error $`e`$ and the unitary $`V`$.

### Exact Bell kernels

In the orthonormal unitary Bell basis associated with $`I,U,Y,V`$, direct Kraus expansion gives

```math
J(\mathcal E_1(G))=
\begin{pmatrix}
3t/2&0&0&0\\
0&3o/4&0&\chi\\
0&0&(t+v)/2&0\\
0&\chi^*&0&t+v/4
\end{pmatrix}.
```

The diagonal Bell kernels are

```math
M_I=3T/2,\qquad M_U=3O/4,\qquad
M_Y=(T+E)/2,\qquad M_V=T+E/4.
```

Every cross Bell kernel except $`U,V`$ is zero. The remaining coherence comes only from the $`O,E`$ controller block, and can be specified explicitly as

```math
\chi=\frac18\sum_{\pi,\sigma}L_{\pi\sigma}G_{\pi\sigma},\qquad
L=\begin{pmatrix}
0&1&-1&1&-1&0\\
1&0&-1&0&-1&1\\
-1&1&0&1&0&-1\\
-1&0&1&0&1&-1\\
1&-1&0&-1&0&1\\
0&-1&1&-1&1&0
\end{pmatrix},\qquad OLE=L.
```

These finite identities cover complex $`G`$ and complex $`\chi`$. They can be reproduced without floating-point arithmetic: use integer Kraus matrices $`\sqrt2K_0`$ and $`\sqrt2K_1`$, form all eight environmental histories and six word products, and contract their vectorizations. The diagonal kernels have denominators $`32,64,32,64`$ for $`I,U,Y,V`$; the unscaled integer $`U,V`$ kernel is $`8L`$ with denominator $`64`$. Each target-coordinate kernel commutes with the six left slot-relabeling actions. Its blocks between the trivial, sign and standard isotypic sectors vanish, and its sign-sector rows and columns vanish. The exact [invariance verifier](../verification/native_benchmark/verify_exact_invariance.py) checks these identities, including the sole cross-kernel support.

### All-unitary upper bound

Every definite order gives $`\mathcal M^3=\Delta`$. Global dephasing therefore adds an isotropic Choi term:

```math
J(\mathcal E_\eta(G))=
\eta J(\mathcal E_1(G))+(1-\eta)qI_4/4.
```

Decompose the controller into its positive trivial, sign and standard compressions. On the standard sector, the positive $`U,V`$ block has trace $`3o/4+v/4`$, and the separate $`Y`$ entry is $`v/2`$. Its largest noisy Choi eigenvalue satisfies

```math
\lambda_{\max}(J_{\mathrm{std},\eta})
\le\frac{1-\eta}{4}(o+v)
+\eta\max\{v/2,3o/4+v/4\}
\le\frac{1+2\eta}{4}o+\frac{1+\eta}{4}v.
```

The positive block's trace bounds its largest eigenvalue for every complex coherence. This is the step that rules out improvement from continuously rotated $`U,V`$ measurements. The trivial sector has largest eigenvalue $`(1+5\eta)t/4`$; the sign sector contributes $`(1-\eta)a/4`$. Final target unitaries preserve Choi eigenvalues. Subadditivity of the largest eigenvalue on positive sums gives

```math
f\le\sum_y\left[
\frac{1+5\eta}{4}t_y+
\frac{1-\eta}{4}a_y+
\frac{1+2\eta}{4}o_y+
\frac{1+\eta}{4}v_y\right].
```

The Bell-input accepted trace is exactly

```math
p=\sum_y\left[(1+2\eta)t_y+(1-\eta)a_y+
\frac{4-\eta}{4}(o_y+v_y)\right].
```

It equals the input-independent total acceptance assumed in the task, even when the fine probabilities depend on the input. Constant projector diagonals and the diagonal resource budget give

```math
\sum_yt_y\le1/6,\qquad\sum_ya_y\le1/6,
\qquad\sum_yo_y\le1/3,\qquad\sum_yv_y\le1/3.
```

Rescaling these four masses gives the flag variables $`x_T,x_A,x_O,x_E\in[0,1]`$ and exactly the theorem's linear-program Bell bound. Testing a normalized Bell input against its ideal output yields $`\frac12\lVert\Lambda-\mathrm{id}\rVert_\diamond\ge1-F_\Phi`$. This proves the full-controller converse.

### Native attainment and exact diamond distance

Before noise and corrections, the uniform-controller native operations are

```math
\mathcal E_T=\tfrac14\mathrm{id}+\tfrac1{12}\mathrm{Ad}_Y+\tfrac16\mathrm{Ad}_V,
\qquad \mathcal E_O=\tfrac14\mathrm{Ad}_U,
```

```math
\mathcal E_E=\tfrac16\mathrm{Ad}_Y+\tfrac1{12}\mathrm{Ad}_V,
\qquad \mathcal E_A=0.
```

For projector $`P`$, global noise replaces this by $`\eta\mathcal E_P+(1-\eta)\mathrm{rank}(P)\Delta/6`$. The native corrections attain the Bell bound separately on every flag. Each corrected native operation is a positive mixture of the common unitary channels $`I,U,Y,V`$, so every accepted native mixture is $`p`$ times a trace-preserving channel. Its normalized identity weight is $`w_0=f/p`$.

The three other unitary Bell states are orthogonal to $`\Phi`$. A Bell input therefore gives half trace distance $`1-w_0`$, while the channel triangle inequality gives

```math
\lVert\Lambda-\mathrm{id}\rVert_\diamond
\le(1-w_0)+\sum_{j\ne0}w_j=2(1-w_0).
```

The half-diamond error is exactly $`1-f/p`$, proving attainment without target twirling or an assumed covariance reduction.

## 4. Closed forms and boundary cases

The decreasing flag-fidelity order is:

| Visibility interval | Order |
|---|---|
| $`0\lt\eta\lt1/7`$ | $`T,O,E,A`$ |
| $`1/7\lt\eta\lt7/13`$ | $`O,T,E,A`$ |
| $`7/13\lt\eta\le1`$ | $`O,E,T,A`$ |

Merge tied ratios at the boundaries. At $`\eta=0`$ all ratios are $`1/4`$; at $`\eta=1`$ omit the zero-probability sign flag. For sorted nonzero flags $`z_1,\ldots,z_k`$, set

```math
P_j=\sum_{i=1}^jp_{z_i},\qquad F_j=\sum_{i=1}^jf_{z_i},
\qquad P_0=F_0=0,\qquad r_j=f_{z_j}/p_{z_j}.
```

For $`P_{j-1}\le p\le P_j`$,

```math
f_{3,*}(p)=F_{j-1}+r_j(p-P_{j-1}),
\qquad
e_{3,*}(p)=1-r_j-\frac{F_{j-1}-r_jP_{j-1}}p.
```

At error budget $`e`$, maximum success is the largest $`p`$ satisfying $`f_{3,*}(p)\ge(1-e)p`$. It is zero below $`1-r_1`$, and one for $`e\ge(9-5\eta)/12`$. On the crossed segment between those cases,

```math
p_{3,*}(e,\eta)=\frac{F_{j-1}-r_jP_{j-1}}{1-e-r_j},
\qquad C_{3,*}=3/p_{3,*}.
```

At a plateau budget choose its largest feasible endpoint. At $`\eta=0`$, success is zero for $`e\lt3/4`$ and one otherwise. The point $`p=0`$ has no conditional channel.

At ideal coherence,

```math
e_{3,*}(p,1)=
\begin{cases}
0,&0\lt p\le1/4,\\
1/3-1/(12p),&1/4\lt p\le1/2,\\
1/2-1/(6p),&1/2\lt p\le1,
\end{cases}
```

```math
C_{3,*}(e,1)=
\begin{cases}
12-36e,&0\le e\le1/6,\\
9-18e,&1/6\le e\le1/3,\\
3,&e\ge1/3.
\end{cases}
```

These costs count compulsory noisy-slot calls per accepted fresh trial. Six coherent order labels, routing implementation, final controller processing and correction costs are separate resources. The [randomized call-cost closure](randomized_call_cost.md) charges unsuccessful attempts and compares classical mixtures of one-, two- and three-call trials at a fixed external error tolerance. Equal logical visibility at different call counts is a model choice, not a consequence of equal elapsed time.
