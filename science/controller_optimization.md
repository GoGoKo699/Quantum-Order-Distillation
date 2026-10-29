# All-controller reduction for noisy recovery

An accepted controller outcome weights interference between the occupied channel orders. Its coefficient matrix gives an exact description of the allowed preparation and final instrument, including mixed controllers and arbitrary finite ancillary controllers. This reduction supplies the upper bounds behind the [orthogonal](recovery_frontier.md) and [tilted](tilted_readout.md) recovery frontiers.

Throughout this page, use the [raw-slot access contract](noisy_routing_model.md#shared-access). Write $r=|\Pi|$ and let $A_{\boldsymbol a,\pi}$ denote its word Kraus products. Define, for a coefficient matrix $H$,

$$
\mathcal E_H(X)
=\sum_{\pi,\sigma}H_{\pi\sigma}
\sum_{\boldsymbol a}A_{\boldsymbol a,\pi}X
A_{\boldsymbol a,\sigma}^\dagger.
\tag{1}
$$

If $H\succeq0$, this map is completely positive. Indeed, decompose $H=\sum_j c_jc_j^\dagger$; its Kraus operators are $\sum_\pi(c_j)_\pi A_{\boldsymbol a,\pi}$.

<a id="controller-feasibility"></a>
## 1. Exact feasible set for every accepted outcome

**Proposition.** The effective matrices of an arbitrary allowed controller preparation and accepted instrument are exactly the families

$$
H_x\succeq0,
\qquad \sum_xH_x\preceq\operatorname{diag}(t),
\qquad t_\pi\ge0,
\qquad \sum_\pi t_\pi=1.
\tag{2}
$$

Mixed preparation and an ancillary controller do not enlarge this feasible set. The unlisted final effect is rejection.

**Proof.** Let $F_x$ be the accepted effects on $CA$, with $F_x\succeq0$ and $\sum_xF_x\preceq I$. For routing-register blocks $\omega_{\pi\sigma}=\langle\pi|\omega|\sigma\rangle$ and similarly for $F_x$, contraction of the controller gives

$$
(H_x)_{\pi\sigma}
=\operatorname{Tr}_A\bigl[(F_x)_{\sigma\pi}\omega_{\pi\sigma}\bigr].
\tag{3}
$$

Purify $\omega$ on another register $B$ as $|\psi\rangle=\sum_\pi|\pi\rangle|a_\pi\rangle_{AB}$. Define $V|\pi\rangle=|\pi\rangle|a_\pi\rangle$. Its columns are orthogonal and

$$
V^\dagger V=\operatorname{diag}(t),
\qquad t_\pi=\langle a_\pi|a_\pi\rangle,
\qquad
H_x^T=V^\dagger(F_x\otimes I_B)V.
$$

Compression and transposition give (2).

Conversely, set $D=\operatorname{diag}(\sqrt t)$ and prepare $|\alpha\rangle=\sum_\pi\sqrt{t_\pi}|\pi\rangle$. On the support of $t$, use effects

$$
F_x=D^{-1}H_x^T D^{-1}.
$$

They are positive and sum to at most the identity. A zero entry of $t$ forces the corresponding row and column of every positive $H_x$ to vanish, so the inverse is needed only on its support. The resulting coefficients are $DF_x^TD=H_x$. This realizes every family in (2). $\square$

A controller Schur channel sends these physical coefficients to

$$
\widehat H_x=\Gamma\circ H_x,
\tag{4}
$$

where $\circ$ denotes entrywise multiplication. Constraint (2) applies to $H_x$ **before** this noise, not to an independently chosen noisy coefficient matrix.

<a id="uniform-preparation"></a>
## 2. When uniform preparation loses nothing

For identical independent slot channels, a common relabeling of their slots leaves every cross map unchanged:

$$
\sum_{\boldsymbol a}A_{\boldsymbol a,g\pi}X
A_{\boldsymbol a,g\sigma}^\dagger
=\sum_{\boldsymbol a}A_{\boldsymbol a,\pi}X
A_{\boldsymbol a,\sigma}^\dagger.
\tag{5}
$$

This follows by relabeling the independent Kraus indices in the sum. Suppose a group of these relabelings acts transitively on $\Pi$ and the fixed noise also obeys $\Gamma_{g\pi,g\sigma}=\Gamma_{\pi\sigma}$. Average each $H_x$ and $t$ over that group. The entire noisy outcome map is preserved, while the averaged diagonal budget becomes uniform. Therefore every outcome-map family has a representative satisfying

$$
H_x\succeq0,
\qquad \sum_xH_x\preceq I_r/r.
\tag{6}
$$

Conversely, uniform pure preparation and effects $F_x=rH_x^T$ realize every family in (6). This is a statement about complete maps, stronger than preservation of one objective.

For the [physical eight-order code](noisy_routing_model.md#physical-codebook), put $a=(12)$, $b=(34)$ and $s_0=(13)(24)$. The occupied orders form

$$
\langle a,b,s_0\rangle
\cong(C_2\times C_2)\rtimes C_2\cong D_8,
\qquad \pi(c,x,y)=a^xb^ys_0^c.
$$

Left multiplication by $g=\pi(c_g,x_g,y_g)$ acts on the physical code by

$$
c\longmapsto c_g\oplus c,
\qquad
(x,y)\longmapsto(x_g,y_g)\oplus
\operatorname{swap}^{c_g}(x,y).
$$

This is a bit permutation followed by bit flips, hence a Hamming isometry. Both global dephasing and equal independent local dephasing are invariant. Constraint (6), with $r=8$, is consequently exact for the stated architecture.

Uniform preparation cannot be inferred from identical slots alone for an arbitrary encoding, unequal local-noise strengths, or an arbitrary fixed $\Gamma$. Without the required noise symmetry, retain (2).

<a id="bell-gram-conventions"></a>
## 3. Bell-input objective matrices and the transpose convention

Use $|\Phi_d\rangle=d^{-1/2}\sum_i|i\rangle_R|i\rangle_T$ and define

$$
|v_{\boldsymbol a,\pi}\rangle
=(I_R\otimes A_{\boldsymbol a,\pi})|\Phi_d\rangle,
\qquad
f^U_{\boldsymbol a,\pi}
=\operatorname{Tr}(UA_{\boldsymbol a,\pi})/d.
$$

The conventional Gram matrices are

$$
B_{\pi\sigma}
=\sum_{\boldsymbol a}\langle v_{\boldsymbol a,\pi}|v_{\boldsymbol a,\sigma}\rangle
=\frac1d\sum_{\boldsymbol a}
\operatorname{Tr}(A_{\boldsymbol a,\pi}^\dagger A_{\boldsymbol a,\sigma}),
$$

$$
A^U_{\pi\sigma}
=\sum_{\boldsymbol a}
\overline{f^U_{\boldsymbol a,\pi}}f^U_{\boldsymbol a,\sigma}.
$$

With a final correction $U_x$ on outcome $x$, acceptance $p$ and unnormalized corrected Bell overlap $u$ are

$$
p=\sum_x\operatorname{Tr}(B_\Gamma H_x),
\qquad
u=\sum_x\operatorname{Tr}(A^{U_x}_\Gamma H_x),
\tag{7}
$$

$$
B_\Gamma=B\circ\Gamma^T,
\qquad A^U_\Gamma=A^U\circ\Gamma^T.
\tag{8}
$$

The transpose in (8) matters for complex correlation matrices: the physical coefficients in (1) are $H_{\pi\sigma}$, whereas the trace objective contracts the opposite matrix indices. All named noise masks here are real symmetric, so their transposes can be omitted after this convention is fixed. For $p>0$, $F=u/p$.

<a id="real-controller-reduction"></a>
## 4. Correction sectors and constant-success representatives

Assume the qubit slot channel has real Kraus matrices and is covariant under conjugation by Pauli $Y$. The [measure–prepare family](noisy_routing_model.md#channel-family) has these properties: conjugation by $Y$ exchanges its two pure Kraus operators up to signs, and its depolarizing Kraus family can be chosen real.

**Proposition.** For these channels and real symmetric $\Gamma$, every achievable Bell pair $(p,u)$ has a realization with real symmetric $H_x$ and corrections on the two circles

$$
U_\parallel(\varphi)
=\cos\varphi\,I+\sin\varphi\,iY,
\qquad
U_\perp(\varphi)
=i(\cos\varphi\,X+\sin\varphi\,Z).
\tag{9}
$$

Each nonzero accepted corrected map in this realization is a scalar multiple of a trace-preserving unital channel.

**Proof.** Conjugating all slot Kraus operators by $Y$ gives another representation of the same channel. Its induced change of environmental indices is the same in every occupied order. Thus the mixed Gram

$$
T(U,V)_{\pi\sigma}
=\frac1{d^2}\sum_{\boldsymbol a}
\overline{\operatorname{Tr}(UA_{\boldsymbol a,\pi})}
\operatorname{Tr}(VA_{\boldsymbol a,\sigma})
$$

is invariant under $(U,V)\mapsto(YUY,YVY)$. It vanishes when one correction commutes with $Y$ and the other anticommutes.

Ignoring a global phase, every qubit unitary is

$$
U=aI+b(iY)+i(cX+dZ),
\qquad a,b,c,d\in\mathbb R,
\qquad a^2+b^2+c^2+d^2=1.
$$

Set $w=a^2+b^2$. The vanishing mixed Grams imply

$$
A^U_\Gamma
=wA^{U_\parallel}_\Gamma
+(1-w)A^{U_\perp}_\Gamma,
$$

using normalized sector unitaries and omitting a zero-weight sector. Split an outcome $H_x$ into $wH_x$ and $(1-w)H_x$ with these corrections. This preserves its $p,u$ contributions and its total coefficient budget.

On the parallel circle the trace amplitudes are real; on the perpendicular circle they are purely imaginary. Its Gram matrices and $B_\Gamma$ are therefore real symmetric. Replacing each split matrix by

$$
\operatorname{Re}H_x=(H_x+\overline H_x)/2
$$

preserves positivity, feasibility and both objective values.

For real $H_x$, the uncorrected noisy map preserves real matrices. The same common Kraus change proves $Y$ covariance of each cross map. Hence both its success effect $\mathcal E^\dagger(I)$ and its output $\mathcal E(I)$ are real symmetric matrices commuting with $Y$. A qubit matrix with both properties is scalar. Equality of their traces identifies the same scalar $p_x$:

$$
\mathcal E^\dagger(I)=p_xI,
\qquad \mathcal E(I)=p_xI.
$$

A final target unitary preserves these properties. The accepted map is $p_x$ times a trace-preserving unital channel for $p_x>0$ and is zero for $p_x=0$. $\square$

The unrestricted Bell-input optimization is thus tight even when recovery is required to have input-independent acceptance. This is an existence statement for representatives with the same objectives. It does **not** say that every permitted complex controller has constant success, or that realification preserves every entanglement property of its accepted state.

<a id="controller-primal-dual"></a>
## 5. Exact optimization and a finite dual certificate

Under the transitive symmetry, the complete frontier at prescribed acceptance $p$ is

$$
\sup_{\{H_x,U_x\}}
\left\{
\sum_x\operatorname{Tr}(A^{U_x}_\Gamma H_x):
H_x\succeq0,\quad\sum_xH_x\preceq I_r/r,\quad
\sum_x\operatorname{Tr}(B_\Gamma H_x)=p
\right\}.
\tag{10}
$$

The two correction circles suffice. General outcome spaces add no objective pairs: finite-dimensional convex decomposition preserves the matrix sum $\sum_xH_x$ and the two real numbers $p,u$ with finitely many outcome types.

A dual description is

$$
\inf_{\lambda\in\mathbb R,\,Q\succeq0}
\left\{
\lambda p+\frac{\operatorname{Tr}Q}{r}:
Q+\lambda B_\Gamma\succeq A^U_\Gamma
\text{ for every }U\text{ on both circles}
\right\}.
\tag{11}
$$

Every feasible dual pair bounds (10): multiply each matrix inequality by its positive $H_x$, take traces and sum, then use the budget $\sum_xH_x\preceq I_r/r$. Finite-dimensional convex separation of the compact attainable set proves equality, allowing an infimum at a boundary value of $p$.

On either circle,

$$
A^U_\Gamma=C_0+C_c\cos(2\varphi)+C_s\sin(2\varphi).
$$

The first-degree matrix Fejér–Riesz factorization turns the universal angle constraint into a finite semidefinite constraint. Put

$$
Q_0=Q+\lambda B_\Gamma-C_0,
\qquad Q_1=-(C_c-iC_s)/2.
$$

For each circle, require a Hermitian $Z$ such that

$$
\begin{pmatrix}
Z&Q_1\\Q_1^\dagger&Q_0-Z
\end{pmatrix}\succeq0.
\tag{12}
$$

This uses the matrix factorization ingredient described by [Chu (2016)](https://doi.org/10.1016/j.laa.2015.09.051). The [frontier](recovery_frontier.md#frontier-theorem) and [tilted](tilted_readout.md#tilted-optimality) proofs below supply explicit operator majorants rather than relying on numerical optimization.

Without the relabeling symmetry, the primal retains (2). Its dual replaces $\operatorname{Tr}Q/r$ by a scalar $\nu$ with additional constraints $Q_{\pi\pi}\le\nu$ for every routing label. Neither the uniform preparation nor the equal-noise frontier is asserted beyond its symmetry assumptions.
