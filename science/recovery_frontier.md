# The complete orthogonal-family recovery frontier

For the [raw measure–prepare family](noisy_routing_model.md#channel-family) at $\theta=\pi/2$, a simple controller readout attains the optimum over every allowed controller preparation, ancillary controller, accepted effect, and final target unitary. The result is the complete tradeoff between acceptance probability and conditional Bell fidelity on the specified physical eight-order architecture.

Write $s=1-\epsilon$, use the [physical codebook](noisy_routing_model.md#physical-codebook), and define

$$
U_-=(X-Z)/\sqrt2,
\qquad U_+=(X+Z)/\sqrt2.
$$

Global phases of target corrections do not affect their channels. The [controller theorem](controller_optimization.md) fixes the Gram conventions, the exact coefficient budget, and the input-independent attaining representatives used below.

<a id="finite-coefficient-identities"></a>
## 1. Exact coefficient identities

Let

$$
V_{\pi k}=2^{-3/2}(-1)^{b_\pi\cdot k},
\qquad k\in\{0,1\}^3,
$$

be the physical Walsh basis, where $b_\pi$ is the actual routing code. Measuring the three physical controller qubits individually in the $X$ basis gives these outcomes.

**Lemma.** At $\theta=\pi/2$, both $V^TBV$ and $V^TA^UV$ are diagonal, coefficient by coefficient in $s$, for

$$
U\in\{I,iY,iU_+,iU_-\}.
$$

The sine harmonic on the $I/iY$ correction circle vanishes; the cosine harmonic on the $iX/iZ$ circle vanishes. The same statements hold after global dephasing and equal independent local dephasing.

**Finite exact proof.** Use six integer matrices

$$
L_0=\begin{pmatrix}1&0\\1&0\end{pmatrix},
\qquad
L_1=\begin{pmatrix}0&1\\0&-1\end{pmatrix},
\qquad
(L_2,L_3,L_4,L_5)=(E_{00},E_{01},E_{10},E_{11}).
$$

The physical Kraus operators are $\sqrt{s/2}\,L_0$, $\sqrt{s/2}\,L_1$, and $\sqrt{(1-s)/2}\,L_j$ for $j\ge2$. The sign choice for $L_1$ differs from a natural pure-state Kraus choice only by a global Kraus sign, which cancels in the channel and its cross maps.

For all $6^4=1296$ tuples $\boldsymbol a$, form the eight integer products $P_{\boldsymbol a,\pi}$. Let $k(\boldsymbol a)$ count tuple labels below 2. Then

$$
B_{\pi\sigma}(s)
=\frac1{32}\sum_{\boldsymbol a}
s^k(1-s)^{4-k}
\operatorname{Tr}(P_{\boldsymbol a,\pi}^TP_{\boldsymbol a,\sigma}),
\tag{1}
$$

$$
A^U_{\pi\sigma}(s)
=\frac1{64}\sum_{\boldsymbol a}
s^k(1-s)^{4-k}
\overline{\operatorname{Tr}(UP_{\boldsymbol a,\pi})}
\operatorname{Tr}(UP_{\boldsymbol a,\sigma}).
\tag{2}
$$

The squared scale of four Kraus factors supplies $1/16$; normalized Choi trace and Bell overlap supply the additional factors $1/2$ and $1/4$. Thus these are degree-four polynomial identities with fixed, explicit normalizations.

Expand the finite sums into integer coefficients and multiply each coefficient matrix by the unnormalized Walsh matrix on both sides. Every off-diagonal coefficient at the four correction endpoints vanishes, as do the stated circle harmonics. The [independent verifier](../verification/audit_verify.py) reconstructs all 1296 tuples and complete coefficient arrays, without importing a producer script. The arrays are retained in [family8_polynomials.npz](../verification/certificates/family8_polynomials.npz).

For global noise, every diagonal order branch is $\mathcal N^4=\Delta$. Consequently

$$
B_\eta=\eta B+(1-\eta)I,
\qquad
A^U_\eta=\eta A^U+(1-\eta)I/4.
\tag{3}
$$

For local noise, $\gamma^{d_H(b_\pi,b_\sigma)}$ is XOR-circulant in the physical codes. Walsh-diagonal matrices form exactly this circulant algebra, which is closed under entrywise products. The mask therefore preserves endpoint diagonality and the vanishing harmonics. Equivalently, split every degree-four coefficient by Hamming distance $h=0,1,2,3$ and verify the identities before multiplying by $\gamma^h$. All 100 masked success/endpoint coefficient matrices are checked exactly, so no interpolation from sampled $s,\gamma$ values is used. $\square$

<a id="frontier-theorem"></a>
## 2. Attaining the full frontier by physical $X$ readout

Let $b_j$ be the eigenvalue of the noisy acceptance Gram in physical Walsh mode $j$, and let

$$
a_j=\max_U(A^U_\Gamma)_{jj}
$$

in that basis. The four endpoint corrections suffice for this maximum.

**Theorem.** For every $0\le\epsilon\le1$ and every $0\le\eta,\gamma\le1$, the exact input-independent recovery upper boundary on the specified eight orders is

$$
u_*(p)=pF_*(p)
=\max\left\{
\sum_j a_jx_j:
0\le x_j\le\frac18,
\quad\sum_jb_jx_j=p
\right\}.
\tag{4}
$$

Uniform preparation, individual controller-qubit $X$ measurements, each outcome's best final target correction, and classical partial acceptance attain it. The optimization includes arbitrary permitted ancillas, mixed preparations and collective controller effects.

**Proof.** The [two correction circles](controller_optimization.md#real-controller-reduction) and the zero-harmonic identities put every unitary objective matrix in the convex hull of the four diagonal endpoint matrices. Therefore

$$
D=V\operatorname{diag}(a_j)V^T\succeq A^U_\Gamma
\qquad\text{for every final }U.
$$

For an arbitrary corrected instrument, set $H_{\mathrm{tot}}=\sum_xH_x$. The uniform coefficient budget gives

$$
u\le\operatorname{Tr}(DH_{\mathrm{tot}}),
\qquad
p=\operatorname{Tr}(B_\Gamma H_{\mathrm{tot}}),
\qquad
0\preceq H_{\mathrm{tot}}\preceq I/8.
$$

Only the Walsh-diagonal entries contribute to the two objectives; call them $x_j$. This is the upper bound (4). Conversely, each allocation $x_j$ is realized by physical $X$ outcome $j$, its best endpoint correction, and acceptance probability $8x_j$. These real branches have input-independent success by the controller proposition. The Bell-input relaxation is therefore attained in the required recovery class. $\square$

To evaluate (4), omit zero-capacity modes, sort $a_j/b_j$ from largest to smallest, and accept modes in that order. A mode has probability capacity $b_j/8$ and reward capacity $a_j/8$. Accept the last needed mode only partially. The sort depends on the parameters; no single order is asserted for the entire parameter square. Fidelity is defined for $p>0$, and the boundary at $p=0$ is understood through $u_*(0)=0$.

<a id="global-mode-table"></a>
## 3. Global-noise mode table and deterministic optimum

The following table gives **ideal** mode data. For a global visibility $\eta$, replace each entry by

$$
b\longmapsto1-\eta+\eta b,
\qquad
a\longmapsto(1-\eta)/4+\eta a
$$

before applying (4). Multiplicities count independent physical outcomes. The mode labels are decimal values of the Walsh characters.

| Class | Physical Walsh modes | Multiplicity | Ideal acceptance $b$ | Best ideal overlap $a$ | Correction |
|---|---|---:|---|---|---|
| $S$ | 0 | 1 | $(19+8s^2+s^4)/8$ | $(7+4s^2+s^4)/8$ | $I$ |
| $P$ | 3 | 1 | $(3+s^4)/8$ | $(3+s^4)/8$ | $I$ |
| $Q$ | 4 | 1 | $(3+s^2)/4$ | $(1+s)(1+s^2)/4$ | $U_-$ |
| $R$ | 7 | 1 | $3(1-s^2)/4$ | $(1+s)(1-s^2)/4$ | $U_+$ |
| $E_-$ | 5, 6 | 2 | $(3+s^2)/4$ | $(1+s^2)/4$ | $iY$ |
| $E_+$ | 1, 2 | 2 | $(9-4s^2-s^4)/8$ | $\max\{(3-2s^2-s^4)/8,(1+s)/4\}$ | $I$ or $U_-$ |

These entries follow from (1)–(2). The verifier checks their full polynomial identities and confirms that the omitted endpoints are dominated throughout $s\in[0,1]$, using exact nonnegative Bernstein coefficients.

Accepting all outcomes gives

$$
F_{\mathrm{det}}
=\frac14+\frac\eta{32}
\max\{4+2s+2s^2,\;3+4s+4s^2+s^4\}.
\tag{5}
$$

The reward switch occurs at the unique root in $[0,1]$ of $s^3-s^2+3s-1=0$, approximately $0.36110308$. It need not coincide with a change in the highest-fidelity herald.

<a id="physical-local-point"></a>
## 4. The physical local-noise frontier at one exact point

For $\epsilon=1/100$ and $\gamma=4/5$, the physical local-noise modes sort as follows. Capacities are per individual outcome; displayed fidelities are decimal evaluations of exact rational rewards and capacities.

| Walsh mode | Multiplicity | Probability capacity per mode | Best conditional fidelity | Correction |
|---|---:|---:|---:|---|
| 3 | 1 | 0.06196376409 | 0.758264362729 | $I$ |
| 4 | 1 | 0.14724419601 | 0.618616111659 | $U_-$ |
| 1, 2 | 2 | 0.09903123591 | 0.475904411037 | $U_-$ |
| 5, 6 | 2 | 0.11176080399 | 0.419928976211 | $iY$ |
| 0 | 1 | 0.33997376409 | 0.400293871071 | $I$ |
| 7 | 1 | 0.02923419601 | 0.376694744615 | $U_+$ |

In particular,

$$
F_{\max}=\frac{1566163803}{2065458803},
\qquad
p_{\max\,\mathrm{at}\,F_{\max}}
=\frac{6196376409}{10^{11}},
\qquad
F_{\mathrm{det}}=\frac{23664816609}{50000000000}.
\tag{6}
$$

At $p=0.1245025$, the frontier has $F_*(p)\simeq0.6881177785722$: accept mode 3 fully and part of mode 4. The exact local data are retained in [local_physical_frontier.json](../verification/certificates/local_physical_frontier.json), with redundant local polynomial arrays in [family8_local_polynomials.npz](../verification/certificates/family8_local_polynomials.npz).

The deterministic optimum in (6) is below $1/2$. This bounds the best corrected Bell fidelity; it does not imply that every possible deterministic output is separable. Realification preserves $(p,u)$, not all entanglement properties of arbitrary complex-controller outputs.

<a id="preferred-herald"></a>
## 5. The inherited exact branch and a change of preferred herald

In displayed-word order, the signs

$$
(+,-,-,+,+,-,-,+)
$$

give Walsh mode 3 and realize the central-polynomial branch

$$
P_2/8,
\qquad P_2=\{[X_1,X_2],[X_3,X_4]\}.
$$

Its ideal branch is proportional to the identity for arbitrary qubit slots; on this family its probability is $(3+s^4)/64$. This polynomial construction is inherited from the central-polynomial/resetting and nested-switch literature, including [Navascués (2018)](https://doi.org/10.1103/PhysRevX.8.031008) and [Kechrimparis, Moran, and Kwon (2025), Section IV](https://doi.org/10.1103/mvrg-yj3l). The [qubit proof](qubit.md) explains its algebraic origin and normalization.

Under global control noise, this branch becomes

$$
\mathcal E_P
=\frac{\eta(3+s^4)}{64}\,\mathrm{id}
+\frac{1-\eta}{8}\Delta.
\tag{7}
$$

The alternative signs $(-,-,-,-,+,+,+,+)$ give mode 4. After correction $U_-$, its global-noise performance is

$$
p_Q=\frac{4-\eta(1-s^2)}{32},
\qquad
F_Q=\frac{1+\eta(s+s^2+s^3)}{4-\eta(1-s^2)}.
\tag{8}
$$

For example, at $\epsilon=.01,\eta=.8$, the $Q$ branch has $p_Q=.1245025$ and $F_Q\simeq.84142868617$, while $P$ has $p_P=.074507450125$ and $F_P\simeq.74834731334$. At the physical local point with $\gamma=.8$, the table instead gives higher fidelity to $P$ and higher probability to $Q$.

The distinction survives matching scalar controller-coherence measures. The noisy uniform-controller fidelity is $(1+\gamma)^3/8$ locally and $(1+7\eta)/8$ globally. Matching it, and the normalized $l_1$ coherence, requires

$$
\eta_{\mathrm{match}}
=\frac{(1+\gamma)^3-1}{7}.
$$

At $\gamma=4/5$, this is $604/875$, and both controller fidelities equal $.729$. For $\epsilon=.01$, matched global noise gives $F_P\simeq.643437$ and $F_Q\simeq.760039$, whereas local noise gives $F_P\simeq.758264$ and $F_Q\simeq.618616$. Thus a scalar coherence amount does not determine the best recovery herald. The spatial pattern of surviving order coherences matters.

## 6. Scope and subsequent connections

The theorem optimizes over controllers and final target unitaries on the stated eight-order raw-slot router. It does not optimize over all 24 orders, arbitrary per-slot target circuits, unequal noise strengths, or different physical encodings. The [strict Schur-noise theorem](noisy_recovery_obstruction.md#strict-schur-obstruction) explains why full-rank slots with strict dephasing cannot retain an exact nonzero unitary branch.

The [tilted theorem](tilted_readout.md) shows that physical $X$ readout remains optimal on a certified interval even when simultaneous diagonalization fails. The [sister entanglement comparison](https://github.com/GoGoKo699/Quantum-Order-Entanglement/blob/main/science/noisy_separation.md) uses a fixed real $P_2$ outcome as a sufficient entanglement witness; it does not infer retained controller–target separability from this recovery frontier.
