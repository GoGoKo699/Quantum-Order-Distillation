# Strict Schur noise obstructs exact recovery

Ideal order interference can produce an exact central-polynomial identity branch even from full-rank noisy slots. A strictly positive Schur-noise floor leaves an unavoidable positive diagonal contribution in every nonzero accepted map. That contribution has full Choi rank, so the accepted map cannot be an exact unitary channel.

The argument applies to any finite number of independent slots and any occupied order support. It uses the [raw-slot access model](noisy_routing_model.md#shared-access) and the exact [controller coefficient representation](controller_optimization.md#controller-feasibility), without requiring uniform preparation or the eight-order symmetry.

<a id="strict-schur-obstruction"></a>
## 1. The obstruction theorem

**Theorem.** Let the target dimension be $`d\ge2`$. Suppose each independent supplied slot has full Choi rank, and the routing-noise correlation matrix satisfies

```math
\Gamma\succeq\lambda I_r
\qquad\text{for some }\lambda>0.
```

Every nonzero accepted target map has full Choi rank. Exact nonzero heralded identity or unitary recovery is impossible at every finite call count, including arbitrary allowed controller ancillas, mixed preparations, accepted effects, and herald-dependent final target unitaries.

**Proof.** Let $`H\succeq0`$ be an accepted coefficient matrix before controller noise. By the Schur product theorem,

```math
\Gamma\circ H\succeq\lambda\mathrm{diag}(H).
\qquad\text{(1)}
```

The map from a positive coefficient matrix to its accepted Choi matrix is positive: a positive coefficient decomposition gives a Kraus decomposition, as in [controller optimization, equation (1)](controller_optimization.md). Thus

```math
J(\mathcal E_{\Gamma\circ H})
\succeq
\lambda J(\mathcal E_{\mathrm{diag}(H)}),
\qquad\text{(2)}
```

```math
\mathcal E_{\mathrm{diag}(H)}
=\sum_\pi H_{\pi\pi}\mathcal N_\pi,
\qquad\text{(3)}
```

where $`\mathcal N_\pi`$ is the ordinary definite-order composition. Either normalized or unnormalized Choi conventions give the same rank argument.

A nonzero accepted map requires $`H\ne0`$. Since $`H`$ is positive, $`\mathrm{Tr}H>0`$, so at least one diagonal weight $`H_{\pi\pi}`$ is strictly positive. A successful effect cannot discard the entire diagonal contribution while retaining a nonzero map.

Full Choi rank of each slot means that its Kraus operators span $`M_d(\mathbb C)`$. The products in any definite-order composition also span that algebra. To see this, expand an arbitrary target matrix in one slot's Kraus span and expand the identity in each other slot's span; multiplying the expansions expresses that arbitrary matrix in the span of the word products. Consequently

```math
J(\mathcal N_\pi)\succ0
\qquad\text{for every occupied }\pi.
```

At least one term of (3) has positive weight. Its Choi matrix is positive definite, and every other term is positive semidefinite. The right side of (2), and hence the accepted Choi matrix, is positive definite.

Final target unitaries preserve positive definiteness. Summing nonzero corrected accepted outcomes cannot remove it. A nonzero unitary channel has a rank-one Choi matrix, which contradicts full Choi rank for $`d\ge2`$. $`\square`$

## 2. The two physical dephasing models

Global unobserved dephasing has

```math
\Gamma=\eta\boldsymbol1\boldsymbol1^T+(1-\eta)I_r.
```

For every $`0\le\eta<1`$, the theorem applies with $`\lambda=1-\eta`$. In the [physical three-qubit local model](noisy_routing_model.md#physical-codebook), $`\Gamma`$ is a simultaneous row-and-column permutation of

```math
\begin{pmatrix}1&\gamma\\\gamma&1\end{pmatrix}^{\otimes3}.
```

For every $`0\le\gamma<1`$, its smallest eigenvalue is $`(1-\gamma)^3>0`$. Thus the noisy [recovery frontiers](recovery_frontier.md) and [tilted theorem](tilted_readout.md) describe approximate recovery whenever $`\epsilon>0`$ and either of these strict dephasing models is present.

The ideal boundary $`\eta=1`$ or $`\gamma=1`$ has a singular correlation matrix, so the theorem makes no assertion there. This is consistent with the [exact four-call qubit construction](qubit.md).

## 3. Attribution and operational limits

The general positive-component obstruction to probabilistic resource distillation is established in [Regula and Takagi (2021), around equation (34)](https://doi.org/10.1038/s41467-021-24699-0). The model-specific step here proves that every nonzero permitted herald retains a diagonal full-rank component. Selection cannot remove that component within this routing model.

The hypotheses are sufficient conditions, not an if-and-only-if characterization of exact recovery. The proof does not apply to accessible noise flags or environments, replacement routing coherence supplied after the specified noise, intermediate target operations, or a noise model with a protected sector and no strictly positive Schur lower bound. In particular, it does not say that every form of controller noise destroys every exact conditional transformation. [Molitor et al. (2024)](https://doi.org/10.1038/s42005-024-01843-y) study an open-control model with special conditional-state protection under some dissipative conditions; that is a different noise model and access statement.

The theorem is analytical. The finite polynomial and grid calculations used elsewhere are not needed for its rank conclusion.
