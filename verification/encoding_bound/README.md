# Exact cost bound over all six-label encodings

The [complete proof](../../science/all_encoding_cost_bound.md) fixes the
rank-two measure-Z/prepare-X channel and the local phase-noise benchmark
`x = 9/10`, `q² = 729/1000`. It covers all 20,160 injective assignments of the
six three-call order labels to six computational three-bit words, arbitrary
controller preparation and final instruments, and final target **unitary**
corrections. The 97 certificates prove the common three-call support ceiling
`27/10 < 3`. Together with the one-/two-call support witness, this gives the
optimal mixed one-/two-/three-call cost `101/32` at error `3/8`, with a strict
penalty for a positive raw three-call share.

All commands below use only the Python standard library. Run them from the
repository root; `--output-dir` keeps regenerated records separate from the
committed certificates.

```bash
python3 verification/encoding_bound/encoding_certificate/verify_certificate.py --output-dir /tmp/qod-encoding-bound/positive-matrices
python3 verification/encoding_bound/encoding_orbits/verify_coverage.py --output-dir /tmp/qod-encoding-bound/coverage
python3 verification/encoding_bound/encoding_certificate_audit/verify_reduction.py --output-dir /tmp/qod-encoding-bound/reduction
python3 verification/encoding_bound/encoding_certificate_audit/verify_baseline.py --output-dir /tmp/qod-encoding-bound/baseline
```

| Program | What is recomputed | Expected result |
|---|---|---|
| `encoding_certificate/verify_certificate.py` | Integer Kraus kernels, rational Gram residuals, algebraic visibility enclosure, perturbation margin, coverage and attaining protocol | `PASS`, 97 certificates and 388 positive-matrix checks |
| `encoding_orbits/verify_coverage.py` | Every encoding's explicit symmetry transformation, cube classes, combined classes and an independent Burnside count | 20,160 encodings, 420 cube classes, 97 combined classes; orbit sizes 48, 144 and 288 with multiplicities 3, 49 and 45 |
| `encoding_certificate_audit/verify_reduction.py` | Bell cross-kernel zeros, physical-slot invariance, and support normalization | `PASS`, 180 zero entries and 1,296 exact kernel invariances |
| `encoding_certificate_audit/verify_baseline.py` | One-/two-call support and matching raw mixture | `PASS`, raw weights `21/61`, `40/61`, acceptance `32/61`, cost `101/32` |

The main certificate uses the rational outward enclosure
`170763/200000 - 1/10000000 < q < 170763/200000` and verifies it by squaring.
For each matrix it checks a rational factorization residual by strict diagonal
dominance. The smallest exact diagonal-dominance margin is
`503251611580961/512000000000000000`, larger than the uniform operator
perturbation allowance `3/100000`. No numerical eigensolver, optimization
solver, floating-point Cholesky factor, or solver log is trusted in this check.

`encoding_dual/certificate_candidates.json` contains the rational dual
matrices; `encoding_certificate/rational_gram_factors.json` contains rational
auxiliary factors. Their original relative directory layout is retained so
the verifier resolves both inputs from its own location. The other JSON files
record the independently reproducible checks. Omitting `--output-dir`
regenerates records next to the corresponding program.

The finite certificate proves its exact matrix inequalities. The scientific
page supplies the controller-Gram reduction, continuous-unitary contraction,
cost duality, common-error interval and scope. The encoding class is six
injectively assigned computational words; the result does not optimize
redundant eight-word or noncomputational encodings, general target recovery
channels, or a full noise/accuracy phase diagram.
