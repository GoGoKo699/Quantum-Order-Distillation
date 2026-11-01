# Local phase-noise controller certificates

The [complete scientific account](../../science/local_controller_encoding.md)
derives the physical three-bit controller model, the pairwise-code frontier,
and a lower-noise advantage from an alternative encoding. These are results
for the fixed rank-two measure-Z/prepare-X channel, compulsory independent
call slots, and final target **unitary** corrections. They concern recovered
reference–target entanglement and accepted target-channel error.

Run from the repository root:

```bash
python3 verification/local_controller/verify_spectral.py --output-dir /tmp/qod-local-controller/spectral
python3 verification/local_controller/exact_flags.py --output-dir /tmp/qod-local-controller/flags
python3 verification/local_controller/certify_thresholds.py --output-dir /tmp/qod-local-controller/thresholds
python3 verification/local_controller/physical_proportional_certificate.py --output-dir /tmp/qod-local-controller/alternative
```

The flag reconstruction uses NumPy with exact `Fraction` entries. The other
three programs use only the standard library for their exact checks;
`verify_spectral.py` also performs supplementary random checks when NumPy is
available. SciPy is not required. All outputs follow `--output-dir`; omitting
it regenerates the corresponding committed JSON records next to the scripts.

| Program | Exact check | Expected result |
|---|---|---|
| `verify_spectral.py` | Four projectors, eight polynomial kernel identities, 216 slot/Hamming invariances, symmetric and antisymmetric cross kernels, and the determinant threshold polynomial | `PASS`; complex-controller spectral support for the pairwise code at `q ≥ q₀ ≈ 0.811155656051003` |
| `exact_flags.py` | Direct integer Kraus products and every cubic coefficient of the corrected flag Choi maps | Four flags and four corrections, each coefficient Bell diagonal and scaled trace preserving |
| `certify_thresholds.py` | Rational series bounds for the logarithm/exponential and exact signs at outward threshold brackets | Routing-time ratio in `(1.14647049186791, 1.14647049186792)`; proportional-schedule two-call visibility near `0.9763304783204344` |
| `physical_proportional_certificate.py` | Exact Choi positivity, Bell diagonality, scaled trace preservation, all native flag subsets and the attaining raw mixture | Error `3/8`; cost `36661573412240983/12985210999106966 < 26609/9208` |

In lexicographic matrix-product order, the pairwise code is
`000,001,100,110,011,111`. Its spectral proof uses
`8q⁴ + 2q³ + q² + q − 6 ≥ 0`; below its positive root the native flag maps
remain a construction, without a four-flag optimality claim. The spectral
program checks the full complex cross-kernel identity, so real-controller
tests are not substituted for the complex-controller converse. Its 320
random complex-Gram/unitary trials are supplementary evidence, not the proof.

The alternative code is `000,011,101,010,100,001`. Its exact construction uses
`t = 49/50`, two-call visibility `x = t²`, three-call per-bit visibility
`q = t³`, and routed windows `Tₘ = mτ`. Reject the alternating flag and retain
the other three with their specified unitary corrections; mixing this
three-call protocol with the one-call baseline beats the optimal one-/two-call
cost by about 2.299%. The script optimizes only its finite native subset hull;
it makes no global-optimum claim at this lower-noise point.

The [separate all-encoding certificate](../encoding_bound/README.md) proves a
global cost bound at the fixed noisier benchmark. Global logical visibility,
physical per-bit visibility and the two-/three-call routing windows are
different parameters and are not interchanged by these checks. Preparation,
readout and gate faults are outside the charged routing-window model.
