# Sharp coefficient budgets for certified singularity counts

Companion artifact for **Sharp Coefficient Budgets for Certified Singularity Counts in Periodic Burgers Flow**, Shivam Gupta, 17 October 2026. This is a research manuscript, not a claim of acceptance or peer review.

[Read the original public paper](paper/main.pdf).

**Revision clarification (19 October 2026):** initial positivity and initial regularity are separate. An l1 perturbation can define only a continuous initial potential. Positive-time heat smoothing remains valid; smooth initial velocity additionally requires a smooth perturbation. Finite Fourier perturbations attain the boundary-zero distances, so the distances remain sharp within the smooth class for a smooth nominal potential. The original PDF predates this clarification; see REVISION_NOTES_2026-10-19.txt for the journal revision and related-work disclosure.

The implementation counts zeros of periodic heat potentials using ball arithmetic, includes the infinite Bessel-series remainder, and encloses the distance to a boundary zero under initial Fourier coefficient l1 uncertainty. Conjugate-symmetric errors preserve real initial potentials and yield a different value set from arbitrary complex errors.

## Reproduce

Tested with Python 3.12.14 on macOS arm64. Install Python and create an isolated environment:

```sh
python3 -m venv .venv
. .venv/bin/activate
pip install -r requirements.txt
python -m unittest discover -s tests -v
python experiments/run.py
python experiments/validate.py
python experiments/figures.py
```

Build the manuscript with [Tectonic](https://tectonic-typesetting.github.io/):

```sh
cd paper
tectonic main.tex
```

The numerical commands overwrite generated results with the new run; timings will differ. All direct scientific dependencies are pinned. The environment record describes the original execution; this is not a container image or a guarantee of bitwise equality across library builds.

## Minimal use

```python
import sys
sys.path.insert(0, "src")
from flint import ctx
from pole_budgets import Potential, certified_count, radius_bounds
ctx.prec = 160
p = Potential(tau="0.5", a="2", N=40)
rect = ("0", "2", "-1", "1")
# Horizontal coordinates are multiples of pi; vertical coordinates are unscaled.
print(certified_count(p, rect))
print(radius_bounds(p, rect, real=True))
```

A budget enclosure has `lower`, `upper` (Arb bounds), approximate float summaries, and `converged`. Use a strict inequality against the rigorous lower endpoint for an acceptance rule. `converged` means the requested enclosure-width criterion succeeded, not that a larger user-specified error is safe. A returned `unresolved` count is not a zero count. The finite-support helper `FinitePotential` permits additional fixtures; the production infinite-series evaluator is specialized to exp(a cos x).

## Recorded evidence

- 16 primary contour counts and 32 global two-sided uncertainty-distance enclosures, all certified/converged.
- Twelve unit-test methods, including twelve analytic inner/outer zero-count fixtures and a boundary-zero refusal.
- 48 seeded perturbed-contour checks across both perturbation classes.
- Nine precision/truncation sensitivity runs.
- An exact two-block geometric construction with single-pair cost 1.2352341853... and two-block cost one.

See `experiments/protocol.json`, `results/grid.json`, and `results/validation.json`. The protocol is a recorded study specification, not a public preregistration. Ordinary floating point proposes support directions and renders figures; all accepted bounds and counts use verified ball inequalities. This is tested scientific software, not a formally verified program.

## Interpretation

The norm is on initial **heat-potential coefficients**, not on velocity samples. The table counts heat zeros with multiplicity. The Burgers logarithmic derivative has a simple pole with residue -2*nu*m at a heat zero of multiplicity m. A potential-distance bound alone does not establish initial positivity: for the paper's exp(a cos x) family, smooth real initial data additionally require error < exp(-a). Absolute budgets use the paper's fixed potential normalization.

Distances apply to the entire cut rectangle, including its vertical edges. They can be conservative for a periodic-cylinder count. At the boundary-zero distance a boundary zero can be created; the paper does not assert a changed count for every larger error. No claim of universal novelty, patentability, commercial demand, or performance superiority is made.

The paper situates its contribution relative to Cole–Hopf analysis, structured zero exclusion/pseudospectra, and validated contour root counting. All are substantial prior art; the claim is the explicit Fourier heat specialization, real-symmetric geometry, sparse witnesses, and validated budget study.

## License

Code is MIT licensed. Manuscript and figures are provided for scholarly review and citation; the code license does not imply a separate patent or publication right.
