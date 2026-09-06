# The mixed-base tetration paper series

Preprints (J. Justus, July 2026), permanently archived on Zenodo — the
DOIs below are concept DOIs and always resolve to the latest version.
The series develops the mathematics behind the base atlas in this
repository (`fatou_backend.basechange`); this repository is the
computational companion cited by the papers.

© 2026 Janis Justus, all rights reserved (per the Zenodo records) —
the PDFs are *not* covered by the repository's MIT license.

| # | File | Title | DOI |
|---|------|-------|-----|
| I | [`paper1_phase_law.pdf`](paper1_phase_law.pdf) | A Phase Law for Mixed-Base Tetration: Logarithmic Tails, Cocycle Structure, and Singular Asymptotics | [10.5281/zenodo.21346803](https://doi.org/10.5281/zenodo.21346803) |
| II | [`paper2_circle_diffeomorphism.pdf`](paper2_circle_diffeomorphism.pdf) | The Mixed-Base Tetration Phase Map Is a Circle Diffeomorphism | [10.5281/zenodo.21361967](https://doi.org/10.5281/zenodo.21361967) |
| III | [`paper3_stopped_channels.pdf`](paper3_stopped_channels.pdf) | Stopped Channels and Offset Alignment: De-Conditionalizing the Real-Channel Hypotheses for Mixed-Base Tetration | [10.5281/zenodo.21362206](https://doi.org/10.5281/zenodo.21362206) |
| IV | [`paper4_fine_structure.pdf`](paper4_fine_structure.pdf) | Fine Structure of Mixed-Base Tetration Phases: Certified Evaluation, a Universal Mode-Decay Law, and a Boundary Essential Singularity | [10.5281/zenodo.21363593](https://doi.org/10.5281/zenodo.21363593) |
| V | [`paper5_kernel_selection.pdf`](paper5_kernel_selection.pdf) | The Base-Derivative of the Kneser Family: Selection Principle, Uniqueness, and Existence via a Riemann–Hilbert Problem | [10.5281/zenodo.21363809](https://doi.org/10.5281/zenodo.21363809) |
| VI | [`paper6_anchor_pure_reconstruction.pdf`](paper6_anchor_pure_reconstruction.pdf) | Anchor-Pure Reconstruction of Regular Tetration: Riemann–Hilbert Selection, Analytic Base Flow, and Blind Validation | *not yet archived* |

How the code relates to the papers:

- The **phase function Φ and its Fourier modes** (Papers I–II) are what
  `research/reference/phi_modes.json` tabulates and
  `basechange.sexp_anchor` / `slog_anchor` consume.
- The **certified evaluation and mode-decay law** (Paper IV) back the
  universal-tail model in the mode table and the ball-arithmetic
  certification prototype (`research/tools/m54_cert_proto.py`).
- The **µ-hub** (`research/reference/mu_hub.json`) tabulates the constants
  µ_{e,b} appearing throughout the series.
- The **anchor-pure reconstruction** (Paper VI) is implemented in
  `fatou_backend.reconstruction` (three modes: fast, validated, certified) and
  audited by `fatou_backend.certification`; the segment certificate for `8r1`
  is produced and independently replayed under `research/certification8r1/`.

Paper VI carries no DOI yet — the row above says so rather than leaving the
column to look archived. Add it once the Zenodo record exists.

**Reproducing Paper VI's segment certificate.** The `8r1` trajectory is not
imported: it is computed here, from this repository's own engine, by
`research/certification8r1/produce_8r1.py`. On 2026-07-31 that reproduction
matched the trajectory shipped with the paper *bit for bit*
(`sha256(nodes) = da3898e549bf88743df5a1295093626e09db197a7addca3a87970269540ed791`)
across a different fork revision, a different PARI binary and a different
operating system — the provenance manifest records those differences instead
of hiding them. See `research/certification8r1/README.md`.
