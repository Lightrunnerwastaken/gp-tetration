# The mixed-base tetration paper series

Preprints (J. Justus, July 2026). The series develops the mathematics behind
the base atlas in this repository (`fatou_backend.basechange`); this
repository is the computational companion cited by the papers.

| # | File | Title |
|---|------|-------|
| I | [`paper1_phase_law.pdf`](paper1_phase_law.pdf) | A Phase Law for Mixed-Base Tetration: Logarithmic Tails, Cocycle Structure, and Singular Asymptotics |
| II | [`paper2_circle_diffeomorphism.pdf`](paper2_circle_diffeomorphism.pdf) | The Mixed-Base Tetration Phase Map Is a Circle Diffeomorphism |
| III | [`paper3_stopped_channels.pdf`](paper3_stopped_channels.pdf) | Stopped Channels and Offset Alignment: De-Conditionalizing the Real-Channel Hypotheses for Mixed-Base Tetration |
| IV | [`paper4_fine_structure.pdf`](paper4_fine_structure.pdf) | Fine Structure of Mixed-Base Tetration Phases: Certified Evaluation, a Universal Mode-Decay Law, and a Boundary Essential Singularity |
| V | [`paper5_kernel_selection.pdf`](paper5_kernel_selection.pdf) | The Base-Derivative of the Kneser Family: Selection Principle, Uniqueness, and Existence via a Riemann–Hilbert Problem |

How the code relates to the papers:

- The **phase function Φ and its Fourier modes** (Papers I–II) are what
  `research/reference/phi_modes.json` tabulates and
  `basechange.sexp_anchor` / `slog_anchor` consume.
- The **certified evaluation and mode-decay law** (Paper IV) back the
  universal-tail model in the mode table and the ball-arithmetic
  certification prototype (`research/tools/m54_cert_proto.py`).
- The **µ-hub** (`research/reference/mu_hub.json`) tabulates the constants
  µ_{e,b} appearing throughout the series.
