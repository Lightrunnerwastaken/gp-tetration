"""Scan knob `thmul`: force MORE theta harmonics per iteration than the
accuracy-coupled formula allows.

Why this is the exponent question and not another tuning knob.

E5a (2026-07-13) measured that every Taylor mode k = 1..1700 decays at the same
1.67-1.74 digits/step, and named the physics: the slow part is the *sequential
resolution of the theta harmonics*, mixed uniformly into all Taylor degrees.
That closes two-grid preconditioning in Taylor space -- but it also says
exactly where Theta(p) comes from. Put the measured numbers together:

    ~1.28 digits per theta harmonic          (journal, theta decay)
    ~1.3-1.6 harmonics resolved per iteration (thsamples ~ 0.43*re)
    => I(p) ~ p / (1.28 * 1.45) ~ 0.54 p     vs measured 0.59 p

So the iteration count IS the harmonic cascade. Nothing else in the loop is
Theta(p).

The cascade is enforced by one line in loop():

    thsamples = floor(re*tht_re_mult + 6)

i.e. the number of harmonics resolved is tied to the accuracy *already*
achieved. Harmonic m needs ~1.28m digits to be visible, and those digits come
from the harmonics -- chicken and egg. The open question is whether that
coupling is FUNDAMENTAL (harmonic m genuinely cannot be resolved before the
accuracy exists) or merely CONSERVATIVE (the engine could ask for more
harmonics per step and get them).

If conservative, iterations fall roughly like 1/thmul while only the theta
block gets more expensive -- and the theta block is ~19% of the run (measured,
doubling probe, dps 400). Then choosing thmul growing with p would lower the
exponent itself, not just the constant. If fundamental, forcing it either does
not reduce the iteration count or destroys accuracy, and Theta(p) is real.

Either answer is worth having. Writes research/tools/_thmul.gp.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src" / "fatou_backend" / "vendor" / "fatou_fork.gp"
DST = Path(__file__).resolve().parent / "_thmul.gp"

PATCHES: list[tuple[str, str]] = [
    ("ctrmul = 1;  /* exp-023 scan knob: scales the theta sampling radius */",
     "ctrmul = 1;  /* exp-023 scan knob: scales the theta sampling radius */\n"
     "thmul = 1;   /* scan knob: scales the theta HARMONIC count per iteration */"),
    # loop(): the accuracy-coupled harmonic count
    ("      thsamples=max(6, floor(re*tht_re_mult+6));",
     "      thsamples=max(6, floor(re*tht_re_mult+6));\n"
     "      if (thmul != 1, thsamples = max(6, floor(thsamples*thmul)));"),
]


def main() -> int:
    text = SRC.read_text(encoding="utf-8")
    for i, (pat, rep) in enumerate(PATCHES):
        if text.count(pat) != 1:
            print(f"patch {i} matched {text.count(pat)} times -- abort", file=sys.stderr)
            return 1
        text = text.replace(pat, rep, 1)
    DST.write_bytes(text.encode("utf-8"))
    print(f"wrote {DST}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
