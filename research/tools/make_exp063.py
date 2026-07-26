"""exp-063: two hygiene fixes in the engine. No arithmetic changes.

Both come from the 2026-07-26 review. Neither is reachable in any run measured
so far; both are one line and remove a whole failure class.

1. **staylor leaks five locals into the global namespace.** `exdig`, `exin`,
   `exp0`, `exsc` and `icsc` are assigned in the body but missing from the
   `local(...)` list, so in this GP dialect they become globals. They are
   genuine temporaries -- `exp0` in particular holds the saved working precision
   across the reduced-precision extraction block. Consequences today are
   confined to hygiene, but not zero: `state_cache.discover_state_var_names`
   dumps *every* user variable, so these ride along into every serialized
   engine state, and any future function reusing one of those names would
   silently share storage with the extraction.

2. **thsamples can go negative, and goes straight into `vector()`.**
   `tht_re_mult` is initialized NEGATIVE (line ~1757: log(10)/log(7/2000) =
   -0.407) and is only replaced by a positive value further down, inside
   `if (thlog<>0, ...)` -- i.e. only when the sampled theta coefficient is
   nonzero. On the first pass re = -3, so `floor(re*tht_re_mult+6)` = 7 and all
   is well. But if that coefficient is ever exactly zero, the multiplier keeps
   its negative seed, and from re > ~14.7 onward `thsamples` turns negative and
   is handed to `vector(samples)` unguarded.

   Exactly-zero is not something any measured run has produced, and the review
   rated it minor for that reason. It is still a floor worth having: the failure
   mode is a PARI error deep inside thtaylor with no hint of where it came from.
   Clamped at 6, which is a length the code already uses deliberately
   (`thtaylor(1,6)` in matrix_ir).

The acceptance test for both is that nothing moves: the anchor must stay
bit-identical and the gate must pass, because neither change touches a value.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src" / "fatou_backend" / "vendor" / "fatou_fork.gp"
DST = Path(__file__).resolve().parent / "_exp063.gp"

PATCHES: list[tuple[str, str]] = [
    # 1. declare the five leaked temporaries
    ("  local(rinv,s,t,x1,y,y0,y1,y2,st,z,tot,t_est,tcrc,halfsamples,wtaylor,terms,om,c0,mu,c1,G,coeffs);",
     "  /* exp-063: exdig/exin/exp0/exsc/icsc were assigned but not declared, so\n"
     "     they leaked into the global namespace and rode along in every dumped\n"
     "     engine state. They are temporaries of the extraction block. */\n"
     "  local(rinv,s,t,x1,y,y0,y1,y2,st,z,tot,t_est,tcrc,halfsamples,wtaylor,terms,om,c0,mu,c1,G,coeffs,\n"
     "        exdig,exin,exp0,exsc,icsc);"),
    # 2. clamp thsamples
    ("      thsamples=floor(re*tht_re_mult+6);",
     "      /* exp-063: tht_re_mult starts NEGATIVE (log(10)/log(7/2000)) and only\n"
     "         turns positive once a nonzero theta coefficient is seen. If that\n"
     "         coefficient were ever exactly zero the seed would survive, and past\n"
     "         re ~ 14.7 this expression goes negative straight into vector(). 6 is\n"
     "         a length the code already uses on purpose (thtaylor(1,6)). */\n"
     "      thsamples=max(6, floor(re*tht_re_mult+6));"),
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
