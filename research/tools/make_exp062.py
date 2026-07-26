"""exp-062: put the working precision into the exmap cache key.

`exmap` is the extraction's coefficient mapping -- (rinv^s/samples)*conj(c0)^s*om^(-s)
in the complex branch, conj(c1)^s*mu^(-s) in the real one. It is cached on
[samples, w, r, branch]: grid identity only, no precision.

But the extraction runs inside a REDUCED-precision block
(`default(realprecision, exdig)`, exdig = realprecision - re + 40), and while
the working-precision ladder is still climbing, exdig climbs with it. So a map
built early in a grid stretch is reused later in the same stretch at a higher
precision than it carries -- silently, because a t_REAL simply keeps the
precision it was created with.

This is the same defect class as the mixdft twiddle cache (exp-061), where
dropping the precision clause cost a measured 232 digits. The 2026-07-26 review
found it and its verifier downgraded it to a nit. Counted instead of argued:

    base e,   dps 300 : 223 cache hits, 68 of them above the build precision,
                        worst shortfall 58 digits
    base 1+I, dps 120 :  90 cache hits, 27 above, worst shortfall 58 digits

The fix mirrors bluedft/mixdft: keep the entry only while it was built at a
precision >= the current one, so a map built HIGH is still reusable LOW (that
direction is free), and rebuild otherwise.

Whether it moves the final digit count is a separate question this patch does
not assume -- measure it. The mapping multiplies the *diff* coefficients, whose
own scale is 10^-re, so part of the shortfall may sit below what survives the
addition to excoef.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src" / "fatou_backend" / "vendor" / "fatou_fork.gp"
DST = Path(__file__).resolve().parent / "_exp062.gp"

PATCHES: list[tuple[str, str]] = [
    ("exmap=0; exmapkey=0;",
     "exmap=0; exmapkey=0; exmapprec=0;   /* exp-062: precision the map carries */"),
    ("""      if ((exmapkey != [samples, w, r, 1]) || (type(exmap) != "t_VEC"),
        exmap = vector(terms, s, (rinv^s/samples) * conj(c0)^s * om^(-s));
        exmapkey = [samples, w, r, 1];
      );""",
     """      /* exp-062: the key must carry the working precision. The extraction
         runs at the reduced exdig, which CLIMBS while the precision ladder
         does, so a map built early in a grid stretch was reused later at a
         precision it does not carry (measured at dps 300: 68 of 223 hits,
         worst shortfall 58 digits). Reuse only downwards, as bluedft and
         mixdft already do. */
      if ((exmapkey != [samples, w, r, 1]) || (type(exmap) != "t_VEC")
          || (exmapprec < default(realprecision)),
        exmap = vector(terms, s, (rinv^s/samples) * conj(c0)^s * om^(-s));
        exmapkey = [samples, w, r, 1];
        exmapprec = default(realprecision);
      );"""),
    ("""      if ((exmapkey != [samples, w, r, 2]) || (type(exmap) != "t_VEC"),
        exmap = vector(terms, s, conj(c1)^s * mu^(-s));
        exmapkey = [samples, w, r, 2];
      );""",
     """      /* exp-062: same for the real branch (this is the one base e takes). */
      if ((exmapkey != [samples, w, r, 2]) || (type(exmap) != "t_VEC")
          || (exmapprec < default(realprecision)),
        exmap = vector(terms, s, conj(c1)^s * mu^(-s));
        exmapkey = [samples, w, r, 2];
        exmapprec = default(realprecision);
      );"""),
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
