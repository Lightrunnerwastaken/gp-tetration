"""F6 probe: how many digits does a universal asymptotic form carry for the
theta harmonics?

Why this is THE question for the exponent. With I(p) = Theta(p) Picard
iterations, every iteration must touch the Theta(p^2) digits of ct, so the
wall-clock floor is Theta(p^3) no matter how fast the evaluation primitive
gets. Exponent 2.1 therefore requires I(p) to fall -- and the only live
candidate is: predict the high theta harmonics from a closed form with O(1)
parameters and iterate only the deviation.

loop() resolves ~1.3-1.6 theta harmonics per iteration (thsamples grows like
0.43*re) at a decay of ~1.28 digits per harmonic. If a fitted asymptotic form
reproduces theta_m to d digits at large m, the iteration only has to resolve
the remaining (1.28*m - d) digits per harmonic -- i.e. the harmonic count that
has to be *iterated* shrinks by d/1.28 harmonics, a CONSTANT, unless d grows
with m.

So the measurement is not "does a form fit" but "does the fit residual
IMPROVE with m". Constant residual => F6 is a constant-factor idea and the
Theta(p) iteration count stands. Growing accuracy => the iteration count can
be broken and 2.1 is back on the table.

    python research/tools/theta_modes.py --dps 300
"""
from __future__ import annotations

import argparse
import math
import re as _re
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
FORK = REPO / "src" / "fatou_backend" / "vendor" / "fatou_fork.gp"
GP = r"C:\Program Files (x86)\Pari64-2-17-3\gp.exe"


def dump(dps: int, base: str) -> list[complex]:
    script = "\n".join([
        "default(parisizemax, 8589934592);",
        f"default(realprecision, {dps});",
        f'read("{FORK.as_posix()}");',
        f"default(realprecision, {dps});",
        "quietmode=1;",
        f"tbase = {base};",
        "sexpinit(tbase,0,0,0);",
        "tv = Vecrev(tht);",
        'print("THETA-LEN ", #tv);',
        'for (i=1, #tv, print("THETA ", real(tv[i]), "|", imag(tv[i])));',
        "quit;",
    ]) + "\n"
    p = subprocess.run([GP, "-q", "-f"], input=script, capture_output=True,
                       text=True, encoding="utf-8", errors="replace", cwd=str(REPO))
    out = []
    for line in p.stdout.splitlines():
        if line.startswith("THETA "):
            a, _, b = line[6:].partition("|")
            f = lambda s: float(s.replace(" E", "e").replace(" ", "") or 0)
            try:
                out.append(complex(f(a), f(b)))
            except ValueError:
                out.append(0j)
    if not out:
        raise SystemExit(p.stdout[-1500:] + p.stderr[-500:])
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dps", type=int, default=300)
    ap.add_argument("--base", default="exp(1)")
    args = ap.parse_args()
    th = dump(args.dps, args.base)
    mags = [abs(z) for z in th]
    print(f"theta harmonics: {len(th)}")
    print(f"{'m':>5s} {'|theta_m|':>12s} {'log10':>9s} {'d log10/dm':>11s}")
    prev = None
    rows = []
    for m, a in enumerate(mags):
        if a == 0:
            continue
        l = math.log10(a)
        d = (l - prev) if prev is not None else float("nan")
        prev = l
        rows.append((m, a, l, d))
    for r in rows[:6] + rows[len(rows)//2 - 2: len(rows)//2 + 2] + rows[-6:]:
        print(f"{r[0]:5d} {r[1]:12.4e} {r[2]:9.3f} {r[3]:11.4f}")

    # Fit log10|theta_m| = A - B*m + C*log(m) on the tail, then report how
    # many digits the fit carries at each m. The slope of that curve is the
    # answer: flat => constant-factor idea, rising => the iteration count can
    # be attacked.
    tail = [r for r in rows if r[0] >= max(8, len(rows) // 4) and r[1] > 0]
    if len(tail) >= 8:
        import statistics
        n = len(tail)
        X1 = [r[0] for r in tail]
        X2 = [math.log(r[0]) for r in tail]
        Y = [r[2] for r in tail]
        # normal equations for  Y = A + B*X1 + C*X2
        def dot(u, v): return sum(a * b for a, b in zip(u, v))
        one = [1.0] * n
        M = [[dot(one, one), dot(one, X1), dot(one, X2)],
             [dot(X1, one), dot(X1, X1), dot(X1, X2)],
             [dot(X2, one), dot(X2, X1), dot(X2, X2)]]
        b = [dot(one, Y), dot(X1, Y), dot(X2, Y)]
        # 3x3 solve
        for i in range(3):
            piv = M[i][i]
            for j in range(i + 1, 3):
                f = M[j][i] / piv
                for k in range(3):
                    M[j][k] -= f * M[i][k]
                b[j] -= f * b[i]
        c = [0.0] * 3
        for i in (2, 1, 0):
            c[i] = (b[i] - sum(M[i][k] * c[k] for k in range(i + 1, 3))) / M[i][i]
        print(f"\nfit log10|theta_m| = {c[0]:.4f} + {c[1]:.6f}*m + {c[2]:.4f}*log(m)"
              f"   (tail m >= {tail[0][0]}, {n} points)")
        print(f"{'m':>6s} {'digits carried by fit':>22s}")
        for r in tail[::max(1, n // 12)]:
            pred = c[0] + c[1] * r[0] + c[2] * math.log(r[0])
            print(f"{r[0]:6d} {abs(pred - r[2]):22.4f}")
        res = [abs(c[0] + c[1] * r[0] + c[2] * math.log(r[0]) - r[2]) for r in tail]
        half = n // 2
        print(f"\nmean |residual| first half {statistics.fmean(res[:half]):.4f} dec, "
              f"second half {statistics.fmean(res[half:]):.4f} dec")
        print("VERDICT: residual must FALL with m for the iteration count to be "
              "attackable; flat or rising means F6 is a constant-factor idea.")


if __name__ == "__main__":
    main()
