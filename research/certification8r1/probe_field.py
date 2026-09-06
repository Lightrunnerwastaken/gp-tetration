"""Exploratory probe of the geometric field v = Gamma(T)/T' off the real axis.

This is a *diagnostic*, not part of the proof path: it locates the analytic
strip and the endpoint behaviour so that the certified quadrature radii can
be chosen before the interval machinery is run.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from mpmath import mp


def cheb_eval(coefficients, x):
    s = 2 * x - 1
    b1 = mp.mpf(0)
    b2 = mp.mpf(0)
    for coefficient in coefficients[:0:-1]:
        b0 = 2 * s * b1 - b2 + coefficient
        b2, b1 = b1, b0
    return s * b1 - b2 + coefficients[0]


class Koenigs:
    def __init__(self, beta):
        self.c = mp.exp(beta)
        self.L = -mp.lambertw(-self.c, -1) / self.c
        self.lam = self.c * self.L
        self.log_lam = mp.log(self.lam)
        self.dL = self.c * self.L * self.L / (1 - self.lam)
        self.dlam = self.c * (self.dL + self.L)
        self.dlog_lam = self.dlam / self.lam

    def field(self, y, y_prime, stop="1e-16", max_depth=600):
        w = mp.mpc(y)
        v = mp.mpc(0)
        s = mp.mpc(1)
        previous = None
        for depth in range(1, max_depth + 1):
            w_next = mp.log(w) / self.c
            v_next = v / (w * self.c) - w_next
            s_next = s / (w * self.c)
            w, v, s = w_next, v_next, s_next
            delta = w - self.L
            raw = mp.arg(delta)
            if previous is None:
                argument = raw
            else:
                predicted = previous - mp.arg(self.lam)
                turns = mp.nint((predicted - raw) / (2 * mp.pi))
                argument = raw + 2 * mp.pi * turns
            previous = argument
            log_s = depth * self.log_lam + mp.log(abs(delta)) + 1j * argument
            dlog_s = depth * self.dlog_lam + (v - self.dL) / delta
            slog_derivative = s / delta
            chi = log_s / self.log_lam
            gamma = (dlog_s - chi * self.dlog_lam) / slog_derivative
            if abs(delta) < mp.mpf(stop):
                return gamma / y_prime, depth
        raise RuntimeError("no contraction")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trajectory", type=Path, required=True)
    parser.add_argument("--node", type=int, default=3)
    args = parser.parse_args()

    mp.dps = 60
    trajectory = json.loads(args.trajectory.read_text(encoding="utf-8"))
    node = trajectory["nodes"][args.node]
    beta = mp.mpf(node["beta"])
    coefficients = [mp.mpf(value) for value in node["state"]["coefficients"]]
    derivatives = [
        mp.mpf(value) for value in node["state"]["derivative_coefficients"]
    ]
    koenigs = Koenigs(beta)

    def value_at(t):
        y = cheb_eval(coefficients, t)
        y_prime = cheb_eval(derivatives, t)
        return koenigs.field(y, y_prime)

    print(f"base node {args.node}: beta = {mp.nstr(beta, 20)}")
    print(f"L = {mp.nstr(koenigs.L, 15)}   lambda = {mp.nstr(koenigs.lam, 15)}")
    print()
    print("real axis: |v| and Im v")
    for exponent in range(-18, 0):
        t = mp.mpf(10) ** exponent
        try:
            v, depth = value_at(t)
            print(
                f"  t=1e{exponent:<4d} |v|={mp.nstr(abs(v),8):>14s} "
                f"Im v={mp.nstr(mp.im(v),8):>14s} depth={depth}"
            )
        except Exception as error:  # noqa: BLE001
            print(f"  t=1e{exponent:<4d} FAILED {error}")
    for t in ("0.25", "0.5", "0.75", "0.9", "0.99", "0.999999"):
        v, depth = value_at(mp.mpf(t))
        print(
            f"  t={t:<9s} |v|={mp.nstr(abs(v),8):>14s} "
            f"Im v={mp.nstr(mp.im(v),8):>14s} depth={depth}"
        )

    print()
    print("off-axis growth of |v| (analytic strip probe)")
    header = "   height   " + "".join(
        f"{name:>13s}" for name in ("t=0.05", "t=0.25", "t=0.5", "t=0.75", "t=0.95")
    )
    print(header)
    for height_text in (
        "0.01",
        "0.05",
        "0.1",
        "0.2",
        "0.3",
        "0.5",
        "0.8",
        "1.2",
    ):
        height = mp.mpf(height_text)
        row = f"  {height_text:>8s}   "
        for centre in ("0.05", "0.25", "0.5", "0.75", "0.95"):
            t = mp.mpc(mp.mpf(centre), height)
            try:
                v, _ = value_at(t)
                row += f"{mp.nstr(abs(v), 6):>13s}"
            except Exception:  # noqa: BLE001
                row += f"{'break':>13s}"
        print(row)

    print()
    print("periodicity check  v(t+1) vs v(t) (RH field should be 1-periodic)")
    for centre in ("0.25", "0.5"):
        a, _ = value_at(mp.mpf(centre))
        b, _ = value_at(mp.mpf(centre) + 1)
        print(
            f"  t={centre}: |v(t)|={mp.nstr(abs(a),10)}  "
            f"|v(t+1)|={mp.nstr(abs(b),10)}  "
            f"|diff|={mp.nstr(abs(a-b),6)}"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
