from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import mpmath as mp

from .gp_backend import FatouGP


PAIRS = [("e", "2"), ("2", "e"), ("2", "10"), ("e", "10"), ("10", "e"), ("10", "2")]
REDUCTION_HEIGHT = {"2": 3, "e": 2, "10": 1}
REDUCTION_SHIFT = {"2": 4, "e": 5, "10": 6}


def mp_base(base: str) -> mp.mpf:
    return mp.e if base == "e" else mp.mpf(base)


def log_base(x: mp.mpf, base: mp.mpf) -> mp.mpf:
    return mp.log(x) / mp.log(base)


@dataclass
class ProfileResult:
    outer: str
    inner: str
    samples: int
    dps: int
    theta: list[mp.mpf]
    raw: list[mp.mpf]
    centered: list[mp.mpf]
    mu: mp.mpf
    fourier: list[tuple[mp.mpf, mp.mpf]]
    rms: list[mp.mpf]
    max_roundtrip: mp.mpf


def reduction_argument(outer: str, inner: str, local_values: list[mp.mpf]) -> list[mp.mpf]:
    b = mp_base(outer)
    d = mp_base(inner)
    c = log_base(d, b)
    alpha = log_base(c, b)
    args: list[mp.mpf] = []
    if inner == "2":
        for t in local_values:
            y2 = alpha + c * mp.power(d, t)
            y3 = log_base(y2, b)
            y4 = log_base(y3, b)
            args.append(y4)
        return args
    if inner == "e":
        for t in local_values:
            y2 = alpha + c * mp.exp(mp.exp(t))
            y3 = log_base(y2, b)
            y4 = log_base(y3, b)
            y5 = log_base(y4, b)
            args.append(y5)
        return args
    if inner == "10":
        for t in local_values:
            y5 = alpha + c * t
            y6 = log_base(y5, b)
            args.append(y6)
        return args
    raise ValueError(f"Unsupported inner base: {inner}")


def _require_real(values: list[mp.mpc]) -> list[mp.mpf]:
    result: list[mp.mpf] = []
    for value in values:
        if abs(mp.im(value)) > mp.mpf("1e-20"):
            raise RuntimeError(f"Expected real GP output, got {value!r}")
        result.append(mp.re(value))
    return result


def reconstruct_profile(gp: FatouGP, outer: str, inner: str, samples: int, harmonics: int) -> ProfileResult:
    shift = REDUCTION_SHIFT[inner]
    local_height = REDUCTION_HEIGHT[inner]
    theta = [mp.mpf(j) / samples for j in range(samples)]
    local_x = [mp.mpf(local_height) + t for t in theta]
    local_tetration = _require_real(gp.sexp_batch(inner, local_x))
    reduced_args = reduction_argument(outer, inner, local_tetration)
    local_slog = _require_real(gp.slog_batch(outer, reduced_args))
    roundtrip = gp.roundtrip_residuals(outer, reduced_args[:: max(1, samples // 16)])

    raw = [mp.mpf(shift) + val - (mp.mpf(6) + th) for th, val in zip(theta, local_slog)]
    mu = mp.fsum(raw) / samples
    centered = [val - mu for val in raw]

    fourier: list[tuple[mp.mpf, mp.mpf]] = []
    two_pi = 2 * mp.pi
    for n in range(1, harmonics + 1):
        a_n = (2 / samples) * mp.fsum(centered[j] * mp.cos(two_pi * n * theta[j]) for j in range(samples))
        b_n = (2 / samples) * mp.fsum(centered[j] * mp.sin(two_pi * n * theta[j]) for j in range(samples))
        fourier.append((a_n, b_n))

    rms: list[mp.mpf] = []
    for k in range(1, harmonics + 1):
        approx = []
        for th in theta:
            y = mp.mpf("0")
            for n in range(1, k + 1):
                a_n, b_n = fourier[n - 1]
                y += a_n * mp.cos(two_pi * n * th) + b_n * mp.sin(two_pi * n * th)
            approx.append(y)
        resid = [centered[j] - approx[j] for j in range(samples)]
        rms.append(mp.sqrt(mp.fsum(r * r for r in resid) / samples))

    return ProfileResult(
        outer=outer,
        inner=inner,
        samples=samples,
        dps=gp.dps,
        theta=theta,
        raw=raw,
        centered=centered,
        mu=mu,
        fourier=fourier,
        rms=rms,
        max_roundtrip=max(abs(v) for v in roundtrip) if roundtrip else mp.mpf("0"),
    )


def amplitude_phase(a1: mp.mpf, b1: mp.mpf) -> tuple[mp.mpf, mp.mpf]:
    amplitude = mp.sqrt(a1 * a1 + b1 * b1)
    phase = mp.atan2(a1, b1)
    if phase < 0:
        phase += 2 * mp.pi
    return amplitude, phase


def save_csv(path: Path, result: ProfileResult) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["theta", "raw", "centered"])
        for theta, raw, centered in zip(result.theta, result.raw, result.centered):
            writer.writerow([mp.nstr(theta, 30), mp.nstr(raw, 30), mp.nstr(centered, 30)])


def print_result(result: ProfileResult) -> None:
    a1, b1 = result.fourier[0]
    amp1, phase1 = amplitude_phase(a1, b1)
    ratio21 = (
        mp.sqrt(result.fourier[1][0] ** 2 + result.fourier[1][1] ** 2) / amp1 if len(result.fourier) >= 2 else mp.mpf("0")
    )
    print(f"({result.outer},{result.inner})")
    print(f"  mu      = {mp.nstr(result.mu, 20)}")
    print(f"  A1      = {mp.nstr(amp1, 12)}")
    print(f"  phase   = {mp.nstr(phase1, 12)}")
    print(f"  A2/A1   = {mp.nstr(ratio21, 12)}")
    for idx, (a_n, b_n) in enumerate(result.fourier, start=1):
        print(f"  a{idx},b{idx} = {mp.nstr(a_n, 12)}, {mp.nstr(b_n, 12)}")
    for idx, rms in enumerate(result.rms, start=1):
        print(f"  RMS{idx}    = {mp.nstr(rms, 12)}")
    print(f"  max roundtrip residual = {mp.nstr(result.max_roundtrip, 8)}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Mixed-base phase profiles via fatou.gp")
    parser.add_argument("--outer", choices=["2", "e", "10"])
    parser.add_argument("--inner", choices=["2", "e", "10"])
    parser.add_argument("--all-pairs", action="store_true")
    parser.add_argument("--samples", type=int, default=256)
    parser.add_argument("--harmonics", type=int, default=3)
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--gp-exe")
    parser.add_argument("--fatou-gp")
    parser.add_argument("--csv")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not args.all_pairs and (args.outer is None or args.inner is None):
        raise SystemExit("Use --all-pairs or specify both --outer and --inner.")
    mp.mp.dps = max(args.dps + 30, 100)
    gp = FatouGP(gp_exe=args.gp_exe, fatou_gp=args.fatou_gp, dps=args.dps,
                 looplim=max(35, args.dps - 20), n_workers=args.workers)
    pairs = PAIRS if args.all_pairs else [(args.outer, args.inner)]
    results: dict[tuple[str, str], ProfileResult] = {}
    for outer, inner in pairs:
        result = reconstruct_profile(gp, outer, inner, args.samples, args.harmonics)
        print_result(result)
        results[(outer, inner)] = result
        if args.csv and len(pairs) == 1:
            save_csv(Path(args.csv), result)
    if args.all_pairs:
        for (b, d), result in results.items():
            if (d, b) in results and b < d:
                antisym = result.mu + results[(d, b)].mu
                print(f"mu({b},{d}) + mu({d},{b}) = {mp.nstr(antisym, 12)}")


if __name__ == "__main__":
    main()
