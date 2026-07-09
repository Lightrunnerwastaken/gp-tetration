"""Calibrate the benchmark noise threshold from repeated warm runs."""
from __future__ import annotations

import argparse
import json
import statistics
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
NOISE_PATH = REPO / "research" / "noise.json"


def noise_threshold_pct(samples: list[float]) -> float:
    if len(samples) < 3:
        raise ValueError("need at least 3 samples")
    median = statistics.median(samples)
    spread_pct = 100.0 * (max(samples) - min(samples)) / median
    return max(3.0, 2.0 * spread_pct)


def main() -> None:
    parser = argparse.ArgumentParser(description="calibrate benchmark noise")
    parser.add_argument("--repeat", type=int, default=5)
    args = parser.parse_args()

    label = "noise-calibration"
    subprocess.run(
        [sys.executable, str(REPO / "bench" / "benchmark.py"),
         "--mode", "throughput", "--repeat", str(args.repeat), "--label", label],
        check=True,
    )
    results = sorted((REPO / "bench" / "results").glob(f"*-{label}.json"))
    payload = json.loads(results[-1].read_text(encoding="utf-8"))
    # total warm digits/s per repetition, summed over groups
    per_rep: dict[int, float] = {}
    for run in payload["runs"]:
        if run["mode"] == "warm":
            per_rep[run["rep"]] = per_rep.get(run["rep"], 0.0) + run["digits_per_second"]
    samples = [per_rep[k] for k in sorted(per_rep)]
    threshold = noise_threshold_pct(samples)
    NOISE_PATH.write_text(
        json.dumps({"threshold_pct": threshold, "samples": samples,
                    "source": results[-1].name}, indent=1),
        encoding="utf-8")
    print(f"[calibrate_noise] threshold {threshold:.1f}% from {len(samples)} runs -> {NOISE_PATH}")


if __name__ == "__main__":
    main()
