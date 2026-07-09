"""Persist sexpinit state per (base, dps, knobs, fatou-hash) via writebin/read.

Discovery uses the \\uv metacommand (verified on gp 2.17.3): it lists all user
variables as "name =" lines with indented value lines. We dump every variable
into one binary vector; functions are NOT dumped (re-reading fatou.gp is cheap
and restores them).
"""
from __future__ import annotations

import hashlib
import json
import os
import re
from pathlib import Path

import mpmath as mp

ANCHOR_EXPR = "sexp(0.5)"
_NAME_RE = re.compile(r"^([A-Za-z_][A-Za-z0-9_]*) =$")
# our own temp variables; note GP identifiers must start with a letter
_EXCLUDED = {"vv", "fatoustate"}


def cache_dir() -> Path:
    env = os.getenv("FATOU_CACHE_DIR")
    base = Path(env) if env else Path.home() / ".cache" / "fatou_backend"
    base.mkdir(parents=True, exist_ok=True)
    return base


def cache_key(base_expr: str, dps: int, nlim: int, nskip: int, looplim: int,
              fatou_gp: Path, gp_exe: Path) -> str:
    fatou_hash = hashlib.sha256(Path(fatou_gp).read_bytes()).hexdigest()
    blob = "|".join([base_expr, str(dps), str(nlim), str(nskip), str(looplim),
                     fatou_hash, str(gp_exe)])
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


def discover_state_var_names(worker) -> list[str]:
    output = worker.eval_raw(["\\uv"], "__UVDONE__", timeout=120.0)
    names = []
    for line in output:
        match = _NAME_RE.match(line)
        if match and match.group(1) not in _EXCLUDED:
            names.append(match.group(1))
    return names


def dump_state(worker, names: list[str], bin_path: Path) -> None:
    gp_path = bin_path.as_posix()
    worker.eval_raw([f'writebin("{gp_path}", [{", ".join(names)}]);'],
                    "__DUMPDONE__", timeout=300.0)
    anchor = worker.eval([ANCHOR_EXPR])[0]
    sidecar = {
        "names": names,
        "anchor_real": mp.nstr(mp.re(anchor), 40, min_fixed=0, max_fixed=0),
        "anchor_imag": mp.nstr(mp.im(anchor), 40, min_fixed=0, max_fixed=0),
    }
    bin_path.with_suffix(".json").write_text(json.dumps(sidecar), encoding="utf-8")


def restore_lines(names: list[str], bin_path: Path) -> list[str]:
    gp_path = bin_path.as_posix()
    lines = [f'fatoustate = read("{gp_path}");']
    for idx, name in enumerate(names, start=1):
        lines.append(f"{name} = fatoustate[{idx}];")
    return lines


def load_sidecar(bin_path: Path) -> dict | None:
    sidecar_path = bin_path.with_suffix(".json")
    if not bin_path.exists() or not sidecar_path.exists():
        return None
    try:
        data = json.loads(sidecar_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    if not isinstance(data.get("names"), list) or "anchor_real" not in data:
        return None
    return data


def drop_cache(bin_path: Path) -> None:
    bin_path.unlink(missing_ok=True)
    bin_path.with_suffix(".json").unlink(missing_ok=True)
