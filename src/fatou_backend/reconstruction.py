"""Three-mode anchor-pure tetration reconstruction calculator.

The shipped legacy phase atlas supports fast evaluation without initializing
a target engine at run time.  Its base-2/base-3 entries were originally
sampled from target-base states, so they are labelled as legacy inputs rather
than silently advertised as anchor-pure proof data.  This module turns the
atlas capability into a strict reconstruction pipeline:

* ``fast`` chooses the nearest anchor and uses floating-point phase data;
* ``validated`` adds independent resolution, lift, normalization,
  functional-equation, and roundtrip gates;
* ``certified`` is reserved for a complete Paper VI interval certificate
  with independent array-level replay.  The current auditor deliberately
  rejects every final certified claim until that replay engine exists.

For a target not present in the atlas, a local RH marcher must be supplied.
The calculator never silently initializes the target-base ``fatou.gp``
engine.  That guard is important: doing so would make the result fast after
caching, but it would no longer be anchor-pure reconstruction.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Protocol, Sequence

import mpmath as mp

from . import basechange
from .certification import (
    CertificateAudit,
    CertificateError,
    load_and_audit_certificate,
)


class ReconstructionMode(str, Enum):
    FAST = "fast"
    VALIDATED = "validated"
    CERTIFIED = "certified"


class AtlasCoverageError(RuntimeError):
    """Raised when no exact atlas state or local marcher is available."""


class ValidationError(RuntimeError):
    """Raised when validated mode cannot support the requested digits."""


@dataclass(frozen=True)
class PhaseMap:
    """Fourier representation of H(x)-x for a degree-one circle lift."""

    base: mp.mpf
    beta: mp.mpf
    mu: mp.mpf
    modes: tuple[mp.mpc, ...]
    provenance: str
    metadata: Mapping[str, Any] = field(default_factory=dict)

    @staticmethod
    def _canonical_mpf(value: Any) -> tuple[int, int, int, int]:
        raw = mp.mpf(value)._mpf_
        return tuple(int(part) for part in raw)

    def fingerprint(self) -> str:
        """SHA-256 of the exact in-memory binary phase coefficients."""

        payload = {
            "schema": "paper-vi-phase-map-v1",
            "base": self._canonical_mpf(self.base),
            "beta": self._canonical_mpf(self.beta),
            "mu": self._canonical_mpf(self.mu),
            "modes": [
                {
                    "re": self._canonical_mpf(value.real),
                    "im": self._canonical_mpf(value.imag),
                }
                for value in self.modes
            ],
        }
        encoded = json.dumps(
            payload,
            sort_keys=True,
            separators=(",", ":"),
        ).encode("ascii")
        return hashlib.sha256(encoded).hexdigest()

    def value(self, theta: Any, k_head: int | None = None) -> mp.mpf:
        theta = mp.mpf(theta)
        modes = self.modes if k_head is None else self.modes[:k_head]
        total = mp.mpf(self.mu)
        for k, coefficient in enumerate(modes, start=1):
            total += (
                coefficient
                * mp.exp(mp.mpc(0, 2) * mp.pi * k * theta)
            ).real
        return total

    def derivative(self, theta: Any, k_head: int | None = None) -> mp.mpf:
        theta = mp.mpf(theta)
        modes = self.modes if k_head is None else self.modes[:k_head]
        total = mp.mpf(0)
        for k, coefficient in enumerate(modes, start=1):
            total += (
                mp.mpc(0, 2) * mp.pi * k * coefficient
                * mp.exp(mp.mpc(0, 2) * mp.pi * k * theta)
            ).real
        return total

    def omitted_head_bound(self, k_head: int) -> mp.mpf:
        """Exact l1 bound for modes stored but omitted at a coarse resolution."""

        return mp.fsum(abs(value) for value in self.modes[k_head:])

    def tail_forecast(self) -> mp.mpf:
        """Non-rigorous tail forecast from the last stored mode ratios.

        This is intentionally used only in fast/validated mode.  Certified
        mode obtains its tail radius from the interval certificate.
        """

        magnitudes = [abs(value) for value in self.modes if value]
        if not magnitudes:
            return mp.mpf(0)
        if len(magnitudes) < 3:
            return magnitudes[-1]
        ratios = [
            magnitudes[index] / magnitudes[index - 1]
            for index in range(len(magnitudes) - 2, len(magnitudes))
            if magnitudes[index - 1] != 0
        ]
        ratio = max(ratios) if ratios else mp.mpf("0.5")
        ratio = min(max(ratio, mp.mpf(0)), mp.mpf("0.95"))
        return magnitudes[-1] * ratio / (1 - ratio)


@dataclass(frozen=True)
class AtlasAnchor:
    key: str
    phase: PhaseMap

    @property
    def base(self) -> mp.mpf:
        return self.phase.base

    @property
    def beta(self) -> mp.mpf:
        return self.phase.beta

    @property
    def anchor_pure(self) -> bool:
        return bool(self.phase.metadata.get("anchor_pure", False))


class AtlasIndex:
    """Read-only atlas index; it never invokes a target-base engine."""

    def __init__(self, anchors: Sequence[AtlasAnchor], path: str):
        self.anchors = tuple(sorted(anchors, key=lambda item: item.beta))
        self.path = path

    @classmethod
    def load(cls, path: str | Path | None = None) -> "AtlasIndex":
        atlas_path = Path(path or basechange.PHI_TABLE_PATH)
        with atlas_path.open(encoding="utf-8") as handle:
            payload = json.load(handle)
        entries: list[AtlasAnchor] = []
        for key, item in payload.get("bases", {}).items():
            base = basechange._base_value(key)
            modes = tuple(
                mp.mpc(mp.mpf(mode["re"]), mp.mpf(mode["im"]))
                for mode in item.get("modes", [])
            )
            phase = PhaseMap(
                base=base,
                beta=mp.log(mp.log(base)),
                mu=mp.mpf(item["mu"]),
                modes=modes,
                provenance=f"legacy stored phase atlas entry {key}",
                metadata={
                    "atlas_key": key,
                    "n_grid": item.get("n_grid", payload.get("meta", {}).get("n_grid")),
                    "k_head": len(modes),
                    "measured": item.get("measured"),
                    "anchor_pure": bool(item.get("anchor_pure", False)),
                    "construction_definition": payload.get("meta", {}).get("definition"),
                },
            )
            entries.append(AtlasAnchor(key=key, phase=phase))
        # The e anchor has the identity phase map exactly.
        entries.append(AtlasAnchor(
            key="e",
            phase=PhaseMap(
                base=mp.e,
                beta=mp.mpf(0),
                mu=mp.mpf(0),
                modes=(),
                provenance="exact base-e identity anchor",
                metadata={"anchor_pure": True, "atlas_key": "e"},
            ),
        ))
        return cls(entries, str(atlas_path))

    def exact(
        self,
        base: Any,
        tolerance: Any = "1e-40",
        require_anchor_pure: bool = False,
    ) -> AtlasAnchor | None:
        target = basechange._base_value(base)
        tol = mp.mpf(tolerance)
        for anchor in self.anchors:
            if (
                (not require_anchor_pure or anchor.anchor_pure)
                and abs(anchor.base - target) <= tol * max(1, abs(target))
            ):
                return anchor
        return None

    def nearest(
        self,
        base: Any,
        require_anchor_pure: bool = False,
    ) -> AtlasAnchor:
        candidates = tuple(
            anchor
            for anchor in self.anchors
            if not require_anchor_pure or anchor.anchor_pure
        )
        if not candidates:
            qualifier = " anchor-pure" if require_anchor_pure else ""
            raise AtlasCoverageError(f"the phase atlas contains no{qualifier} anchors")
        target_beta = mp.log(mp.log(basechange._base_value(base)))
        return min(candidates, key=lambda item: abs(item.beta - target_beta))


@dataclass(frozen=True)
class LocalMarchResult:
    """Output contract for a floating-point local RH marcher.

    ``relative_phase`` is H_{c,b}-id, where c is the selected atlas anchor.
    The calculator composes it with H_{e,c} to obtain H_{e,b}.
    """

    relative_phase: PhaseMap
    residual: mp.mpf
    path_split_delta: mp.mpf | None = None
    resolution_delta: mp.mpf | None = None
    contraction_q_estimate: mp.mpf | None = None
    diagnostics: Mapping[str, Any] = field(default_factory=dict)


class LocalRHMarcher(Protocol):
    def march(
        self,
        anchor: AtlasAnchor,
        target_base: mp.mpf,
        digits: int,
    ) -> LocalMarchResult:
        ...


@dataclass(frozen=True)
class ReconstructionState:
    target_base: mp.mpf
    target_beta: mp.mpf
    phase: PhaseMap
    anchor: AtlasAnchor
    local_beta_distance: mp.mpf
    source: str
    diagnostics: Mapping[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class ErrorForecast:
    mode_tail: mp.mpf
    resolution: mp.mpf
    lift: mp.mpf
    normalization: mp.mpf
    functional_equation: mp.mpf
    roundtrip: mp.mpf
    path_split: mp.mpf | None = None
    residual: mp.mpf | None = None
    contraction_q: mp.mpf | None = None
    safety_factor: mp.mpf = mp.mpf(3)

    @property
    def residual_to_error(self) -> mp.mpf | None:
        if self.residual is None:
            return None
        if self.contraction_q is None:
            return self.residual
        if not (0 <= self.contraction_q < 1):
            return mp.inf
        return self.residual / (1 - self.contraction_q)

    @property
    def predicted_absolute(self) -> mp.mpf:
        components = [
            self.mode_tail,
            self.resolution,
            self.lift,
            self.normalization,
            self.functional_equation,
            self.roundtrip,
        ]
        if self.path_split is not None:
            components.append(self.path_split)
        if self.residual_to_error is not None:
            components.append(self.residual_to_error)
        return self.safety_factor * mp.fsum(abs(value) for value in components)

    @property
    def predicted_digits(self) -> float:
        value = self.predicted_absolute
        if value == 0:
            return float("inf")
        return float(-mp.log10(value))

    def as_dict(self) -> dict[str, Any]:
        def show(value: mp.mpf | None) -> str | None:
            return None if value is None else mp.nstr(value, 20)

        return {
            "mode_tail": show(self.mode_tail),
            "resolution": show(self.resolution),
            "lift": show(self.lift),
            "normalization": show(self.normalization),
            "functional_equation": show(self.functional_equation),
            "roundtrip": show(self.roundtrip),
            "path_split": show(self.path_split),
            "residual": show(self.residual),
            "contraction_q": show(self.contraction_q),
            "residual_to_error": show(self.residual_to_error),
            "safety_factor": show(self.safety_factor),
            "predicted_absolute": show(self.predicted_absolute),
            "predicted_digits": self.predicted_digits,
        }


@dataclass(frozen=True)
class ValidationReport:
    passed: bool
    requested_digits: int
    forecast: ErrorForecast
    failed_gates: tuple[str, ...] = ()
    notes: tuple[str, ...] = ()

    def require_pass(self) -> "ValidationReport":
        if not self.passed:
            gates = (
                ""
                if not self.failed_gates
                else f"; missing/failed gates: {', '.join(self.failed_gates)}"
            )
            raise ValidationError(
                f"validated reconstruction forecast="
                f"{self.forecast.predicted_digits:.2f} digits, "
                f"requested={self.requested_digits}{gates}"
            )
        return self


@dataclass(frozen=True)
class ReconstructionValue:
    height: mp.mpf
    center: mp.mpf
    radius: mp.mpf | None = None


@dataclass(frozen=True)
class ReconstructionResult:
    mode: ReconstructionMode
    state: ReconstructionState
    values: tuple[ReconstructionValue, ...]
    validation: ValidationReport | None
    certificate: CertificateAudit | None


def compose_phase_maps(
    outer: PhaseMap,
    inner: PhaseMap,
    target_base: Any,
    n_grid: int = 128,
    k_head: int | None = None,
) -> PhaseMap:
    """Compose H_outer o H_inner and return its periodic correction.

    If ``outer`` is H_{e,c} and ``inner`` is H_{c,b}, the result is
    H_{e,b}.  Composition is evaluated on a uniform circle grid and
    transformed back to real Fourier modes.
    """

    if n_grid < 8 or n_grid % 2:
        raise ValueError("n_grid must be an even integer >= 8")
    if k_head is None:
        k_head = min(n_grid // 2 - 1, max(len(outer.modes), len(inner.modes), 20))
    values: list[mp.mpf] = []
    for j in range(n_grid):
        theta = mp.mpf(j) / n_grid
        inner_value = inner.value(theta)
        values.append(inner_value + outer.value(theta + inner_value))
    mu = mp.fsum(values) / n_grid
    modes: list[mp.mpc] = []
    for k in range(1, k_head + 1):
        coefficient = mp.fsum(
            value * mp.exp(mp.mpc(0, -2) * mp.pi * j * k / n_grid)
            for j, value in enumerate(values)
        )
        modes.append(2 * coefficient / n_grid)
    base = basechange._base_value(target_base)
    return PhaseMap(
        base=base,
        beta=mp.log(mp.log(base)),
        mu=mu,
        modes=tuple(modes),
        provenance=f"composition of [{outer.provenance}] and [{inner.provenance}]",
        metadata={"n_grid": n_grid, "k_head": k_head},
    )


def _sexp_from_phase(
    gp_anchor: Any,
    phase: PhaseMap,
    y: Any,
    n_lift: int = 8,
    k_head: int | None = None,
) -> mp.mpf:
    y = mp.mpf(y)
    theta = y - mp.floor(y)
    h = n_lift + y + phase.value(theta, k_head=k_head)
    ln_b = mp.log(phase.base)
    lnln_b = mp.log(ln_b)
    threshold = mp.mpf(10) ** (mp.mp.dps // 2 + 10)
    j = max(int(mp.ceil(h - 1 - mp.mpf("3.2"))), 1)
    anchor_height = h - 1 - j
    tower = mp.mpf(gp_anchor.sexp("exp(1)", anchor_height).real)
    while tower <= threshold and j > 1:
        tower = mp.exp(tower)
        j -= 1
    value = tower - lnln_b
    while j < n_lift:
        value = mp.log(value) - lnln_b
        j += 1
    return mp.exp(value)


def _slog_from_phase(
    gp_anchor: Any,
    phase: PhaseMap,
    w: Any,
    k_head: int | None = None,
) -> mp.mpf:
    w = mp.mpf(w)
    ln_b = mp.log(phase.base)
    value = w
    levels = 0
    while value <= basechange.TOWER_CUT:
        value = mp.exp(ln_b * value)
        levels += 1
    transformed = mp.log(ln_b) + value * ln_b
    peels = 0
    while transformed > basechange.PEEL_CUT:
        transformed = mp.log(transformed)
        peels += 1
    anchor_height = mp.mpf(gp_anchor.slog("exp(1)", transformed).real) + 2 + peels
    target = anchor_height - (levels + 2)
    x = target - phase.mu
    for _ in range(5):
        theta = x - mp.floor(x)
        residual = x + phase.value(theta, k_head=k_head) - target
        derivative = 1 + phase.derivative(theta, k_head=k_head)
        x -= residual / derivative
    return x


class ReconstructionCalculator:
    def __init__(
        self,
        gp_anchor: Any,
        atlas: AtlasIndex | None = None,
        local_marcher: LocalRHMarcher | None = None,
    ):
        self.gp_anchor = gp_anchor
        self.atlas = atlas or AtlasIndex.load()
        self.local_marcher = local_marcher

    def reconstruct_state(
        self,
        base: Any,
        digits: int = 12,
        require_anchor_pure: bool = False,
    ) -> ReconstructionState:
        target = basechange._base_value(base)
        if target <= mp.e ** (1 / mp.e):
            raise ValueError("Paper VI reconstruction requires base > e^(1/e)")
        beta = mp.log(mp.log(target))
        exact = self.atlas.exact(
            target,
            require_anchor_pure=require_anchor_pure,
        )
        if exact is not None:
            return ReconstructionState(
                target_base=target,
                target_beta=beta,
                phase=exact.phase,
                anchor=exact,
                local_beta_distance=mp.mpf(0),
                source=(
                    "stored anchor-pure atlas state"
                    if exact.anchor_pure
                    else "legacy stored target-sampled atlas state"
                ),
            )

        anchor = self.atlas.nearest(
            target,
            require_anchor_pure=require_anchor_pure,
        )
        distance = abs(beta - anchor.beta)
        if self.local_marcher is None:
            raise AtlasCoverageError(
                "target base is not stored in the atlas and no local RH "
                f"marcher was supplied; nearest anchor={mp.nstr(anchor.base, 16)}, "
                f"|delta beta|={mp.nstr(distance, 8)}. Refusing to initialize "
                "the target-base engine because that would violate anchor purity."
            )
        marched = self.local_marcher.march(anchor, target, digits)
        composed = compose_phase_maps(
            anchor.phase,
            marched.relative_phase,
            target,
            n_grid=max(128, 4 * len(marched.relative_phase.modes)),
        )
        diagnostics = dict(marched.diagnostics)
        diagnostics.update({
            "residual": marched.residual,
            "path_split_delta": marched.path_split_delta,
            "resolution_delta": marched.resolution_delta,
            "contraction_q_estimate": marched.contraction_q_estimate,
        })
        return ReconstructionState(
            target_base=target,
            target_beta=beta,
            phase=composed,
            anchor=anchor,
            local_beta_distance=distance,
            source="nearest atlas anchor plus local floating-point RH march",
            diagnostics=diagnostics,
        )

    def sexp(
        self,
        state: ReconstructionState,
        height: Any,
        n_lift: int = 8,
        k_head: int | None = None,
    ) -> mp.mpf:
        return _sexp_from_phase(
            self.gp_anchor, state.phase, height, n_lift=n_lift, k_head=k_head
        )

    def slog(
        self,
        state: ReconstructionState,
        value: Any,
        k_head: int | None = None,
    ) -> mp.mpf:
        return _slog_from_phase(self.gp_anchor, state.phase, value, k_head=k_head)

    def validate(
        self,
        state: ReconstructionState,
        heights: Sequence[Any],
        requested_digits: int,
        require_independent_gates: bool = True,
    ) -> ValidationReport:
        probes = tuple(mp.mpf(value) for value in heights) or (mp.mpf("0.5"),)
        probes = tuple(dict.fromkeys((mp.mpf(0), mp.mpf("0.25"), *probes)))
        fine = [self.sexp(state, x) for x in probes]
        coarse_k = max(0, len(state.phase.modes) - 4)
        coarse = [
            self.sexp(state, x, k_head=coarse_k)
            if coarse_k < len(state.phase.modes)
            else value
            for x, value in zip(probes, fine)
        ]
        lift = [self.sexp(state, x, n_lift=10) for x in probes]
        resolution_delta = max(abs(a - b) for a, b in zip(fine, coarse))
        lift_delta = max(abs(a - b) for a, b in zip(fine, lift))
        normalization = abs(self.sexp(state, 0) - 1)
        functional = max(
            abs(self.sexp(state, x + 1) - mp.power(state.target_base, value))
            for x, value in zip(probes, fine)
        )
        roundtrip = max(
            abs(self.slog(state, value) - x)
            for x, value in zip(probes, fine)
        )
        diagnostics = state.diagnostics
        forecast = ErrorForecast(
            mode_tail=state.phase.tail_forecast(),
            resolution=max(
                resolution_delta,
                mp.mpf(diagnostics.get("resolution_delta", 0) or 0),
            ),
            lift=lift_delta,
            normalization=normalization,
            functional_equation=functional,
            roundtrip=roundtrip,
            path_split=(
                None
                if diagnostics.get("path_split_delta") is None
                else mp.mpf(diagnostics["path_split_delta"])
            ),
            residual=(
                None
                if diagnostics.get("residual") is None
                else mp.mpf(diagnostics["residual"])
            ),
            contraction_q=(
                None
                if diagnostics.get("contraction_q_estimate") is None
                else mp.mpf(diagnostics["contraction_q_estimate"])
            ),
        )
        failed_gates: list[str] = []
        if require_independent_gates:
            if diagnostics.get("path_split_delta") is None:
                failed_gates.append("path_splitting")
            if diagnostics.get("resolution_delta") is None:
                failed_gates.append("independent_resolution")
            contraction_q = diagnostics.get("contraction_q_estimate")
            if contraction_q is None or not (0 <= mp.mpf(contraction_q) < 1):
                failed_gates.append("residual_to_error")
        passed = (
            forecast.predicted_digits >= requested_digits
            and not failed_gates
        )
        notes = (
            (
                "Validated gates were required."
                if require_independent_gates
                else "Fast mode used the available internal forecast only."
            ),
            "The forecast is an internal floating-point reliability estimate, not a proof.",
            "Certified mode replaces this forecast by an outward-rounded "
            "certificate radius.",
        )
        return ValidationReport(
            passed=passed,
            requested_digits=requested_digits,
            forecast=forecast,
            failed_gates=tuple(failed_gates),
            notes=notes,
        )

    def evaluate(
        self,
        base: Any,
        heights: Sequence[Any],
        mode: ReconstructionMode | str = ReconstructionMode.FAST,
        digits: int = 12,
        certificate_path: str | Path | None = None,
    ) -> ReconstructionResult:
        mode = ReconstructionMode(mode)
        with mp.workdps(max(50, digits + 35)):
            state = self.reconstruct_state(
                base,
                digits=digits,
                require_anchor_pure=mode is ReconstructionMode.CERTIFIED,
            )
            centers = tuple(self.sexp(state, height) for height in heights)
            validation: ValidationReport | None = self.validate(
                state,
                heights,
                digits,
                require_independent_gates=(
                    mode is not ReconstructionMode.FAST
                ),
            )
            certificate: CertificateAudit | None = None
            radius: mp.mpf | None = None

            if mode in (ReconstructionMode.VALIDATED, ReconstructionMode.CERTIFIED):
                validation.require_pass()

            if mode is ReconstructionMode.CERTIFIED:
                if certificate_path is None:
                    raise CertificateError(
                        "certified mode requires a machine-readable interval certificate"
                    )
                certificate = load_and_audit_certificate(certificate_path).require_pass()
                target_text = mp.mpf(certificate.target_base)
                if target_text != state.target_base:
                    raise CertificateError(
                        "certificate target_base does not match the reconstruction request"
                    )
                with Path(certificate_path).open(encoding="utf-8") as handle:
                    payload = json.load(handle)
                evaluation = payload["evaluation_enclosure"]
                height_min = mp.mpf(evaluation["height_min"])
                height_max = mp.mpf(evaluation["height_max"])
                outside = [
                    height for height in heights
                    if not (height_min <= mp.mpf(height) <= height_max)
                ]
                if outside:
                    raise CertificateError(
                        "requested height lies outside the certificate's "
                        f"evaluation domain [{height_min}, {height_max}]"
                    )
                if evaluation["phase_map_sha256"] != state.phase.fingerprint():
                    raise CertificateError(
                        "certificate phase_map_sha256 does not match the "
                        "phase map used for evaluation"
                    )
                radius = mp.mpf(payload["evaluation_enclosure"]["radius"])

            values = tuple(
                ReconstructionValue(mp.mpf(height), center, radius)
                for height, center in zip(heights, centers)
            )
            return ReconstructionResult(
                mode=mode,
                state=state,
                values=values,
                validation=validation,
                certificate=certificate,
            )


def _result_json(result: ReconstructionResult, digits: int) -> dict[str, Any]:
    return {
        "mode": result.mode.value,
        "base": mp.nstr(result.state.target_base, digits),
        "beta": mp.nstr(result.state.target_beta, digits),
        "anchor_base": mp.nstr(result.state.anchor.base, digits),
        "local_beta_distance": mp.nstr(result.state.local_beta_distance, digits),
        "source": result.state.source,
        "values": [
            {
                "height": mp.nstr(value.height, digits),
                "center": mp.nstr(value.center, digits),
                "radius": None if value.radius is None else mp.nstr(value.radius, digits),
            }
            for value in result.values
        ],
        "validation": (
            None
            if result.validation is None
            else {
                "passed": result.validation.passed,
                "requested_digits": result.validation.requested_digits,
                "forecast": result.validation.forecast.as_dict(),
                "failed_gates": list(result.validation.failed_gates),
                "notes": list(result.validation.notes),
            }
        ),
        "certificate": (
            None
            if result.certificate is None
            else {
                "passed": result.certificate.passed,
                "segment_count": result.certificate.segment_count,
                "global_error_bound": result.certificate.global_error_bound,
            }
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Anchor-pure tetration reconstruction calculator"
    )
    parser.add_argument("command", choices=["sexp"])
    parser.add_argument("--base", required=True)
    parser.add_argument("--values", nargs="+", required=True)
    parser.add_argument(
        "--mode",
        choices=[mode.value for mode in ReconstructionMode],
        default=ReconstructionMode.FAST.value,
    )
    parser.add_argument("--digits", type=int, default=12)
    parser.add_argument("--dps", type=int, default=60)
    parser.add_argument("--atlas-table")
    parser.add_argument("--certificate")
    parser.add_argument("--gp-exe")
    parser.add_argument("--fatou-gp", default="fork")
    args = parser.parse_args()

    from .gp_backend import FatouGP

    mp.mp.dps = max(args.dps, args.digits + 35)
    gp = FatouGP(
        dps=args.dps,
        gp_exe=args.gp_exe,
        fatou_gp=args.fatou_gp,
    )
    atlas = AtlasIndex.load(args.atlas_table)
    calculator = ReconstructionCalculator(gp_anchor=gp, atlas=atlas)
    result = calculator.evaluate(
        args.base,
        args.values,
        mode=args.mode,
        digits=args.digits,
        certificate_path=args.certificate,
    )
    print(json.dumps(_result_json(result, max(20, args.digits + 5)), indent=2))


if __name__ == "__main__":
    main()
