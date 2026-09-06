"""Machine-checkable *precertificate* audit for Paper VI reconstruction.

This module is deliberately a scalar/structural auditor, not the interval
arithmetic producer or the independent array-level replay engine.  The
producer must evaluate the RH map, admissibility margins, residuals, endpoint
gluing, and cut-plane templates with outward-rounded interval arithmetic.
The auditor verifies the resulting finite record's dependency structure,
scalar inequalities, and bound-file integrity.

No record can currently pass as a complete computer-assisted certificate:
``audit_certificate`` has an explicit fail-closed array-replay gate.  This is
intentional.  A collection of self-declared ``pass=true`` flags and hashes is
not a proof, even when every scalar inequality is internally consistent.

Only the conservative Ovsyannikov contraction profile is accepted for a
``certified`` result.  Its segment radii polynomial is

    p(r) = Y + (q - 1) r,

where Y encloses the Volterra residual and q < 1 is an interval upper bound
for the full Volterra map on the whole admissible ball.  A quadratic
Newton--Kantorovich profile may be added later, but it requires a separate
proof that D F is Lipschitz in the *same* fixed path norm; the analytic-scale
estimate alone does not imply that fact.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from decimal import (
    Inexact,
    Decimal,
    InvalidOperation,
    ROUND_CEILING,
    ROUND_FLOOR,
    localcontext,
)
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
from typing import Any, Mapping


SCHEMA_VERSION = "paper-vi-cert-v2"
THEOREM_PROFILE = "ovsyannikov-contraction-v1"
ARRAY_REPLAY_PROFILE = "paper-vi-array-replay-v1"
PAULSEN_COWGILL_DOI = "10.1007/s10444-017-9524-1"
BASE_THRESHOLD = Decimal(
    "1.4446678610097661336583391085964302230585954532422531658205226643038549"
)
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
_RH_COMPONENTS = (
    "qh_models",
    "koenigs_models",
    "geometric_quotients",
    "periodic_trace",
    "hardy_projection",
    "normalization",
    "lipschitz_and_tail",
)


class CertificateError(ValueError):
    """Raised when a certificate cannot support a certified claim."""


@dataclass(frozen=True)
class CertificateAudit:
    """Result of a complete, non-short-circuiting certificate audit."""

    passed: bool
    issues: tuple[str, ...]
    target_base: str | None = None
    target_beta: str | None = None
    global_error_bound: str | None = None
    segment_count: int = 0
    checks: Mapping[str, bool] = field(default_factory=dict)

    def require_pass(self) -> "CertificateAudit":
        if not self.passed:
            joined = "\n".join(f"- {issue}" for issue in self.issues)
            raise CertificateError(f"certificate audit failed:\n{joined}")
        return self


def _decimal(value: Any, path: str, issues: list[str]) -> Decimal | None:
    if isinstance(value, bool) or value is None:
        issues.append(f"{path}: expected a finite decimal string")
        return None
    try:
        out = Decimal(str(value))
    except (InvalidOperation, ValueError):
        issues.append(f"{path}: invalid decimal value {value!r}")
        return None
    if not out.is_finite():
        issues.append(f"{path}: value must be finite")
        return None
    return out


def _pass_flag(obj: Any, path: str, issues: list[str]) -> bool:
    if not isinstance(obj, Mapping) or obj.get("pass") is not True:
        issues.append(f"{path}: explicit pass=true is required")
        return False
    return True


def _positive_margins(obj: Any, path: str, issues: list[str]) -> bool:
    if not isinstance(obj, Mapping) or not obj:
        issues.append(f"{path}: non-empty interval margins are required")
        return False
    ok = True
    for name, raw in obj.items():
        value = _decimal(raw, f"{path}.{name}", issues)
        if value is None or value <= 0:
            if value is not None:
                issues.append(f"{path}.{name}: margin must be strictly positive")
            ok = False
    return ok


def _sha256(value: Any, path: str, issues: list[str]) -> bool:
    if not isinstance(value, str) or not _SHA256_RE.fullmatch(value):
        issues.append(f"{path}: expected lowercase SHA-256")
        return False
    return True


def _check_source_manifest(
    obj: Any,
    certificate_dir: Path | None,
    issues: list[str],
) -> bool:
    if not isinstance(obj, Mapping) or not obj:
        issues.append("source_files: at least one source file is required")
        return False
    ok = True
    root = None if certificate_dir is None else certificate_dir.resolve()
    for name, record in obj.items():
        path = f"source_files.{name}"
        if not isinstance(record, Mapping):
            issues.append(f"{path}: expected an object")
            ok = False
            continue
        relative = record.get("path")
        digest = record.get("sha256")
        if (
            not isinstance(relative, str)
            or not relative
            or PurePosixPath(relative).is_absolute()
            or ".." in PurePosixPath(relative).parts
        ):
            issues.append(f"{path}.path: expected a safe relative POSIX path")
            ok = False
        if not _sha256(digest, f"{path}.sha256", issues):
            ok = False
        if "phase_map_sha256" in record and not _sha256(
            record.get("phase_map_sha256"),
            f"{path}.phase_map_sha256",
            issues,
        ):
            ok = False
        if root is not None and isinstance(relative, str) and relative:
            candidate = (root / relative).resolve()
            try:
                candidate.relative_to(root)
            except ValueError:
                issues.append(f"{path}.path: resolves outside certificate directory")
                ok = False
                continue
            if not candidate.is_file():
                issues.append(f"{path}.path: source file is missing")
                ok = False
            elif isinstance(digest, str) and sha256_file(candidate) != digest:
                issues.append(f"{path}.sha256: source file hash mismatch")
                ok = False
    return ok


def _check_data_hash_bindings(
    payload: Mapping[str, Any],
    issues: list[str],
) -> bool:
    source_files = payload.get("source_files")
    source_digests: set[Any] = set()
    if isinstance(source_files, Mapping):
        source_digests = {
            record.get("sha256")
            for record in source_files.values()
            if isinstance(record, Mapping)
        }
    ok = True

    def visit(value: Any, path: str) -> None:
        nonlocal ok
        if isinstance(value, Mapping):
            for key, child in value.items():
                child_path = f"{path}.{key}" if path else str(key)
                if key == "data_sha256" and child not in source_digests:
                    issues.append(
                        f"{child_path}: hash is not bound to a source_files entry"
                    )
                    ok = False
                visit(child, child_path)
        elif isinstance(value, list):
            for index, child in enumerate(value):
                visit(child, f"{path}[{index}]")

    visit(payload, "")
    return ok


def _outward_record(obj: Any, path: str, issues: list[str]) -> bool:
    if not isinstance(obj, Mapping):
        issues.append(f"{path}: object is required")
        return False
    ok = _pass_flag(obj, path, issues)
    if obj.get("rounding") != "outward":
        issues.append(f"{path}.rounding: must be 'outward'")
        ok = False
    if not _sha256(obj.get("data_sha256"), f"{path}.data_sha256", issues):
        ok = False
    return ok


def _check_interval_rh(obj: Any, path: str, issues: list[str]) -> bool:
    ok = _outward_record(obj, path, issues)
    if not isinstance(obj, Mapping):
        return False
    components = obj.get("components")
    if not isinstance(components, Mapping):
        issues.append(f"{path}.components: component map is required")
        ok = False
    else:
        for name in _RH_COMPONENTS:
            if components.get(name) is not True:
                issues.append(f"{path}.components.{name}: must be true")
                ok = False
    if not _positive_margins(obj.get("margins"), f"{path}.margins", issues):
        ok = False
    for key in ("tail_upper", "lipschitz_upper"):
        value = _decimal(obj.get(key), f"{path}.{key}", issues)
        if value is None or value < 0:
            if value is not None:
                issues.append(f"{path}.{key}: must be non-negative")
            ok = False
    return ok


def _check_end_gluing(obj: Any, path: str, issues: list[str]) -> bool:
    ok = _outward_record(obj, path, issues)
    if not isinstance(obj, Mapping):
        return False
    eta = _decimal(obj.get("eta_upper"), f"{path}.eta_upper", issues)
    decay = _decimal(obj.get("koenigs_decay_rate_lower"), f"{path}.koenigs_decay_rate_lower", issues)
    margin = _decimal(obj.get("interval_newton_margin_lower"), f"{path}.interval_newton_margin_lower", issues)
    boxes = obj.get("overlap_box_count")
    if eta is None or not (Decimal(0) <= eta < Decimal(1)):
        if eta is not None:
            issues.append(f"{path}.eta_upper: must satisfy 0 <= eta < 1")
        ok = False
    for value, key in (
        (decay, "koenigs_decay_rate_lower"),
        (margin, "interval_newton_margin_lower"),
    ):
        if value is None or value <= 0:
            if value is not None:
                issues.append(f"{path}.{key}: must be strictly positive")
            ok = False
    if isinstance(boxes, bool) or not isinstance(boxes, int) or boxes < 1:
        issues.append(f"{path}.overlap_box_count: positive integer required")
        ok = False
    if obj.get("base_uniform") is not True:
        issues.append(f"{path}.base_uniform: must be true")
        ok = False
    return ok


def _check_cut_plane(obj: Any, issues: list[str]) -> bool:
    path = "cut_plane_continuation"
    ok = _outward_record(obj, path, issues)
    if not isinstance(obj, Mapping):
        return False
    if obj.get("domain") != "C minus (-infinity,-2]":
        issues.append(f"{path}.domain: exact cut-plane domain is required")
        ok = False
    for key in ("transition_box_count", "invariant_template_count", "overlap_check_count"):
        value = obj.get(key)
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            issues.append(f"{path}.{key}: positive integer required")
            ok = False
    if obj.get("induction_verified") is not True:
        issues.append(f"{path}.induction_verified: must be true")
        ok = False
    if obj.get("upper_lower_compatible") is not True:
        issues.append(f"{path}.upper_lower_compatible: must be true")
        ok = False
    if not _positive_margins(obj.get("margins"), f"{path}.margins", issues):
        ok = False
    return ok


def _decimal_exp_upper(value: Decimal) -> Decimal:
    """Echte obere Schranke fuer exp(value).

    Decimal.exp() ist laut Sprachdefinition *correctly rounded* mit
    ROUND_HALF_EVEN und ignoriert den Rundungsmodus des Kontexts -- gemessen:
    ROUND_CEILING und ROUND_FLOOR liefern Ziffer fuer Ziffer dasselbe. Die
    frueher hier gesetzte Zeile `context.rounding = ROUND_CEILING` hatte also
    keine Wirkung, und der Name versprach eine Schranke, die nur zufaellig
    stimmte: faellt das Runden nach unten, ist das Ergebnis KLEINER als
    exp(value).

    Deshalb mit Schutzstellen rechnen, dann bewusst nach oben runden und einen
    ulp aufschlagen. Der Aufschlag bei 80 Stellen deckt den Rundungsfehler bei
    100 Stellen um zwanzig Dekaden.
    """
    with localcontext() as context:
        # exp() ist auf 100 Stellen korrekt gerundet, also hoechstens einen
        # halben ulp daneben; ein ganzer ulp deckt das sicher nach oben ab.
        # Bewusst NICHT auf 80 Stellen vorrunden: der Aufrufer rechnet
        # ohnehin in einem ROUND_CEILING-Kontext, und dort verschwindet der
        # Aufschlag. Rundet man ihn hier selbst hoch, wird aus einem
        # unsichtbaren ulp bei 100 Stellen ein sichtbarer bei 80 -- fuer
        # exp(0) = 1, das exakt ist, sogar der einzige Unterschied.
        context.prec = 100
        context.clear_flags()
        raw = value.exp()
        # NUR aufschlagen, wenn wirklich gerundet wurde. exp(0) = 1 ist exakt,
        # und ein ulp darauf ist nicht bloss unnoetig, sondern schaedlich: bei
        # einem einzigen Segment ist der Gronwall-Exponent die leere Summe,
        # also exp(0), und der Aufschlag laesst das Budget die Schranke um
        # einen ulp uebersteigen, obwohl beide gleich sind.
        return raw.next_plus() if context.flags[Inexact] else raw


def _condition_set(obj: Any, issues: list[str]) -> bool:
    """Check the exact hypotheses of Paulsen--Cowgill Proposition 2."""

    required = (
        "analytic_on_half_plane_Re_gt_minus_2",
        "functional_equation_on_half_plane",
        "normalization_F0_eq_1",
        "upper_fixed_point_limit_for_every_x_gt_minus_2",
        "lower_fixed_point_limit_for_every_x_gt_minus_2",
    )
    if not isinstance(obj, Mapping):
        issues.append("static_uniqueness.conditions: condition map is required")
        return False
    ok = True
    for key in required:
        if obj.get(key) is not True:
            issues.append(
                "static_uniqueness.conditions."
                f"{key}: exact Paulsen--Cowgill hypothesis must be true"
            )
            ok = False
    return ok


def _check_cut_plane_bridge(obj: Any, issues: list[str]) -> bool:
    """Check the explicit bridge from Proposition 2 to the cut plane.

    Proposition 2 itself is a comparison theorem on ``Re(z) > -2``.  The
    reconstructed and Kneser branches are first identified there.  Equality
    on the connected cut plane is a separate identity-theorem step.
    """

    path = "static_uniqueness.cut_plane_bridge"
    if not isinstance(obj, Mapping):
        issues.append(f"{path}: object is required")
        return False
    ok = True
    expected = {
        "source_domain": "C minus (-infinity,-2]",
        "comparison_domain": "Re(z)>-2",
        "restriction_verified": True,
        "cut_plane_connected": True,
        "identity_theorem_extension": True,
    }
    for key, value in expected.items():
        if obj.get(key) != value:
            issues.append(f"{path}.{key}: expected {value!r}")
            ok = False
    return ok


def _check_array_replay_gate(obj: Any, issues: list[str]) -> bool:
    """Fail closed until the independent interval-array replay exists."""

    path = "array_replay"
    if not isinstance(obj, Mapping):
        issues.append(f"{path}: object is required")
    else:
        if obj.get("profile") != ARRAY_REPLAY_PROFILE:
            issues.append(
                f"{path}.profile: expected {ARRAY_REPLAY_PROFILE!r}"
            )
        if obj.get("status") != "not-implemented":
            issues.append(
                f"{path}.status: this release only accepts the explicit "
                "'not-implemented' marker"
            )
        source_name = obj.get("source_file")
        if source_name not in (None, "REQUIRED_REPLAY_TRACE_SOURCE"):
            issues.append(
                f"{path}.source_file: a replay trace cannot be trusted "
                "before the independent replay implementation exists"
            )
    issues.append(
        "array_replay: independent array-level interval replay is not "
        "implemented in this release; certified status cannot be issued"
    )
    return False


def audit_certificate(
    payload: Mapping[str, Any],
    certificate_dir: Path | None = None,
) -> CertificateAudit:
    """Audit a precertificate and refuse a final certified claim.

    Component pass flags are checked together with their required margins,
    metadata, scalar reductions, and bound-file hashes.  The final
    ``array_replay`` check is deliberately false until an independent module
    regenerates every component claim from the interval arrays.
    """

    issues: list[str] = []
    checks: dict[str, bool] = {}

    if payload.get("schema_version") != SCHEMA_VERSION:
        issues.append(
            f"schema_version: expected {SCHEMA_VERSION!r}, "
            f"got {payload.get('schema_version')!r}"
        )
    if payload.get("theorem_profile") != THEOREM_PROFILE:
        issues.append(
            "theorem_profile: certified mode accepts only "
            f"{THEOREM_PROFILE!r}"
        )
    if payload.get("status") != "complete":
        issues.append("status: must be 'complete'")

    target_base = payload.get("target_base")
    target_beta = payload.get("target_beta")
    base_dec = _decimal(target_base, "target_base", issues)
    beta_dec = _decimal(target_beta, "target_beta", issues)
    beta_tolerance = _decimal(
        payload.get("target_beta_tolerance"),
        "target_beta_tolerance",
        issues,
    )
    if base_dec is not None and base_dec <= BASE_THRESHOLD:
        issues.append("target_base: Paper VI requires base > e^(1/e)")
    parameter_ok = (
        beta_tolerance is not None
        and Decimal(0) < beta_tolerance <= Decimal("1e-30")
    )
    if beta_tolerance is not None and not parameter_ok:
        issues.append(
            "target_beta_tolerance: must satisfy 0 < tolerance <= 1e-30"
        )
    if (
        base_dec is not None
        and base_dec > 0
        and beta_dec is not None
        and parameter_ok
    ):
        # base_dec > 0 reicht NICHT: fuer 0 < base <= 1 ist ln(base) <= 0 und
        # das zweite ln() wirft InvalidOperation. Das entkam der
        # CertificateError-Konvention und brach ein Audit, das sich als
        # vollstaendig ausgibt, mitten in der Pruefliste ab.
        if base_dec <= 1:
            issues.append(
                "target_base: log(log(target_base)) is not real for "
                "0 < base <= 1; this regime is outside the certificate"
            )
            parameter_ok = False
            expected_beta = None
        else:
            with localcontext() as context:
                context.prec = 90
                expected_beta = base_dec.ln().ln()
        if expected_beta is not None and abs(beta_dec - expected_beta) > beta_tolerance:
            issues.append(
                "target_beta: inconsistent with log(log(target_base))"
            )
            parameter_ok = False
    checks["target_parameter"] = parameter_ok

    backend = payload.get("interval_backend")
    backend_ok = isinstance(backend, Mapping)
    if not backend_ok:
        issues.append("interval_backend: object is required")
    else:
        assert isinstance(backend, Mapping)
        if backend.get("directed_rounding") is not True:
            issues.append("interval_backend.directed_rounding: must be true")
            backend_ok = False
        for key in ("name", "version"):
            if not isinstance(backend.get(key), str) or not backend.get(key):
                issues.append(f"interval_backend.{key}: non-empty string required")
                backend_ok = False
    checks["interval_backend"] = backend_ok

    provenance = payload.get("anchor_pure_provenance")
    provenance_ok = isinstance(provenance, Mapping)
    if not provenance_ok:
        issues.append("anchor_pure_provenance: object is required")
    else:
        assert isinstance(provenance, Mapping)
        for key, expected in (
            ("pass", True),
            ("initial_anchor", "e"),
            ("target_engine_data_used", False),
        ):
            if provenance.get(key) != expected:
                issues.append(
                    f"anchor_pure_provenance.{key}: expected {expected!r}"
                )
                provenance_ok = False
        if not _sha256(
            provenance.get("lineage_sha256"),
            "anchor_pure_provenance.lineage_sha256",
            issues,
        ):
            provenance_ok = False
    checks["anchor_pure_provenance"] = provenance_ok

    checks["source_files"] = _check_source_manifest(
        payload.get("source_files"), certificate_dir, issues
    )
    checks["data_hash_bindings"] = _check_data_hash_bindings(payload, issues)
    checks["anchor"] = _outward_record(
        payload.get("anchor_certificate"), "anchor_certificate", issues
    )
    if isinstance(payload.get("anchor_certificate"), Mapping):
        anchor = payload["anchor_certificate"]
        if not _positive_margins(
            anchor.get("margins"), "anchor_certificate.margins", issues
        ):
            checks["anchor"] = False
        remainder = _decimal(
            anchor.get("remainder_upper"),
            "anchor_certificate.remainder_upper",
            issues,
        )
        if remainder is None or remainder < 0:
            if remainder is not None:
                issues.append(
                    "anchor_certificate.remainder_upper: must be non-negative"
                )
            checks["anchor"] = False
    checks["interval_rh"] = _check_interval_rh(
        payload.get("interval_rh"), "interval_rh", issues
    )

    segments = payload.get("segments")
    if not isinstance(segments, list) or not segments:
        issues.append("segments: a non-empty ordered segment list is required")
        segments = []

    previous_beta_b: Decimal | None = None
    budget_rows: list[tuple[Decimal, Decimal, Decimal]] = []
    segment_checks = True
    for index, segment in enumerate(segments):
        path = f"segments[{index}]"
        if not isinstance(segment, Mapping):
            issues.append(f"{path}: expected an object")
            segment_checks = False
            continue

        beta_a = _decimal(segment.get("beta_a"), f"{path}.beta_a", issues)
        beta_b = _decimal(segment.get("beta_b"), f"{path}.beta_b", issues)
        if beta_a is not None and beta_b is not None and beta_a == beta_b:
            issues.append(f"{path}: zero-length segment")
            segment_checks = False
        if index == 0 and beta_a is not None and beta_a != 0:
            issues.append(f"{path}.beta_a: first segment must start at beta=0")
            segment_checks = False
        if previous_beta_b is not None and beta_a is not None and beta_a != previous_beta_b:
            issues.append(
                f"{path}.beta_a: does not glue exactly to the previous beta_b"
            )
            segment_checks = False
        if beta_b is not None:
            previous_beta_b = beta_b

        delta = None if beta_a is None or beta_b is None else beta_b - beta_a
        if delta is not None and beta_dec is not None and delta * beta_dec <= 0:
            issues.append(f"{path}: segment direction is inconsistent with target_beta")
            segment_checks = False

        a_upper = _decimal(segment.get("a_upper"), f"{path}.a_upper", issues)
        rho0 = _decimal(segment.get("rho0_lower"), f"{path}.rho0_lower", issues)
        kappa_lo = _decimal(
            segment.get("kappa_lower"), f"{path}.kappa_lower", issues
        )
        kappa_hi = _decimal(
            segment.get("kappa_upper"), f"{path}.kappa_upper", issues
        )
        rho_end = _decimal(
            segment.get("rho_end_upper"), f"{path}.rho_end_upper", issues
        )
        l_upper = _decimal(
            segment.get("L_upper"), f"{path}.L_upper", issues
        )
        for value, key in (
            (a_upper, "a_upper"),
            (rho0, "rho0_lower"),
            (kappa_lo, "kappa_lower"),
            (kappa_hi, "kappa_upper"),
            (rho_end, "rho_end_upper"),
        ):
            if value is None or value <= 0:
                if value is not None:
                    issues.append(f"{path}.{key}: must be strictly positive")
                segment_checks = False
        if (
            delta is not None
            and a_upper is not None
            and a_upper < abs(delta)
        ):
            issues.append(f"{path}.a_upper: must enclose |beta_b-beta_a|")
            segment_checks = False
        if (
            kappa_lo is not None
            and kappa_hi is not None
            and kappa_lo > kappa_hi
        ):
            issues.append(f"{path}: kappa_lower must not exceed kappa_upper")
            segment_checks = False

        y_head = _decimal(
            segment.get("Y_head_upper"), f"{path}.Y_head_upper", issues
        )
        y_tail = _decimal(
            segment.get("Y_tail_upper"), f"{path}.Y_tail_upper", issues
        )
        y = _decimal(segment.get("Y_upper"), f"{path}.Y_upper", issues)
        q = _decimal(segment.get("q_upper"), f"{path}.q_upper", issues)
        radius = _decimal(segment.get("radius"), f"{path}.radius", issues)
        path_endpoint_error = _decimal(
            segment.get("path_endpoint_error_upper"),
            f"{path}.path_endpoint_error_upper",
            issues,
        )
        q_error = _decimal(
            segment.get("q_error_upper"),
            f"{path}.q_error_upper",
            issues,
        )
        representation_error = _decimal(
            segment.get("representation_error_upper"),
            f"{path}.representation_error_upper",
            issues,
        )
        endpoint_error = _decimal(
            segment.get("endpoint_error_upper"),
            f"{path}.endpoint_error_upper",
            issues,
        )
        stored_p = _decimal(
            segment.get("radii_polynomial_upper"),
            f"{path}.radii_polynomial_upper",
            issues,
        )
        for value, key in (
            (y_head, "Y_head_upper"),
            (y_tail, "Y_tail_upper"),
            (y, "Y_upper"),
            (l_upper, "L_upper"),
        ):
            if value is None or value < 0:
                if value is not None:
                    issues.append(f"{path}.{key}: must be non-negative")
                segment_checks = False
        if (
            y_head is not None
            and y_tail is not None
            and y is not None
            and y < y_head + y_tail
        ):
            issues.append(
                f"{path}.Y_upper: must enclose Y_head_upper + Y_tail_upper"
            )
            segment_checks = False
        if q is not None and not (Decimal(0) <= q < Decimal(1)):
            issues.append(f"{path}.q_upper: must satisfy 0 <= q < 1")
            segment_checks = False
        if radius is not None and radius <= 0:
            issues.append(f"{path}.radius: must be strictly positive")
            segment_checks = False
        if path_endpoint_error is not None and path_endpoint_error <= 0:
            issues.append(
                f"{path}.path_endpoint_error_upper: must be positive"
            )
            segment_checks = False
        for value, key in (
            (q_error, "q_error_upper"),
            (representation_error, "representation_error_upper"),
        ):
            if value is None or value < 0:
                if value is not None:
                    issues.append(f"{path}.{key}: must be non-negative")
                segment_checks = False
        if endpoint_error is not None and endpoint_error <= 0:
            issues.append(f"{path}.endpoint_error_upper: must be positive")
            segment_checks = False

        if None not in (q, l_upper, kappa_lo):
            assert q is not None and l_upper is not None and kappa_lo is not None
            with localcontext() as context:
                # sqrt() ignoriert den Rundungsmodus genauso wie exp(); der
                # Faktor wird daher mit Schutzstellen gebildet und erst danach
                # nach oben gerundet, damit required_q wirklich eine obere
                # Schranke ist -- sonst faellt das Gate zu grosszuegig aus.
                context.prec = 100
                context.clear_flags()
                root_two = Decimal(2).sqrt()
                if context.flags[Inexact]:
                    root_two = root_two.next_plus()
            with localcontext() as context:
                context.prec = 80
                context.rounding = ROUND_CEILING
                required_q = Decimal(4) * root_two * l_upper / kappa_lo
            if q < required_q:
                issues.append(
                    f"{path}.q_upper: smaller than 4*sqrt(2)*L_upper/kappa_lower"
                )
                segment_checks = False

        strip_slack: Decimal | None = None
        if None not in (rho0, kappa_hi, a_upper, rho_end):
            assert rho0 is not None and kappa_hi is not None
            assert a_upper is not None and rho_end is not None
            with localcontext() as context:
                context.prec = 80
                context.rounding = ROUND_FLOOR
                strip_slack = rho0 - kappa_hi * a_upper - rho_end
            if strip_slack <= 0:
                issues.append(
                    f"{path}: endpoint strip slack rho0-kappa*a-rho_end "
                    "must be positive"
                )
                segment_checks = False
        if (
            strip_slack is not None
            and strip_slack > 0
            and radius is not None
            and path_endpoint_error is not None
        ):
            with localcontext() as context:
                # untere Schranke: sqrt() rundet half-even, also mit
                # Schutzstellen rechnen und einen ulp abziehen.
                context.prec = 100
                context.clear_flags()
                root_raw = strip_slack.sqrt()
                if context.flags[Inexact]:
                    root_raw = root_raw.next_minus()
            with localcontext() as context:
                context.prec = 80
                context.rounding = ROUND_FLOOR
                root_slack = +root_raw
            with localcontext() as context:
                context.prec = 80
                context.rounding = ROUND_CEILING
                required_endpoint_error = radius / root_slack
            if path_endpoint_error < required_endpoint_error:
                issues.append(
                    f"{path}.path_endpoint_error_upper: smaller than "
                    "radius/sqrt(endpoint strip slack)"
                )
                segment_checks = False
        if None not in (
            path_endpoint_error,
            q_error,
            representation_error,
            endpoint_error,
        ):
            assert path_endpoint_error is not None
            assert q_error is not None and representation_error is not None
            assert endpoint_error is not None
            with localcontext() as context:
                context.prec = 80
                context.rounding = ROUND_CEILING
                required_complete_endpoint = (
                    path_endpoint_error + q_error + representation_error
                )
            if endpoint_error < required_complete_endpoint:
                issues.append(
                    f"{path}.endpoint_error_upper: must enclose path, "
                    "raw-normalization, and representation errors"
                )
                segment_checks = False

        if None not in (y, q, radius, stored_p):
            assert y is not None and q is not None
            assert radius is not None and stored_p is not None
            with localcontext() as context:
                context.prec = 80
                context.rounding = ROUND_CEILING
                recomputed = y + (q - Decimal(1)) * radius
            if stored_p < recomputed:
                issues.append(
                    f"{path}.radii_polynomial_upper: is smaller than the "
                    "value recomputed from its declared upper bounds"
                )
                segment_checks = False
            if stored_p >= 0:
                issues.append(
                    f"{path}.radii_polynomial_upper: must be strictly negative"
                )
                segment_checks = False

        if not _positive_margins(
            segment.get("admissibility_margins"),
            f"{path}.admissibility_margins",
            issues,
        ):
            segment_checks = False
        if not _check_interval_rh(
            segment.get("interval_rh"), f"{path}.interval_rh", issues
        ):
            segment_checks = False

        endpoint = segment.get("endpoint_ball")
        endpoint_ok = _outward_record(endpoint, f"{path}.endpoint_ball", issues)
        if isinstance(endpoint, Mapping):
            declared_endpoint = _decimal(
                endpoint.get("error_upper"),
                f"{path}.endpoint_ball.error_upper",
                issues,
            )
            if (
                declared_endpoint is None
                or endpoint_error is None
                or declared_endpoint < endpoint_error
            ):
                if declared_endpoint is not None and endpoint_error is not None:
                    issues.append(
                        f"{path}.endpoint_ball.error_upper: must enclose "
                        "endpoint_error_upper"
                    )
                endpoint_ok = False
            if index < len(segments) - 1:
                center_delta = _decimal(
                    endpoint.get("center_distance_upper"),
                    f"{path}.endpoint_ball.center_distance_upper",
                    issues,
                )
                next_radius = _decimal(
                    endpoint.get("next_start_radius_lower"),
                    f"{path}.endpoint_ball.next_start_radius_lower",
                    issues,
                )
                inclusion_margin = _decimal(
                    endpoint.get("inclusion_margin_lower"),
                    f"{path}.endpoint_ball.inclusion_margin_lower",
                    issues,
                )
                if None not in (
                    declared_endpoint,
                    center_delta,
                    next_radius,
                    inclusion_margin,
                ):
                    assert declared_endpoint is not None
                    assert center_delta is not None
                    assert next_radius is not None
                    assert inclusion_margin is not None
                    if (
                        center_delta < 0
                        or next_radius <= 0
                        or inclusion_margin <= 0
                        or declared_endpoint + center_delta + inclusion_margin
                        > next_radius
                    ):
                        issues.append(
                            f"{path}.endpoint_ball: strict full-ball handoff "
                            "inclusion failed"
                        )
                        endpoint_ok = False
            elif endpoint.get("is_target") is not True:
                issues.append(f"{path}.endpoint_ball.is_target: must be true")
                endpoint_ok = False
        if not endpoint_ok:
            segment_checks = False

        if not _check_end_gluing(
            segment.get("segment_end_gluing"),
            f"{path}.segment_end_gluing",
            issues,
        ):
            segment_checks = False

        m1 = _decimal(segment.get("M1_upper"), f"{path}.M1_upper", issues)
        if m1 is not None and m1 < 0:
            issues.append(f"{path}.M1_upper: must be non-negative")
            segment_checks = False
        if (
            endpoint_error is not None
            and endpoint_error > 0
            and m1 is not None
            and m1 >= 0
            and delta is not None
        ):
            budget_rows.append((endpoint_error, m1, abs(delta)))

    if segments and beta_dec is not None and previous_beta_b != beta_dec:
        issues.append("segments: final beta_b must equal target_beta exactly")
        segment_checks = False
    checks["segments"] = segment_checks and bool(segments)

    global_end = payload.get("global_end_gluing")
    global_end_ok = _outward_record(
        global_end, "global_end_gluing", issues
    )
    if isinstance(global_end, Mapping):
        if global_end.get("segment_count") != len(segments):
            issues.append(
                "global_end_gluing.segment_count: must equal segment count"
            )
            global_end_ok = False
        margin = _decimal(
            global_end.get("compatibility_margin_lower"),
            "global_end_gluing.compatibility_margin_lower",
            issues,
        )
        if margin is None or margin <= 0:
            if margin is not None:
                issues.append(
                    "global_end_gluing.compatibility_margin_lower: "
                    "must be strictly positive"
                )
            global_end_ok = False
    checks["global_end_gluing"] = global_end_ok
    checks["cut_plane_continuation"] = _check_cut_plane(
        payload.get("cut_plane_continuation"), issues
    )

    uniqueness = payload.get("static_uniqueness")
    uniqueness_ok = isinstance(uniqueness, Mapping)
    if not uniqueness_ok:
        issues.append("static_uniqueness: object is required")
    else:
        assert isinstance(uniqueness, Mapping)
        doi = str(uniqueness.get("doi", "")).lower()
        if doi != PAULSEN_COWGILL_DOI:
            issues.append(
                "static_uniqueness.doi: must cite Paulsen--Cowgill 2017 "
                f"({PAULSEN_COWGILL_DOI})"
            )
            uniqueness_ok = False
        if uniqueness.get("proposition") != 2:
            issues.append("static_uniqueness.proposition: must equal 2")
            uniqueness_ok = False
        if uniqueness.get("comparison_domain") != "Re(z)>-2":
            issues.append(
                "static_uniqueness.comparison_domain: Proposition 2 "
                "requires 'Re(z)>-2'"
            )
            uniqueness_ok = False
        fixed_points = uniqueness.get("fixed_points")
        if not isinstance(fixed_points, Mapping):
            issues.append("static_uniqueness.fixed_points: object is required")
            uniqueness_ok = False
        else:
            for key, expected in (
                ("upper", "L1"),
                ("lower", "L2"),
                ("conjugate_pair_for_real_base", True),
            ):
                if fixed_points.get(key) != expected:
                    issues.append(
                        f"static_uniqueness.fixed_points.{key}: "
                        f"expected {expected!r}"
                    )
                    uniqueness_ok = False
        if not _condition_set(uniqueness.get("conditions"), issues):
            uniqueness_ok = False
        if not _check_cut_plane_bridge(
            uniqueness.get("cut_plane_bridge"), issues
        ):
            uniqueness_ok = False
    checks["static_uniqueness"] = uniqueness_ok

    global_bound = _decimal(
        payload.get("global_error_bound"), "global_error_bound", issues
    )
    if global_bound is not None and global_bound <= 0:
        issues.append("global_error_bound: must be strictly positive")
    computed_global: Decimal | None = None
    if len(budget_rows) == len(segments) and budget_rows:
        with localcontext() as context:
            context.prec = 80
            context.rounding = ROUND_CEILING
            computed_global = Decimal(0)
            for index, (endpoint_error, _, _) in enumerate(budget_rows):
                exponent = sum(
                    (
                        m1 * delta
                        for _, m1, delta in budget_rows[index + 1:]
                    ),
                    Decimal(0),
                )
                computed_global += endpoint_error * _decimal_exp_upper(exponent)
        if global_bound is not None and global_bound < computed_global:
            issues.append(
                "global_error_bound: smaller than recomputed Gronwall budget"
            )
    checks["global_error_bound"] = (
        global_bound is not None
        and global_bound > 0
        and computed_global is not None
        and global_bound >= computed_global
    )

    evaluation = payload.get("evaluation_enclosure")
    eval_ok = _outward_record(
        evaluation, "evaluation_enclosure", issues
    )
    if isinstance(evaluation, Mapping):
        phase_fingerprint = evaluation.get("phase_map_sha256")
        if not _sha256(
            phase_fingerprint,
            "evaluation_enclosure.phase_map_sha256",
            issues,
        ):
            eval_ok = False
        phase_source = evaluation.get("phase_map_source")
        source_files = payload.get("source_files")
        source_record = (
            source_files.get(phase_source)
            if isinstance(source_files, Mapping) and isinstance(phase_source, str)
            else None
        )
        if (
            not isinstance(phase_source, str)
            or not isinstance(source_record, Mapping)
            or source_record.get("phase_map_sha256") != phase_fingerprint
        ):
            issues.append(
                "evaluation_enclosure.phase_map_source: must name a "
                "source_files record carrying the same phase_map_sha256"
            )
            eval_ok = False
        radius = _decimal(
            evaluation.get("radius"), "evaluation_enclosure.radius", issues
        )
        height_min = _decimal(
            evaluation.get("height_min"),
            "evaluation_enclosure.height_min",
            issues,
        )
        height_max = _decimal(
            evaluation.get("height_max"),
            "evaluation_enclosure.height_max",
            issues,
        )
        factor = _decimal(
            evaluation.get("phase_to_value_lipschitz_upper"),
            "evaluation_enclosure.phase_to_value_lipschitz_upper",
            issues,
        )
        raw_error = _decimal(
            evaluation.get("raw_branch_error_upper"),
            "evaluation_enclosure.raw_branch_error_upper",
            issues,
        )
        floating_error = _decimal(
            evaluation.get("floating_evaluation_error_upper"),
            "evaluation_enclosure.floating_evaluation_error_upper",
            issues,
        )
        if radius is None or radius <= 0:
            if radius is not None:
                issues.append("evaluation_enclosure.radius: must be positive")
            eval_ok = False
        if (
            height_min is None
            or height_max is None
            or height_min > height_max
        ):
            if height_min is not None and height_max is not None:
                issues.append(
                    "evaluation_enclosure: height_min must not exceed height_max"
                )
            eval_ok = False
        if (
            factor is None
            or factor < 0
            or raw_error is None
            or raw_error < 0
            or floating_error is None
            or floating_error < 0
        ):
            if factor is not None and factor < 0:
                issues.append(
                    "evaluation_enclosure.phase_to_value_lipschitz_upper: "
                    "must be non-negative"
                )
            if raw_error is not None and raw_error < 0:
                issues.append(
                    "evaluation_enclosure.raw_branch_error_upper: "
                    "must be non-negative"
                )
            if floating_error is not None and floating_error < 0:
                issues.append(
                    "evaluation_enclosure.floating_evaluation_error_upper: "
                    "must be non-negative"
                )
            eval_ok = False
        if None not in (
            radius,
            factor,
            raw_error,
            floating_error,
            global_bound,
        ):
            assert radius is not None and factor is not None
            assert raw_error is not None and floating_error is not None
            assert global_bound is not None
            with localcontext() as context:
                context.prec = 80
                context.rounding = ROUND_CEILING
                required_radius = (
                    factor * global_bound + raw_error + floating_error
                )
            if radius < required_radius:
                issues.append(
                    "evaluation_enclosure.radius: smaller than propagated "
                    "phase and raw-branch error"
                )
                eval_ok = False
    checks["evaluation_enclosure"] = eval_ok

    checks["array_replay"] = _check_array_replay_gate(
        payload.get("array_replay"), issues
    )

    passed = not issues and all(checks.values())
    return CertificateAudit(
        passed=passed,
        issues=tuple(issues),
        target_base=None if target_base is None else str(target_base),
        target_beta=None if target_beta is None else str(target_beta),
        global_error_bound=None if global_bound is None else str(global_bound),
        segment_count=len(segments),
        checks=checks,
    )


def load_and_audit_certificate(path: str | Path) -> CertificateAudit:
    certificate_path = Path(path)
    with certificate_path.open(encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, Mapping):
        raise CertificateError("certificate root must be a JSON object")
    return audit_certificate(payload, certificate_dir=certificate_path.parent)


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Audit a Paper VI precertificate; final certified status is "
            "blocked until independent array replay is implemented"
        )
    )
    parser.add_argument("certificate")
    args = parser.parse_args()
    audit = load_and_audit_certificate(args.certificate)
    print(json.dumps({
        "claim_level": "precertificate-only",
        "certified_available": False,
        "passed": audit.passed,
        "issues": list(audit.issues),
        "target_base": audit.target_base,
        "target_beta": audit.target_beta,
        "global_error_bound": audit.global_error_bound,
        "segment_count": audit.segment_count,
        "checks": dict(audit.checks),
    }, indent=2))
    return 0 if audit.passed else 2


if __name__ == "__main__":
    raise SystemExit(main())
