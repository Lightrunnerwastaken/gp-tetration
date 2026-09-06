"""Pure-Python tests for the Paper VI reconstruction control layer."""
from __future__ import annotations

import json
import hashlib
from pathlib import Path
import sys
import tempfile
import unittest

import mpmath as mp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from fatou_backend.certification import audit_certificate, load_and_audit_certificate
from fatou_backend.reconstruction import (
    AtlasIndex,
    AtlasCoverageError,
    PhaseMap,
    compose_phase_maps,
)


class AtlasTests(unittest.TestCase):
    def test_shipped_atlas_is_read_only_and_has_exact_anchors(self):
        atlas = AtlasIndex.load()
        self.assertIsNotNone(atlas.exact(2))
        self.assertIsNotNone(atlas.exact(3))
        self.assertIsNotNone(atlas.exact(mp.e))
        self.assertIsNone(atlas.exact(2, require_anchor_pure=True))
        self.assertEqual(atlas.nearest(2, require_anchor_pure=True).key, "e")
        nearest = atlas.nearest(mp.mpf("2.05"))
        self.assertEqual(nearest.key, "2")

    def test_phase_composition_identity(self):
        with mp.workdps(70):
            phase = PhaseMap(
                base=mp.mpf(2),
                beta=mp.log(mp.log(2)),
                mu=mp.mpf("0.1"),
                modes=(mp.mpc("1e-3", "-2e-4"),),
                provenance="test",
            )
            identity = PhaseMap(
                base=mp.e,
                beta=mp.mpf(0),
                mu=mp.mpf(0),
                modes=(),
                provenance="identity",
            )
            composed = compose_phase_maps(identity, phase, 2, n_grid=64, k_head=8)
            for x in ("0.03", "0.41", "0.88"):
                self.assertLess(
                    abs(composed.value(x) - phase.value(x)),
                    mp.mpf("1e-60"),
                )

    def test_phase_fingerprint_binds_coefficients(self):
        first = PhaseMap(
            base=mp.mpf(2),
            beta=mp.log(mp.log(2)),
            mu=mp.mpf("0.1"),
            modes=(mp.mpc("1e-3", "2e-4"),),
            provenance="test",
        )
        same = PhaseMap(
            base=first.base,
            beta=first.beta,
            mu=first.mu,
            modes=first.modes,
            provenance="different non-mathematical label",
        )
        changed = PhaseMap(
            base=first.base,
            beta=first.beta,
            mu=first.mu,
            modes=(mp.mpc("1.0001e-3", "2e-4"),),
            provenance="test",
        )
        self.assertEqual(first.fingerprint(), same.fingerprint())
        self.assertNotEqual(first.fingerprint(), changed.fingerprint())


def _interval_rh_record() -> dict:
    return {
        "pass": True,
        "rounding": "outward",
        "data_sha256": "a" * 64,
        "components": {
            "qh_models": True,
            "koenigs_models": True,
            "geometric_quotients": True,
            "periodic_trace": True,
            "hardy_projection": True,
            "normalization": True,
            "lipschitz_and_tail": True,
        },
        "margins": {
            "T_prime": "0.8",
            "koenigs_overlap": "0.1",
            "fourier_tail": "1e-30",
        },
        "tail_upper": "1e-30",
        "lipschitz_upper": "0.01",
    }


def _end_gluing_record() -> dict:
    return {
        "pass": True,
        "rounding": "outward",
        "data_sha256": "a" * 64,
        "eta_upper": "0.01",
        "koenigs_decay_rate_lower": "0.5",
        "interval_newton_margin_lower": "0.02",
        "overlap_box_count": 8,
        "base_uniform": True,
    }


def _valid_certificate() -> dict:
    with mp.workdps(80):
        beta = str(mp.log(mp.log(2)))
    return {
        "schema_version": "paper-vi-cert-v2",
        "theorem_profile": "ovsyannikov-contraction-v1",
        "status": "complete",
        "target_base": "2",
        "target_beta": beta,
        "target_beta_tolerance": "1e-70",
        "interval_backend": {
            "name": "test-interval",
            "version": "1",
            "directed_rounding": True,
        },
        "anchor_pure_provenance": {
            "pass": True,
            "initial_anchor": "e",
            "target_engine_data_used": False,
            "lineage_sha256": "d" * 64,
        },
        "source_files": {
            "candidate": {
                "path": "data/candidate.bin",
                "sha256": "a" * 64,
                "phase_map_sha256": "4" * 64,
            },
        },
        "anchor_certificate": {
            "pass": True,
            "rounding": "outward",
            "data_sha256": "a" * 64,
            "remainder_upper": "1e-40",
            "margins": {
                "T_prime": "1",
                "transfer_cut": "0.1",
            },
        },
        "interval_rh": _interval_rh_record(),
        "segments": [{
            "beta_a": "0",
            "beta_b": beta,
            "a_upper": "0.37",
            "rho0_lower": "1",
            "kappa_lower": "1",
            "kappa_upper": "1",
            "rho_end_upper": "0.5",
            "L_upper": "0.01",
            "Y_head_upper": "5e-15",
            "Y_tail_upper": "5e-15",
            "Y_upper": "1e-14",
            "q_upper": "0.1",
            "radius": "2e-14",
            "path_endpoint_error_upper": "6e-14",
            "q_error_upper": "1e-15",
            "representation_error_upper": "1e-15",
            "endpoint_error_upper": "6.2e-14",
            "radii_polynomial_upper": "-8e-15",
            "admissibility_margins": {
                "strip": "0.2",
                "derivative": "0.8",
                "cut": "0.1",
            },
            "interval_rh": _interval_rh_record(),
            "endpoint_ball": {
                "pass": True,
                "rounding": "outward",
                "data_sha256": "a" * 64,
                "error_upper": "6.2e-14",
                "is_target": True,
            },
            "segment_end_gluing": _end_gluing_record(),
            "M1_upper": "0.01",
        }],
        "global_end_gluing": {
            "pass": True,
            "rounding": "outward",
            "data_sha256": "a" * 64,
            "segment_count": 1,
            "compatibility_margin_lower": "0.01",
        },
        "cut_plane_continuation": {
            "pass": True,
            "rounding": "outward",
            "data_sha256": "a" * 64,
            "domain": "C minus (-infinity,-2]",
            "transition_box_count": 10,
            "invariant_template_count": 3,
            "overlap_check_count": 12,
            "induction_verified": True,
            "upper_lower_compatible": True,
            "margins": {
                "cut": "0.01",
                "zero": "0.02",
            },
        },
        "static_uniqueness": {
            "doi": "10.1007/s10444-017-9524-1",
            "proposition": 2,
            "comparison_domain": "Re(z)>-2",
            "fixed_points": {
                "upper": "L1",
                "lower": "L2",
                "conjugate_pair_for_real_base": True,
            },
            "conditions": {
                "analytic_on_half_plane_Re_gt_minus_2": True,
                "functional_equation_on_half_plane": True,
                "normalization_F0_eq_1": True,
                "upper_fixed_point_limit_for_every_x_gt_minus_2": True,
                "lower_fixed_point_limit_for_every_x_gt_minus_2": True,
            },
            "cut_plane_bridge": {
                "source_domain": "C minus (-infinity,-2]",
                "comparison_domain": "Re(z)>-2",
                "restriction_verified": True,
                "cut_plane_connected": True,
                "identity_theorem_extension": True,
            },
        },
        "array_replay": {
            "profile": "paper-vi-array-replay-v1",
            "status": "not-implemented",
            "source_file": "REQUIRED_REPLAY_TRACE_SOURCE",
        },
        "global_error_bound": "6.2e-14",
        "evaluation_enclosure": {
            "pass": True,
            "rounding": "outward",
            "data_sha256": "a" * 64,
            "phase_map_sha256": "4" * 64,
            "phase_map_source": "candidate",
            "height_min": "-1",
            "height_max": "2",
            "phase_to_value_lipschitz_upper": "2",
            "raw_branch_error_upper": "1e-14",
            "floating_evaluation_error_upper": "1e-15",
            "radius": "1.4e-13",
        },
    }


def _key_tree(value):
    """Nested key structure of a JSON object; leaves collapse to None.

    Lists are compared element-wise, so a template with one segment and a
    fixture with one segment have the same tree.
    """

    if isinstance(value, dict):
        return {key: _key_tree(child) for key, child in value.items()}
    if isinstance(value, list):
        return [_key_tree(child) for child in value]
    return None


class CertificateTests(unittest.TestCase):
    def test_template_matches_auditor_layout(self):
        # research/certificates/paper_vi_certificate_template.json is written
        # from the auditor, and this is what binds it there: its key tree must
        # equal the tree of the certificate fixture above, which the auditor
        # accepts on every check but the fail-closed array-replay gate.  The
        # template's own "_template" note is the one key not read by the
        # auditor.
        template_path = (
            Path(__file__).resolve().parents[1]
            / "research"
            / "certificates"
            / "paper_vi_certificate_template.json"
        )
        with template_path.open(encoding="utf-8") as handle:
            template = json.load(handle)
        note = template.pop("_template")
        self.assertIn("REQUIRED", note["how_to_read"])
        self.assertEqual(_key_tree(template), _key_tree(_valid_certificate()))

        # As it stands the template is not a certificate, and it must not
        # slip through: every pass flag is false and every number is a
        # placeholder, so the audit fails on each of them, not only on the
        # array-replay gate.
        audit = audit_certificate(template)
        self.assertFalse(audit.passed)
        self.assertEqual(template["status"], "template")
        self.assertFalse(any(audit.checks.values()), audit.checks)
        self.assertTrue(
            any("invalid decimal value 'REQUIRED" in issue for issue in audit.issues),
            audit.issues,
        )

    def test_scalar_precertificate_is_blocked_without_array_replay(self):
        audit = audit_certificate(_valid_certificate())
        self.assertFalse(audit.passed)
        self.assertEqual(audit.segment_count, 1)
        self.assertTrue(
            all(
                passed
                for name, passed in audit.checks.items()
                if name != "array_replay"
            ),
            audit.issues,
        )
        self.assertFalse(audit.checks["array_replay"])
        self.assertTrue(
            any("certified status cannot be issued" in issue for issue in audit.issues)
        )

    def test_nonnegative_radii_polynomial_fails(self):
        payload = _valid_certificate()
        payload["segments"][0]["q_upper"] = "0.75"
        payload["segments"][0]["radii_polynomial_upper"] = "5e-15"
        audit = audit_certificate(payload)
        self.assertFalse(audit.passed)
        self.assertTrue(any("strictly negative" in issue for issue in audit.issues))

    def test_weaker_uniqueness_hypothesis_fails(self):
        payload = _valid_certificate()
        del payload["static_uniqueness"]["conditions"][
            "lower_fixed_point_limit_for_every_x_gt_minus_2"
        ]
        audit = audit_certificate(payload)
        self.assertFalse(audit.passed)
        self.assertTrue(
            any("lower_fixed_point_limit" in issue for issue in audit.issues)
        )

    def test_wrong_uniqueness_domain_fails(self):
        payload = _valid_certificate()
        payload["static_uniqueness"]["comparison_domain"] = (
            "C minus (-infinity,-2]"
        )
        audit = audit_certificate(payload)
        self.assertFalse(audit.passed)
        self.assertTrue(
            any("comparison_domain" in issue for issue in audit.issues)
        )

    def test_declared_complete_replay_cannot_bypass_missing_engine(self):
        payload = _valid_certificate()
        payload["array_replay"] = {
            "profile": "paper-vi-array-replay-v1",
            "status": "complete",
            "source_file": "candidate",
        }
        audit = audit_certificate(payload)
        self.assertFalse(audit.passed)
        self.assertFalse(audit.checks["array_replay"])
        self.assertTrue(
            any("cannot be trusted" in issue for issue in audit.issues)
        )

    def test_endpoint_path_weight_is_recomputed(self):
        payload = _valid_certificate()
        payload["segments"][0]["path_endpoint_error_upper"] = "2e-14"
        audit = audit_certificate(payload)
        self.assertFalse(audit.passed)
        self.assertTrue(any("radius/sqrt" in issue for issue in audit.issues))

    def test_component_hash_must_be_bound_to_source_file(self):
        payload = _valid_certificate()
        payload["segments"][0]["interval_rh"]["data_sha256"] = "9" * 64
        audit = audit_certificate(payload)
        self.assertFalse(audit.passed)
        self.assertTrue(any("not bound" in issue for issue in audit.issues))

    def test_file_loader_rejects_hash_mismatch(self):
        payload = _valid_certificate()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            data = root / "data"
            data.mkdir()
            source = data / "candidate.bin"
            source.write_bytes(b"certified candidate")
            digest = hashlib.sha256(
                source.read_bytes()
            ).hexdigest()
            payload["source_files"]["candidate"]["sha256"] = digest

            def bind_data_hashes(value):
                if isinstance(value, dict):
                    for key, child in value.items():
                        if key == "data_sha256":
                            value[key] = digest
                        else:
                            bind_data_hashes(child)
                elif isinstance(value, list):
                    for child in value:
                        bind_data_hashes(child)

            bind_data_hashes(payload)
            certificate = root / "certificate.json"
            certificate.write_text(json.dumps(payload), encoding="utf-8")
            before = load_and_audit_certificate(certificate)
            self.assertFalse(before.passed)
            self.assertTrue(before.checks["source_files"])
            self.assertFalse(before.checks["array_replay"])
            source.write_bytes(b"tampered")
            audit = load_and_audit_certificate(certificate)
            self.assertFalse(audit.passed)
            self.assertTrue(any("hash mismatch" in issue for issue in audit.issues))


if __name__ == "__main__":
    unittest.main()
