# Paper VI certificate layout

`paper_vi_certificate_template.json` is the schema-v2 layout that
`fatou_backend.certification.audit_certificate` reads (console entry point
`tetration-cert-audit`). The implementation specification lists it among the
artefacts; the shipped packages did not contain it, so it is written here
**from the auditor**, not from any computed result.

It is a template, not a certificate:

- every numeric field, integer count and hash is a string beginning with
  `REQUIRED`, stating the inequality the auditor will check;
- every `pass` flag and every boolean claim is `false`. The auditor requires
  `true`, so a producer has to assert each claim itself rather than inherit
  it from the layout. The one document the auditor warns against is a
  collection of self-declared `pass=true` flags;
- the fields the auditor pins to a fixed value (`schema_version`,
  `theorem_profile`, the cut-plane domain string, the Paulsen–Cowgill DOI and
  proposition number, `array_replay.profile`) carry that value.

Running the auditor on the template therefore lists what is missing:

```bash
tetration-cert-audit research/certificates/paper_vi_certificate_template.json
```

The structure is bound to the auditor by
`tests/test_reconstruction_modes.py::CertificateTests::test_template_matches_auditor_layout`,
which compares the template's key tree with the certificate fixture the
auditor's own tests use.

## What can and cannot be filled today

The only segment with interval arrays and an independent replay is `8r1`
(`research/certification8r1/`, β ≈ −0.321 … −0.336). Its replay report carries
`Y`, `L`, `q`, `r`, `p(r)` and the `Z1…Z7` margins, which is the material for
one `segments[]` entry. It cannot produce a passing certificate on its own,
for two reasons the auditor enforces by design:

1. the segment chain must start at β = 0 (the anchor `e`) and be gapless up to
   `target_beta`; `8r1` is one interior segment of that chain;
2. `array_replay` only accepts the explicit `not-implemented` marker in this
   release and always adds an issue, so **no** certificate can pass the
   auditor until the independent array replay is wired to it for the whole
   chain. `replay_8r1.py` is that replay for one segment.

Both are stated in `research/certification8r1/README.md` under *What is not
closed*.
