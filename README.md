# Fatou Backend

Dieses Projekt kapselt `fatou.gp` als duennes Python-Backend um PARI/GP. Ziel ist
nicht ein neuer Tetrationsport, sondern eine reproduzierbare und kleine API auf
dem bestehenden GP-Rechenkern.

## Ziel

- `gp.exe` + `fatou.gp` als primaeren Rechenkern benutzen
- eine kleine Python-API fuer `sexp`, `slog`, Batch-Auswertung und Roundtrips anbieten
- Paper- und Analysewerkzeuge wie `mixed_phase` auf denselben Wrapper setzen

## Struktur

- `fatou_backend.gp_backend.FatouGP`: Haupt-Wrapper fuer PARI/GP
- `fatou_backend.gp_backend.FatouGPSession`: basisgebundene Session-Hilfe
- `fatou_backend.cli`: allgemeine CLI fuer `sexp`, `slog`, `roundtrip`, `eval`
- `fatou_backend.mixed_phase`: Mixed-Base-Phase-Profile auf demselben Backend
- `src/fatou_backend/vendor/fatou.gp`: vendorte GP-Datei

## Pfadauflosung

`fatou.gp` wird in dieser Reihenfolge gesucht:

1. expliziter Konstruktorparameter `fatou_gp=...`
2. Umgebungsvariable `FATOU_GP_FILE`
3. vendorte Repo-Datei `src/fatou_backend/vendor/fatou.gp`

Es gibt keinen stillen Fallback mehr auf `Downloads`.

`gp.exe` wird in dieser Reihenfolge gesucht:

1. expliziter Konstruktorparameter `gp_exe=...`
2. Umgebungsvariable `FATOU_GP_EXE`
3. Windows-Standardpfade fuer PARI/GP

## Python-API

```python
from fatou_backend import FatouGP

gp = FatouGP(dps=80)
value = gp.sexp("e", 0.5)

session = gp.session("1+I")
vals = session.sexp_batch([0.25, 0.4 + 0.2j])
logs = session.slog_batch(vals)
```

## CLI

```powershell
python -m pip install -e .
python -m fatou_backend.cli sexp --base e --values 0.5
python -m fatou_backend.cli slog --base 2 --values 1.5
python -m fatou_backend.cli eval --base 1+I --expressions "sexp(0.5)"
python -m fatou_backend.mixed_phase --outer e --inner 2
```

## Tests

Schneller Standardlauf:

```powershell
python -m unittest discover -s tests -p "test_gp_backend.py" -v
```

Mit Slow-Block:

```powershell
$env:FATOU_BACKEND_RUN_SLOW = "1"
python -m unittest discover -s tests -v
```
