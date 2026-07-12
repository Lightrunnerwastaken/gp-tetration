# M5 — Base-Change-Fastpath (Konsument des basechange-modes-Research)

Stand: 2026-07-12. Quelle: `Tetration/Research/basechange-modes/RESEARCH_LOG.md`
(§ "Empfehlung Datenformat (für den Rechner)"), dort AKTIVES Research —
wir lesen nur; alle Implementierung hier im Repo.

## Ziel

`sexp_b` / `slog_b` für beliebige Basen b über den **Anker e** plus
Base-Change-Funktion Φ_{b,e}, statt pro Basis einen tiefen Kontur-Init:

    slog_b(w) = (n+θ) + Φ_{b,e}(θ),  θ aus slog_e(w) (Turm-Peel)
    sexp_b     = Inversion (H̃' ∈ [0.997, 1.003] ⇒ Newton perfekt konditioniert)

Ökonomie (deren Messung): D(k)-Ausbeute 47@k=60 … ~148@k=300;
~500 digits ⇔ k ≈ 1300–1700 Moden ⇔ nlim ≈ 600–750, Init 1–2 h original
(≈ 20–40 min mit 21-Keep-Fork), EINMALIG + disk-cached pro Basis.

## Bausteine (nach deren Empfehlungs-Liste)

1. **µ-Hub-Tabelle** `research/reference/mu_hub.json`: m(b) = µ_{e,b}
   hochpräzise (N Basen statt N² Paare); Paar-µ via
   µ_{b,d} = m(d) − m(b) + ∆⁽²⁾ (Formel (42), validiert ~1e-13 absolut).
   Initialdaten: deren results/mb_raster.csv (24 Basen, 40 digits, READ-ONLY
   übernehmen) + m(2) auf 150+ digits.
2. **Φ-Fourier-Tabelle** für Arbeitsbereich ≤ ~150 digits (k ≤ 300):
   pro Paar low-k-Kopf (k ≤ 10–20) + universelle Kurve
   A_k ≈ C₁·|ln α|·G(k)·e^(δk) (3 Parameter/Paar). C₁(e) =
   0.00038859729244516473 (deren R2-2, ~17 digits).
3. **Punktweise σ-Leiter** (> ~150 digits): Turm explizit hochziehen bis
   T > 1e4, Zwei-Level-Einstieg W₂ = α̂ + c·T_d(m+θ), numerisches Peeling
   bis W ≤ 50, EIN slog-Call. Zertifizierbar (Trunkierung ≤ C₀·d^(−10⁴);
   Gate: Extra-Turmstufe ändert Φ um < 1e-120).
4. **nlim-Regel** (Moden-Läufe): Band ≈ 2.3·nlim, Floor ≈ e^(−2π·0.47·nlim).
   ACHTUNG: auf dem Fork neu kalibrieren (Theta-Quantisierung/Radius-Keeps
   verschieben das Trunkierungs-Artefakt).
5. **Gates als Beweis-Residuen**: µ-Antisymmetrie (≤1e-38), Groupoid-Residuum
   Φ_{d,b}(H̃(θ)) + Φ_{b,d}(θ) = 0, Trunkierungs-Gate.

## Implementierungs-Reihenfolge

- M5.1: `fatou_backend/basechange.py` — Peel-Leiter (Baustein 3) + Tests
  gegen deren publizierte Werte (µ_{e,2} = −1.1284037766289240…,
  µ_{3,5} = 0.5788350021355712…).
- M5.2: Φ-Moden-Lauf-Treiber (nlim-skaliert, Fork) + Floor-Rekalibrierung.
- M5.3: Fourier-Tabelle + universelle Kurve; slog_b/sexp_b-Fastpath ≤150 digits.
- M5.4: Deep-Pfad (Moden bis k~1500, einmalig/Basis) + Beweis-Gates.
- Abhängigkeit: deren σ∞/C_k(b,d)-Front ist AKTIV — Datenformate mit dem
  Nutzer abstimmen, bevor M5.3 festzurrt.

## Voraussetzungen aus M4 (erfüllt/laufend)

- Anker-Referenzen: e|500 ✓, e|700 ✓, e|1000 (Partner-Lauf läuft).
- Fork-Speed für Moden-Inits: 21 Keeps ✓.
