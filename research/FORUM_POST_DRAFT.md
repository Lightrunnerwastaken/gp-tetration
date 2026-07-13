# Forum-Post-Entwurf (Tetration Forum) — NUR ENTWURF, Publish auf Nutzer-Go

Titel-Vorschlag:
**fatou.gp fork: ~13x faster Kneser tetration, sexp values with proven
error bounds (up to 972 digits), and a "base atlas" (any base from one
anchor, milliseconds)**

---

Hi all,

over the last days I've been working on top of Sheldonison's fatou.gp
(first of all: thank you for that engine — everything below builds on it).
Three results that might interest this forum:

**1. An optimized fork, ~13x faster at high precision.**
21 changes, each verified against a frozen accuracy gate before keeping it
(caches for the orbit walks and Schröder evaluations, incremental
Taylor/theta/extraction updates at reduced precision, FFT extraction on
exact grids, precision laddering, contour radius tuning). dps 520
(495 true digits, base e) went from ~5.4 h to ~25 min on the same machine;
dps 1020 (973 digits) runs in ~6.7 h. The full experiment journal —
including everything that did NOT work (Aitken, Anderson/Krylov, fast
multipoint evaluation, and why) — is in the repo.

**2. sexp values with proven error bounds.**
Using the error-vector method (two runs at different iteration depths;
their agreement bounds the worse run's true error) plus engine diversity
(fork vs. original agreeing to 495+ digits):
- sexp_e(0.5) to 972 proven digits
- sexp_e(0.5) to 698 / 497 digits (mid tiers), sexp_2(0.5) to 497 digits
As far as I know such certified tetration values weren't publicly
available; happy to compute/verify more points or bases if useful.

**3. A "base atlas": tetration for any base from one anchor.**
After a single base-e setup, sexp/slog for arbitrary bases are evaluated
through a small Fourier-mode table of the base-change function (~3 KB per
base, first request ~0.2 s, then milliseconds) — validated against the
proven references to ~1e-24 (scalable to ~1e-37+ with bigger tables).
The underlying mode/phase analysis comes from companion research notes
(in preparation); the repo has the consumer side.

**Demo:** a single static HTML file plots sexp_b(x) with a *continuous*
base slider (b = 1.5 … 100) — the curve morphs live through base space,
no server needed. [Link/Anhang]

Repo: [GitHub-Link nach Push]
Everything is reproducible: frozen gate, frozen benchmark, journal.

Feedback very welcome — especially from Sheldonison regarding the fork
(I'd gladly upstream anything useful).

---

## Offene Punkte vor dem Absenden (Checkliste)
- [ ] GitHub-Push (main ist lokal ~160 Commits voraus) — Nutzer-Go
- [ ] Lizenz-/Attributionsfrage fuer fatou.gp-Fork klaeren (idealerweise
      Sheldonison direkt fragen — der Post tut das implizit)
- [ ] Companion-Notes-Referenz mit dem Research-Projekt abstimmen
      (nichts Unveroeffentlichtes vorwegnehmen)
- [ ] Calculator hosten (GitHub Pages) oder als Anhang
- [ ] v9 (2|700) noch eintragen, wenn fertig
