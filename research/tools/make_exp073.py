"""exp-073: Schalter fuer die exp-048-Gradabschneidung — in einer TESTKOPIE.

Befund, der dazu fuehrt (Journal 2026-07-30, Versions-Sweep ueber alle 34
Fork-Versionen mit b = 1.4494 bei dps 150, einheitliches 240-s-Budget):

    Original                   70 Iter,  60.670 Stellen, Gitter  184, haelt
    exp-018 ... exp-032       ~195 Iter, ~122    Stellen, Gitter 3308, Timeout
    exp-033 ... exp-042b       233 Iter, 137.626 Stellen, Gitter 7120, Timeout
    af28712 (exp-048/051)       77 Iter,   1.837 Stellen, Gitter  292, haelt
    exp-056 ... exp-063         77 Iter,   1.837 Stellen, Gitter  292, haelt

Ein Commit, ein Einbruch, danach acht Keeps bitgleich. af28712 enthaelt ZWEI
Keeps, und dieser Schalter trennt sie:

  exp-048 -- Gradabschneidung pro Punkt (icbuild). Behaelt Term k, solange
     exponent(c_k) + k*log2(rho) >= exponent(max|c|) - dig*log2(10) - 16.
     Die Zuordnung eines Punktes zu seiner Radius-Stufe klemmt in sfunc:
     `lv = 1 + floor(icnlev*abs(zz-circc)/icrmax); if (lv > icnlev, lv = icnlev)`
     -- ein Punkt AUSSERHALB von icrmax landet damit still in der obersten
     Stufe, deren Radius kleiner ist als sein eigener, und wird zu stark
     abgeschnitten. Bei b nahe eta ist die Geometrie gestreckt, also ist genau
     das der Verdacht.
  exp-051 -- Quantisierung des Theta-Gitters unter 64. Fuer diese Basis ist
     thetamode = 0, das Theta-Gitter also gar nicht in Gebrauch -- deshalb der
     zweite, unwahrscheinlichere Kandidat.

ic48on = 1 ist bitgleich zur Vorlage. Mit ic48on = 0 kehrt icbuild sofort
zurueck, icdlon bleibt 0, und sfunc nimmt den vollen Grad wie vor exp-048.
Kommen damit die ~137 Stellen zurueck, ist exp-048 der Verursacher; bleiben
sie weg, ist es exp-051.
"""
from __future__ import annotations

from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SRC = REPO / "src" / "fatou_backend" / "vendor" / "fatou_fork.gp"
DST = REPO / "research" / "tools" / "_ic48.gp"

DECL_ANCHOR = "icnlev=8; swlev=0; icdl=0; icrmax=1; icdlon=0;"
DECL_NEW = (
    "icnlev=8; swlev=0; icdl=0; icrmax=1; icdlon=0;\n"
    "ic48on=1;  /* exp-073 scan knob: 0 disables the exp-048 per-point degree\n"
    "   truncation entirely (icbuild returns early, icdlon stays 0). */"
)

# Anker bewusst OHNE die local()-Zeile: die aendert sich bei jedem Umbau von
# icbuild -- exp-074 hat bv/icsq ergaenzt, exp-074b mxi wieder entfernt -- und
# hat diesen Patcher schon einmal totgelegt.
BODY_ANCHOR = (
    '  icdlon = 0;\n'
    '  if (type(pp) != "t_POL", return(0));'
)
BODY_NEW = (
    '  icdlon = 0;\n'
    '  if (ic48on == 0, return(0));  /* exp-073 */\n'
    '  if (type(pp) != "t_POL", return(0));'
)


def main() -> None:
    text = SRC.read_text(encoding="utf-8", errors="surrogateescape")
    for anchor in (DECL_ANCHOR, BODY_ANCHOR):
        n = text.count(anchor)
        if n != 1:
            raise SystemExit(f"Anker {n}x gefunden, erwartet 1x: {anchor[:50]!r}")
    text = text.replace(DECL_ANCHOR, DECL_NEW).replace(BODY_ANCHOR, BODY_NEW)
    DST.write_text(text, encoding="utf-8", errors="surrogateescape")
    print(f"geschrieben: {DST}  ({len(text.splitlines())} Zeilen)")
    print("Fork UNBERUEHRT. Test: ic48on=1 (Vorlage) gegen ic48on=0 (exp-048 aus).")


if __name__ == "__main__":
    main()
