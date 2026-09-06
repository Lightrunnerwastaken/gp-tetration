"""exp-072: macht die Theta-Abschaltschwelle zu einem Knopf — in einer TESTKOPIE.

Hintergrund
-----------
In fatou.gp (Z.22) und im Fork steht wörtlich:

    /* I added || (real(Period)>47) to handle speed for sexpinit(1.4494);
       takes the place of theta0lim=0.224 */

und die Bedingung lautet

    if ((ctr==1) || (arp<theta0lim) || (real(Period)>47), thetamode=0, thetamode=1);

Gemessen 2026-07-29 bei b = 1.4494 (Re(Period) = 47.227, also knapp über der
Schwelle, thetamode wird 0):

    dps 150 -> 1.8369421892 Konturstellen
    dps 300 -> 1.8369000816 Konturstellen

Eine harte Mauer: mehr Präzision anfordern ändert nichts. Die Engine gibt
trotzdem 38 Stellen aus, ohne Fehler. Die Schwelle wurde ausgerechnet für
diese Basis eingebaut und nimmt ihr die Genauigkeit vollständig.

Hypothese: mit erzwungenem Theta kommen die Stellen zurück — langsam statt
falsch. Dieses Skript erzeugt die Testkopie, um das zu prüfen, OHNE den
ausgelieferten Fork anzufassen (gleiche Vorgehensweise wie make_exp070.py).

Knopf: perlim, Vorgabe 47 -> bitgleiches Verhalten zur Vorlage.

VORBEHALT (2026-07-30): die hier erzeugte Testkopie enthaelt WEDER exp-071
(safefs/sub-eta) NOCH exp-074 (Abschneide-Schwelle). Sie ist als A/B in sich
geschlossen und die damit gemessenen Zahlen bleiben gueltig, aber sie ist
nicht der ausgelieferte Fork. Wer sie fuer neue Messungen benutzt, misst einen
Stand von vor 0.1.1.
"""
from __future__ import annotations

from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SRC = REPO / "src" / "fatou_backend" / "vendor" / "fatou_fork.gp"
DST = REPO / "research" / "tools" / "_perlim.gp"

# einzeiliger, eindeutiger Anker (Patcher-Regel aus dem Journal)
ANCHOR = "  if ((ctr==1) || (arp<theta0lim) || (real(Period)>47), thetamode=0,thetamode=1);"
REPLACEMENT = "  if ((ctr==1) || (arp<theta0lim) || (real(Period)>perlim), thetamode=0,thetamode=1);"

# perlim global einführen, direkt neben den anderen Knöpfen
DECL_ANCHOR = "ctrmul = 1;  /* exp-023 scan knob: scales the theta sampling radius */"
DECL_NEW = (
    "ctrmul = 1;  /* exp-023 scan knob: scales the theta sampling radius */\n"
    "perlim = 47; /* exp-072 scan knob: real(Period) above which theta is\n"
    "   switched OFF. Upstream hardcodes 47 for speed at b~1.4494 and thereby\n"
    "   caps that base at 1.84 contour digits regardless of the request. */"
)


def main() -> None:
    text = SRC.read_text(encoding="utf-8", errors="surrogateescape")
    for anchor in (ANCHOR, DECL_ANCHOR):
        n = text.count(anchor)
        if n != 1:
            raise SystemExit(f"Anker {n}x gefunden, erwartet 1x: {anchor[:60]!r}")
    text = text.replace(DECL_ANCHOR, DECL_NEW)
    text = text.replace(ANCHOR, REPLACEMENT)
    DST.write_text(text, encoding="utf-8", errors="surrogateescape")
    print(f"geschrieben: {DST}  ({len(text.splitlines())} Zeilen)")
    print("Fork UNBERUEHRT. Test: perlim=47 (Vorlage) gegen perlim=1e9 (Theta erzwungen).")


if __name__ == "__main__":
    main()
