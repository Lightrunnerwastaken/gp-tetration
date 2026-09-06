"""Erzeugt _pre074.gp: der aktuelle Fork MIT exp-071, aber OHNE exp-074.

HISTORISCH -- LAEUFT NICHT MEHR
-------------------------------
Der Anker unten beschreibt den Fork-Stand VOR exp-074b. Seit die Schwelle als
Vektoroperation (`bv = vecmax(lg2 + lr2*icsq)`) im Fork steht, findet er nichts
mehr und das Skript bricht mit einer Ankermeldung ab. Das ist Absicht: es wird
als Beleg des Experiments aufbewahrt, nicht als Werkzeug. Wer den Vorher-Stand
wieder braucht, holt ihn aus der Versionsgeschichte statt hier einen Anker
nachzuziehen, der beim naechsten Umbau erneut bricht.

Wozu: die Zeitmessung fuer exp-074 braucht ein A/B, das sich in GENAU EINER
Sache unterscheidet. `git show HEAD:...fatou_fork.gp` waere dafuer untauglich,
weil HEAD auch die exp-071-Fixes noch nicht hat -- der Vergleich haette dann
zwei Unterschiede und der gemessene Faktor waere nicht zuordenbar.

Also wird exp-074 aus dem aktuellen Stand zurueckgedreht: die Schwelle haengt
wieder am groessten Koeffizienten mx statt am groessten Beitrag der Stufe.
Der Rest der Datei ist byteweise der Fork.

Erwartung fuer Basis e: kein messbarer Unterschied. icbuild laeuft EINMAL pro
Gitterabschnitt, und der Zusatzaufwand ist 8n Vergleiche gegen N^2
Polynomauswertungen pro Durchgang -- bei n ~ N ~ 7000 also ~0.1 %. Dass die
Gate-Werte bitgleich sind, zeigt zudem, dass dieselbe Arbeit geleistet wird.
Aber ~0.1 % ist eine Rechnung, keine Messung, und die Projektregel verlangt
einen Exponenten oder Faktor nur mit Drift-Klammer.
"""
from __future__ import annotations

from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SRC = REPO / "src" / "fatou_backend" / "vendor" / "fatou_fork.gp"
DST = REPO / "research" / "tools" / "_pre074.gp"

# exp-074 rueckwaerts: nur der Anker der Schwelle, nichts sonst.
NEW_BLOCK = """    if ((mxi == 1) && (lr2 < 0),
      bv = mx
    ,
      bv = -1000000;
      for (i3 = 1, n, cv = lg2[i3] + (i3-1)*lr2; if (cv > bv, bv = cv))
    );
    tol2 = bv - dig*3.3219280948873623 - 16;"""
OLD_BLOCK = """    tol2 = mx - dig*3.3219280948873623 - 16;  /* pre-exp-074: Anker am groessten KOEFFIZIENTEN */"""


def main() -> None:
    text = SRC.read_text(encoding="utf-8", errors="surrogateescape")
    n = text.count(NEW_BLOCK)
    if n != 1:
        raise SystemExit(
            f"exp-074-Block {n}x gefunden, erwartet 1x -- dieses Skript ist "
            f"HISTORISCH und beschreibt den Stand vor exp-074b. Siehe Dateikopf.")
    DST.write_text(text.replace(NEW_BLOCK, OLD_BLOCK),
                   encoding="utf-8", errors="surrogateescape")
    print(f"geschrieben: {DST}")
    print("Unterschied zum Fork: EIN Anker (mx statt bv). Sonst byteweise gleich.")


if __name__ == "__main__":
    main()
