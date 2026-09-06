"""exp-074: exp-048s Abschneide-Schwelle an den groessten BEITRAG haengen.

HISTORISCH -- LAEUFT NICHT MEHR
-------------------------------
Der Anker unten beschreibt den Fork-Stand VOR exp-074b. Seit die Schwelle als
Vektoroperation (`bv = vecmax(lg2 + lr2*icsq)`) im Fork steht, findet er nichts
mehr und das Skript bricht mit einer Ankermeldung ab. Das ist Absicht: es wird
als Beleg des Experiments aufbewahrt, nicht als Werkzeug. Wer den Vorher-Stand
wieder braucht, holt ihn aus der Versionsgeschichte statt hier einen Anker
nachzuziehen, der beim naechsten Umbau erneut bricht.

Befund (Journal 2026-07-30). exp-048 behaelt Term k, solange

    exponent(c_k) + k*log2(rho)  >=  mx - dig*log2(10) - 16,     mx = max_k exponent(c_k)

Die Schwelle haengt am groessten KOEFFIZIENTEN. Gemeint ist aber "dig Stellen
unter dem groessten BEITRAG |c_k|*rho^k". Beides ist dasselbe, solange die
Koeffizienten FALLEN -- dann liegt das Maximum von lg2 bei k=1, und dort ist
(k-1)*log2(rho) = 0. Der Kommentar der Funktion nennt die Annahme selbst:
"ct's coefficients decay like circr^-k", und das gilt genau fuer circr > 1.

Bei Basis e ist circr = 1.3372 > 1, die Annahme haelt.
Bei b = 1.4494 ist circr ~ 0.133 (gemessen) < 1, die Koeffizienten WACHSEN, mx sitzt am
hoechsten Index -- gemessen 160 von 161 -- und der groesste Beitrag bei
Index 2. Die Schwelle liegt damit 284 Bit (~85 Dezimalstellen) zu hoch, kk
faellt von 161 auf 45, sfunc liefert Muell, die Iteration stagniert sofort.
Ergebnis: 1.837 statt 140.3 Stellen bei dps 150, ohne Fehlermeldung.
A/B mit dem Schalter ic48on (exp-073) hat exp-048 als Ursache bestaetigt:
    ic48on=1 -> 77 Iter, 1.837 Stellen     ic48on=0 -> 240+ Iter, 140.3 Stellen

FIX: die Schwelle pro Stufe aus dem groessten Beitrag DIESER Stufe bilden.
Kosten: ein zusaetzlicher O(n)-Durchgang je Stufe, achtmal pro Aufbau.

WICHTIGE EIGENSCHAFT: wo die alte Annahme gilt, ist bv = mx und der Fix ist
ein No-op -- bei fallenden Koeffizienten liegt max(lg2) bei k=1 und
(k-1)*lr2 = 0. Wo sie kippt, korrigiert er in Richtung MEHR behaltener Terme,
also nie in Richtung falsch. Deshalb sollten die Gate-Werte bitgleich bleiben;
das ist zu PRUEFEN, nicht zu behaupten.
"""
from __future__ import annotations

from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SRC = REPO / "src" / "fatou_backend" / "vendor" / "fatou_fork.gp"
DST = REPO / "research" / "tools" / "_ic48fix.gp"

OLD = """  local(v, n, lg2, tol2, j, lr2, kk, ii, mx);
  icdlon = 0;
  if (type(pp) != "t_POL", return(0));
  v = Vecrev(pp);
  n = #v;
  if (n < 64, return(0));
  lg2 = vector(n, kk, if (v[kk]==0, -1000000, exponent(v[kk])));
  mx = vecmax(lg2);
  if (mx <= -1000000, return(0));
  tol2 = mx - dig*3.3219280948873623 - 16;
  icdl = vector(icnlev);
  for (j=1, icnlev,
    lr2 = log(icrmax*j/icnlev)/log(2);
    kk = n;
    while ((kk > 1) && (lg2[kk] + (kk-1)*lr2 < tol2), kk--);"""

NEW = """  local(v, n, lg2, tol2, j, lr2, kk, ii, mx, bv, cv, i3);
  icdlon = 0;
  if (type(pp) != "t_POL", return(0));
  v = Vecrev(pp);
  n = #v;
  if (n < 64, return(0));
  lg2 = vector(n, kk, if (v[kk]==0, -1000000, exponent(v[kk])));
  mx = vecmax(lg2);
  if (mx <= -1000000, return(0));
  icdl = vector(icnlev);
  for (j=1, icnlev,
    lr2 = log(icrmax*j/icnlev)/log(2);
    /* exp-074: Schwelle am groessten BEITRAG |c_k|*rho^k dieser Stufe, nicht
       am groessten Koeffizienten. Identisch, solange die Koeffizienten fallen
       (Maximum bei k=1, dort ist (k-1)*lr2 = 0); korrigiert den Fall
       wachsender Koeffizienten (circr < 1, Basen nahe eta), wo mx am hoechsten
       Index sitzt und die alte Schwelle ~85 Dezimalstellen zu hoch lag. */
    bv = -1000000;
    for (i3 = 1, n, cv = lg2[i3] + (i3-1)*lr2; if (cv > bv, bv = cv));
    tol2 = bv - dig*3.3219280948873623 - 16;
    kk = n;
    while ((kk > 1) && (lg2[kk] + (kk-1)*lr2 < tol2), kk--);"""


def main() -> None:
    text = SRC.read_text(encoding="utf-8", errors="surrogateescape")
    n = text.count(OLD)
    if n != 1:
        raise SystemExit(
            f"Anker {n}x gefunden, erwartet 1x -- dieses Skript ist HISTORISCH "
            f"und beschreibt den Fork-Stand vor exp-074b. Siehe Kopf der Datei.")
    DST.write_text(text.replace(OLD, NEW), encoding="utf-8", errors="surrogateescape")
    print(f"geschrieben: {DST}")
    print("Fork UNBERUEHRT. Zu pruefen: b=1.4494 erholt sich, Gate-Werte bitgleich.")


if __name__ == "__main__":
    main()
