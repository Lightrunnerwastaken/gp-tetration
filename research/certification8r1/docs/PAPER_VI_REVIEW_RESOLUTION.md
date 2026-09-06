# Paper VI – Review-Auflösung und verbleibende Beweisbarrieren

Stand: 29. Juli 2026

Bezug: Review „Paper VI – Was noch fehlt“

## Ergebnis in einem Satz

Die falsche Instanziierung von Paulsen–Cowgill Proposition 2 ist korrigiert,
und `8r1` liegt jetzt als erstes echtes outward-rounded Segment mit
unabhängigem PASS-Replay vor. Der vollständige Marsch \(T_e\to T_2\) ist
weiterhin **kein computerassistierter Beweis**, weil die übrigen Segmente,
ihre Endpunktinklusionen, das segmentuniforme End-Gluing und die
Schnittebenenzertifikate noch fehlen.

Diese Grenze ist nun sowohl im Paper als auch im Rechner technisch
fail-closed: Ein skalar konsistentes JSON mit `pass=true` kann den Modus
`certified` nicht freischalten.

## Punktgenaue Auflösung des Reviews

| Nr. | Review-Forderung | Bearbeitung | Belastbarer Status |
|---:|---|---|---|
| 1 | vollständige Intervall-RH-Daten | QH-, Koenigs-, Quotienten-, graded-Trace-, Fourier-/Chebyshev-, Hardy-, Tail-, Cut- und Lipschitzarrays für `8r1` erzeugt und hashgebunden | für `8r1` erfüllt; übrige Segmente offen |
| 2 | mindestens ein echtes Segmentzertifikat | `8r1` auf dem exakten \(\beta\)-Intervall neu gerechnet und vom getrennten Array-Replay akzeptiert | **erfüllt**: \(q<1\), \(p(r)<0\), Z1–Z7 positiv |
| 3 | alle Marschsegmente zertifizieren | Reihenfolge „erst `8r1`, dann identische Pipeline für alle Segmente“ festgeschrieben; volle Kugelinklusion und echte Grönwall-Summe sind Pflicht | offen |
| 4 | segmentuniformes End-Gluing | obere und untere Endmodelle, Basisboxuniformität, \(\sup|\vartheta_\beta'|<1\), Intervall-Newton und Re-Anker-Overlaps getrennt verlangt | analytischer Vertrag vorhanden, Zahlen offen |
| 5 | Cut-plane-Fortsetzung | endliche Übergangsboxen, Links-/Rechtsrekursion, Cut-Sektoren, Induktion, Overlaps und Monodromie als Replay-Pflichten festgeschrieben | Templates fehlen |
| 6 | Kneser-Eindeutigkeit exakt | Proposition 2 auf \(\Re z>-2\) mit **beiden** Grenzwerten \(L_1,L_2\) korrigiert; Identitätssatz zur Schnittebene als separater Schritt ergänzt | logische Form korrigiert, Intervallvoraussetzungen offen |
| 7 | maschinenlesbares Gesamtzertifikat | unabhängiger Arb-Array-Replay für das reale `8r1` implementiert; deterministische Wiederholung und Manipulationsablehnung geprüft | lokales Segmentbundle erfüllt; globales Bundle offen |
| 8 | Rechner fertigstellen | Garantiegrenzen vereinheitlicht; 8r1-Produzent und -Replay implementiert | Controller plus erster Intervallkern vorhanden; produktiver allgemeiner Float-/Intervallmarcher noch offen |

## Reales Ergebnis des ersten Segment-Replays

Die historische Rohdatei wurde in keiner Quelle gefunden. Deshalb wurde
nicht aus gerundeten Logwerten rekonstruiert, sondern das exakte Segment

\[
\beta_a=\frac78\log\log2,\qquad
\beta_b=\frac{11}{12}\log\log2
\]

mit sieben Basis-Lobattozuständen neu gerechnet. Die outward-rounded
Primitivdateien enthalten 576 uniforme Taylor-/Arrayzeilen und 5184
graded-Hilbert-Zeilen. Der unabhängige Replay liefert:

| Größe | Arb-Intervall |
|---|---|
| \(Y\) | \([2.89398623506807\cdot10^{-10}\pm8.15\cdot10^{-25}]\) |
| \(L\) | \([13613.30588130968\pm8.70\cdot10^{-12}]\) |
| \(q\) | \([0.1495393411962912\pm9.19\cdot10^{-17}]\) |
| \(r\) | \([4.2535568887144\cdot10^{-10}\pm3.01\cdot10^{-24}]\) |
| \(p(r)\) | \([-7.234965587670\cdot10^{-11}\pm2.51\cdot10^{-24}]\) |

Alle Z1–Z7-Margen sind strikt positiv. Der kleinste Rand ist die
superexponentielle Tail-Marge aus Z5,
\(9.64335466598625\cdot10^{-51}\). Der konservative Lipschitzmajorant
erzwingt eine uniforme mathematische Unterteilung in \(262144\)
Mikrosegmente. Ein zweiter Replay war byte-identisch; eine absichtlich
veränderte Arraydatei wurde wegen des kanonischen Payload-Hashs abgelehnt.

Dieser Abschluss ist lokal und setzt eine zertifizierte Eingangskugel
voraus. Er zertifiziert noch nicht den vorangehenden Marschabschnitt oder
die globale Identifikation bei Basis \(2\). Der Zielbasis-Engine wurde zur
Erzeugung der beliebigen Kandidatenzentren und nur für den blinden
Stencilvergleich verwendet; Intervallproduzenten und Replay rufen ihn nicht
auf. Damit ist der lokale a-posteriori Segmenttest echt, seine Eingangslinie
aber noch nicht als anchor-pure vom \(e\)-Anker her bewiesen.

## Wesentliche mathematische Korrektur

Die frühere Revision bezeichnete folgende Kombination als „exakt
Paulsen–Cowgill Proposition 2“:

- Holomorphie auf
  \(\mathbb C\setminus(-\infty,-2]\);
- Funktionalgleichung;
- Normalisierung und Konjugationssymmetrie;
- nur der obere Fixpunktgrenzwert.

Das ist nicht die exakte Aussage der Proposition. Der korrekte Schluss hat
zwei Stufen.

### Stufe 1 – Vergleich auf der Halbebene

Für zwei Kandidaten \(F_1,F_2\), analytisch auf

\[
H_{-2}=\{z\in\mathbb C:\Re z>-2\},
\]

sind nach Paulsen–Cowgill Proposition 2 zu prüfen:

1. \(F_j(z+1)=b^{F_j(z)}\) auf dem zulässigen Teil von \(H_{-2}\);
2. \(F_j(0)=1\);
3. für jedes \(x>-2\)
   \[
   \lim_{y\to+\infty}F_j(x+iy)=L_1;
   \]
4. für jedes \(x>-2\)
   \[
   \lim_{y\to-\infty}F_j(x+iy)=L_2.
   \]

Dann gilt \(F_1=F_2\) auf \(H_{-2}\).

### Stufe 2 – Fortsetzung auf die Schnittebene

Sind der rekonstruierte Zweig und der Kneser-Zweig zusätzlich holomorph auf
der zusammenhängenden Menge

\[
\mathbb C\setminus(-\infty,-2],
\]

so erweitert der Identitätssatz die bereits auf \(H_{-2}\) bewiesene
Gleichheit auf die ganze Schnittebene. Die Schnittebenenholomorphie ist also
nicht Ersatz für einen der beiden Fixpunktgrenzwerte.

Referenz: William Paulsen und Samuel Cowgill, *Solving
\(F(z+1)=b^{F(z)}\) in the complex plane*, Advances in Computational
Mathematics 43 (2017), Proposition 2,
DOI `10.1007/s10444-017-9524-1`.

## Erreichter erster Meilenstein: `8r1`

Die verlangte Abnahme durch einen vom Produzenten unabhängigen Checker ist
für das im Forschungslog genannte Kantensegment

\[
b: 2.0661\ldots\longrightarrow 2.0435\ldots .
\]

erfolgt. Der folgende Katalog ist damit keine Vorlage mehr, sondern die
Inhaltsbeschreibung des ausgelieferten und hashgebundenen Datensatzes.

### Exportierte Eingaben

- vollständige Lobatto-Knoten und Gewichte;
- alle sieben Kandidatenzustände;
- HI-\(\bar G\)-Werte und RH-Fouriermoden je Knoten;
- Base-\(\beta\)-Intervall und Branch-Konventionen;
- Arbeitspräzision und Intervallbackend;
- validierte Transformationsmatrizen oder ein reproduzierbarer
  Intervall-FFT-/Chebyshev-Pfad;
- endliche RH-Köpfe und mathematisch begründete analytische Tails;
- Hashes der Rohdaten, des Produzentencodes und der Konventionen.

### Vom Checker neu berechnete Größen

\[
\begin{aligned}
Y &=Y_{\mathrm{head}}+Y_{\mathrm{tail}},\\
L&=\sup_{\mathfrak B}\|D\mathscr V\|,\\
q&\ge \frac{4\sqrt2\,L}{\kappa},\\
p(r)&=Y+(q-1)r,\\
\varepsilon_{\mathrm{path}}
&\ge \frac{r}{
  \sqrt{\rho_0-\kappa a-\rho_{\mathrm{end}}}}.
\end{aligned}
\]

Die lokale Segmentakzeptanz erfordert gleichzeitig:

- \(0\le q<1\);
- \(p(r)<0\);
- positive Strip-, Ableitungs-, Cut-, Koenigs- und Quotientenränder auf der
  **ganzen** Radii-Kugel;
- eine outward-rounded Endpunktkugel.

Diese Bedingungen werden vom 8r1-Replay erfüllt. Für die Aufnahme des
Segments in die globale Marschkette sind zusätzlich die zertifizierte
Eingangskugel und bei einem Nicht-Endsegment die strikte Inklusion
  \[
  E_{\mathrm{out}}+D_{\mathrm{centers}}+m_{\mathrm{strict}}
  \le R_{\mathrm{next}}.
  \]

zu prüfen. Diese globale Handoff-Inklusion ist für 8r1 noch offen. Der
veröffentlichte Report nennt \(Y,L,q,r,p(r)\) als Intervalle und führt die
Neuauswertung auf konkrete Replay-Knoten und Quelldateihashes zurück.

## Suchergebnis und Neuberechnung von `8r1`

Der übergebene Forschungslog nennt:

- `save_traj_8r1.py`;
- `trajectory_8r1.json`;
- `harvest_r16.json`;
- weitere R16-Marschdaten.

Diese Dateien waren weder im aktuellen Git-Repository, in erreichbaren
Git-Objekten und Archiven noch in den übergebenen Anlagen vorhanden.
Gerundete Resultate wurden nicht zu angeblichen Intervallarrays
hochgerechnet. Stattdessen wurde der zweite saubere Weg vollständig
ausgeführt: Das exakte Segment wurde mit dokumentierten Konventionen neu
gerechnet; sämtliche Zwischenarrays, Tails, Transformationen und Hashes
wurden exportiert; danach wurde das lokale Segment von einem getrennten
Arb-Array-Replay akzeptiert.

Die neu erzeugten Dateien tragen aus Kompatibilitätsgründen wieder die Namen
`save_traj_8r1.py` und `trajectory_8r1.json`. Ihre Herkunft ist im
Zertifikatsbericht ausdrücklich als Neuberechnung und nicht als
Wiederauffinden des historischen Rechenstands gekennzeichnet.

## Härtegrad des Rechners nach dieser Bearbeitung

### Fast

Vorhanden:

- Modussteuerung;
- read-only Atlas und exakter \(e\)-Anker;
- Phasenabbildung und Komposition;
- Fehlerprognose und fehlertolerante Ausgabe für vorhandene Atlaswerte.

Fehlt:

- produktiver Float-RH-Marcher für eine nicht gespeicherte Zielbasis.

### Validated

Vorhanden:

- Verträge und Gates für Pfadsplitting, unabhängige Auflösung,
  Normierung, Funktionalgleichung, Roundtrip und Residuum-zu-Fehler.

Fehlt:

- der echte lokale Marcher, der diese voneinander unabhängigen Daten
  erzeugt.

### Certified

Vorhanden:

- Schema v2;
- Hash-, Parameter-, Segmentketten-, Radii-, Endpunkt-, Grönwall-,
  Auswertungs- und statische Eindeutigkeitsprüfungen;
- vollständige outward-rounded RH-/Tail-/Transformationsarrays für `8r1`;
- producer-unabhängiger Arb-Replay für `8r1`, einschließlich
  \(Y,L,q,r,p(r)\), Z1–Z7, Normierung, Funktionalgleichung und
  Manipulationsablehnung;
- expliziter Fail-closed-Blocker für eine unvollständige globale
  Segmentkette.

Fehlt:

- die Verallgemeinerung des Intervall-RH-Produzenten und Replays auf alle
  übrigen Marschsegmente;
- die zertifizierte Eingangskugel und globale Handoff-Inklusion für `8r1`;
- die verbleibenden segmentweisen Intervallzertifikate;
- Intervall-End-Gluing auf der gesamten Kette;
- Cut-plane-Templates;
- unabhängiger Replay des vollständigen globalen Bundles.

Damit ist das lokale Profil `8r1` ein nutzbares Beweisartefakt. Der
öffentliche globale Modus `certified` bleibt bis zur vollständigen Kette
eine reservierte, fail-closed Garantieklasse.

## Freigaberegel

Die Formulierung „Paper VI ist computerassistiert zertifiziert“ darf erst
verwendet werden, wenn alle folgenden Dateien aus einem sauberen Replay
hervorgegangen sind:

1. ein vollständiger, hashgebundener Base-\(e\)-Anker;
2. ein akzeptiertes reales `8r1`-Zertifikat;
3. entsprechende Zertifikate für sämtliche Segmente von \(e\) bis \(2\);
4. alle Endpunktkugel-Inklusionen;
5. beidseitige, segmentuniforme End-Gluing-Daten;
6. Cut-plane-Templates und Monodromie-/Overlap-Checks;
7. ein aus echten Segmentradien aufgebautes globales Grönwall-Budget;
8. ein unabhängiger Replay-Trace mit erfolgreichem Endstatus.

Vor diesem Zeitpunkt lauten die korrekten Bezeichnungen:

- **analytisch bewiesener Zertifizierungsvertrag**;
- **blind numerisch validierte Rekonstruktion**;
- **noch nicht computerassistiert zertifizierter Basis-2-Marsch**.
