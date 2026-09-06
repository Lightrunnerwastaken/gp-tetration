# Paper VI – Outward-rounded Segmentzertifikat `8r1`

Stand: 29. Juli 2026

## 1. Ergebnis und Reichweite

Das Kantensegment

\[
\beta_a=\frac78\log\log2,\qquad
\beta_b=\frac{11}{12}\log\log2
\]

beziehungsweise

\[
b:\
2.0660566229920658700359362681437378\ldots
\longrightarrow
2.0434599358867413651231620349778924\ldots
\]

wurde vollständig neu gerechnet. Der vom Produzenten getrennte
Arb-Array-Replay akzeptiert das Segment:

| Größe | outward-rounded Arb-Intervall |
|---|---|
| \(Y_{\rm head}\) | \([1.034581190661640\cdot10^{-14}\pm5.10\cdot10^{-30}]\) |
| \(Y_{\rm tail}\) | \([2.89388277694900\cdot10^{-10}\pm4.31\cdot10^{-25}]\) |
| \(Y\) | \([2.89398623506807\cdot10^{-10}\pm8.15\cdot10^{-25}]\) |
| \(L\) | \([13613.30588130968\pm8.70\cdot10^{-12}]\) |
| \(q\) | \([0.1495393411962912\pm9.19\cdot10^{-17}]\) |
| \(r\) | \([4.2535568887144\cdot10^{-10}\pm3.01\cdot10^{-24}]\) |
| \(p(r)\) | \([-7.234965587670\cdot10^{-11}\pm2.51\cdot10^{-24}]\) |

Damit sind \(0\le q<1\) und \(p(r)<0\) als strikte
Intervallaussagen bewiesen. Z1–Z7, Strip-Slack, Normalisierung,
Funktionalgleichung und Rohtransformation bestehen ebenfalls.

Dies ist ein lokales a-posteriori Segmentzertifikat. Es ist noch kein
globaler Beweis des Marsches \(T_e\to T_2\): Die zertifizierte
Eingangskugel aus dem vorangehenden Segment, die übrigen Segmente,
Endpunktinklusionen, End-Gluing und Cut-Plane-Templates fehlen weiterhin.

## 2. Provenienz

Die im alten Forschungslog genannten Dateien `trajectory_8r1.json` und
`save_traj_8r1.py` wurden im Repository, in Git-Objekten, Archiven und im
Dateispeicher nicht gefunden. Es wurden keine Zahlen aus gerundeten
Berichtstabellen zurückgerechnet.

Der neue Kandidatenproduzent verwendet die Zielbasis-Engine ausschließlich
zur Erzeugung eines beliebigen guten Zentrums und für einen blinden
\(\beta\)-Stencilvergleich. Das steht im Datensatz als
`candidate_target_engine_used=true`. Weder die beiden Intervallproduzenten
noch der unabhängige Replay starten PARI/GP oder importieren
Zielbasiswerte. Der lokale a-posteriori Beweis ist von der Herkunft des
Zentrums unabhängig; die globale anchor-pure Elternkette ist damit jedoch
noch nicht hergestellt.

## 3. Dateien und Aufgaben

| Datei | Aufgabe |
|---|---|
| `produce_8r1.py` | sieben Basis-Lobattoknoten, 48 Zustandsknoten, Rohtransformation, Koenigs-RH, graded Hilbert, Tails und Floatdiagnostik |
| `save_traj_8r1.py` | kompatibler Einstieg unter dem im Forschungslog genannten Namen |
| `trajectory_8r1.json` | vollständige Kandidatenbahn, Transformationen, Quadratur, Zustände, RH-Werte, Moden, Tails und Quellhashes |
| `certify_8r1.py` | 448-Bit-Arb-Taylorboxen für \(T,T',T''\), \(h,h',h''\), \(\Gamma,\Gamma'\), \(S,\partial_T S,\partial_{T'}S\) |
| `segment_8r1_interval_primitives.json` | 576 uniforme outward-rounded Primitive Rows |
| `certify_hilbert_8r1.py` | outward-rounded Koenigs-Spur an allen graded Gauss-Legendre-Knoten |
| `segment_8r1_hilbert_primitives.json` | 5184 graded-Hilbert-Taylorzeilen |
| `replay_8r1.py` | unabhängiger RH-, Hardy-, Volterra-, Lipschitz-, Transfer- und Radii-Replay |
| `segment_8r1_replay_report.json` | akzeptierter Report mit allen 95 aufgelösten Moden je Basispanel, Intervallen und Hashbindungen |

## 4. Kandidatenrechnung

Die Basisvariable ist \(\beta=\log\log b\). Auf dem Segment werden sieben
Chebyshev-Lobattoknoten verwendet. Jeder Zustand besitzt 48
Chebyshev-Lobattowerte auf \(x\in[0,1]\); die RH-Spur wird an 96 uniformen
Punkten und auf 36 symmetrisch graded Panels mit jeweils 24
Gauss-Legendre-Punkten berechnet.

Der RH-Wert wird nicht aus der Zielbasis-Bahn differenziert. Für
\(c=e^\beta\) werden

\[
L=-\frac{W_{-1}(-c)}c,\qquad \lambda=cL
\]

und der prinzipielle Rückwärtsorbit

\[
w_{n+1}=\frac{\log w_n}{c}
\]

verwendet. Gleichzeitig laufen die exakten \(\beta\)- und
\(y\)-Ableitungsrekursionen. Daraus entstehen
\(\Gamma=\partial_\beta\chi/\chi'\), \(\Gamma'\) und
\(v_{\rm geo}=\Gamma(T)/T'\). Der imaginäre Randwert wird mit dem
subtrahierten periodischen Cotangenskern transformiert; die Normierung bei
\(x=0\) wird anschließend exakt abgezogen.

Der unabhängige Zielbasis-Stencil dient nur als Diagnose. Über alle sieben
Basisknoten beträgt seine maximale Abweichung vom engine-freien RH-Feld

\[
1.253930445303475\cdot10^{-15}.
\]

## 5. Outward-rounded Primitive

Die beiden Intervallproduzenten verwenden `python-flint`/Arb mit 448 Bit.
Jedes der sechs aufeinanderfolgenden Basis-Lobattointervalle wird als
Taylorbox der Ordnung 22 behandelt. Dadurch bleiben die starken
Korrelationen in \(\beta\) bis nach DFT, Hardy-Projektion und
Normalisierung erhalten.

Der gestoppte Koenigs-Lauf weist folgende globale Gates auf:

| Gate | Intervall/Schranke |
|---|---|
| kleinste Phasenmarge, uniforme Rows | \(>1.9057\) rad |
| kleinste Phasenmarge, graded Rows | \(>1.6071\) rad |
| \(1/|\lambda|\) | \(<0.806839\) |
| terminale Fixpunktdistanz | \(<9.99978\cdot10^{-17}\) |
| \(\Gamma\)-Rest | \(<1.77496\cdot10^{-12}\) |
| \(\Gamma'\)-Rest | \(<3.80530\cdot10^{-10}\) |
| Tail, letzter Term | \(<3.56646\cdot10^{-52}\) |

Die seriellen Intervalle werden beim Export zusätzlich auf nach außen
erweiterte IEEE-754-Endpunkte gerundet. Der Replay liest diese Endpunkte
wieder als Arb-Kugeln ein.

## 6. Unabhängiger Replay

`replay_8r1.py` importiert keinen Produzentencode. Er führt in dieser
Reihenfolge neu aus:

1. kanonische Payload- und Dateihashprüfung;
2. Rekonstruktion aller \(\beta\)-Taylorboxen;
3. Quotienten \(v_{\rm geo}=\Gamma/T'\);
4. graded Cotangentintegral und Hardy-Normalisierung;
5. komplexen Tail \(S\) und \(G=\mathcal R-S\);
6. normalisiertes Feld
   \[
   \Psi=-(1+h')G-\bigl(-(1+h')G\bigr)(0);
   \]
7. Volterraresiduum und alle 95 aufgelösten Fouriermoden;
8. QH-, Trace-, RH- und Tail-Lipschitzmajoranten;
9. uniforme mathematische Segmentunterteilung;
10. \(Y,L,q,r,p(r)\), Endpunktgewicht und Z1–Z7;
11. Normalisierung, Funktionalgleichung samt Ableitungsrandgleichung und
    Rohtransformationsidentitäten.

Der konservative Full-Ball-Wert \(L\) führt zu \(262144\)
Mikrosegmenten. Pro Zelle gelten

\[
\rho_0=0.05,\qquad \rho_{\rm end}=0.01,\qquad
\rho_0-\kappa a-\rho_{\rm end}=0.01.
\]

## 7. Z1–Z7

| Margin | akzeptiertes Intervall |
|---|---|
| Z1 Hub-Kugel | \([0.09886741207065093\pm4.39\cdot10^{-18}]\) |
| Z2 \(U'\) | \([0.9776218658730432\pm5.55\cdot10^{-17}]\) |
| Z2 \(Q'\) | \([0.8826871485040080\pm4.20\cdot10^{-17}]\) |
| Z3 \(T'\) | \([0.8329225120637840\pm5.58\cdot10^{-17}]\) |
| Z3 Reziprokal | \([0.7994081255863313\pm7.15\cdot10^{-17}]\) |
| Z4 Cut/Fixpunkt | \([0.99888478371099921\pm2.84\cdot10^{-18}]\) |
| Z5 Transfer/Tail | \([9.64335466598625\cdot10^{-51}\pm4.01\cdot10^{-69}]\) |
| Z6 \(|1-\lambda|\) | \([1.16784210316836834\pm4.62\cdot10^{-19}]\) |
| Z6 \(|\lambda|-1\) | \([0.206316559281329926\pm4.57\cdot10^{-19}]\) |
| Z7 Branch-Overlap | \([0.193161183046091034\pm1.25\cdot10^{-19}]\) |

Der endliche Transferteil verwendet \(J_0=2\). Danach wird nicht mit einer
explodierenden globalen Turmbox weitergerechnet; der superexponentielle
Reziprokalprodukttail wird separat replayt. Sein letzter Term liegt bei
\(3.57\cdot10^{-52}\), unter dem Stoppwert \(10^{-50}\).

## 8. Reproduktionsbefehle

Voraussetzungen:

- PARI/GP 2.17.4;
- Python-Umgebung des Repositories;
- `python-flint==0.9.0`;
- `mpmath` und `numpy`.

```bash
export FATOU_GP_EXE=/absolute/path/to/gp

.venv/bin/python research/certification8r1/save_traj_8r1.py \
  --gp-exe "$FATOU_GP_EXE" \
  --output research/certification8r1/trajectory_8r1.json

.venv/bin/python research/certification8r1/certify_8r1.py \
  --trajectory research/certification8r1/trajectory_8r1.json \
  --output research/certification8r1/segment_8r1_interval_primitives.json

.venv/bin/python research/certification8r1/certify_hilbert_8r1.py \
  --trajectory research/certification8r1/trajectory_8r1.json \
  --output research/certification8r1/segment_8r1_hilbert_primitives.json

.venv/bin/python research/certification8r1/replay_8r1.py \
  --trajectory research/certification8r1/trajectory_8r1.json \
  --primitives research/certification8r1/segment_8r1_interval_primitives.json \
  --hilbert research/certification8r1/segment_8r1_hilbert_primitives.json \
  --report research/certification8r1/segment_8r1_replay_report.json
```

## 9. Hashes

| Artefakt | SHA-256 |
|---|---|
| `trajectory_8r1.json` | `bb6dfe2fdc17946eb485abe3e89def6c89d18982d853b650910d0fcb6782fd67` |
| `segment_8r1_interval_primitives.json` | `e8a0493baad0c742a6220d1a76385b8b24f08e0c6d0e86418cc732a97c669765` |
| `segment_8r1_hilbert_primitives.json` | `321efb1a21f92049c69ea5bdbba718a509704b339ba071e83f9f109e199137d3` |
| `segment_8r1_replay_report.json` | `d5ee3f697bb27fc2459ffa17fb383550b469bd6f547115656b0a896136a6e905` |
| `produce_8r1.py` | `3acc4e9c9cf96bcc4f19731a1b2c224826627c41636a2897a7308c27007862f3` |
| `certify_8r1.py` | `0e7e608222ad5b2ec8cbd354837f700e67fc23f9b33fcf5a7366408265c2f0e2` |
| `certify_hilbert_8r1.py` | `74997b67838776ed9080d7631db5715f4567a59d5a23acc4c0778d1c117aca1d` |
| `replay_8r1.py` | `965bd832130645d32ffd540e6e041489bfbd37d88134b88fde8a3d3f7530ccb8` |

Die JSON-Dateien enthalten zusätzlich kanonische Payload-Hashes. Der
akzeptierte Replayreport trägt den Payload-Hash
`b92ccde12f5891058706bb50165afb96bccf9162c0e8732654f1b373e65b1be1`.

## 10. Negativ- und Determinismustests

- Zwei unveränderte Replayläufe erzeugten byte-identische Reports.
- Nach Änderung eines einzelnen `x_index`-Feldes im Primitive-JSON endete
  der Checker mit Exitcode 1 und
  `primitive canonical payload hash mismatch`.
- Ein enger vorläufiger Z3-Reziprokalwert wurde nicht stillschweigend
  toleriert: Der erste Replay lehnte ab. Der endgültige Checker konstruiert
  \(\mu_T\) aus dem replayten \(m_{T'}\) und prüft den positiven Rest.
- Eine ungetrennte vierte Transferbox wickelte und wurde abgelehnt. Der
  endgültige Beweis verwendet den mathematisch getrennten endlichen Teil
  \(J_0=2\) und den independently replayten Reziprokalprodukttail.

## 11. Nächste verbindliche Arbeiten

1. Zertifizierte Eingangskugel aus dem vorhergehenden Segment einschließen.
2. Dasselbe Produzenten-/Replayprofil auf alle übrigen Marschsegmente
   anwenden.
3. Strikte Endpunktkugel-Inklusionen und die reale Grönwall-Summe bilden.
4. Segmentuniforme obere und untere Endmodelle prüfen.
5. Cut-Plane-Templates und anschließend den
   Paulsen–Cowgill-/Identitätssatz-Schluss replayen.
6. Erst dann den globalen `Certified`-Modus freischalten.
