# Failed approaches — gp-tetration

_Append-only. **Vor jedem neuen Ansatz prüfen** (Verfassungsregel).
Neueste Einträge zuerst. Wurzelursache, nicht Symptom._

## „Alle Basen"-Keeps ohne Shell-Thron-Probe (2026-07-15, Commits exp-032…034)

- Was versucht wurde: Fork-Optimierungen (exp-032 ff.) wurden als für
  „alle Basen" gültig behalten, verifiziert aber nur auf Gate-Basen
  (reell > η, komplex).
- Warum es scheiterte (Wurzelursache): Shell-Thron-Interior (1 < b < η,
  kc ≤ 0) war blinder Fleck — exp-032 (6c63a09) brach dort die reelle
  reguläre Tetration (b=1.2 lieferte 7.49−2.79i statt 1.13626…);
  gefunden per Commit-Bisektion, gefixt in 10e2350 (subeta-Flag + vier
  Guards).
- Nicht wiederholen, außer: b=1.2 läuft als Diversitäts-Probe mit
  (jetzt Slow-Regressionstest „Engine-Diversität b=1.2").
- Quelle: journal.md 2026-07-15 („Sub-eta-Regression gefunden und gefixt").

<!-- Vorlage:

## <Kurzname> (<Datum>, exp-<NNN> / Run <id>)
- Was versucht wurde:
- Warum es scheiterte (Wurzelursache):
- Nicht wiederholen, außer:
-->

## Warmstart via State-Dump über dps-Grenzen (2026-07-25, exp-039-Pilot, Run 20260725-005345_exp039-pilot)

- Was versucht wurde: einen bei niedrigem dps konvergierten sexpinit-Zustand als
  Startpunkt für einen höheren dps wiederverwenden (writebin/read, wie der
  bestehende State-Cache), um Init-Iterationen zu sparen.
- Warum es scheiterte (Wurzelursache): `sexpinit` → `loop` → `initsch` berechnet
  **unbedingt** jede Zustandsgröße bei der *aktuellen* `realprecision` neu
  (`precis=precision(z)`, `L=fixedk(kd)`, `sfunczero`, `z0h/z0l`, die
  formal(i)schroder-Reihen); danach setzt `loop()` `ct=0; re=-3; n=0`. Ein
  geladener Zustand wird also stillschweigend überschrieben — es existiert kein
  Resume-Einstiegspunkt. Empirisch im *günstigsten* Fall bestätigt (dps-60-State
  live im selben Prozess, ohne writebin-Roundtrip): 11469 ms kalt vs 10672 ms
  warm = **6.9 % Ersparnis** (H2 verlangte ≥ 30 %), Anker bit-gleich. Da schon
  der Bestfall scheitert, kann kein Dump-Schema funktionieren.
- Nicht wiederholen, außer: der Fork bekommt einen echten `loopresume()`-Pfad,
  der die Taylor-/Theta-Arrays aus dem geladenen Zustand *erweitert* statt sie
  neu zu bauen. Das ist Leiter-Stufe 3 (algorithmisch) und eine eigene Mutation
  mit Gate — nicht der hier ausgeschlossene mutationsfreie Pilot.
- Kosten der Erkenntnis: 27 s Rechenzeit (statischer Fund zuerst, dann
  Bestätigung), nicht ein voller Warmstart-Aufbau.
- Werkzeug: `research/tools/warmstart_probe.py` (neu, keine Fork-Mutation).

## ctsamples-Cap zur Init-Beschleunigung (2026-07-10, exp-009 — hier nachgetragen 2026-07-25)

- Was versucht wurde: das Sample-Budget deckeln (`ctsamples`-Cap 5*re+120), weil
  die späten Iterationen mit grossem Grid teuer sind.
- Warum es scheiterte (Wurzelursache): Gate FAIL — e-Gruppen fielen auf ~49
  Stellen. Die ~1988 Terme sind ECHT nötig (~9.8 Terme/Digit für Basis e); die
  4.3er-Schätzung galt für ein anderes ctr-Regime. `stopterms` overshootet nicht.
- Nicht wiederholen, außer: es gibt einen Nachweis, dass Terme über einem
  bestimmten Index nichts zur Genauigkeit beitragen — nicht nur eine Schätzung.
- Nachgetragen, weil exp-040 (H3 supported: späte Iterationen 3.60× teurer) genau
  diesen naheliegenden Fix wieder attraktiv macht. Die Kostenkurve rechtfertigt
  KEIN Kürzen des Sample-Budgets; sie zeigt nur, wo die Zeit liegt.
- Quelle: research/journal.md, exp-009.

## Produktbaum-Multipoint (P0) fuer die ct-Auswertung (2026-07-25, F1-Verdikt)

- Was versucht wurde: die N Horner-Auswertungen pro Iteration durch schnelle
  Mehrpunktauswertung (FLINT/Arb `_acb_poly_evaluate_vec_fast`) zu ersetzen,
  mit ehrlichem Guard-Budget und Radius-Schalen — also genau das, was beide
  Strategiepapiere als „nicht widerlegt, nur unterbudgetiert" fuehren.
- Warum es scheiterte (Wurzelursache): die Fehlerverstaerkung des Produktbaums
  haengt an der Punkt-GEOMETRIE, nicht am Praezisionsbudget. Gemessen auf einem
  echten dps-300-Zustandsdump: ~0.85 Stellen (~2.8 Bit) Verlust PRO PUNKT im
  Block, und ein Guard-Sweep von 0 bis 20480 Bit aendert daran GAR NICHTS
  (1024 von 1280 Werten nicht-endlich, bei jedem Budget). Kontrolle mit
  Zufallspunkten gleicher Groesse: exakt. Ursache: die Auswertungspunkte sind
  die Bilder EINES Kreises unter drei analytischen Karten, liegen also auf
  glatten Boegen — der Worst Case fuer Teilprodukt-Polynome. Radius-Schalen
  normieren den Betrag, nicht die Winkel-Clusterung, und helfen nicht.
- Folge fuer den Exponenten: der Guard-Bedarf waechst linear in N, also kostet
  der Baum c*N*log^2(N)*M(q+2.8N) statt 0.5*N^2*M(q) — bei p=1020 ist er 7.7x
  LANGSAMER als Horner, Crossover erst bei p ~ 12000.
- Nicht wiederholen, ausser: mit Moroz (FOCS 2021, arXiv:2106.02505), der
  quasi-lineare BIT-Komplexitaet ohne den O(d)-Guard erreicht — eigene
  Forschungsimplementierung, nicht in FLINT.
- Werkzeuge: research/tools/dump_state.py, research/tools/multipoint_bench.py.

## Praezisions-Leiter mit Faktor 1 statt 2 (2026-07-25, exp-044)

- Was versucht wurde: `realprecision = m*re + g` mit m=1 statt 2 fuer die
  e-Familie. Motivation war stark: damit wird icdig = realprecision - re + 40
  KONSTANT statt linear in re, was den Arithmetik-Exponenten 1.29 auf ~0 und
  den Gesamtexponenten im Modell von 3.9 auf 3.05 senken wuerde.
- Warum es scheiterte (Wurzelursache): erreichbare echte Stellen saettigen bei
  ~ 115*m + g + 38 — ein FESTER Offset (213.0 Stellen sowohl bei dps 300 als
  auch dps 400 mit m=1/g=60), nicht proportional zu p. Die noetige Marge
  waechst deshalb linear mit p, womit icdig nie konstant wird; bei der
  ausreichenden Marge g=200 ist die Variante bei dps 400 genauso schnell wie
  die Baseline. Mechanismus: die Arbeitspraezision muss `precis` RECHTZEITIG
  erreichen (mit m=2 bei re=(p-60)/2, auf halber Strecke), sonst wird der
  Zustand nie voll aufgeloest.
- Wichtig fuer die Methodik: der Verlust war NUR gegen die bewiesene Referenz
  sichtbar. Engine-gegen-Engine stimmten die Anker auf 213 Stellen ueberein,
  und das Gate haette es wegen der niedrigen e-Schwellen (58.9/59.2)
  durchgelassen. Seitdem: research/tools/digits_vs_reference.py bei jeder
  Mutation, die an Praezision oder Konvergenz ruehrt.
- Nicht wiederholen, ausser: es gibt einen Nachweis, dass der Zustand ohne das
  fruehe Erreichen von precis vollstaendig aufgeloest werden kann.

## Zweigradius ircircr senken (2026-07-25, irmul-Scan)

- Was versucht wurde: `ircircr` (die direct-vs-theta-Grenze in sfunc) senken,
  damit mehr Samples auf den Theta-Zweig fallen — der seit exp-041b gar kein ct
  mehr auswertet und dessen Reihe Grad ~0.43*re statt ~6.2*re hat.
- Warum es scheiterte (Wurzelursache): der Gewinn dreht sich mit der Tiefe um.
  irmul=0.9 gibt bei dps 200 1.143x (und -22 % Iterationen, echter Effekt), bei
  dps 300 aber 0.937x. Naeher am Zentrum konvergiert die Theta-Fourier-Reihe
  schlechter, thsamples waechst, und dieser Term skaliert mit der
  Iterationszahl. Bei irmul=0.85 bricht die Genauigkeit hart ein (92 statt 193
  Stellen) — das ist die Gueltigkeitsgrenze der Theta-Darstellung, und sie
  erklaert, warum das Original `ir` nur nach oben variiert (7/10, 15/16, 24/25).
- Nicht wiederholen, ausser: die Theta-Reihe wird so umgebaut, dass sie naeher
  am Zentrum konvergiert (das waere eine andere Darstellung, keine Knopfdrehung).
- Werkzeug: research/tools/make_irscan.py + ctr_scan.py --var irmul.

## Universelle Theta-Asymptotik als Iterations-Hebel (2026-07-25, F6)

- Was versucht wurde: die hohen Theta-Harmonischen aus einer geschlossenen Form
  mit O(1) Parametern vorhersagen und nur die Abweichung iterieren — der einzige
  Kandidat gegen I(p) = Theta(p) und damit die einzige Route Richtung 2.1.
- Warum es scheiterte (Wurzelursache): die beste 3-Parameter-Form
  log10|theta_m| = A + B*m + C*log(m) trifft out-of-sample auf 1.0e-3 Dezimalen,
  also ~3 RELATIVE Stellen, und das Residuum STEIGT mit m (3.6e-4 bei m=61 auf
  1.5e-3 bei m=109). Mehr Parameter (D/m, E/m^2) machen es out-of-sample
  messbar schlechter. Bei 2.41 Stellen pro Harmonische entspricht das ~einer
  eingesparten Iteration von ~180.
- Methodik-Warnung: IN-SAMPLE faellt das Residuum ueber die Quartile
  (0.0035 -> 0.0001 dez) und sieht damit genau wie der erhoffte Ausgang aus.
  Nur der out-of-sample-Test (trainieren auf m<=60, vorhersagen m=61..115) zeigt
  das Gegenteil. Fits an Asymptotiken IMMER out-of-sample pruefen.
- Nicht wiederholen, ausser: es gibt eine hergeleitete (nicht gefittete) Form
  aus der Singularitaetsstruktur — Darboux/Resurgenz, also F7-Territorium.
- Werkzeug: research/tools/theta_modes.py.

## Reihen-Komposition fuer den exp-Arc (2026-07-25, Praekondition #1)

- Was versucht wurde: die ct-Auswertung auf dem n=-1-Arc (54 % der verbliebenen
  Arbeit) durch eine Reihenkomposition F(u) = ct(A*e^u + B) plus eine FFT
  ersetzen — die Karte ist ganz, Binaersplitting gibt O(M(N) log N), und PARIs
  `fft` genuegt (kein C/FLINT noetig).
- Warum es scheiterte (Wurzelursache): auf dem Sampling-Kreis |u| = r erreicht F
  bereits **10^1141** (dps-300-Zustand), waehrend die auf dem Arc gebrauchten
  Werte O(1) sind. exp dehnt den Bildbereich weit ueber den Konvergenzradius von
  ct hinaus. Eine Darstellung ueber den vollen Kreis muss diese 1141 Dekaden
  fuehren und wieder wegkuerzen — Guard-Bedarf linear in N, also derselbe
  Killer wie beim Produktbaum. Termbedarf zusaetzlich 4.22*N statt der
  geschaetzten 2.9*N.
- **Verallgemeinerung (wichtig fuer kuenftige Ansaetze):** die Punkte liegen auf
  BOEGEN. Jede globale Darstellung — Produktbaum, Komposition ueber den vollen
  Kreis, Niedrigrang-Kompression — zahlt Guard ~ N. Ein Kandidat muss von
  vornherein nur auf dem Arc arbeiten UND O~(N) sein.
- Nicht wiederholen, ausser: Moroz (FOCS 2021) laesst sich auf Bogenpunkte
  uebertragen — das ist offen und waere eine eigene Forschungsfrage.

## Radius-Trunkierung fuer thtaylors Horner (2026-07-25, exp-049)

- Was versucht wurde: die exp-048-Trunkierung auch auf thfunc anwenden. Die
  Theta-Punkte liegen sogar in einem schmalen Band (|p-circc| in [0.567,0.747]),
  eine einzige Stufe genuegt.
- Warum es scheiterte: 1.011x / 0.994x. Die Trunkierung war VERIFIZIERT aktiv
  (deg 1792 -> 885), also ist der Horner dort schlicht nicht der Kostenblock —
  obwohl die Schleife 21.4 % des Laufs ausmacht.
- Nicht wiederholen. Der Kostenblock war das Cache-Verhalten (exp-051).

## thtaylors exp/log-Rundweg entfernen (2026-07-25, exp-050)

- Was versucht wurde: `tcrc[s]=exp(Pi*I*x1)` und das sofortige
  `y=log(z)/(2*Pi*I)` in thfunc sind fuer x1 in (-1,1) exakt invers, also
  y = x1/2 direkt uebergeben (spart pro Sample ein exp und ein log bei voller
  Praezision und ist sogar genauer, da x1 exakt rational ist).
- Warum es scheiterte: 1.015x / 0.997x. Die Transzendenten sind dort nicht der
  Kostenblock.
- Nicht wiederholen — ausser als Aufraeumung, denn die Identitaet ist echt.

## Schaerferer Stall-Exit (2026-07-25, exp-055, H4 aus exp-040)

- Was versucht wurde: exp-010 gibt pro marginaler Iteration EINEN Stall-Kredit
  aus, also verbrennt ein fertig konvergierter Lauf weiterhin ~3 volle spaete
  Iterationen — die teuersten. Die Strafe an den tatsaechlichen Gewinn koppeln
  (`(re-relast) < 0.005` kostet zwei Kredite statt einem).
- Warum es scheiterte: der Ertrag ist Rauschen. dps 300 1.018x, dps 400 1.071x,
  dps 520 1.021x, Benchmark (dps 80/200) ~1.00 — die 400er-Zahl ist nicht
  reproduzierbar. Gate PASS und Stellen ueberall identisch, aber ~2 % liegt
  unter der Rauschschwelle 6.79 %.
- **Prozess-Lektion (wichtiger als das Ergebnis):** ich hatte nach EINEM
  Tiefenpunkt (dps 400) plus plausiblem Mechanismus auf KEEP entschieden und
  damit die eingefrorene Keep-Regel uebergangen, weil der Benchmark bei dps 200
  endet und den Effekt gar nicht sehen kann. Zwei weitere Tiefenpunkte haben
  den Keep kassiert. Regel daraus: wenn der Benchmark blind ist, braucht ein
  Keep MEHRERE Tiefenpunkte — ein einzelner ist nicht besser als ein
  Benchmark-Ausreisser.
- Nicht wiederholen, ausser: der Stall-Exit wird zusammen mit einer echten
  Konvergenz-Prognose umgebaut (dann ist es eine andere Mutation).
