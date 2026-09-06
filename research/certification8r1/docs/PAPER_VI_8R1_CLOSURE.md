# Paper VI – `8r1` vom Pre-Certificate zum hergeleiteten Segmentzertifikat

Stand: 29. Juli 2026

Bezug: Auftrag „alle gesetzten Sicherheitsradien durch hergeleitete
Intervallabschätzungen ersetzen“, Punkte 1–5.

## 0. Ergebnis in drei Sätzen

Die β-Taylorreste, die Koenigs-Endrekursion und die Tailregel sind jetzt
hergeleitet, implementiert und von einem produzentenunabhängigen Programm
neu erzeugt; die hergeleiteten Schranken sind zwischen zwei und dreihundert
Größenordnungen schärfer als die zuvor gesetzten Konstanten. Beim Aufräumen
des β-Modells ist ein Implementierungsdefekt des Pre-Certificates
aufgefallen: die deklarierte Taylorordnung 22 wurde nie erreicht, gerechnet
wurde mit Ordnung 10. Zwei der fünf Punkte sind damit noch nicht
abgeschlossen — die graded-Gauss–Legendre-Schranke ist hergeleitet, aber der
zugehörige komplexe Tubus-Sup ist noch nicht durchgerechnet, und der
Mikrosegment-Punkt hat sich als inhaltliche und nicht bloß numerische Lücke
herausgestellt (Abschnitt 8).

## 1. Befund: die deklarierte β-Taylorordnung wurde nie erreicht

`certify_8r1.py` setzt `BETA_SERIES_ORDER = 22` und schreibt
`"beta_series_order": 22` in die Payload; `replay_8r1.py` liest daraus
`order = 22` und rekonstruiert 22-gliedrige Modelle.

`python-flint` schneidet jedoch **jede** Serienoperation bei der globalen
Kappungslänge `ctx.cap` ab. Deren Default ist `10`, und keines der drei
Programme setzt sie. Damit sind sämtliche β-Modelle Polynome vom Grad 9.
Nachweis am ausgelieferten Datensatz:

```
segment_8r1_interval_primitives.json, panel 0, row 5
  Gamma[9]  = [-2.52477822350329141e-25, -2.52477822350329049e-25]
  Gamma[10] = [-4.94065645841246544e-324, 4.94065645841246544e-324]
  Gamma[21] = [-4.94065645841246544e-324, 4.94065645841246544e-324]
```

Die Koeffizienten 10…21 sind exakt null (±4.94e-324 ist die
outward-Rundung der Null); dasselbe gilt für `tail`. Die freie Konstante
`BETA_SERIES_REMAINDER = 2e-16` sollte also eine Grad-9-Abschneidung
abdecken, während Payload und Replay eine Grad-21-Abschneidung behaupten.

Numerisch war 2e-16 dabei vermutlich ausreichend — die tatsächlichen
Koeffizienten fallen schnell —, aber die deklarierte Ordnung stimmte nicht
mit der gerechneten überein, und genau diese Diskrepanz macht eine gesetzte
Restkonstante gefährlich: sie verdeckt, welche Größe sie eigentlich
abdecken muss.

Die Closure setzt `ctx.cap` explizit (`tmodel.configure`) und prüft nach
jeder Serienoperation die erreichte Länge (`tmodel.assert_full_length`).

## 2. Inventar der freien Konstanten im Beweispfad

| Konstante | Ort | Rolle | Status nach dieser Arbeit |
|---|---|---|---|
| `BETA_SERIES_REMAINDER = 2e-16` | `certify_8r1.py` | Abschneidefehler aller β-Taylormodelle | **hergeleitet**, Abschnitt 3 |
| `KOENIGS_GAMMA_FACTOR = 128` | `certify_8r1.py` | Γ-Rest der gestoppten Koenigs-Rekursion | **hergeleitet**, Abschnitt 4 |
| `KOENIGS_GAMMA_PRIME_FACTOR = 1024` | `certify_8r1.py` | Γ′-Rest | **hergeleitet**, Abschnitt 4 |
| Tailregel `R ≤ 2·|a_N|`, Ableitungen `4R` | `certify_8r1.py` | Abschneidefehler des Reziprokalprodukt-Tails | **bewiesen**, Abschnitt 5 |
| `HILBERT_REMAINDER = 2e-10` | `certify_hilbert_8r1.py` | graded-Gauss–Legendre-Fehler | **hergeleitet, Rechnung offen**, Abschnitt 6 |
| `ENDPOINT_COLLAPSE_CONSTANT = 128` | `certify_hilbert_8r1.py` | Endpunktanteil derselben Schranke | **hergeleitet, Rechnung offen**, Abschnitt 6 |
| `COMPLEX_DERIVATIVE_INFLATION = 16` | `replay_8r1.py` | Tubusinflation von M_Γ′ | **hergeleitet, Rechnung offen**, Abschnitt 6.4 |
| `M1_UPPER = 0.01` | `replay_8r1.py` | globaler Grönwall-Exponent | **offen — inhaltliche Lücke**, Abschnitt 8 |
| `RADIUS_SAFETY = 1.25` | `replay_8r1.py` | Wahl von r | **keine Schranke, sondern Wahl** — zulässig, Abschnitt 9 |
| `RHO_*`, `SIGMA_GEOMETRIC`, `TARGET_Q_UPPER`, `0.75` in κ | `replay_8r1.py` | Skalenwahl | **Entwurfsparameter**, zulässig, Abschnitt 9 |

Die letzten beiden Zeilen sind bewusst nicht als Mangel geführt: `r` ist
frei wählbar, solange `p(r) < 0` nachgewiesen wird, und die Streifenbreiten
sind Parameter des Ansatzes, deren Zulässigkeit über die Z-Margen geprüft
wird. Alles darüber ist ein Restbudget und muss hergeleitet sein.

## 3. Satz 1 — zertifizierte β-Taylormodelle

**Aufbau.** Ein Taylormodell der Ordnung `N` auf der abgeschlossenen
Einheitsscheibe des skalierten Panelparameters `s` ist ein Paar `(P, r)` mit

```
    f(s) ∈ P(s) + D(0, r)    für alle |s| ≤ 1,
```

wobei `P` Koeffizienten als Arb-Kugeln trägt. `β = β_mid + β_half·s`.

**Satz 1.1 (algebraische Operationen).** Für Modelle `(P₁,r₁), (P₂,r₂)` ist

```
    (P₁P₂)_{<N},   r = Σ_{k≥N} |(P₁P₂)_k| + ‖P₁‖₁ r₂ + ‖P₂‖₁ r₁ + r₁r₂
```

ein gültiges Modell des Produkts, wobei `‖·‖₁` die Koeffizientensumme ist
(eine gültige Schranke auf `|s| ≤ 1`). Der hohe Teil der Faltung wird
vollständig berechnet, nicht abgeschätzt.

**Satz 1.2 (Kehrwert).** Sei `h` die abgeschnittene Reihenumkehr von `P` und
`D = 1 − P·h` der vollständig berechnete Defekt. Ist `‖D‖₁ < 1` und
`m = min_{|s|≤1} |P| > r`, so gilt

```
    1/f ∈ h + D(0, ‖h‖₁‖D‖₁/(1−‖D‖₁) + r/(m(m−r))) .
```

Der erste Term ist der Abschneidefehler, der zweite der Mittelwertsatz für
die Störung `r`. Beide sind exakt, keiner ist geschätzt.

**Satz 1.3 (exp, log; Cauchytail).** Ist `g` auf einer Umgebung von
`P(D̄_R)` holomorph und `M_R = sup_{|s|=R} |g(P(s))|`, so erfüllen die
Taylorkoeffizienten von `g∘P` die Cauchyabschätzung `|a_k| ≤ M_R R^{−k}` und
damit

```
    Σ_{k≥N} |a_k| ≤ M_R R^{−N}/(1 − R^{−1}) .
```

`R` wird aus einer Leiter gewählt, die die Schranke minimiert; die
Zulässigkeit jedes `R` wird aus den Koeffizienten selbst zertifiziert
(für `log`: `Σ_{k≥1}|c_k|R^k < |c_0|`, was `P ≠ 0` auf `D̄_R` beweist und
zugleich `M_R ≤ |log c_0| − log(1 − ρ/|c_0|)` liefert).

**Satz 1.4 (Zweigtreue).** `log` wird stets über die Drehung
`log f = log(f e^{−iθ}) + iθ` mit `θ = arg(Mittelpunkt von c₀)` ausgewertet.
Das ist eine Identität, verschiebt die Einschließung aber vom
Hauptzweigschnitt weg. Ohne diese Drehung bricht die Rekursion zusammen:
der Rückwärts-Koenigs-Orbit **landet konstruktionsbedingt auf der negativen
reellen Achse**, bevor er zum Fixpunkt spiralt, und die Hauptzweig-Auswertung
einer Kugel um diesen Punkt ist unbrauchbar weit. Die Drehung wählt
gleichzeitig die analytische Fortsetzung durch den Modellmittelpunkt, also
den Zweig des reellen Kanals.

**Satz 1.5 (Fixpunkt).** `L` mit `L = exp(cL)` wird nicht aus `W_{−1}`
übernommen, sondern a posteriori zertifiziert: ist `φ(z) = log(z)/c` die
kontrahierende Umkehrung, `d = ‖φ(L̃) − L̃‖` und `κ = 1/min|cL̃| < 1` auf der
Scheibe vom Radius `d/(1−κ)` um `L̃`, so liegt der eindeutige Fixpunkt in
dieser Scheibe. Die Kontraktion wird auf der aufgeblasenen Scheibe erneut
geprüft.

**Gemessene Reste** (Panel 0, Ordnung 22, 448 Bit), gegen die bisherige
freie Konstante 2e-16:

| Größe | hergeleitet | gesetzt |
|---|---|---|
| `exp(β)` | 5.29e-80 | 2e-16 |
| `1/exp(β)` | 1.06e-79 | 2e-16 |
| `L(β)` | 1.29e-45 | 2e-16 |

Alle Modellinvarianten sind gegen unabhängige Punktauswertung auf neun
Stützstellen der Einheitsscheibe geprüft (`test_tmodel.py`), inklusive
Negativkontrollen.

## 4. Satz 2 — Koenigs-Endrekursion

**Aufbau.** `φ(w) = log(w)/c`, Fixpunkt `L`, `φ′(L) = 1/λ`, `|λ| > 1`,
`δ_n = w_n − L`, und

```
    ψ(δ)      = φ(L+δ) − L − δ/λ            (doppelte Nullstelle in 0)
    Φ_β(δ)    = ∂_β [ φ_β(L_β+δ) − L_β ]    (einfache Nullstelle in 0)
```

Beide werden als Taylormodelle ausgewertet, deren konstanter Term auf die
abgeschlossene Scheibe `|δ| ≤ ρ` aufgeblasen ist. Die Schranken gelten damit
gleichzeitig uniform in `δ` und im Basisparameter. Mit
`K = sup|ψ|/ρ²` und `B = sup|Φ_β|/ρ` liefert das Schwarzsche Lemma

```
    |ψ(δ)| ≤ K|δ|²,   |ψ′(δ)| ≤ 2K|δ|(1 + 2|δ|/ρ),
    |ξ(δ)| ≤ B,       |ξ(δ) − ξ(0)| ≤ 4B|δ|/ρ,      ξ(δ) = Φ_β(δ)/δ .
```

**Satz 2.1 (Wert).** Aus `δ_{m+1} = δ_m/λ + ψ(δ_m)` folgt
`λ^{m+1}δ_{m+1} = λ^m δ_m (1 + λ η(δ_m) δ_m)` mit `|η| ≤ K`. Mit
`θ = 1/|λ| + K|δ_n| < 1` und `t = |λ| K |δ_n|/(1−θ) < 1/2` gilt

```
    | log κ − log S^{(n)} | ≤ t/(1−t) =: E_κ .
```

**Satz 2.2 (Zustandsableitung).** Für `g_m = s_m/δ_m` ist

```
    log g_{m+1} − log g_m
      = log(1 + λψ′(δ_m)) − log(1 + λη(δ_m)δ_m),
```

also mit `u = |λ| K |δ_n| (3 + 4|δ_n|/ρ)/(1−θ) < 1/2`

```
    | log g − log g_n | ≤ u/(1−u),
    | g − g_n | ≤ |g_n| (exp(u/(1−u)) − 1) .
```

**Satz 2.3 (Basisableitung).** Für `G_m = (v_m − dL)/δ_m` und
`dlogS^{(m)} = m·dlogλ + G_m` ist

```
    dlogS^{(m+1)} − dlogS^{(m)}
      = [ (ξ(δ_m) − ξ(0)) + dlogλ·η δ_m + G_m(ψ′(δ_m) − η δ_m) ]
        / (1/λ + η δ_m),
```

weil `ξ(0) = −dlogλ/λ` den führenden Term exakt aufhebt. Summation mit
`|G_m| ≤ |G_n| + (m−n)·2|dlogλ|` und `|δ_m| ≤ |δ_n|θ^{m−n}` gibt die in
`KoenigsBase.stopping_bounds` implementierte Schranke `E_β`.

**Satz 2.4 (Weitergabe auf Γ).** Mit `Γ = N/P`,
`N = ∂_β log κ − χ ∂_β log λ`, `P = ∂_y log κ`:

```
    |Γ − Γ_n| ≤ ( E_β + |dlogλ| E_κ/|log λ| + |Γ_n| E_P ) / ( |P_n| − E_P ) .
```

`Γ′` wird nicht separat teleskopiert: derselbe Quotientenfehler wird auf
einer abgeschlossenen Zustandsscheibe vom Radius `ρ_y` erneut ausgewertet,
und die Cauchyabschätzung `|∂_y f| ≤ sup_{Scheibe}|f|/ρ_y` überführt ihn.

**Gemessen** über alle 6 Panels × 96 Knoten (`stop = 1e-18`, `ρ = 0.5`,
`ρ_y = 1e-3`, Ordnung 22, 448 Bit, Laufzeit 1031 s):

| Größe | hergeleitet (Worst Case) | Pre-Certificate (gesetzt) |
|---|---|---|
| `E_κ` (log κ) | 2.50e-17 | — |
| `E_β` (Basisableitung) | 1.31e-14 | — |
| Γ-Rest | 1.70e-14 | < 1.77e-12 |
| Γ′-Rest | 1.78e-11 | < 3.81e-10 |
| Phasenmarge | 1.844 rad | 1.906 rad (nur Mittelpunkt) |

Die Phasenmarge ist bei uns **über das ganze β-Panel** zertifiziert, im
Pre-Certificate nur am Mittelpunktskoeffizienten
(`delta.coeffs()[0].arg()`). Der niedrigere Wert ist also die stärkere
Aussage. Zusätzlich ersetzt die Closure die geratene Windungszahl durch eine
stetige Logarithmusfortsetzung: statt `turn` zu runden, wird
`log δ_{n+1} = log δ_n + log(δ_{n+1}/δ_n)` fortgeschrieben, wobei das
Verhältnis nahe `1/λ` und damit weit vom Schnitt bleibt.

## 5. Satz 3 — Tailregel

**Satz 3.1.** Für `a_j = 1/(c y_j′)` mit `y_{j+1} = e^{c y_j}`,
`y_{j+1}′ = c y_{j+1} y_j′` gilt die **exakte** Identität

```
    a_{j+1}/a_j = y_j′/y_{j+1}′ = 1/(c y_{j+1}) .
```

**Satz 3.2 (Vorwärtsinvarianz).** Sei `c > 0` reell, `c ≤ 2`, und
`c y_N ≥ log(2/c)`. Dann ist `y_{N+1} = exp(c y_N) ≥ 2/c`, also
`c y_{N+1} ≥ 2 ≥ log(2/c)`, und die Voraussetzung reproduziert sich. Folglich
ist `(y_j)_{j≥N}` reell, wachsend und `≥ 2/c`, also `|a_{j+1}/a_j| ≤ 1/2`
und

```
    | Σ_{j>N} a_j | ≤ |a_N| Σ_{k≥1} 2^{−k} = |a_N| ,
```

die Hälfte der zuvor angenommenen Regel. Die Prämissen `c ≤ 2` und
`c y_N ≥ log(2/c)` werden als Intervallaussagen geprüft; ohne sie
verweigert die Implementierung die Abschneidung.

**Satz 3.3 (Abschneideindex).** Das Verhältnis zu *jedem* späteren Term ist
`1/(c y_k)` mit wachsendem `y_k`, also genügt eine untere Schranke für
`|c y_{j}| = |c| exp(Re(c y_{j−1}))`. Diese ist ein Skalar und braucht das
nächste Turmniveau **nicht** als Taylormodell.

Das ist nicht bloß eine Optimierung, sondern notwendig: über einem
β-Panel der Halbbreite 5e-4 wächst schon `d(c y_4)/dβ ≈ 1.8e1`, beim
nächsten Niveau `≈ 3.2e7`, sodass `y_5` über dem Panel um den Faktor
`exp(1.6e4)` variiert. Ein Taylormodell dieser Funktion existiert praktisch
nicht. Das Pre-Certificate rechnet an dieser Stelle sieben Turmniveaus in
`acb_series` über dem β-Panel und deckt den Fehler mit `2e-16` ab; die
Modelle der oberen Niveaus tragen dort keine belastbare Information mehr.

**Gemessen** über alle 576 Zeilen: Verhältnisschranke im Worst Case
`7.12e-9` (weit unter 1/2), Tailrest im Worst Case `1.14e-10`, an den
meisten Knoten astronomisch kleiner (z. B. `1.46e-32366875` bei `x = 0`).
Der Worst Case entsteht dort, wo die Repräsentierbarkeitsschranke aus
Satz 3.3 den Abschluss ein Turmniveau früher erzwingt. Er ist unschädlich:
`S` geht über `G = R − S` in das Residuum ein, mit Wienerfaktor ≈ 591 also
`6.7e-8` in `residual_head` und `8.7e-16` in `Y_head` — gegen
`Y_tail = 2.89e-10` vernachlässigbar. In `M_tail ≈ 3.9` ist `1.14e-10`
ebenfalls belanglos.

**Erratum (2026-09-06), Ableitungstails.** Satz 3.1–3.3 betreffen den
*Wert* `S`. Die Implementierung `tails.py` behauptete zusätzlich, die beiden
Ableitungstails `∂_T S` und `∂_{T′} S` seien durch *dieselbe* geometrische
Reihe wie der Wert dominiert. Für `∂_{T′}` stimmt das (der Turm hängt nicht
von `T′` ab, das Verhältnis ist exakt `1/(c y_{j+1})`). Für `∂_T` ist es
falsch: aus `∂_y y_{j+1} = c y_{j+1} ∂_y y_j` folgt

```
    ∂_y a_{j+1} / ∂_y a_j = (1 + c y_j θ_{j−1}) / (c y_{j+1}),   θ ∈ [0,1],
```

also ein um bis zu `1 + c y_j` — am Abschneideindex dieses Segments etwa 20 —
größeres Verhältnis. Gefunden durch Wertevergleich der Closure-Einschließungen
mit denen des Konstanten-Produzenten, der ein Turmniveau mehr explizit
summiert: 254 von 2880 Stichproben von `tail_partial_T` lagen disjunkt; an
der schlimmsten Zeile (Panel 3, Zeile 5) übertraf der erste vernachlässigte
Term `8.29e-9` das behauptete Restglied `6.16e-10` um den Faktor 3. Korrigiert
in `tails.py` (Schranke `(1 + u_N)/(c e^{u_N})` mit der zertifizierten unteren
Schranke `u_N` von `c y_N`, Monotonie des Turms oberhalb `log(2/c)`, geprüfte
Prämisse `log(2/c) < 2`), Regressionstest `test_tails.py`. Die Zeile 2c der
Statusübersicht gilt seither auch für die Ableitungen; vorher galt sie nur
für den Wert.

## 6. Satz 4 — graded Gauss–Legendre (hergeleitet, Rechnung offen)

### 6.1 Was approximiert wird

Der Hilbert-Anteil der RH-Selektion ist

```
    H[u](x) = ∫₀¹ (u(t) − u(x)) cot(π(x−t)) dt ,    u = Im v,  v = Γ(T)/T′,
```

genähert durch 36 dekadengraduierte Panels mit je 24 Gauss–Legendre-Punkten
auf `[t_min, 1−t_min]`, `t_min = 1e-18`. Der Gesamtfehler zerfällt in

```
    |H − Q| ≤ | ∫₀^{t_min} F | + | ∫_{1−t_min}¹ F | + Σ_p | ∫_{P_p} F − GL₂₄(F,P_p) | .
```

### 6.2 Feldprobe: die Struktur des Integranden

Eine hochgenaue Probe (`probe_field.py`, mpmath, 60 Stellen, Basisknoten 3)
ergibt:

* auf der reellen Achse ist `|v| ≤ 3.83`, Maximum bei `t = 0.5`;
* am Rand gilt `v(t) ≈ t·(A + B log t)` — die Dekadenverhältnisse fallen von
  9.27 auf 9.18 pro Dekade, konsistent mit `A ≈ 9.9 B`, `B ≈ 30`; es liegen
  also **logarithmische Verzweigungspunkte bei `t = 0` und `t = 1`** vor,
  was die Dekadengraduierung rechtfertigt;
* außerhalb der Achse bleibt `|v|` bis Höhe 1.2 beschränkt (Werte 1.3–3.7),
  es gibt also **keine bindende Streifenbreite**;
* `Im v(0) = Im v(1) = 0`, der Integrand ist an beiden Enden hebbar bzw.
  höchstens logarithmisch singulär.

Bindend für die Bernsteinellipsen sind damit die Verzweigungspunkte bei 0
und 1, nicht die Analytizitätsbreite — anders als die Wahl
`SIGMA_GEOMETRIC = 0.10` nahelegt.

### 6.3 Die Schranke

**Satz 4.1 (Panelfehler).** Ist `F` auf der abgeschlossenen Bernsteinellipse
`E_ρ` des auf `[−1,1]` transformierten Panels holomorph mit `|F| ≤ M_p`, so
gilt für die `n`-Punkt-Gauss–Legendre-Regel (Trefethen, ATAP Thm. 19.3)

```
    | ∫_{P_p} F − GL_n(F,P_p) | ≤ h_p · (64/15) · M_p ρ^{−2n} / (1 − ρ^{−2}) .
```

Für die dekadengraduierten Panels ist `d_p/h_p = 5.5/4.5 = 1.2̄`, also
`(ρ+ρ^{−1})/2 < 1.2̄` und `ρ ≲ 1.92`; mit `n = 24` ist
`ρ^{−48} ≈ 2.8e-14`. Für das Panel `[0.1,0.5]` ist `d/h = 1.5`, `ρ ≤ 2.0`,
`ρ^{−48} ≈ 3.6e-15`.

**Satz 4.2 (Faktorisierung).** Für Panels, deren Ellipse den Pol `t = x`
nicht trifft, ist

```
    M_p ≤ ( sup_{E_p} |u| + |u(x)| ) · sup_{E_p} | cot(π(x−t)) | ,
```

und der Kotangensfaktor ist elementar. Für das Panel, das `x` enthält, wird
`F = D·K` mit `D(t) = (u(t)−u(x))/(t−x)` und
`K(t) = −(1/π)·ζ cot ζ`, `ζ = π(x−t)`, zerlegt; `K` ist auf `|ζ| < π`
holomorph und beschränkt, `D` wird über die Cauchy-Darstellung der
dividierten Differenz beschränkt.

**Satz 4.3 (Endstücke).** `|∫₀^{t_min} F| ≤ t_min · sup_{[0,t_min]}|F|`, und
wegen `v(t) = O(t log t)` ist dieser Sup durch
`t_min·(A + B|log t_min|)·sup|cot|` beschränkt, also von der Größenordnung
`1e-17·|cot(πx)|`. Das ersetzt `ENDPOINT_COLLAPSE_CONSTANT = 128`, die im
Pre-Certificate ohnehin nur gegen `HILBERT_REMAINDER` geprüft und dann von
ihr abgezogen wurde — die Zeile
`panel_and_rounding_remainder = HILBERT_REMAINDER − endpoint_remainder`
definiert den Panelanteil als Rest einer gesetzten Zahl und ist zirkulär.

### 6.4 Was zur Ausführung fehlt

`M_p` verlangt `sup_{E_p} |u|`, also eine zertifizierte Schranke des
Koenigs-Feldes auf einer **komplexen Region**, gleichmäßig über das
β-Panel. Drei Wege wurden geprüft:

1. **β als gewöhnliches Intervall.** Scheitert vollständig. Schon Radius
   `1e-6` bricht die Rekursion („Koenigs step ratio leaves the principal
   sector“). Der Wrapping-Effekt über ~200 Schritte ist der Grund, warum das
   Pre-Certificate zu Recht `acb_series` in β verwendet.
2. **`t` als Taylormodell über ein ganzes Panel.** Scheitert an der
   Koeffizientenkonvergenz: bei `d/h = 1.2̄` ist der Cauchyradius kleiner als
   2, die Tailschranke also nicht kontraktiv.
3. **Zweivariables Taylormodell (β × z).** Funktioniert. Implementiert in
   `tmodel2.py` als Modell in `τ` mit einvariablen Modellen in `s` als
   Koeffizienten; alle Skalaroperationen delegieren an den getesteten
   1-D-Code. Der Orbit überlebt für z-Scheiben vom Radius 0.12 um `z = 1.5`
   bis Tiefe 102 in 4–8 s. Nahe den Enden (`z ≈ 1`, `z ≈ 2.07`) müssen die
   Scheiben kleiner gewählt werden.

Weg 3 ist damit gangbar; ausstehend sind die Ellipsengeometrie pro Panel,
die Überdeckung des Tubus, der Kotangensfaktor und die Summation. Der
Aufwand ist rechnerisch etwa eine Stunde, die Implementierung deutlich mehr.

**Sensitivität.** Der Beitrag zum Zertifikat ist gering: `HILBERT_REMAINDER`
geht als reeller Ball in jeden der 96 Hilbertwerte und über die
Wienernorm mit Faktor ≈ 591 in `residual_head` ein. Mit
`Y_head = cell_length·√ρ₀·residual_head = 1.303e-8·residual_head` und
`Y_tail = 2.894e-10` gilt

| ε | Beitrag zu `residual_head` | `Y_head` | Anteil an `Y` |
|---|---|---|---|
| 2e-10 (gesetzt) | 1.18e-7 | 1.5e-15 | 5e-6 |
| 1e-13 (Bernstein, erwartet) | 5.9e-11 | 7.7e-19 | 3e-9 |
| 1e-7 (grobe Schranke) | 5.9e-5 | 7.7e-13 | 0.3 % |

`Y` wird von `Y_tail` dominiert; selbst eine um drei Größenordnungen
schlechtere Schranke als die gesetzte würde das Zertifikat nicht kippen.
Die Herleitung ist also eine Frage der Strenge, nicht der Schärfe — was sie
nicht weniger notwendig macht.

Dieselbe Tubusrechnung liefert `M_Γ` und `M_Γ′` auf dem komplexen Tubus und
ersetzt damit `COMPLEX_DERIVATIVE_INFLATION = 16`, das im Replay ohne
Begründung als Inflationsfaktor von `M_Γ′` auftritt.

## 7. Punkt 3 — Q_β aus dem Anker, SP/SP′

**Stand im Pre-Certificate.** Der Replay prüft die Rohtransformation nur als
Identität an den 48 Zustandsknoten:
`U(x_i) = x_i + q + h(x_i)`, `α = log c_anchor/log b`,
`c₀ = log α/log b`. Er rekonstruiert `Q_β` **nicht** durch Peeling und prüft
`T = Q_β ∘ U` **nicht** auf der Kandidatenkugel. Punkt 3 ist damit offen.

**Was zu zeigen ist.** Nach Definition SP/SP′ (Paper VI, Def. 9.1/8.1) und
Satz „Exact base derivative of the raw branch“:

```
    z_N = α_β T_c(N+x) + c_{0,β},     y_{k+1} = log_{b_β} y_k,   y_0 = z_N,
    ∂_β Q_β = −Q_β′ S_{Q,β},
    S_{Q,β}(x) = Q_β(x)/Q_β′(x) + Σ_{j≥1} 1/( c^j Q_β′(x) Π_{r=1}^{j−1} Q_β(x+r) ) .
```

Ein SP-Zertifikat verlangt Intervall-Untergrenzen aller Zwischenmoduln, die
Abstände zu den logarithmischen Schnitten, Ableitungsschranken und eine
summierbare Tailmajorante `ε_N`; SP′ zusätzlich `|∂_β z_N| ≤ Z_N` mit
`Z_N/|∂_x z_N| → 0` und `Σ ε_N^{(β)} < ∞`.

Die Bausteine dafür liegen mit `tmodel`/`tails` vor — die Peeling-Kette ist
strukturell dieselbe Rekursion wie der Tail, und Satz 3.2 liefert die
Vorwärtsinvarianz. Die Ausführung fehlt.

## 8. Punkt 4 — Mikrosegment-Induktion: eine inhaltliche Lücke

Dies ist der wichtigste Befund dieses Abschnitts, und er ist kein
Rechenproblem.

**Was der Replay tut.** Er bestimmt eine Zellzahl `N` durch

```
    κ = 0.75(ρ₀ − ρ_end)/a_cell ,    q = 4√2 L/κ < 0.24
    ⇒ a_cell < 9.35e-8 ,  N = 262144 ,
```

berechnet `Y`, `r`, `ε_cell = r/√(ρ₀ − κ a_cell − ρ_end)` **für eine Zelle**
und setzt dann

```
    endpoint_total = N · ε_cell · exp(M₁ · |β_b − β_a|) ,   M₁ = 0.01 .
```

Das ist eine lineare Aufsummierung mit einem einzigen globalen
Grönwall-Faktor. Es wird **keine** Inklusion „Ausgangskugel von Zelle `i`
⊂ zulässige Eingangskugel von Zelle `i+1`“ geprüft.

**Warum die naheliegende Induktion nicht durchgeht.** Der Satz
(thm:radiiF) liefert die Ausgangsschranke im Streifen `ρ_end = 0.01`. Die
nächste Zelle braucht ihre Eingangskugel im Streifen `ρ₀ = 0.05`. Eine
Schranke in `A_{0.01}` liefert keine in `A_{0.05}`. Der Streifen wird pro
Zelle um `κ·a_cell = 0.75(ρ₀−ρ_end) = 0.03` verbraucht und im Ansatz nicht
wiederhergestellt.

**Warum der fixe-Streifen-Grönwall die Lücke nicht schließt.** `M₁` müsste
eine Lipschitzkonstante von `Ψ_β` auf einem *festen* Streifen sein. Das
Vektorfeld `Ψ_β[h] = −(1+h′)G_β[q+h] − q′` verliert aber eine Ableitung:
auf festem Streifen ist es nicht Lipschitz. Genau deshalb ist die Theorie in
einer Banachraum-Skala (Ovsyannikov) formuliert. Setzt man stattdessen die
Skalen-Lipschitzkonstante ein, so ist der Ein-Zell-Faktor

```
    exp( L a_cell/(ρ₀ − ρ_end) ) = exp( 13613 · 5.83e-8 / 0.04 ) = 1.0198 ,
```

und über `N = 262144` Zellen `exp(5192)` — das Zertifikat wäre wertlos.

**Warum die globale Pfadnorm allein es auch nicht tut.** Die Alternative des
Auftrags — „eine direkt zertifizierte globale Pfadnorm“ — verlangt
`a < ρ₀/(4√2 L)` (Kor. „A priori segment length“). Mit `L = 13613` und
`ρ₀ = 0.05` ergibt das `a < 6.5e-7`, während das Segment `1.53e-2` lang ist,
ein Faktor 23500. Ein einziger globaler Ovsyannikov-Schritt ist mit diesem
`L` also ausgeschlossen.

**Was bleibt.** Genau eine von zwei Möglichkeiten:

1. **`L` um mindestens Faktor 2.4e4 senken.** `L = 2(M_G + (1+‖h′‖)C_G)` ist
   konservativ aufgebaut (`C_RH = 2(C_geo + 2C_trace)` usw.); ein
   schärferer Intervall-Lipschitzschätzer könnte hier viel gewinnen, aber
   vier Größenordnungen sind nicht plausibel.
2. **Echte Re-Ankerung pro Zelle mit Streifenwiederherstellung.** Das ist
   der im Paper vorgesehene Weg (thm:cont): am Zellende wird re-ankert, der
   Zustand liegt wieder in `A_{ρ₀}`, und die Induktion braucht pro Zelle
   eine zertifizierte Inklusion `E_out + D_centers + m_strict ≤ R_next`.
   Diese Inklusion ist derzeit weder formuliert noch für 8r1 gerechnet — der
   Review-Auflösungsbericht führt sie selbst als offen („Diese globale
   Handoff-Inklusion ist für 8r1 noch offen“), allerdings nur für die
   Segment-, nicht für die Zellebene.

Solange keiner der beiden Wege ausgeführt ist, trägt der Schritt von der
Zell-Radienpolynom-Aussage zur Segment-Endpunktschranke `endpoint_total`
nicht. `M1_UPPER = 0.01` ist dabei nicht nur ungerechtfertigt, sondern
steht für eine Konstante, die in dieser Form nicht existiert.

## 9. Was zu Recht keine Herleitung braucht

* `RADIUS_SAFETY = 1.25`: `r = 1.25 Y/(1−q)` ergibt
  `p(r) = Y + (q−1)r = −0.25 Y < 0`. `r` ist frei wählbar; die Schranke ist
  `p(r) < 0`, und die wird geprüft. Das ist eine Wahl, kein Restbudget.
* `RHO_STATE`, `RHO_ZERO`, `RHO_TRACE`, `RHO_WORK`, `RHO_END`,
  `SIGMA_GEOMETRIC`, `TARGET_Q_UPPER` und der Faktor `0.75` in `κ`:
  Entwurfsparameter der analytischen Skala. Ihre Zulässigkeit ist genau der
  Inhalt der Z1–Z7-Margen. Sie dürfen gesetzt sein, solange die Margen
  positiv nachgewiesen werden — was der Replay tut.
* `KOENIGS_STOP`, `TAIL_STOP`, `KOENIGS_MAX_DEPTH`, `RESIDUAL_HEAD`,
  `HUB_CUTOFF`: Abbruch- und Ordnungsparameter. Sie sind zulässig, sobald
  der zugehörige Rest hergeleitet ist — was für Koenigs und Tail jetzt der
  Fall ist.

## 10. Unabhängiger Primitivverifier

`verify_primitives_8r1.py` teilt keinen Code mit `produce_8r1.py`,
`certify_8r1.py` oder `certify_hilbert_8r1.py`. Er liest ausschließlich die
Rohtrajektorie — die sieben Lobatto-Zustände und die Segmentkonventionen —
und erzeugt die Koenigs- und Tailprimitive neu, mit hergeleiteten Resten.
Anschließend vergleicht er seine Einschließungen an fünf Stützstellen des
β-Panels mit den gespeicherten und schreibt einen eigenen, hashgebundenen
Primitivdatensatz.

Damit ist die Forderung „der Replay erzeugt die primitiven Koenigs-, Tail-
und Hilbertintervalle aus den Rohdaten selbst neu **oder** liefert einen
zweiten, unabhängig implementierten Primitive-Verifier“ für Koenigs und Tail
erfüllt; für Hilbert fehlt der Quadraturrest aus Abschnitt 6.

### 10.1 Einschließungsvergleich

`compare_enclosures.py` prüft anschließend die eigentliche Frage: Der
Vergleich ist nicht „sind die Koeffizienten gleich“ — das können sie nicht
sein, weil das Pre-Certificate mit effektiver Ordnung 10 und einer flachen
Restkonstanten rechnete —, sondern „überlappen die beiden Einschließungen
derselben Funktion“. Für jedes Modell und jede der fünf β-Stützstellen wird

```
    gap(s) = min | P_hergeleitet(s) − P_gespeichert(s) |
```

gegen den hergeleiteten Rest gehalten. Ergebnis über alle
6 × 96 × 5 = 2880 Auswertungen:

| Größe | Worst Gap | hergeleiteter Rest dort | enthalten |
|---|---|---|---|
| `T`, `T′`, `T″` | 0 | 0 | ja |
| `Γ` | 0 | 0 | ja |
| `Γ′` | 0 | 0 | ja |
| `tail` | 1.1228e-10 | 1.1423e-10 | ja |

Null Verletzungen. Die Primitive des Pre-Certificates liegen also innerhalb
der unabhängig hergeleiteten, zertifizierten Einschließungen — die gesetzten
Konstanten waren numerisch nicht falsch, sie waren unbewiesen. Genau diese
Unterscheidung war der Auftrag.

## 11. Statusübersicht

| Punkt | Inhalt | Hergeleitet | Implementiert | Neu gerechnet |
|---|---|---|---|---|
| 1 | graded Gauss–Legendre, `HILBERT_REMAINDER` | ja (Abschn. 6) | Werkzeuge ja (`tmodel2`), Pipeline nein | nein |
| 2a | β-Taylorreste | ja | ja | ja |
| 2b | Koenigs-Endrekursion | ja | ja | ja |
| 2c | Tailregel `R ≤ 2\|a_N\|` | ja | ja | ja |
| 3 | `Q_β` aus dem Anker, SP/SP′, `T = Q_β∘U` auf der Kugel | teilweise (Abschn. 7) | nein | nein |
| 4 | Mikrosegment-Induktion | **Lücke identifiziert** (Abschn. 8) | nein | nein |
| 5 | unabhängiger Primitivverifier | ja | ja (Koenigs, Tail) | ja, 576 Zeilen, PASS |

**Konsequenz für die Sprachregelung.** Die Freigaberegel des
Review-Auflösungsberichts bleibt in Kraft. `8r1` ist nach dieser Arbeit ein
Segmentzertifikat, dessen Koenigs-, Tail- und β-Restbudgets hergeleitet und
unabhängig nachgerechnet sind; es ist weiterhin **kein** abgeschlossener
computerassistierter Beweis, und zwar jetzt aus einem präziser benannten
Grund als vorher: neben den fehlenden Marschsegmenten fehlt die
Zell-Induktion aus Abschnitt 8.

## 12. Reproduktion

```bash
pip install --break-system-packages python-flint==0.9.0 mpmath

cd certificate
python3 test_tmodel.py          # Taylormodell-Invarianten, Negativkontrollen
python3 test_tmodel2.py         # zweivariable Modelle, 36 Bidisk-Stützstellen
python3 probe_field.py          # Feldprobe (Diagnose, nicht im Beweispfad)

python3 verify_primitives_8r1.py \
  --trajectory   .../trajectory_8r1.json \
  --stored-primitives .../segment_8r1_interval_primitives.json \
  --output       segment_8r1_derived_primitives.json \
  --report       segment_8r1_independent_verification.json

python3 compare_enclosures.py \
  --derived segment_8r1_derived_primitives.json \
  --stored  .../segment_8r1_interval_primitives.json \
  --report  segment_8r1_enclosure_containment.json
```

Laufzeit des Verifiers: 1031 s für 6 Panels × 96 Knoten, ein Kern.
Ausgelieferte Artefakte:

| Datei | SHA-256 der Payload |
|---|---|
| `segment_8r1_derived_primitives.json` | `e1dfcd91804800f2f8650f9929079da74673ce28c4909c25f0f575a405ceba37` |
