"""Unabhaengige Referenz fuer reelle Basen 1 < b < eta = e^(1/e).

Warum es das gibt
-----------------
Fuer das Shell-Thron-Innere hatte dieses Projekt KEINE Engine-Diversitaet:
b = 1.2 wurde nur Fork-gegen-Original geprueft (2026-07-15, Journal), und
beide teilen dieselbe Konstruktion, also denselben Fehler. Gemessen am
hier erzeugten Wert liefern BEIDE bei angeforderten 60 Stellen nur ~15
korrekte (Imaginaerteil 7.26e-14 auf einer beweisbar reellen Groesse), und
ab dps ~100 gar nichts Brauchbares mehr. Siehe Journal 2026-07-29.

Methode
-------
Regulaere Iteration am reellen ANZIEHENDEN Fixpunkt — kein Kneser, keine
Kontur, keine Theta-Funktion. Fuer 1 < b < eta existiert dieser Fixpunkt und
der Multiplikator ist reell in (0,1), damit ist die Loesung reell und eindeutig:

    f(z)        = b^z,   Fixpunkt L: b^L = L,   lam = L*ln b = ln L
    sigma(z)    = lim_n (f^n(z) - L) / lam^n                   (Koenigs)
    sigma^-1(u) = lim_n f^-n(L + lam^n u),  f^-1(z) = ln z / ln b
    S(x)        = sigma^-1( sigma(1) * lam^x )

    =>  S(0) = 1,  S(x+1) = b^S(x)

Praezisions-Falle (kostete mich einen Fehlversuch)
--------------------------------------------------
n und Arbeitspraezision sind NICHT unabhaengig waehlbar:
  * sigma^-1 braucht lam^n * u oberhalb der Aufloesung von L, sonst geht die
    Information beim Addieren verloren;
  * f^-1 expandiert mit 1/lam, jeder Restfehler wird mit lam^-n verstaerkt.
Mit dps=160 und n=400 (lam^400 ~ 1e-256) kam S(0) = 14.77 statt 1 heraus.
Daher: dps ~ n*log10(1/lam) + Ziel + Reserve, und die Selbsttests sind
PFLICHT — eine still falsche Referenz ist schlimmer als keine.
Und sie muessen etwas PRUEFEN: S(0)=1 taugt nicht, das reduziert sich
algebraisch auf f^-n(f^n(1)) und zeigt nur, dass log und exp invers sind.
Geprueft werden S(1)=b und S(2)=b^b.

Kosten (gemessen 2026-07-29, b=1.2)
-----------------------------------
    Ziel   60 Stellen:     11 ms     Ziel  300:    394 ms
    Ziel  120 Stellen:     37 ms     Ziel 1000:   15.9 s
    Ziel 2000 Stellen:  100.5 s
Exponent ~2.6-3.0 in der Stellenzahl (n ~ p mal Multiplikation M(p)).
Zum Vergleich: fatou.gp braucht fuer Basis e bei 2000 Stellen 26.69 h.
"""
from mpmath import mp, mpf, log, lambertw


def sexp_regular(b_str, x_str, target, _check=True):
    """sexp_b(x) fuer reelles 1 < b < e^(1/e) auf `target` Stellen.

    Setzt mp.dps selbst. Wirft ValueError, wenn die Selbsttests
    S(1)=b / S(2)=b^b das Ziel nicht halten.
    """
    b0 = mpf(b_str)
    lb0 = log(b0)
    if not 1 < b0 < mp.e ** (1 / mp.e):
        raise ValueError(f"b={b_str} liegt nicht in (1, eta)")
    L0 = -lambertw(-lb0) / lb0
    lam0 = L0 * lb0
    shrink = -mp.log10(lam0)                  # Stellen, die lam^n frisst
    n = int((target + 25) / shrink) + 1       # Koenigs-Konvergenz ~ lam^n
    mp.dps = int(n * shrink) + target + 60

    b = mpf(b_str)
    lb = log(b)
    L = -lambertw(-lb) / lb
    lam = L * lb

    def sigma(z):
        v = z
        for _ in range(n):
            v = b ** v
        return (v - L) / lam ** n

    def sigma_inv(u):
        v = L + lam ** n * u
        for _ in range(n):
            v = log(v) / lb
        return v

    s1 = sigma(mpf(1))
    val = sigma_inv(s1 * lam ** mpf(x_str))

    if _check:
        # S(0)=1 waere TAUTOLOGISCH: sigma^-1(sigma(1)) reduziert sich
        # algebraisch auf f^-n(f^n(1)) = 1 und prueft nur, dass log und exp
        # invers sind. Geprueft werden zwei echte Schritte der
        # Funktionalgleichung: S(1) = b und S(2) = b^b.
        e0 = abs(sigma_inv(s1 * lam ** 2) - b**b)
        e1 = abs(sigma_inv(s1 * lam ** 1) - b)
        held = int(-mp.log10(max(e0, e1, mpf(10) ** (-10 * target))))
        if held < target:
            raise ValueError(
                f"Selbsttest haelt nur {held} Stellen, Ziel {target} "
                f"— Referenz unbrauchbar")
    return val


def regime(b_str, target):
    """Alles Basis-Abhaengige EINMAL berechnen.

    sexp_regular() baut sigma(1) bei jedem Aufruf neu, was fuer ein einzelnes
    Ergebnis egal ist, aber fuer ein Gitter aus hunderten Punkten nicht: jeder
    Aufbau kostet regn ~ 2*target Exponentialfunktionen. Der Rueckgabewert
    haelt L, lam, regn, die Arbeitspraezision und sigma(1).
    """
    b0 = mpf(b_str)
    lb0 = log(b0)
    if not 1 < b0 < mp.e ** (1 / mp.e):
        raise ValueError(f"b={b_str} liegt nicht in (1, eta)")
    L0 = -lambertw(-lb0) / lb0
    shrink = -mp.log10(L0 * lb0)
    n = int((target + 25) / shrink) + 1
    wp = int(n * shrink) + target + 60

    saved = mp.dps
    mp.dps = wp
    b = mpf(b_str)
    lb = log(b)
    L = -lambertw(-lb) / lb
    lam = L * lb
    reg = {"b": b, "lb": lb, "L": L, "lam": lam, "n": n, "wp": wp,
           "target": target, "s1": None}
    reg["s1"] = _sigma(reg, mpf(1))
    mp.dps = saved
    return reg


def _sigma(reg, z):
    v = +z
    for _ in range(reg["n"]):
        v = reg["b"] ** v
    return (v - reg["L"]) / reg["lam"] ** reg["n"]


def _sigma_inv(reg, u):
    v = reg["L"] + reg["lam"] ** reg["n"] * u
    for _ in range(reg["n"]):
        v = log(v) / reg["lb"]
    return v


def sexp_from(reg, x):
    """T_b(x) = sigma^-1( sigma(1) * lam^x )."""
    saved = mp.dps
    mp.dps = reg["wp"]
    y = _sigma_inv(reg, reg["s1"] * reg["lam"] ** mpf(x))
    # ZUERST auf die Zielstellenzahl runden, DANN die Praezision des Aufrufers
    # wiederherstellen. Andersherum rundet das unaere Plus auf mp.dps des
    # Aufrufers -- bei der Vorgabe 15 kaemen 16 Stellen zurueck statt der
    # gerechneten, und zwar still.
    mp.dps = reg["target"]
    y = +y
    mp.dps = saved
    return y


def slog_from(reg, y):
    """A_b(y) = log(sigma(y)/sigma(1)) / log(lam).

    Definiert fuer y unterhalb des anziehenden Fixpunkts L; bei y = L sitzt
    die logarithmische Verzweigung (A_b -> +unendlich).
    """
    saved = mp.dps
    mp.dps = reg["wp"]
    z = log(_sigma(reg, mpf(y)) / reg["s1"]) / log(reg["lam"])
    # ZUERST auf die Zielstellenzahl runden, DANN die Praezision des Aufrufers
    # wiederherstellen. Andersherum rundet das unaere Plus auf mp.dps des
    # Aufrufers -- bei der Vorgabe 15 kaemen 16 Stellen zurueck statt der
    # gerechneten, und zwar still.
    mp.dps = reg["target"]
    z = +z
    mp.dps = saved
    return z


if __name__ == "__main__":
    from mpmath import nstr
    v = sexp_regular("1.2", "0.5", 120)
    print("sexp_1.2(0.5) =", nstr(v, 120))
