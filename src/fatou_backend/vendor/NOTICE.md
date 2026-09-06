# Third-party notice: fatou.gp

`fatou.gp` is by **Sheldon Levenstein** ("sheldonison"). It was published
on the Tetration Forum in the thread
["new fatou.gp program"](https://tetrationforum.org/showthread.php?tid=1017)
(first release July 10, 2015; later revisions in the same thread). It
computes Kneser's real-analytic tetration — sexp(z) and slog(z) to
arbitrary precision for a wide range of real and complex bases — and is,
to our knowledge, the reference implementation of that construction.

No explicit license accompanies the original publication. The file is
vendored here **unmodified**, with attribution, in good faith, for
reproducibility of the results in this repository. Copyright remains with
the author. If the author states license terms, this repository will
follow them; on request we will relicense our derivative or remove the
files.

`fatou_fork.gp` is a **derivative work** of `fatou.gp`: the same engine
with this project's 30 gate-verified performance optimizations, three
correctness fixes, and one construction the original does not carry —
regular iteration for real bases 1 < b < e^(1/e) (see
`research/METHODS.md` §2).

Two of those three defects are present in the unmodified original as
well, and are described in METHODS §2 in enough detail to be fixed there:
the saturation sentinel `safefs()` returns for large towers (base 2 above
x ≈ 4.9 gives `1E400/ln 2` instead of `2^65536`, with a `random()`
imaginary part), and the sub-eta regime, where the Kneser construction has
no complex fixed-point pair to work with and yields about 15 correct
digits regardless of the requested precision. Both are offered upstream
unconditionally, as are the optimizations. The third defect was introduced
by this fork and never affected the original.

The underlying code remains © Sheldon Levenstein; the modifications are
© 2026 Janis Justus and are offered under the repository's MIT license to
the extent they are separable.

Everything else in this repository builds on this engine. Thank you,
Sheldon.
