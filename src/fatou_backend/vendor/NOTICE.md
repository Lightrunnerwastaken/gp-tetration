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
with this project's 30 gate-verified performance optimizations (see
`research/METHODS.md` §2). The underlying code remains © Sheldon
Levenstein; the modifications are © 2026 Janis Justus and are offered
under the repository's MIT license to the extent they are separable —
and offered upstream to the original author unconditionally.

Everything else in this repository builds on this engine. Thank you,
Sheldon.
