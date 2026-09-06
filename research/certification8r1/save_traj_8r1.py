"""Compatibility entry point for rebuilding the missing 8r1 trajectory.

The historical research log names ``save_traj_8r1.py``.  The complete
implementation lives in :mod:`produce_8r1`; keeping this thin entry point
restores the recorded command name without duplicating any numerical code.
"""

from produce_8r1 import main


if __name__ == "__main__":
    raise SystemExit(main())
