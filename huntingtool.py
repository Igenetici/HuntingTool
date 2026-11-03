#!/usr/bin/env python3
"""Run HuntingTool without installing it: `python huntingtool.py` opens the graphical
interface; `python huntingtool.py run -i INPUT -o OUTDIR` runs the command line."""
import multiprocessing
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

if __name__ == "__main__":
    multiprocessing.freeze_support()
    from huntingtool.cli import main
    sys.exit(main())
