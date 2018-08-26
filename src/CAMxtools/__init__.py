from __future__ import print_function
__doc__ = r"""
REpyprimer provides a collection of pre- and post-
processing tools for CAMx, CMAQ, and GEOS-Chem
"""

__all__ = []
import sys
import os
import warnings
warn = warnings.warn
def clean_showwarning(message, category, filename, lineno, file=None, line=None):
    print('**REpp:%s:%s:%s:\n  %s' % ((filename), lineno, category.__name__, message), file = sys.stderr)
    return
warnings.showwarning = clean_showwarning
