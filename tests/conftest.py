"""Shared pytest config.

Adds the project root to sys.path so test modules can ``import router``,
``import qiskit_pass`` etc. without an editable install.
"""
import os
import sys

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)
