"""
SNF: Sliding Lomb–Scargle Periodogram analysis package

This package provides tools for:
- Sliding Lomb–Scargle periodogram (sLSP) analysis
- Super-Nyquist Frequency (SNF) identification
- Publication-quality visualization

Main public interfaces:
- SLSP : core analysis class
- Plot : visualization helper
"""
from importlib.resources import files
from .version import __version__

# ------------------------------------------------------------------
# Core analysis
# ------------------------------------------------------------------
from .main import SLSP

# ------------------------------------------------------------------
# Plotting interface
# ------------------------------------------------------------------
from .plotting import Plot

# ------------------------------------------------------------------
# Configuration helpers
# ------------------------------------------------------------------
from .config import (
    LABEL_SIZE,
    set_mpl_style,
    set_default_colormap,
)

# ------------------------------------------------------------------
# Public API
# ------------------------------------------------------------------
__all__ = [
    "SLSP",
    "Plot",
    "LABEL_SIZE",
    "set_mpl_style",
    "set_default_colormap",
    "__version__"
]