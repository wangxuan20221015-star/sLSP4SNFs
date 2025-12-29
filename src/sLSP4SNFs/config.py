"""
Global matplotlib configuration for SNF analysis.

This module centralizes:
- matplotlib rcParams
- default colormap
- shared plotting constants

Recommended usage:
    from snf.config import setMplStyle, setDefaultColormap, LABEL_SIZE
"""

import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap


# Global plotting constants
LABEL_SIZE = 24


def set_mpl_style():
    """
    Apply publication-quality matplotlib style.
    Should be called ONCE at the beginning of a script.
    """
    
    mpl.rcParams.update({
        "font.sans-serif": ["Arial"],
        "font.family": "sans-serif",

        "xtick.direction": "in",
        "ytick.direction": "in",

        "xtick.minor.visible": True,
        "ytick.minor.visible": True,

        "xtick.labelsize": 16,
        "ytick.labelsize": 16,

        "xtick.major.width": 1.5,
        "ytick.major.width": 1.5,

        "xtick.major.size": 8,
        "ytick.major.size": 8,

        "xtick.minor.size": 4,
        "ytick.minor.size": 4,

        "xtick.top": True,
        "ytick.right": True,
    })
    

def set_default_colormap(colors=None):
    """
    Create and return the default SNF colormap.

    Parameters
    ----------
    colors : list of str, optional
        List of color names defining the gradient.

    Returns
    -------
    matplotlib.colors.LinearSegmentedColormap
    """

    if colors is None:
        colors = ["blue", "royalblue", "red", "yellow"]
        
    custom_cmap = LinearSegmentedColormap.from_list(
        "custom_gradient", colors
    )
    
    return custom_cmap