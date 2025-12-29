# sLSP4SNFs -- sliding Lomb–Scargle periodogram (sLSP) for super-Nyquist frequencies (SNFs)

**sLSP4SNF** is a Python package for identifying **super-Nyquist frequencies (SNFs)** in stellar photometric time series using the **sliding Lomb–Scargle periodogram (sLSP)** technique.

The package is designed for asteroseismic studies based on space photometry from **Kepler**, **TESS**, and future missions such as **PLATO** **ET2.0**.

---

## Features

- Sliding Lomb–Scargle periodogram (sLSP)
- Identification of super-Nyquist frequencies
- Publication-ready visualizations
- Modular, pip-installable research code
- Optimized for large stellar samples

---

## Installation

install from source:

```bash
git clone https://github.com/wangxuan20221015-star/sLSP4SNFs.git
cd sLSP4SNFs
pip install .
```

---

## Quick Start

```python
from sLSP4SNFs import *

time, flux = load_lightcurve()

# Initialize the SLSP
slsp_ob = SLSP(time,flux)

# sLSP calculation
slsp_ob.compute_slsp()

# Determine whether this frequency is an SNF.
slsp_ob.analyze_snf(189.33)

# plot
ploter = slsp_ob.plot()
ploter.SNF()
```

---

## Scientific Background

Frequencies above the Nyquist limit can appear as aliases in discretely sampled photometric data.  
By tracking frequency variations with a sliding Lomb–Scargle periodogram, genuine **super-Nyquist frequencies** can be distinguished from real pulsation signals.

---

## Requirements

- Python ≥ 3.8
- numpy
- scipy
- matplotlib
- pandas

---

## License

MIT License

---

## Citation

If you use **sLSP4SNF** in your research, please cite:

> Wang, Xuan, Zong, Weikai, et al. 2025.
> Mo Yanqi, Zong, Weikai, et al. 2026.
