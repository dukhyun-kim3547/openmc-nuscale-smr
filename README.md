# OpenMC NuScale VOYGR SMR Study

Monte Carlo neutron transport simulations of NuScale VOYGR SMR using OpenMC 0.15.3.  
Conducted as part of a computational nuclear engineering portfolio for PhD applications.

## Environment

- OpenMC 0.15.3 (conda-forge)
- Nuclear data: ENDF/B-VIII.0 HDF5
- Python 3.11 / WSL2 Ubuntu

## Simulations

### 1. NuScale vs PWR Comparison (`nuscale_vs_pwr.py`)
Pin cell k-infinity comparison between Large PWR and NuScale VOYGR.

![NuScale vs PWR](nuscale_vs_pwr.png)

| Design | Enrichment | Coolant Density | k-infinity |
|--------|------------|-----------------|------------|
| Large PWR (AP1000) | 3.10% | 0.71 g/cc | 1.366 ± 0.001 |
| NuScale VOYGR | 4.95% | 0.74 g/cc | 1.448 ± 0.001 |

**Finding:** Higher enrichment in NuScale provides +0.082 excess reactivity,  
enabling sufficient fuel life in a compact core.

## References

- NuScale FSAR (Publicly available NRC Design Certification Document)
- OpenMC Documentation: https://docs.openmc.org
- ENDF/B-VIII.0 Nuclear Data Library
