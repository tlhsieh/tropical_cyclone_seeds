[![DOI](https://zenodo.org/badge/436679362.svg)](https://zenodo.org/badge/latestdoi/436679362)

### tracker.py

Tracks convective clusters and tropical cyclone seeds given six hourly climate model output. Detailed documentation in [Hsieh, T. L., Vecchi, G. A., Yang, W., Held, I. M., & Garner, S. T. (2020). Large-scale control on the frequency of tropical cyclones and seeds: a consistent relationship across a hierarchy of global atmospheric models. Climate Dynamics, 55(11), 3177-3196.](https://link.springer.com/article/10.1007/s00382-020-05446-5)

Example: 
python tracker.py path_to_6hrly_output 111 150

### spi.py

Computes the seed propensity index given monthly mean omega_500 and vort_850 defined in [Hsieh, T. L., Yang, W., Vecchi, G. A., & Zhao, M. (2022). Model spread in the tropical cyclone frequency and seed propensity index across global warming and ENSO‐like perturbations. Geophysical Research Letters, 49(7), e2021GL097157.](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2021GL097157)
