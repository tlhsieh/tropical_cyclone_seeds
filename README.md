### tracker.py

Tracker for convective clusters and tropical cyclone seeds, using six hourly climate model output. Detailed documentation in [Hsieh, T. L., Vecchi, G. A., Yang, W., Held, I. M., & Garner, S. T. (2020). Large-scale control on the frequency of tropical cyclones and seeds: a consistent relationship across a hierarchy of global atmospheric models. Climate Dynamics, 55(11), 3177-3196.](https://link.springer.com/article/10.1007/s00382-020-05446-5)

Example: 
python tracker.py path_to_6hrly_output 111 150

### spi.py

Computes the seed propensity index (Hsieh, Yang, Vecchi, and Zhao, 2022) given the monthly mean omega_500 and vort_850. Requires helper functions in [gfd.py](https://github.com/tlhsieh/geophysical_fluid_analysis/blob/master/gfd.py)