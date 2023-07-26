#!/usr/bin/python
## Tsung-Lin Hsieh (hsiehtl@princeton.edu)

import numpy as np
import xarray as xr

def get_fcor(lat):
    return 2*2*np.pi/86164*np.sin(lat/180*np.pi)

def get_beta(lat):
    return 2*2*np.pi/86164/6371e3*np.cos(lat/180*np.pi)
    
def Zparam(f, b, U=20): 
    """the non-dimensional Z parameter defined in Hsieh et al. (2020): Z = L_beta/L_f
    
    Args:
        f: the cyclonic absolute vorticity
        b: meridional gradient of the absolute vorticity
        U: the characteristic velocity scale for pre-seed disturbances
    """
    
    return np.maximum(f, 0)/abs(b)**0.5/U**0.5

def Z2prob(Z, sigma=0.69):
    """convert Z to the probability of transition from a convective cluster to a seed"""
    
    return 1/(1 + Z**(-1/sigma))

def spi(omega500, vort850):
    """the seed propensity index defined in Hsieh et al. (2022)"""
    
    fcor = get_fcor(vort850.lat)
    beta = get_beta(vort850.lat)

    absvort = (fcor + vort850)*np.sign(fcor) # such that positive is cyclonic
    effbeta = beta + vort850.differentiate('lat')*(1/111e3) # 1 deg lat \approx 111 km
    
    seed_propensity_index = np.maximum(-omega500, 0)*Z2prob(Zparam(absvort, abs(effbeta)))
    
    return seed_propensity_index