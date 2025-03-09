import pandas as pd
import numpy as np
from scipy.stats import truncnorm
from scipy.integrate import quad
import matplotlib.pyplot as plt

dat = pd.read_csv('../Data/pseudo_unbinned_cat3.csv')
obs = np.array(dat['x'])

l = 110; u = 160
mean_sig = 124.82; sd_sig = 2.86
aL = 1.53; aR = 1.59;
nL = 4.62; nR = 10.5
rate_gb = 1.2
eps = 5*1e-2

def fs_unnormed(x):
    l_std = (l-mean_sig)/sd_sig
    u_std = (u-mean_sig)/sd_sig
    t = (x - mean_sig)/sd_sig
    if t<aR and t>-aL:
        val = np.exp(-(t**2)/2)
    elif t >= aR and t<=u_std:
        val = np.exp(-(aR**2)/2) * ((aR/nR)*(np.abs(nR/aR - aR + t)))**(-nR)
    elif t<=-aL and t >= l_std:
        val = 
    return truncnorm.pdf(x, a = (l-mean_sig)/sd_sig, b = (u-mean_sig)/sd_sig,
                         loc = mean_sig, scale = sd_sig)

def qb: 