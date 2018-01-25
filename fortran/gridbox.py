import os
import time
import math
import time
import locale
import netCDF4
import analysis
import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import decorator

def load_data():
    """
    Load the data in pctau.txt.
    This data is returned as taupc to match the shape of variable 'pctau_dist'
    in netCDF files (pctau_dist data is arranged more like taupc)
    
    RETURN
    ------
    taupc_histograms: float32 numpy array of shape (12, 6, 7)
    """
    pctaulist = []
    tmp = []
    count = 0
    with open('./pctau.txt') as pctau:
        for line in pctau:
            for value in line.split():
                tmp.append(float(value))
                count += 1
            if count == 42:
                count = 0 
                pctaulist.append(np.asarray(tmp).reshape(7, 6))
                tmp = []

    pctau_histograms = np.asarray(pctaulist)
    
    # Transform pctau_histograms to taupc to match pctaudist (matrix transpose)
    taupc_histograms = np.zeros((12, 6, 7))#np.copy(pctau_histograms).reshape((12, 6, 7))
    for i in xrange(pctau_histograms.shape[0]):
        taupc_histograms[i,:,:] = pctau_histograms[i,:,:].T

    taupc_histograms = taupc_histograms.astype(np.float32)
    return taupc_histograms