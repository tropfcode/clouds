import os
import time
import math
import time
import locale
import netCDF4
import analysis
import decorators
import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors

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

def gridboxes(datapath):
    """
    Finda Frequency of Occurance for each Weather State and total.
    Currently only for one year at a time.
    
    RETURN
    ------
    wsgrids: np.array, shape=(360, 180, 12), type=int
        Gridboxes for each weather state
        
    grid: np.array, shape=(360, 180), type=int
        Total gridbox
    """
    wsgrids = np.zeros((360, 180, 12))
    grid = np.zeros((360, 180))
    masks = np.zeros((6, 7, 41252), dtype=bool)
    histograms = load_data()
    for f1 in sorted(os.listdir(datapath)):
        t1 = time.time()
        print('this is f1', f1)
        for f2 in sorted(os.listdir(datapath+f1+'/')):
            fullpath = datapath+f1+'/'+f2
            data = netCDF4.Dataset(fullpath)
            pctaudist = data['n_pctaudist']
            sqlonbeg = data['sqlon_beg'][:]
            sqlonend = data['sqlon_end'][:]
            eqlat = data['eqlat_index'][:]
            ntotal = data['n_total'][:]
            masks[:,:,:] = np.ma.getmaskarray(pctaudist[:,:,:])
            pctaudistnorm = analysis.fnorm(pctaudist[:,:,:], ntotal[:], masks[:,:,:])
            wsgrids, grid = analysis.findws(pctaudistnorm[:,:,:], histograms[:,:,:]/100.0, masks[:,:,:], 
                                                wsgrids[:,:,:], grid[:,:], sqlonbeg[:], sqlonend[:], eqlat[:])
        t2 = time.time()
        print('run time for entire month ', t2-t1)
        return wsgrids, grid