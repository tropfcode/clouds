import numpy as np
import math

def weighted_rfo(maps):
    """
    Computes Weighted RFO for a Equal-Area gridbox map
    
    Parameters
    ----------
    maps: 1D np.ndarray, dtype=float
        Gridbox map of shape (64800), (14400)
        
    Returns
    -------
    rfo: float
        Weighted RFO of maps
        
    TODO
    ----
    Implement RFO calcultion for D1 data (degree of 2.5)
    """
    if maps.shape[0] == 64800:
        latbox = np.arange(.5, 90, 1.0)
        deg = 90
        longitude = 360
        latitude = 180
    elif maps.shape[0] == 16200:
        latbox = np.arange(1, 90, 2)
        deg = 45
        longitude = 180
        latitude = 90
    else:
        print('Incorrect value given for degree: must be integer 1 or 2')
        return False
    
    total = 0
    weightsum = 0
    tmp_map = maps.reshape(longitude, latitude)
    for j, lat in enumerate(latbox[::-1]):
        if np.all(tmp_map[:,j].mask):
            continue
        avg = np.mean(tmp_map[:,j])
        total += avg*math.cos(math.radians(lat))
        weightsum += math.cos(math.radians(lat))
    for j, lat in enumerate(latbox):
        if np.all(tmp_map[:,(j+deg)].mask):
            continue
        avg = np.mean(tmp_map[:,(j+deg)])
        total += avg*math.cos(math.radians(lat))
        weightsum += math.cos(math.radians(lat))
    rfo = total/(weightsum)
    return rfo


# Get rid of this after adding gridsize of 2.5 to above code
"""
def compute_rfo(maps, gridsize=1):
        latbox = np.arange(1.25, 90, 2.5)
        for i in range(12):
            total = 0
            weightsum = 0
            for j, lat in enumerate(latbox[::-1]):
                if np.all(maps[j,:,i].mask):
                    continue
                avg = np.mean(maps[j,:,i])
                total += avg*math.cos(math.radians(lat))
                weightsum += math.cos(math.radians(lat))
            for j, lat in enumerate(latbox):
                if np.all(maps[(j+36),:,i].mask):
                    continue
                avg = np.mean(maps[(j+36),:,i])
                total += avg*math.cos(math.radians(lat))
                weightsum += math.cos(math.radians(lat))
            totallist.append(total/weightsum)
    #print(totallist)
    #print('')
    return totallist
"""