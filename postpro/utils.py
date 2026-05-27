import numpy as np
import xarray as xr


def fftfit(x):
    while not (x & 1):
        x = x >> 1
    return x == 1 or x == 3


def mesh2k(phi, mode_size):
    """1D wavenumbers: 0..mode_size, then negative modes."""
    k = np.zeros(2*mode_size + 1)
    k[:mode_size + 1] = np.arange(0, mode_size + 1)*phi
    k[mode_size + 1:] = -np.arange(mode_size, 0, -1)*phi
    return k

def foldAlongDim(data: xr.DataArray, oldVar="y", inverse = False, return_half=False):
    newVar = oldVar + '_folded'
    data =  data.rename({oldVar: 'dic'})
    if oldVar == "y":
        inv = data.isel(dic=slice(None, None,-1)).assign_coords(coords=data.coords)
    elif oldVar == "kz":
        inv = data.assign_coords(dic = -data.dic.values)

    returnData = (data + (-1 if inverse else 1) * inv)/2
    if return_half:
        return returnData.isel(dic = slice(0,(data["dic"].size+1) // 2)).rename({'dic': newVar})
    else:
        return returnData.rename({'dic': oldVar})