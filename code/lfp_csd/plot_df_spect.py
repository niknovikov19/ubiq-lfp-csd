from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr


# =============================================================================
# def plot_xarray_2d(X):
#     coord_x = X.dims[1]
#     coord_y = X.dims[0]
#     cx = X.coords[coord_x].values
#     cy = X.coords[coord_y].values
#     ext = [cx[0], cx[-1], cy[-1], cy[0]]
#     plt.imshow(X.values, aspect='auto', extent=ext, origin='upper')
#     plt.xlabel(coord_x)
#     plt.ylabel(coord_y)
# =============================================================================

def plot_xarray_2d(X, need_log=False, cmap=None):
    coord_x = X.dims[1]
    coord_y = X.dims[0]
    cx = X.coords[coord_x].values
    cy = X.coords[coord_y].values
    Cx, Cy = np.meshgrid(cx, cy)
    ax = plt.gca()
    X_ = X.values
    if need_log:
        X_ = np.log(X_)
    if cmap is not None:
        cmap = plt.get_cmap(cmap)
    m = ax.pcolormesh(Cx, Cy, X_, shading='auto', cmap=cmap)
    #ax.set_xscale('log')
    ax.invert_yaxis()
    plt.xlabel(coord_x)
    plt.ylabel(coord_y)
    plt.colorbar(m)
    plt.title(X.name)


dirpath_data = r'D:\WORK\Salvador\repo\ubiquitous-main\lfp_csd\data'

#fpath_in = Path(dirpath_data) / f'lfp_csd_avgtrial_3.nc'
fpath_in = Path(dirpath_data) / 'lfp_csd_avgtrial_(n=100000_s=2).nc'
Xavg = xr.open_dataset(fpath_in)

#data_types_vis = ['lfp', 'csd', 'csd2']
#data_types_vis = ['lfp', 'csd2']
#data_types_vis = ['lfp', 'csd']
data_types_vis = ['lfp_interp', 'csd2_interp']

nx = 2
ny = len(data_types_vis)

plt.figure(figsize=(12, 10))
for n, data_type in enumerate(data_types_vis):
    plt.subplot(ny, nx, n * nx + 1)
    W = Xavg[f'{data_type}_pow']
    W /= W.max(dim='chan')
    plot_xarray_2d(W.T, need_log=True)
    plt.subplot(ny, nx, n * nx + 2)
    Z = xr.apply_ufunc(np.angle, Xavg[f'{data_type}_complex'])
    plot_xarray_2d(Z.T, cmap='hsv')

Xavg.close()
