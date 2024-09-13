from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from numpy import concatenate as concat
from numpy import diff, sqrt
import xarray as xr


dirpath_data = r'D:\WORK\Salvador\repo\ubiquitous-main\lfp_csd\data'

#fpath_in = Path(dirpath_data) / f'lfp_csd_avgtrial_5.nc'
fpath_in = Path(dirpath_data) / 'lfp_csd_avgtrial_(n=100000_s=3).nc'
Xavg = xr.open_dataset(fpath_in)

fband = (70, 110)
#fband = (90, 95)
#fband = (110, 145)
#fband = (20, 22)
#fband = (10, 30)

fband= fband[::-1]

def norm(x):
    return x / x.max()

Xavg['csd2_interp_complex'] *= -1

W, Z = {}, {}
for data_type in ['lfp_interp', 'csd2_interp']:
    W_ = Xavg[f'{data_type}_pow']
    Z_ = Xavg[f'{data_type}_complex']
    #W_ = norm(W_).sel(freq=slice(*fband)).mean(dim='freq')
    W_ = W_.sel(freq=slice(*fband)).mean(dim='freq')
    Z_ = Z_.sel(freq=slice(*fband)).mean(dim='freq')
    data_type_new = data_type[:3]
    W[data_type_new] = W_
    Z[data_type_new] = xr.apply_ufunc(np.angle, Z_)

def diff_kd(x, n=1):
    """Derivative with preserved array size. """
    npre = int(np.ceil(n / 2))
    npost = int(np.floor(n / 2))
    #x = concat((x[: npre], x, x[len(x) - npost :]))
    dx = diff(x, n)
    dx = concat((dx[: npre], dx, dx[len(dx) - npost :]))
    return dx

# LFP amplitude and phase profiles
r2 = W['lfp']
r = np.sqrt(r2)
phi = Z['lfp'] * 0

# Derivatives of LFP amplitude and phase profiles
dr = diff_kd(r)
d2r = diff_kd(r, 2)
dphi = diff_kd(phi)
d2phi = diff_kd(phi, 2)

# Analytic calculation of CSD amplitude profile
r2_csd = (d2phi**2 * r**2 + 4 * d2phi * dphi * dr * r +
          d2r**2 - 2 * d2r * dphi**2 * r +
          dphi**4 * r**2 + 4 * dphi**2 * dr**2)
r_csd = sqrt(r2_csd)
What = {'csd': r2_csd}

# Plot LFP and CSD profilles
plt.figure(figsize=(8, 8))
nx, ny = 2, 2
chans = Xavg.chan.values
fband_str = f'{fband[1]}-{fband[0]} Hz'
for n, data_type in enumerate(['lfp', 'csd']):
    plt.subplot(ny, nx, n * nx + 1)
    plt.plot(norm(W[data_type]), chans)
    #plt.plot(W[data_type], chans)
    #if data_type in What:
        #plt.plot(norm(What[data_type]), chans, 'r--')
    plt.xlim(0, 1)
    plt.gca().invert_yaxis()
    plt.ylabel('Channel')
    plt.title(f'{data_type} ({fband_str}), power')
    plt.subplot(ny, nx, n * nx + 2)
    plt.plot(Z[data_type], chans)
    plt.gca().invert_yaxis()
    #plt.ylabel('Channel')
    plt.title(f'{data_type} ({fband_str}), phase')
#plt.subplot(2, 2, 3)
#plt.xlabel()
