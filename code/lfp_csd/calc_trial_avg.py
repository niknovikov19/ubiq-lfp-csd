from pathlib import Path

import dask.array as da
import h5py
import numpy as np
import scipy.ndimage
import xarray as xr


fpath_in = r'D:\WORK\Salvador\repo\ubiquitous-main\A1_SpontWaveComplex.mat'

dirpath_data = r'D:\WORK\Salvador\repo\ubiquitous-main\lfp_csd\data'

dt = 0.001

nsamp_used = 100000
#nsamp_used = 5000

# Total: (sample: 5000, freq: 73, chan: 21)
#chunks = (5000, -1, -1)
chunks = (-1, 13, -1)

ref_chan = 0

# Vertical smoothing
s = 3

need_compute = False


# Open mat-file
fid = h5py.File(fpath_in, 'r')

X, Xavg = {}, {}

for data_type in ['lfp', 'csd']:    
    print(f'Load: {data_type}')
    
    # Get data from mat-file
    X_ = da.from_array(fid[f'data/{data_type}'], chunks=chunks)
    
    # Select time interval and convert to complex format
    X_ = X_[:nsamp_used]
    X_ = X_['real'] + 1j * X_['imag']
    if need_compute:
        X_ = X_.compute()
    
    nsamples = X_.shape[0]
    nfreqs = X_.shape[1]
    nchans = X_.shape[2]

    ff = np.array(fid['data/freqs']).ravel()    
    tt = np.arange(nsamples) * dt
    
    # Convert to xarray
    dims = ['sample', 'freq', 'chan']
    coords = {
        'sample': np.arange(nsamples),
        'time': ('sample', tt),
        'freq': ff,
        'chan': np.arange(nchans)
        }
    X_ = xr.DataArray(X_, dims=dims, coords=coords, name=data_type)
    
    # Re-reference
    if ref_chan is not None:
        X_ -= X_.sel(chan=ref_chan)
        
    X[data_type] = X_

def diff_keepdims(X, n, axis):
    X = np.concatenate(
        [np.take(X, [0], axis), X, np.take(X, [-1], axis)],
        axis=axis
        )
    return np.diff(X, n ,axis)

# Smooth over depth
X_ = X['lfp'].copy()
#Y_ = X_.isel(chan=slice(None, None, 2))
#Y_ = Y_.interp(chan=X_['chan'], method='cubic')
if s > 0:
    Y_ = xr.apply_ufunc(
        scipy.ndimage.gaussian_filter1d,
        X_,
        kwargs={'sigma': s, 'axis': X_.get_axis_num('chan')},
        dask='parallelized',
        output_dtypes=[X_.dtype]
    )
else:
    Y_ = X_
X['lfp_interp'] = Y_

# Re-calculate CSD
for data_type in ['lfp', 'lfp_interp']:
    X_ = X[data_type]
    Y_ = xr.apply_ufunc(
        diff_keepdims, X_,
        input_core_dims=[['chan']], output_core_dims=[['chan']],
        kwargs={'n': 2, 'axis': X_.get_axis_num('chan')},
        dask='parallelized', output_dtypes=[np.complex128]
    )
    Y_ = Y_.assign_coords(chan=X_.chan)
    #Y_ = Y_.compute()
    data_type_new = data_type.replace('lfp' ,'csd2')
    Y_.name = data_type_new
    X[data_type_new] = Y_

# Average power and complex amplitudes over time
#data_types = ['lfp', 'csd', 'csd2', 'lfp_interp', 'csd2_interp']
data_types = ['lfp_interp', 'csd2_interp']
#data_types = ['lfp', 'csd']
for data_type in data_types:
    print(f'Average over time: {data_type}')
    
    X_ = X[data_type]
    
    W_ = (X_ * X_.conj()).real.mean(dim='sample')
    W_ = W_.compute()
    
    #Z_ = (X_ / X_.sel(chan=0)).mean(dim='sample')
    Z_ = (X_ * X_.sel(chan=0).conj()).mean(dim='sample')
    Z_ = Z_.compute()
    
    Xavg[data_type + '_pow'] = W_
    Xavg[data_type + '_complex'] = Z_

X = xr.Dataset(X)
Xavg = xr.Dataset(Xavg)

fpath_out = (Path(dirpath_data) / 
             f'lfp_csd_avgtrial_(ref={ref_chan}_n={nsamp_used}_s={s}).nc')
#fpath_out = Path(dirpath_data) / f'lfp_csd_avgtrial_5.nc'
Xavg.to_netcdf(fpath_out, engine='h5netcdf', invalid_netcdf=True)

fid.close()


