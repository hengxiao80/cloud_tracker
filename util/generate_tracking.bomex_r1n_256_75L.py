import os, sys, glob
import numpy as np
import ujson as json

import xarray as xr
from netCDF4 import Dataset as nc

cp = 1004.    # Heat capacity at constant pressure for dry air [J kg^-1 K^-1]
Rv = 461.      # Gas constant of water vapor [J kg^-1 K^-1]
Rd = 287.      # Gas constant of dry air [J kg^-1 K^-1]
p_0 = 100000.
lam = Rv/Rd - 1.

def T_v(T, qv, qn, qp):
    return T*(1. + lam*qv - qn - qp)

def theta(p, T): return T*(p_0/p)**(Rd/cp)

def theta_v(p, T, qv, qn, qp):
    return theta(p, T_v(T, qv, qn, qp))

def generate_tracking(time, file):
    with xr.open_dataset(file) as f:
        # try:
        #     f = f.squeeze('time')
        # except:
        #     pass
        print(f)

        print("\t Calculating velocity fields...")
        w_field = f['W'][:]
        w_field = (w_field + np.roll(w_field, -1, axis=1)) / 2

        u_field = f['U'][:]
        u_field = (u_field + np.roll(u_field, -1, axis=2)) / 2

        v_field = f['V'][:]
        v_field = (v_field + np.roll(v_field, -1, axis=3)) / 2

        print("\t Calculating buoynacy fields...")
        qn_field = f['QN'][:] / 1e3
        thetav_field = theta_v(f['p']*100, f['TABS'][:], 
                               f['QV'][:] / 1e3, qn_field, 0)
        buoy_field = (thetav_field > 
                      np.mean(thetav_field, axis=(2,3)))

        print("\t Calculating tracer fields...")
        tr_field = np.array(f['TR01'][:])
        tr_mean = np.mean(tr_field, axis=(2,3))
        tr_stdev = np.std(tr_field, axis=(2,3))
        tr_min = .05 * np.cumsum(tr_stdev, axis=1)/(np.arange(len(tr_stdev[0,:]))+1)
        tr_limit = np.maximum(tr_mean+tr_stdev, tr_min)

        #---- Dataset for storage 
        print("\t Saving DataArrays...")
        ds = xr.Dataset(coords= {'time': f.time, 'z': f.z, 'y':f.y, 'x':f.x})
        
        ds['u'] = u_field
        ds['v'] = v_field
        ds['w'] = w_field
    
        # ds['core'] = (w_field > 0.) & (buoy_field > 0.) & (qn_field > 0)
        ds['core'] = (w_field > 0.) & buoy_field & (qn_field > 0)
        ds['condensed'] = (qn_field > 0)
        ds['plume'] = xr.DataArray(tr_field > tr_limit[:, :, None, None], 
            dims=['time', 'z', 'y', 'x'])

        # Flag for netCDF compression (very effective for large, sparse arrays)
        encoding = {}
        if False:
            encoding = {var: dict(zlib=True) for var in ds.data_vars}
        # ds.to_netcdf('data/cloudtracker_input_%08g.nc' % time, encoding=encoding)
        ds.to_netcdf('data/cloudtracker_input_%08g.nc' % time)

def main():
    global model_config
    with open('model_config.json', 'r') as json_file:
        model_config = json.load(json_file)
        nt = model_config['config']['nt']

    filelist = sorted(glob.glob('%s/*.nc' % model_config['variables']))
    for time in range(nt):
        print('\t Working...%s/%s                        ' % (time, nt), end='\r')
        generate_tracking(time, filelist[time])

if __name__ == "__main__":
    main()
