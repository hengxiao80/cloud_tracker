#!/usr/bin/env python
"""
This program generates a pkl file containing a list of dictionaries.
Each dictionary in the list represents a condensedlet.
The dictionaries have the structure:
{'core': array of ints of core points,
'condensed': array of ints of condensed points,
'plume': array of ints of plume points,
'u_condensed': ,
'v_condensed': ,
'w_condensed': ,
'u_plume': ,
'v_plume': ,
'w_plume': }
pkl files are saved in pkl/ subdirectory indexed by time
"""

from netCDF4 import Dataset as nc

import xarray, dask, h5py
import numpy as np

from .utility_functions import index_to_zyx, expand_indexes
from .load_config import c

def expand_cloudlet(cloudlet, indexes):
    """Given an array of indexes composing a cloudlet and a boolean mask
    array indicating if each model index may be expanded into (True) or
    not (False), expand the cloudlet into the permissable indicies that
    are find all indicies adjacent to the cloudlet.

    Returns an array of the indicies composing the expanded cloudlet, and
    an array of the remaining indicies that may be expanded into.
    """

    # Expand the cloudlet indexes into their nearest neighbours
    expanded_cloudlet = expand_indexes(cloudlet)

    # Find the mask values of the expanded indexes
    mask = indexes[expanded_cloudlet]

    # Select the expanded cloudlet indexes that may be expanded into
    new_points = expanded_cloudlet[mask]

    # Remove the indicies that have been added to the cloudlet
    indexes[new_points] = False

    return new_points, indexes

#---------------------

def expand_current_cloudlets(key, cloudlets, mask):
    cloudlet_points = []
    for cloudlet in cloudlets:
        cloudlet_points.append( [cloudlet[key]] )

    cloudlet_expand_indexes = range(len(cloudlet_points))

    while cloudlet_expand_indexes:
        next_loop_cloudlet_list = []

        # Go through the current list of cloudlets
        for n in cloudlet_expand_indexes:
            expanded_points, mask = expand_cloudlet(cloudlet_points[n][-1], mask)

            if len(expanded_points) > 0:
                cloudlet_points[n].append(expanded_points)
                next_loop_cloudlet_list.append(n)

        cloudlet_expand_indexes = next_loop_cloudlet_list

    for n, cloudlet in enumerate(cloudlet_points):
        cloudlets[n][key] = np.hstack(cloudlet)

    return cloudlets, mask


def make_new_cloudlets(key, mask):
    indexes = np.arange(c.nz * c.ny * c.nx)[mask]
    cloudlets = []

    for n in indexes:
        if mask[n]:
            mask[n] = False
            cloudlet_indexes = [np.array((n,))]

            # add_new_cloudlet
            done = False
            while not done:
                new_indexes, mask = expand_cloudlet(cloudlet_indexes[-1], mask)

                if len(new_indexes) > 0:
                    cloudlet_indexes.append( new_indexes )
                else:
                    # If the number of points in the cloudlet has not changed, we are done
                    done = True
        
            cloudlet = {}
            cloudlet[key] = np.hstack(cloudlet_indexes)
            cloudlets.append( cloudlet )
            
    return cloudlets

#-----------------

def find_mean_cloudlet_velocity(cloudlets, u, v, w):
    dx, dy, dz = c.dx, c.dy, c.dz

    for cloudlet in cloudlets:
        if len(cloudlet['condensed']) > 0:
            K, J, I = index_to_zyx(cloudlet['condensed'])
            # find the mean motion of the cloudlet
            u_mean = np.mean(u[K, J, I]) - c.ug
            v_mean = np.mean(v[K, J, I]) - c.vg
            w_mean = np.mean(w[K, J, I])
        
            cloudlet['u_condensed'] = (u_mean * c.dt / dx)
            cloudlet['v_condensed'] = (v_mean * c.dt / dy)
            cloudlet['w_condensed'] = (w_mean * c.dt / dz)
        else:
            cloudlet['u_condensed'] = 0.
            cloudlet['v_condensed'] = 0.
            cloudlet['w_condensed'] = 0.

        K, J, I = index_to_zyx(cloudlet['plume'])
        # find the mean motion of the cloudlet
        u_mean = np.mean(u[K, J, I]) - c.ug
        v_mean = np.mean(v[K, J, I]) - c.vg
        w_mean = np.mean(w[K, J, I])
        
        cloudlet['u_plume'] = (u_mean * c.dt / dx)
        cloudlet['v_plume'] = (v_mean * c.dt / dy)
        cloudlet['w_plume'] = (w_mean * c.dt / dz)

    return cloudlets

#----------------------------

def find_cloudlets(time, core, condensed, plume, u, v, w): 
    # find the indexes of all the core and plume points
    core = np.ravel(core)
    condensed = np.ravel(condensed)
    plume = np.ravel(plume)

    plume[condensed] = False
    condensed[core] = False

    # Create the list that will hold the cloudlets
    cloudlets = make_new_cloudlets('core', core)
                    
    for cloudlet in cloudlets:
        cloudlet['condensed'] = cloudlet['core'][:]
            
    ncore = len(cloudlets)
    print(" \t%d core cloudlets" % ncore) 

    cloudlets, condensed = expand_current_cloudlets('condensed', 
                                                    cloudlets,
                                                    condensed)

    # Add any remaining points that have not been added to cloudlets 
    # as new cloudlets.
    condensed_cloudlets = make_new_cloudlets('condensed', condensed)

    for cloudlet in condensed_cloudlets:
        cloudlet['core'] = np.array([], dtype=np.int)
        cloudlets.append(cloudlet)

    for cloudlet in cloudlets:
        cloudlet['plume'] = cloudlet['condensed'][:]

    ncondensed = len(cloudlets)
    print(" \t%d condensed cloudlets" % (ncondensed-ncore))


    cloudlets, plume = expand_current_cloudlets('plume', 
                                                 cloudlets,
                                                 plume)

    # Add any remaining points that have not been added to cloudlets 
    # as new cloudlets.
    plume_cloudlets = make_new_cloudlets('plume', plume)

    for cloudlet in plume_cloudlets:
        cloudlet['core'] = np.array([], dtype=np.int)
        cloudlet['condensed'] = np.array([], dtype=np.int)
        cloudlets.append(cloudlet)

    nplume = len(cloudlets)
    
    print(" \t%d plume cloudlets" % (nplume-ncondensed))

    u, v, w = u.values, v.values, w.values
    cloudlets = find_mean_cloudlet_velocity(cloudlets, u, v, w)

    cloudlet_items = ['core', 'condensed', 'plume', 'u_condensed', \
        'v_condensed', 'w_condensed', 'u_plume', 'v_plume', 'w_plume']
    filename = './hdf5/cloudlets_%08g.h5' % time
    print(" \t%s\n " % filename)
    with h5py.File(filename, 'w') as f:
        # Filter out noisy cloud region
        for n in range(len(cloudlets)):
            # if (cloudlets[n]['plume'].size > 7
            #     or (cloudlets[n]['condensed'].size > 1)
            #     or (cloudlets[n]['core'].size > 0)):

            grp = f.require_group(str(n))
            for var in cloudlet_items:
                if(var in ['core', 'condensed', 'plume']):
                    dset = grp.create_dataset(var, data=cloudlets[n][var][...])
                else:
                    dset = grp.create_dataset(var, data=cloudlets[n][var])

    return 

def load_data(time, ds):
    core = np.array(ds.core.values, dtype=np.bool)
    condensed = np.array(ds.condensed.values, dtype=np.bool)
    plume = np.array(ds.plume.values, dtype=np.bool)
    u = ds.u
    v = ds.v
    w = ds.w

    return core, condensed, plume, u, v, w

def generate_cloudlets():
    # ds = xarray.open_mfdataset((("./data/*.nc")), concat_dim="time")
    ds = xarray.open_mfdataset("data/*.nc")

    for time in range(c.nt):
        core, condensed, plume, u, v, w = load_data(time, ds.isel(time=time))
        find_cloudlets(time, core, condensed, plume, u, v, w)
    
if __name__ == "__main__":
    import doctest
    doctest.testmod()
