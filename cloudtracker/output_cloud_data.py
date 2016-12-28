#!/usr/bin/env python
# Runtime (690, 130, 128, 128): 1 hour 20 minutes

import pickle
import h5py
import networkx
import numpy as np
from .utility_functions import zyx_to_index, index_to_zyx, calc_radii, \
    expand_indexes
import sys, gc
#import scipy.io

def calc_shell(index, MC):
    # Expand the cloud points outward
    maskindex = expand_indexes(index, MC['nz'], MC['ny'], MC['nx'])
    
    # From the expanded mask, select the points outside the cloud
    shellindex = np.setdiff1d(maskindex, index, assume_unique=True)

    return shellindex

def calc_edge(index, shellindex, MC):
    # Find all the points just inside the clouds
    maskindex = expand_indexes(shellindex, MC['nz'], MC['ny'], MC['nx'])
    
    # From the expanded mask, select the points inside the cloud
    edgeindex = np.intersect1d(maskindex, index, assume_unique=True)

    return edgeindex

def calc_env(index, shellindex, edgeindex, MC):
    if len(shellindex) > 0:
        K_J_I = index_to_zyx(shellindex, MC['nz'], MC['ny'], MC['nx'])

        n = 6
        for i in range(n):
            # Find all the points just outside the clouds
            stacklist = [K_J_I, ]
            for item in ((0, -1, 0), (0, 1, 0),
                         (0, 0, -1), (0, 0, 1)):
                stacklist.append( K_J_I + np.array(item)[:, np.newaxis] )

            maskindex = np.hstack(stacklist)
            maskindex[1, :] = maskindex[1, :] % MC['ny']
            maskindex[2, :] = maskindex[2, :] % MC['nx']            
            maskindex = np.unique( zyx_to_index(maskindex[0, :],
                                                   maskindex[1, :],
                                                   maskindex[2, :],
                                                   MC['nz'], MC['ny'], MC['nx']) )

            # From the expanded mask, select the points outside the cloud
            envindex = np.setdiff1d(maskindex, index, assume_unique = True)

            K_J_I = index_to_zyx(envindex, MC['nz'], MC['ny'], MC['nx'])
            
        # Select the points within 4 grid cells of cloud
        r = calc_radii(envindex, edgeindex, MC)
        mask = r < 4.5
        envindex = envindex[mask]
    else:
        envindex = []

    return envindex

def calculate_data(cluster, MC):
    result = {}

    result['plume'] = cluster['plume']

    condensed = cluster['condensed']
    result['condensed'] = condensed
    condensed_shell = calc_shell(condensed, MC)
    result['condensed_shell'] = condensed_shell
    condensed_edge = calc_edge(condensed, condensed_shell, MC)
    result['condensed_edge'] = condensed_edge
    result['condensed_env'] = calc_env(condensed, condensed_shell, condensed_edge, MC)

    core = cluster['core']
    result['core'] = core
    core_shell = calc_shell(core, MC)
    result['core_shell'] = core_shell
    core_edge = calc_edge(core, core_shell, MC)
    result['core_edge'] = core_edge
    result['core_env'] = calc_env(core, core_shell, core_edge, MC)

    return result

# def save_text_file(clouds, t, MC):
#     count = 0

#     for id in clouds:
#         for point_type in clouds[id]:
#             count = count + len(clouds[id][point_type])

#     recarray = np.zeros(count, dtype=[('id', 'i4'),('type', 'a14'),('x','i4'),('y', 'i4'), ('z', 'i4')])

#     count = 0
#     for id in clouds:
#         for point_type in clouds[id]:
#             data = clouds[id][point_type]
#             n = len(data)
#             if n == 0: continue
#             recarray['id'][count:n + count] = id
#             recarray['type'][count:n + count] = point_type
#             z, y, x = index_to_zyx(data, MC)
#             recarray['x'][count:n + count] = x
#             recarray['y'][count:n + count] = y
#             recarray['z'][count:n + count] = z
#             count = count + n
            
#     recarray.tofile(open('output/clouds_at_time_%08g.txt' % t, 'w'), '\r\n')

def output_cloud_data(cloud_graphs, cloud_noise, t):
    print('Timestep:', t)

    cluster = {}
    clusters = {}
    items = ['core', 'condensed', 'plume']
    attribute_items = ['time', 'nx', 'ny', 'nz', 'dx', 'dy', 'dz', 'dt', 'ug', 'vg']
    
    with h5py.File('cloudtracker/hdf5/clusters_%08g.h5' % t, 'r') as cluster_dict:
        # Read model parameters
        MC = {}
        for item in attribute_items:
            MC[item] = cluster_dict.attrs[item] 

        keys = np.array(list(cluster_dict.keys()), dtype=int)
        keys.sort()
        for id in keys:
            key = "%08g|%08g" % (t, id)

            clusters[key] = dict(zip(items, np.array([cluster_dict['%s/%s' % (id, 'core')][...], \
                cluster_dict['%s/%s' % (id, 'condensed')][...], cluster_dict['%s/%s' % (id, 'plume')][...]])))
            # for var in items:
            #     cluster[var] = cluster_dict['%s/%s' % (id, var)][...]
            # clusters[key] = cluster

    clouds = {}
    id = 0
    for subgraph in cloud_graphs:
        # Grab the nodes at the current time 
        # that all belong to subgraph 'id'
        nodes = [item for item in subgraph.nodes() 
                      if item[:8] == ('%08g' % t)]
                      
        if nodes:
            # Pack them into a single cloud object
            core = []
            condensed = []
            plume = []
            for node in nodes:
                core.append(clusters[node]['core'])
                condensed.append(clusters[node]['condensed'])
                plume.append(clusters[node]['plume'])
                
            cloud = {'core': np.hstack(core),
                     'condensed': np.hstack(condensed),
                     'plume': np.hstack(plume)}

            # Calculate core/cloud, env, shell and edge
            clouds[id] = calculate_data(cloud, MC)
        id += 1

    # Add all the noise to a noise cluster
    noise_clust = {'core': [], 'condensed': [], 'plume': []}
    for subgraph in cloud_noise:
        nodes = [item for item in subgraph.nodes() 
                       if item[:8] == ('%08g' % t)]
        if nodes:
            for node in nodes:
                noise_clust['core'].append(clusters[node]['core'])
                noise_clust['condensed'].append(clusters[node]['condensed'])
                noise_clust['plume'].append(clusters[node]['plume'])
                    
    if noise_clust['core']:                    
        noise_clust['core'] = np.hstack(noise_clust['core'])         
    if noise_clust['condensed']: 
        noise_clust['condensed'] = np.hstack(noise_clust['condensed'])
    if noise_clust['plume']:
        noise_clust['plume'] = np.hstack(noise_clust['plume'])

    # Only save the noise if it contains cloud core
    clouds[-1] = calculate_data(noise_clust, MC)
            
    print("Number of Clouds at Current Timestep: ", len(clouds.keys()) + 1)

    items = ['core', 'condensed', 'plume', 'core_shell', 'condensed_shell']#, \
        # 'core_edge', 'condensed_edge', 'core_env', 'condensed_env']
    with h5py.File('cloudtracker/hdf5/clouds_%08g.h5' % t, 'w') as f:
        for id in clouds:
            grp = f.create_group(str(id))
            for point_type in clouds[id]:
                dset = grp.create_dataset(point_type, data=clouds[id][point_type])

    # save_text_file(clouds, t, MC)

#   save .mat file for matlab
#    new_dict = {}
#    for key in clouds:
#        new_dict['cloud_%08g' % key ] = clouds[key]
    
#    savedict = {'clouds': new_dict} 
#    scipy.io.savemat('mat/cloud_data_%08g.mat' % t, savedict)
    

