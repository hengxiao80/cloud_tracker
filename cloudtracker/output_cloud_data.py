import h5py
import networkx
import numpy as np

from .utility_functions import zyx_to_index, index_to_zyx, \
                                calc_radii, expand_indexes

from .load_config import c

def calc_shell(index):
    # Expand the cloud points outward
    maskindex = expand_indexes(index)
    
    # From the expanded mask, select the points outside the cloud
    shellindex = np.setdiff1d(maskindex, index, assume_unique=True)

    return shellindex

def calc_edge(index, shellindex):
    # Find all the points just inside the clouds
    maskindex = expand_indexes(shellindex)
    
    # From the expanded mask, select the points inside the cloud
    edgeindex = np.intersect1d(maskindex, index, assume_unique=True)

    return edgeindex

def calc_env(index, shellindex, edgeindex):
    if len(shellindex) > 0:
        K_J_I = index_to_zyx(shellindex)

        n = 6
        for i in range(n):
            # Find all the points just outside the clouds
            stacklist = [K_J_I, ]
            for item in ((0, -1, 0), (0, 1, 0),
                         (0, 0, -1), (0, 0, 1)):
                stacklist.append( K_J_I + np.array(item)[:, np.newaxis] )

            maskindex = np.hstack(stacklist)
            maskindex[1, :] = maskindex[1, :] % c.ny
            maskindex[2, :] = maskindex[2, :] % c.nx            
            maskindex = np.unique( zyx_to_index(maskindex[0, :],
                                                maskindex[1, :],
                                                maskindex[2, :]) )

            # From the expanded mask, select the points outside the cloud
            envindex = np.setdiff1d(maskindex, index, assume_unique = True)

            K_J_I = index_to_zyx(envindex)
            
        # Select the points within 4 grid cells of cloud
        r = calc_radii(envindex, edgeindex)
        mask = r < 4.5
        envindex = envindex[mask]
    else:
        envindex = []

    return envindex

def calc_regions(cluster):
    result = {}

    result['plume'] = cluster['plume']

    condensed = cluster['condensed']
    result['condensed'] = condensed
    condensed_shell = calc_shell(condensed)
    result['condensed_shell'] = condensed_shell
    condensed_edge = calc_edge(condensed, condensed_shell)
    result['condensed_edge'] = condensed_edge
    result['condensed_env'] = calc_env(condensed, condensed_shell, condensed_edge)

    core = cluster['core']
    result['core'] = core
    core_shell = calc_shell(core)
    result['core_shell'] = core_shell
    core_edge = calc_edge(core, core_shell)
    result['core_edge'] = core_edge
    result['core_env'] = calc_env(core, core_shell, core_edge)

    return result

def output_clouds_at_time(cloud_graphs, cloud_noise, t):
    print('Timestep:', t)

    # cluster = {}
    clusters = {}
    items = ['core', 'condensed', 'plume']
    
    with h5py.File('hdf5/clusters_%08g.h5' % t, 'r') as cluster_dict:
        keys = np.array(list(cluster_dict.keys()), dtype=int)
        keys.sort()
        for id in keys:
            key = "%08g|%08g" % (t, id)

            # clusters[key] = dict(zip(items, np.array([cluster_dict['%s/%s' % (id, 'core')][...], \
            #     cluster_dict['%s/%s' % (id, 'condensed')][...], cluster_dict['%s/%s' % (id, 'plume')][...]])))
            clusters[key] = dict(zip(items, [cluster_dict['%s/%s' % (id, 'core')][...], \
                cluster_dict['%s/%s' % (id, 'condensed')][...], cluster_dict['%s/%s' % (id, 'plume')][...]]))
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
            clouds[id] = calc_regions(cloud)
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
    clouds[-1] = calc_regions(noise_clust)
            
    print("Number of Clouds at Current Timestep: ", len(clouds.keys()) + 1)

    items = ['core', 'condensed', 'plume', 'core_shell', 'condensed_shell']#, \
        # 'core_edge', 'condensed_edge', 'core_env', 'condensed_env']
    with h5py.File('hdf5/clouds_%08g.h5' % t, 'w') as f:
        for id in clouds:
            grp = f.create_group(str(id))
            for point_type in clouds[id]:
                dset = grp.create_dataset(point_type, data=clouds[id][point_type])

def output_cloud_data(cloud_graphs, cloud_noise):
    for time in range(c.nt):
        print("\n Outputting cloud data, time step: %d" % time)
        output_clouds_at_time(cloud_graphs, cloud_noise, time)
