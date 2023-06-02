#!/usr/bin/env python
#Runtime (690,130,128,128): ?

import pickle

import networkx

import numpy as np
import glob, h5py, json, dask
from .load_config import c

full_output=True

# def full_output(cloud_times, cloud_graphs, merges, splits):
def full_output(cloud_volumes, cloud_times, cloud_graphs, merges, splits):
    # cloud_times = tuple(cloud_times)

    n = 0
    clouds = {}
    cloud_nodes = {}
    cloud_nodes_w_cloudlets = {}
    # for subgraph in cloud_graphs:
    for cloud_volume, cloud_time, subgraph in zip(cloud_volumes, cloud_times, cloud_graphs):
        # events = {'has_condensed': False, 'has_core': False}
        nodes = {}
        nodes_w_cloudlets = {}
        events = {
                  'plume_volume': cloud_volume[0],
                  'cloud_volume': cloud_volume[1],
                  'core_volume': cloud_volume[2],
                  'plume_time': len(cloud_time[0]),
                  'condense_time': len(cloud_time[1]),
                  'core_time': len(cloud_time[2]),
                  'has_condensed': False,
                  'has_core': False,
                 }

        for node in subgraph:
            node_events = []
            t = int(node[:8])
            if t == 0: node_events.append('break_start')
            if t == (c.nt-1): node_events.append('break_end')

            if node in merges:
                node_events.append('merge_end')
            if node in splits:
                node_events.append('split')
                
            info = subgraph.nodes[node]
            if info['merge']:
                node_events.append('merge')
            if info['split']:
                node_events.append('split_start')
                
            events['has_condensed'] = (info['condensed'] > 0) or events['has_condensed']
            events['has_core'] = (info['core'] > 0) or events['has_core']

            if t in events:
                events[t].extend(node_events)
            else:
                events[t] = node_events[:]

            cloudlet_list = [int(cldt) for cldt in info['cloudlets']]
            if t in nodes_w_cloudlets:
                nodes_w_cloudlets[t][int(node[9:])] = cloudlet_list
            else:
                nodes_w_cloudlets[t] = {int(node[9:]):cloudlet_list}

            if t in nodes:
                nodes[t].append(int(node[9:]))
            else:
                nodes[t] = [int(node[9:])]

        cloud_nodes_w_cloudlets[n] = nodes_w_cloudlets
        cloud_nodes[n] = nodes 
        clouds[n] = events
        n = n + 1

    with open('hdf5/events.json', 'w') as f:
        json.dump(clouds, f, indent=4)
    with open('hdf5/cloud_nodes.json', 'w') as f:
        json.dump(cloud_nodes, f, indent=4)
    with open('hdf5/cloud_nodes_w_cloudlets.json', 'w') as f:
        json.dump(cloud_nodes_w_cloudlets, f, indent=4)


#---------------------

def make_graph():
    graph = networkx.Graph()

    merges = {}
    splits = {}

    # calculate the minimum time duration (5 minutes) in timesteps allowed for a cluster that 
    # has splitting or merging events in its history.
    t_min = 5.*60. / c.dt

    for t in range(c.nt):
        with h5py.File('hdf5/clusters_%08g.h5' % t, 'r') as f:
            keys = np.array(list(f.keys()), dtype=int)
            keys.sort()
            for id in keys:
                m_conns = set(f['%s/merge_connections' % id][...])
                s_conns = set(f['%s/split_connections' % id][...])
                core = len(f['%s/core' % id])
                condensed = len(f['%s/condensed' % id])
                plume = len(f['%s/plume' % id])
                cloudlets = set(f['%s/cloudlet_ids' % id][...])
                # attr_dict = {'merge': m_conns,
                #              'split': s_conns,
                #              'core': core,
                #              'condensed': condensed,
                #              'plume': plume}

                # for item in m_conns:
                #     node1 = '%08g|%08g' % (t, id)
                #     node2 = '%08g|%08g' % (t-1, item)
                #     merges[node2] = node1
                # for item in s_conns:
                #     node1 = '%08g|%08g' % (t, id)
                #     node2 = '%08g|%08g' % (t, item)
                #     splits[node2] = node1            
            
                # Construct a graph of the cloudlet connections
                graph.add_node('%08g|%08g' % (t, id), merge = m_conns,
                                                     split = s_conns,
                                                     core = core,
                                                     condensed = condensed,
                                                     plume = plume,
                                                     cloudlets = cloudlets,
                                                     )
                if f['%s/past_connections' % id]:
                    for item in f['%s/past_connections' % id]:
                        graph.add_edge('%08g|%08g' % (t-1, item),
                                       '%08g|%08g' % (t, id))

    # Iterate over every cloud in the graph
    s = [graph.subgraph(c) for c in networkx.connected_components(graph)]
    for subgraph in s:
        # Find the duration over which the cloud_graph has cloudy points.
        condensed_times = set()
        for node in subgraph:
            if subgraph.nodes[node]['condensed'] > 0:
                condensed_times.add(int(node[:8]))

        plume_times = set()
        for node in subgraph:
            if subgraph.nodes[node]['plume'] > 0:
                plume_times.add(int(node[:8]))

        # If cloud exists for less than 5 minutes, check if it has split events
        # If it has split events, remove them and reconnect the cloud
        if (len(condensed_times) < t_min):
            for node in subgraph:
                if subgraph.nodes[node]['split']:
                    item = subgraph.nodes[node]['split'].pop()
                    t = int(node[:8]) 
                    graph.add_edge(node, '%08g|%08g' % (t, item))

        if (len(plume_times) < t_min):
            for node in subgraph:
                if subgraph.nodes[node]['split']:
                    item = subgraph.nodes[node]['split'].pop()
                    t = int(node[:8]) 
                    graph.add_edge(node, '%08g|%08g' % (t, item))

    s = [graph.subgraph(c) for c in networkx.connected_components(graph)]
    for subgraph in s:
        # Find the duration over which the cloud_graph has cloudy points.
        condensed_times = set()
        for node in subgraph:
            if subgraph.nodes[node]['condensed'] > 0:
                condensed_times.add(int(node[:8]))

        plume_times = set()
        for node in subgraph:
            if subgraph.nodes[node]['plume'] > 0:
                plume_times.add(int(node[:8]))

        # If a cloud exists less than 5 minutes, check for merge events
        if (len(condensed_times) < t_min):
            for node in subgraph:
                if subgraph.nodes[node]['merge']:
                    item = subgraph.nodes[node]['merge'].pop()
                    t = int(node[:8])
                    graph.add_edge(node, '%08g|%08g' % (t-1, item))

        if (len(plume_times) < t_min):
            for node in subgraph:
                if subgraph.nodes[node]['merge']:
                    item = subgraph.nodes[node]['merge'].pop()
                    t = int(node[:8])
                    graph.add_edge(node, '%08g|%08g' % (t-1, item))

    for node in graph:
        t = int(node[:8])
        for item in graph.nodes[node]['merge']:
            node2 = '%08g|%08g' % (t-1, item)
            merges[node2] = node
        for item in graph.nodes[node]['split']:
            node2 = '%08g|%08g' % (t, item)
            splits[node2] = node     

    cloud_graphs = []
    cloud_noise = []
    s = [graph.subgraph(c) for c in networkx.connected_components(graph)]
    for subgraph in s:
        plume_time = set()
        condensed_time = set()
        core_time = set()

        plume_volume = 0
        condensed_volume = 0
        core_volume = 0

        for node in subgraph:
            condensed_vol = subgraph.nodes[node]['condensed'] 
            if condensed_vol > 0:
                condensed_volume = condensed_volume + condensed_vol
                condensed_time.add(int(node[:8]))

            core_vol = subgraph.nodes[node]['core'] 
            if core_vol > 0:
                core_volume = core_volume + core_vol
                core_time.add(int(node[:8]))

            plume_vol = subgraph.nodes[node]['plume']
            if plume_vol > 0:
                plume_volume = plume_volume + plume_vol
                plume_time.add(int(node[:8]))

        if (len(condensed_time) < 2) and (len(core_time) == 0):
            cloud_noise.append(subgraph)
        else:
            plume_time = list(plume_time)
            plume_time.sort()
            condensed_time = list(condensed_time)
            condensed_time.sort()
            core_time = list(core_time)
            core_time.sort()
            cloud_graphs.append((plume_volume, condensed_volume, core_volume, subgraph, \
                                plume_time, condensed_time, core_time))
            
    # cloud_graphs.sort(key=lambda x:x[0])
    cloud_graphs.sort(key=lambda x:x[1])
    cloud_graphs.reverse()
    cloud_volumes = [item[0:3] for item in cloud_graphs] 
    cloud_times = [item[4:] for item in cloud_graphs] 
    cloud_graphs = [item[3] for item in cloud_graphs]
    
    if full_output: full_output(cloud_volumes, cloud_times, cloud_graphs, merges, splits)
    # if full_output: full_output(cloud_times, cloud_graphs, merges, splits)
    return cloud_graphs, cloud_noise
