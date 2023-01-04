import sys, os

from .generate_cloudlets import generate_cloudlets
from .cluster_cloudlets import cluster_cloudlets
from .make_graph import make_graph
from .output_cloud_data import output_cloud_data

def main():
    sys.setrecursionlimit(1000000)

    if not os.path.exists('./hdf5'):
        os.mkdir('./hdf5')

    # Cloudlets
    print( " Gathering cloudlets... " )
    generate_cloudlets()

    # Clusters
    print( "\n Creating clusters... " )
    cluster_cloudlets()

    # Cloud graph
    print( "\n Make graph... " )
    cloud_graphs, cloud_noise = make_graph()
    print( "\tFound %d clouds" % len(cloud_graphs))
            
    # Cloud output
    output_cloud_data(cloud_graphs, cloud_noise)
        