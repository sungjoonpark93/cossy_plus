__author__ = 'SungJoonPark'
__author__ = 'SungJoonPark'

import pandas as pd
import numpy as np

def get_network(network_file):
    #get sif file formatted network file and return to adjacency matrix

    #network file format should be sif
    network_dict = {}
    #construct dict
    with open(network_file, 'r') as r:
        for line in r:
            node1= line.strip().split("\t")[0]
            node2= line.strip().split("\t")[2]
            edge_type = line.strip().split("\t")[1]
            if node1 not in network_dict:
                network_dict[node1] = {node2:1}
            else:
                if node2 not in network_dict[node1]:
                    network_dict[node1][node2]=1

            if node2 not in network_dict:
                network_dict[node2] = {node1:1}
            else:
                if node1 not in network_dict[node2]:
                    network_dict[node2][node1]=1

    #append the key itself
    for key in network_dict.keys():
        if key not in network_dict[key].keys():
            network_dict[key][key]=1


    #append ajacency matrix
    network_df = pd.DataFrame(network_dict).fillna(0)
    #filter where gene symbol is like a strange integer
    #filter_network = network_df.iloc[18:,18:]
    #after string excel problem
    filter_network = network_df.iloc[0:,0:]
    #make diagonal matrix where element is rowsum of the filter network
    diagnal_matrix = pd.DataFrame(np.diag(np.sum(filter_network,axis=0)),index=filter_network.index, columns=filter_network.columns)
    #get non zero index from filter network
    non_zero_index = list(filter_network[filter_network != 0].stack().index)
    #normalize filter network
    for non_zero_tuple in non_zero_index:
        filter_network.loc[non_zero_tuple[0],non_zero_tuple[1]] = float(filter_network.loc[non_zero_tuple[0],non_zero_tuple[1]]) / np.sqrt(diagnal_matrix.loc[non_zero_tuple[0],non_zero_tuple[0]]*diagnal_matrix.loc[non_zero_tuple[1],non_zero_tuple[1]])

    return filter_network


