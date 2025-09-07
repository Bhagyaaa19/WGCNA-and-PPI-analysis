# 12/05/2024
# Author - Bhagya WIjeratne

# Input: STRING datafile of global PPI links for Zea mays
# Output: PPI global network

# Creates the global PPI network for maize

import networkx as nx
from collections import OrderedDict
from networkx import neighbors

# creating an empty graph for the imported protein network
protein_graph = nx.Graph()

# opening the network file for all proteins
with open(
        "C:\\Users\\BhagyaWijeratne\\Desktop\\Level IV\\Research\\methodolgy\\PPI network\\4577.protein.links.v12.0.txt") as network_file:
    # skip first line as it contains headers
    next(network_file)

    # Make a network with node 1, node2
    for line in network_file:  # iterating through the file line by line

        columns = line.split(' ')  # split the lines in to columns based on the tab
        # extracting the required information from the file
        # turning all the data in to upper case to avoid mismatches
        # print(columns)
        node1 = columns[0].upper()
        # print(node1)
        node2 = columns[1].upper()
        combined_score = columns[2]
        # adding edges to the graph
        protein_graph.add_edge(node1, node2, weight = combined_score)

    # creating a set to store all proteins without duplicates
    # all_proteins_set = set(list(protein_graph.nodes))

# converting to a file such that cytoscape can open
nx.write_gml(protein_graph, "maize_global.gml")
# Layout >> prefuse force directed OpenCL Layout

print("Program has successfully run")
