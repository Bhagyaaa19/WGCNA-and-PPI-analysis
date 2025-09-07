"""
02. 03. 2025
Author - Bhagya Wijeratne

Input - PPI network module
Output - List of submodules, enrichment analysis and hub analysis

Identify submodules, enrichment analysis and hub analysis
"""
import os
import pandas as pd
import community.community_louvain as community_louvain
from collections import defaultdict
import networkx as nx
from gprofiler import GProfiler
import numpy as np
import re
import traceback

# reading the dataset as a gml file
# NOTE: Change here when the module is changed
module_color = "black_mature_foliar"  # NOTE: Change module name
# graphml_file = r"C:\Users\BhagyaWijeratne\Desktop\Level IV\Research\methodolgy\PPI network\cytoscape_network_data\turquoise_network.graphml"
# graphml_file = r"C:\Users\BhagyaWijeratne\Desktop\Level IV\Research\methodolgy\PPI network\cytoscape_network_data\blue_network.graphml"
# graphml_file = r"C:\Users\BhagyaWijeratne\Desktop\Level IV\Research\methodolgy\PPI network\cytoscape_network_data\brown_network.graphml"
# graphml_file =5 r"C:\Users\BhagyaWijeratne\Desktop\Level IV\Research\methodolgy\PPI network\cytoscape_network_data\Yellow network.graphml"
# graphml_file = r"C:\Users\BhagyaWijeratne\Desktop\Level IV\Research\methodolgy\PPI network\cytoscape_network_data\green_network.graphml"
# graphml_file = r"C:\Users\BhagyaWijeratne\Desktop\Level IV\Research\methodolgy\PPI network\cytoscape_network_data\black network.graphml"
# graphml_file = r"C:\Users\BhagyaWijeratne\Desktop\Level IV\Research\methodolgy\PPI network\cytoscape_network_data\dark_red_FP.graphml"
# graphml_file = r"C:\Users\BhagyaWijeratne\Desktop\Level IV\Research\methodolgy\PPI network\cytoscape_network_data\skyblue and paleturqiouse - FP34.graphml"
# graphml_file = r"C:\Users\BhagyaWijeratne\Desktop\Level IV\Research\methodolgy\PPI network\cytoscape_network_data\HP34_darkmag_mag_greenyellow_darkturq.graphml"
# graphml_file = r"C:\Users\BhagyaWijeratne\Desktop\Level IV\Research\methodolgy\PPI network\cytoscape_network_data\FI_black.graphml"
# graphml_file = r"C:\Users\BhagyaWijeratne\Desktop\Level IV\Research\methodolgy\PPI network\cytoscape_network_data\FI_combined.graphml"
# graphml_file = r"C:\Users\BhagyaWijeratne\Desktop\Level IV\Research\methodolgy\PPI network\cytoscape_network_data\HE_yellow.graphml"
# graphml_file = r"C:\Users\BhagyaWijeratne\Desktop\Level IV\Research\methodolgy\PPI network\cytoscape_network_data\HE combined.graphml"
# graphml_file = r"C:\Users\BhagyaWijeratne\Desktop\Level IV\Research\methodolgy\PPI network\cytoscape_network_data\maize_husk_red.graphml"
# graphml_file = r"C:\Users\BhagyaWijeratne\Desktop\Level IV\Research\methodolgy\PPI network\cytoscape_network_data\blue_primordial_foliar_new.graphml"
# graphml_file = r"C:\Users\BhagyaWijeratne\Desktop\Level IV\Research\methodolgy\PPI network\cytoscape_network_data\primordial_husk_brown.graphml"
# graphml_file = r"C:\Users\BhagyaWijeratne\Desktop\Level IV\Research\methodolgy\PPI network\cytoscape_network_data\Merged Network_mature_foliar_down.graphml"
# graphml_file = r"C:\Users\BhagyaWijeratne\Desktop\Level IV\Research\methodolgy\PPI network\cytoscape_network_data\Merged Network Mature_Husk.graphml"
graphml_file = r"C:\Users\BhagyaWijeratne\Desktop\Level IV\Research\methodolgy\PPI network\cytoscape_network_data\mature_foliar_black.graphml"
G = nx.read_graphml(graphml_file)

# Define base directory for results
# base_dir = os.path.join("C:/Users/BhagyaWijeratne/Desktop/Level IV/Research/methodolgy/PPI network",
#                        "er_analysis_results")
# os.makedirs(base_dir, exist_ok=True)

# Rename nodes using their 'name' attribute
mapping = {node: G.nodes[node].get("name", node) for node in G.nodes}
G = nx.relabel_nodes(G, mapping)


# Get node attributes
# node_attributes = set()
# for _, data in G.nodes(data=True):
#    node_attributes.update(data.keys())

# print("Node attributes in the GraphML file:")
# print(node_attributes)
# # Rename problematic attributes - issues with having spaces in between the column headers
# for node, data in G.nodes(data=True):
#     for key in list(data.keys()):  # Use list() to avoid runtime error while modifying dict
#         new_key = key.replace(" ", "_")  # Replace spaces with underscores
#         if new_key != key:
#             data[new_key] = data.pop(key)  # Rename the attribute
#
# for edge, data in G.edges(data=True):
#     for key in list(data.keys()):
#         new_key = key.replace(" ", "_")
#         if new_key != key:
#             data[new_key] = data.pop(key)
#
# # Save as GML
# nx.write_gml(G, "black_network.gml")
# print("Conversion complete: .gml file saved.")
#
# network_file = "C:\\Users\\BhagyaWijeratne\\Desktop\\Level IV\\Research\\methodolgy\\PPI network\\cytoscape_network_data\\black network.gml"
#
# selected_network = nx.read_gml(network_file)
#
# # remove Isolated Nodes (Proteins Without Interactions)
# isolated_nodes = list(nx.isolates(selected_network))
# selected_network.remove_nodes_from(isolated_nodes)
# print(selected_network)

def clean_attribute_keys(data_dict):
    """Cleans attribute keys to be GML-compliant."""
    # print("clean_attribute_keys called!")
    cleaned = {}
    # print("dictionary", data_dict)
    for key1, value in data_dict.items():
        if key1 is None:
            # print("Invalid key:", key1)
            continue  # Skip None keys

        original_key = str(key1)
        # print(f"Original key bytes: {original_key.encode('utf-8')}")

        # Replace problematic characters
        new_key1 = original_key.replace("::", "_")
        new_key1 = re.sub(r'[^a-zA-Z0-9_]', '_', new_key1)  # Replace non-alphanumeric with "_"
        new_key1 = new_key1.strip('_')  # Remove leading/trailing underscores

        # print(f"Original: {original_key}, Cleaned Key: {new_key1}")

        # If key is empty after cleaning, skip it
        if new_key1:
            cleaned[new_key1] = value
            # print(f"Updated key: {original_key} → {new_key1}")
        # else:
            # print(f"Skipping invalid key: {original_key}")

    return cleaned


def clean_node_attributes(graph):
    """Cleans node attributes in a NetworkX graph."""
    # print("Cleaning the nodes")
    for node, data in list(graph.nodes(data=True)):
        # print(f"Node ({node}) before cleaning: {data}")
        try:
            cleaned_data = clean_attribute_keys(data)
            # print(f"Node ({node}) after cleaning: {cleaned_data}")
            if cleaned_data:
                graph.nodes[node].clear()  # Remove all existing attributes
                graph.nodes[node].update(cleaned_data)  # Apply cleaned attributes
                # print(f"Updated Node {node}: {cleaned_data}")

        except Exception as e:
            # print(f"Error cleaning node {node}: {e}")
            traceback.print_exc()
    return graph


def clean_edge_attributes(graph):
    """Cleans edge attributes in a graph's edge data."""
    for u, v, data in list(graph.edges(data=True)):  # Convert iterator to list to avoid modification issues
        # print(f"Edge ({u}, {v}) data before cleaning: {data}")  # check before cleaning.
        cleaned_data = clean_attribute_keys(data)
        # print(f"Edge ({u}, {v}) data after cleaning: {cleaned_data}")  # check after cleaning.
        if cleaned_data:
            graph.edges[u, v].clear()  # Remove all existing attributes
            graph.edges[u, v].update(cleaned_data)
            # print("Cleaned the edges")

    return graph  # Ensure the function returns the graph


def clean_graph_attributes(graph):
    """Cleans graph attributes in a NetworkX graph."""
    # print("Cleaning the graph attributes")
    cleaned_data = clean_attribute_keys(graph.graph)
    graph.graph.clear()  # clear the old graph attributes.
    graph.graph.update(cleaned_data)
    return graph  # Ensure the function returns the graph


# Clean attribute names (remove spaces)
for node, data in G.nodes(data=True):
    for key in list(data.keys()):
        new_key = key.replace(" ", "_")
        if new_key != key:
            data[new_key] = data.pop(key)

for node1, node2, data in G.edges(data=True):
    for key in list(data.keys()):
        new_key = key.replace(" ", "_")
        if new_key != key:
            data[new_key] = data.pop(key)

# Remove isolated nodes
isolated_nodes = list(nx.isolates(G))
G.remove_nodes_from(isolated_nodes)

print("Graph cleaned and isolated nodes removed.")

if len(G.nodes) == 0:
    raise ValueError("The graph has no nodes left after removing isolated nodes!")

# Convert directed graph to undirected
if G.is_directed():
    G = G.to_undirected()

# debugging the z score issue
# print("edge length:", len(G.edges))
# ("Sample of edges:", list(G.edges)[:10])  # View a sample of edges
#  1931 edges, with the edges represented as pairs of protein IDs.
node_ids = set(G.nodes)
edge_node_ids = {u for u, v in G.edges} | {v for u, v in G.edges}
# print(f"Nodes without edges: {node_ids - edge_node_ids}")

is_connected = nx.is_connected(G)
print("Is the graph connected?", is_connected)

if not is_connected:
    components = list(nx.connected_components(G))
    print(f"The graph has {len(components)} connected components.")

# 4 connected components: Even though no nodes are isolated, the graph as a whole is not connected. It consists of 4 separate groups of nodes, where nodes within each group are connected to each other, but there are no edges between the different groups.
# we have to analyse each component separately

# community detection using Louvain algorithm
print("Starting Louvain clustering...")

for i, component_nodes in enumerate(components, start=1):
    # print(i)
    # print(component_nodes)
    subgraph = G.subgraph(component_nodes).copy()
    print(f"Subgraph {i} has {len(subgraph.nodes)} nodes and {len(subgraph.edges)} edges")

    # Perform Louvain clustering on subgraph
    sub_partitions = community_louvain.best_partition(subgraph, resolution=1.44, random_state=1)
    # Adjust partition numbers (increment by 1) to avoid 0-based indexing
    sub_partitions = {node: community + 1 for node, community in sub_partitions.items()}
    print("Sub-partitions:", sub_partitions)  # Debugging

    print("Done Louvain clustering...")

    # Perform enrichment and hub analysis on each subgraph
    # Perform similar steps for clustering results, enrichment analysis, and hub detection

    # Initial workflow for partitioning ========================================================================================
    # partitions = community_louvain.best_partition(G, resolution=1.44, random_state=1)

    # Clusters are renumbered to avoid zero-based indexing. Adjusts cluster IDs to start from 1 (instead of 0).
    # partitions = {node: community + 1 for node, community in partitions.items()}

    # Proteins are grouped by their clusters.
    # The results are stored in a DataFrame, where columns represent clusters and rows contain proteins.

    results = defaultdict(list)

    for node, cluster in sorted(sub_partitions.items()):
        # Get the database identifier if it exists, otherwise use the node name
        db_id = subgraph.nodes[node].get("name", node)
        results[cluster].append(db_id)  # Store the identifier instead of the node ID

        # Sorting results - so that the cluster numbers appear in order
        sorted_results = dict(sorted(results.items()))

        # Store results in a DataFrame
        df = pd.DataFrame.from_dict(sorted_results, orient='index')
        df = df.transpose()  # Transpose the DataFrame to have clusters as columns

        # Remove "4577." prefix from all values in the DataFrame
        df = df.applymap(lambda x: x.replace("4577.", "") if isinstance(x, str) else x)

        # print(df.head())  # Print first few rows to verify

    # ===================================================================================================
    # # mapping the genes IDs to the string IDs
    # # Load the alias file
    # alias_file = "C:\\Users\\BhagyaWijeratne\\Desktop\\Level IV\\Research\\methodolgy\\PPI network\\4577.protein.aliases\\4577.protein.aliases.v12.0.txt"  # path to the file
    # aliases = pd.read_csv(alias_file, sep="\t", header=None, names=["STRING_ID", "Alias", "Source"])
    # print(aliases.head())  # Show the first 5 rows
    # print(aliases.info())  # Show column types and missing values
    #
    # # Check for missing values
    # if aliases.isnull().sum().sum() > 0:
    #     print("Warning: Missing values detected in alias file!")
    #     aliases = aliases.dropna()
    #
    # # Create dictionary safely
    # # Create a mapping where each gene ID can have multiple STRING IDs
    # gene_to_string = aliases.groupby("Alias")["STRING_ID"].apply(list).to_dict()
    #
    #
    # # Function to map a gene ID to its STRING ID(s)
    # def map_gene_to_string(gene_id):
    #     return gene_to_string.get(gene_id, gene_id)  # Keep original ID if no mapping is found
    #
    #
    # # Apply mapping to the entire DataFrame
    # df = df.applymap(map_gene_to_string)
    #
    # # Flatten lists where multiple mappings exist
    # df = df.applymap(lambda x: x[0] if isinstance(x, list) and len(x) > 0 else x)
    # print(f"Created mapping for {len(gene_to_string)} genes.")
    #
    # # Check a few mappings
    # print(list(gene_to_string.items())[:10])
    #
    # #  replace gene IDs in your DataFrame using this mapping.
    # df.replace(to_replace=gene_to_string, inplace=True)  # This will convert gene IDs into STRING IDs where possible.
    #
    # # Some gene IDs might not have a direct match. To check how many were converted:
    # num_converted = df.isin(gene_to_string.values()).sum().sum()
    # print(f"Converted {num_converted} out of {df.size} gene IDs.")
    # # For missing values, we need to manually check or use alternative identifiers like UniProt.
    # df.to_csv("C:/Users/BhagyaWijeratne/Desktop/Level IV/Research/methodolgy/PPI network/mapped_clusters.csv")

    # Perform Enrichment Analysis (GO Terms)
    gp = GProfiler(return_dataframe=True)  # Initializes GProfiler for Gene Ontology (GO) enrichment analysis.

    # debugging
    # print(f"Nodes in graph: {list(G.nodes)[:10]}")  # Print first 10 nodes
    # print(f"Nodes in graph: {[G.nodes[node].get('name', node) for node in list(G.nodes)[:10]]}")
    # print(f"Nodes in clusters: {list(results.values())[:10]}")  # Print first 10 cluster members

    column_dict = {column: df[column].dropna().tolist() for column in df.columns}
    # NOTE: SAVE POINT
    er_output_file = f"C:/Users/BhagyaWijeratne/Desktop/Level IV/Research/methodolgy/PPI network/results/er analysis results/{module_color}/er_analysis_{module_color}_component_{i}.xlsx"

    # Extract the parent directory
    er_parent_dir = os.path.dirname(er_output_file)

    # Create the directory if it does not exist
    os.makedirs(er_parent_dir, exist_ok=True)

    # Create a workbook and add a default visible sheet
    with pd.ExcelWriter(er_output_file, engine="openpyxl", mode="w") as w:
        workbook = w.book
        default_sheet = workbook.create_sheet("Default")
        workbook.active = default_sheet  # Ensure there's always an active sheet

        valid_sheets = 0  # Track valid clusters written

        for key, val in column_dict.items():
            # print(f"Cluster {key} has {len(val)} proteins")  # Debugging step

            if not val or all(pd.isna(val)):
                print(f"Skipping Cluster {key} (No valid STRING IDs)")
                continue

            try:
                er = gp.profile(query=val, organism='zmays', sources=['GO'])
                if er.empty:
                    print(f"No significant enrichment found for Cluster {key}")
                else:
                    er.to_excel(w, sheet_name=f'cluster_{key}', index=False)
                    valid_sheets += 1  # Increment valid sheet count
            except Exception as e:
                print(f"Error processing Cluster {key}: {e}")

        # If at least one valid sheet was created, delete "Default"
        if valid_sheets > 0 and "Default" in workbook.sheetnames:
            del workbook["Default"]

    # Compute Hub Analysis
    # Computes the average and standard deviation of node degrees within each cluster.
    # Degree = number of direct neighbors within the same cluster.

    average_values = {}  # Stores the average degree of nodes within each cluster
    stdev_values = {}  # Stores the standard deviation of degrees within each cluster
    degrees_within_cluster = {}  # Stores individual node degrees within their respective clusters

    for cluster, proteins in sorted(results.items()):
        # print("Cluster:", cluster, "Proteins:", proteins) - debugging
        degree_list = []
        for node in proteins:
            #     print(node)
            #     for n in G.neighbors(node):
            #         print(G.neighbors(node))
            # print(G.nodes[n].get("name"))
            if node not in subgraph:
                # print(f"Warning: Node {node} not found in graph!")  # Debugging
                continue  # Skip this node
            # neighbors = list(subgraph.neighbors(node))  # Get all neighbors
            # print(f"Node {node} neighbors before filtering: {neighbors}")
            # for n in subgraph.neighbors(node):
            # print(f"Node {n}: Subgraph nodes {subgraph.nodes[n].get('shared_name')}")
            within_degree = len(
                [n for n in subgraph.neighbors(node) if subgraph.nodes[n].get('shared_name') in proteins])
            # print(f"Node {node} in Cluster {cluster} has degree {within_degree}")

            degrees_within_cluster[node] = within_degree
            degree_list.append(within_degree)

        average_values[cluster] = np.average(degree_list)
        stdev_values[cluster] = np.std(degree_list)
    print("Checkpoint1")
    # Calculate Z-Scores and Identify Hubs
    # Z-score measures how different a protein’s degree is compared to others in the same cluster.
    z_scores = []
    intra_hub_data = []
    inter_hub_data = []

    for node in sub_partitions:
        try:
            string_id = df[df.eq(node).any(axis=1)].stack().iloc[0]  # Safely get the first match
        except IndexError:
            string_id = node  # If not found, use the original node name
            # if node not in subgraph.nodes:
            #     print(f"Warning: Node {node} from sub_partitions is missing in the subgraph!")

        cluster = sub_partitions[node]
        # print(f"Assigning cluster {cluster} to node {node}")
        if node in subgraph.nodes:
            subgraph.nodes[node]["main_cluster"] = str(cluster)
            # print(f"Node: {node}, main_cluster: {subgraph.nodes[node].get('main_cluster')}")
        else:
            print(f"Warning: Node {node} from sub_partitions is not in subgraph.nodes!")

        if stdev_values[cluster] != 0:
            z_score = (degrees_within_cluster[node] - average_values[cluster]) / stdev_values[cluster]
            z_scores.append(z_score)
        else:
            z_score = 0
            z_scores.append(z_score)
        subgraph.nodes[node]["z_score"] = z_score

    # debugging
    # z_score_values = [G.nodes[n]["z_score"] for n in G.nodes if "z_score" in G.nodes[n]]
    # print(f"Z-scores: {z_score_values}")

    # Classify as Intra-Hub or Inter-Hub

    # Participation Coefficient (PC) determines intra-hubs vs. inter-hubs:

    # Intra-hubs (PC < 0.5): Highly connected within their own cluster.
    # Inter-hubs (PC ≥ 0.5): Highly connected to multiple clusters.
    # print(f"Node: {node}, Z-score: {z_score}")

    if z_score >= 1.0:
        squared_ratios = []
        degree = subgraph.degree(node)
        neighbors = [m for m in subgraph.neighbors(node)]
        neighbor_clusters = {key: sub_partitions[key] for key in neighbors}

        neighbor_clusters_dict = defaultdict(list)
        for neighbor, partition in sorted(neighbor_clusters.items()):
            neighbor_clusters_dict[partition].append(neighbor)
        # Calculate the participation coefficient (PC)
        for cluster_no, neighbor_list in sorted(neighbor_clusters_dict.items()):
            squared_ratios.append((len(neighbor_list) / degree) ** 2)

        if degree == 0:
            pc = 0  # Avoid division by zero2
        else:
            # squared_ratios = [(len(neighbor_list) / degree) ** 2 for neighbor_list in
            #                  neighbor_clusters_dict.values()]
            pc = 1 - sum(squared_ratios)

        if pc < 0.5:
            subgraph.nodes[node]['hub_type'] = 'intra-hub'
            intra_hub_data.append((node, string_id, pc, z_score, cluster))
            print(f"intra hub identified in component {i}")
        elif pc >= 0.5:
            subgraph.nodes[node]['hub_type'] = 'inter-hub'
            connected_clusters = ', '.join(map(str, neighbor_clusters_dict.keys()))
            inter_hub_data.append((node, string_id, pc, z_score, cluster, connected_clusters))
            print(f"inter hub indentified in component {i}")

    elif z_score < 1.0:
        squared_ratios = []

        degree = subgraph.degree(node)
        neighbors = [m for m in subgraph.neighbors(node)]
        neighbor_clusters = {key: sub_partitions[key] for key in neighbors}

    # Store Hub Data in DataFrames ----------------------------------------------------------------------------

    # Saves the sorted intra-hubs and inter-hubs for further analysis.
    intra_hub_df = pd.DataFrame(intra_hub_data,
                                columns=['Protein', 'String id', 'Participation coefficient', 'Z_score', 'Cluster'])
    inter_hub_df = pd.DataFrame(inter_hub_data,
                                columns=['Protein', 'String id', 'Participation coefficient', 'Z_score', 'Cluster',
                                         'Connected clusters'])
    # print("inter:", inter_hub_df)
    # print("intra:", intra_hub_df)

    intra_hub_df = intra_hub_df.sort_values(by='Z_score', ascending=False)
    inter_hub_df = inter_hub_df.sort_values(by='Participation coefficient', ascending=False)

    # Debug code block
    # print("Checking for 'stringdb::score' before cleaning...")
    for node, data in G.nodes(data=True):
        if "stringdb::score" in data:
            print(f"Found in Node {node}: {data}")

    for u, v, data in G.edges(data=True):
        if "stringdb::score" in data:
            print(f"Found in Edge ({u}, {v}): {data}")

    for key in G.graph.keys():
        if "stringdb::score" in key:
            # print(f"Found in Graph Attributes: {key} → {G.graph[key]}")
            continue


    # print("Checkpoint2")
    # Save Processed Network
    # Stores the processed subgraph as a GML file for future analysis.

    # Apply the cleaning functions before saving
    if len(subgraph.nodes()) > 0:

        print(f"Number of nodes in subgraph: {len(subgraph.nodes())}")

        try:
            subgraph1 = clean_node_attributes(subgraph)
            subgraph1 = clean_edge_attributes(subgraph)
            subgraph1 = clean_graph_attributes(subgraph1)
        except Exception as e:
            print("Caught an exception:", e)
            traceback.print_exc()
    else:
        print(f"Subgraph {i} has no nodes, skipping cleaning.")

    # Save Results
    # NOTE: SAVE POINT
    hub_file_path = f"C:/Users/BhagyaWijeratne/Desktop/Level IV/Research/methodolgy/PPI network/results/hub_analysis/{module_color}/hub_analysis_{module_color}_{i}.xlsx"

    # Extract the directory path
    hub_parent_dir = os.path.dirname(hub_file_path)

    # Ensure the directory exists
    os.makedirs(hub_parent_dir, exist_ok=True)

    with pd.ExcelWriter(hub_file_path, engine="openpyxl", mode="w") as w:

        df.to_excel(w, sheet_name="clustered proteins", index=False)
        intra_hub_df.to_excel(w, sheet_name="intra-hubs", index=False)
        inter_hub_df.to_excel(w, sheet_name="inter-hubs", index=False)

    # Save the graph as GML
    # NOTE: SAVE POINT

    # for node in subgraph.nodes:
    #     print(f"Node: {node}, main_cluster: {subgraph.nodes[node].get('main_cluster')}")

    # Define file path
    subgraph_file_path = f"C:/Users/BhagyaWijeratne/Desktop/Level IV/Research/methodolgy/PPI network/results/processed_subgraphs/{module_color}/{module_color}_module_processed_subgraph_{i}.gml"

    # Extract the parent directory
    subgraph_parent_dir = os.path.dirname(subgraph_file_path)

    # Ensure the directory exists
    os.makedirs(subgraph_parent_dir, exist_ok=True)
    # print("Subgraph before writing:", list(subgraph1.nodes(data=True)))
    nx.write_gml(subgraph1, subgraph_file_path)
    print(f"GML file for {module_color}, component {i} saved successfully.")
