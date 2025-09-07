# 16/04/2025
# Author - Bhagya Wijeratne

# Input: GML network file of the processed subgraph
# Output: GML network file of the processed subgraph with node attribute updated on annotation status for 'leadf development' function

# Identifies the nodes which are known for leaf development function

# ==================================================================================================

import os
import networkx as nx

# === CONFIGURATION ===
base_dir = r"C:\Users\BhagyaWijeratne\Desktop\Level IV\Research\methodolgy\PPI network\results\processed_subgraphs"          # Root directory with subgraphs
gene_list_file = r"C:\Users\BhagyaWijeratne\Desktop\Level IV\Research\methodolgy\PPI network\GO annotations downloading\geneID_all_annotations_updated.txt"           # Text file with gene IDs (one per line)
node_attribute = 'query term'             # Node attribute to match gene IDs

# === LOAD GENE IDS INTO A SET ===
with open(gene_list_file, 'r') as f:
    gene_ids = set(line.strip() for line in f if line.strip())

# === WALK THROUGH DIRECTORY AND PROCESS GML FILES ===
for root, dirs, files in os.walk(base_dir):
    for file in files:
        if file.endswith('.gml'):
            gml_path = os.path.join(root, file)
            print(f"Processing {gml_path}...")

            # Load the graph
            G = nx.read_gml(gml_path)

            # Annotate nodes
            for node in G.nodes():
                node_val = str(G.nodes[node].get(node_attribute, '')).strip()
                G.nodes[node]['annotated'] = 1 if node_val in gene_ids else 0

            # Overwrite original GML
            nx.write_gml(G, gml_path)

print("Annotation completed for all graphs.")
