# 15/12/2024
# Author - Bhagya WIjeratne

# Input: Seed proteins annotated for leaf development in Maize from Amigo, QuickGo and Planteome
# Output: the excel file containing all details without duplicate genes

# Creates an excel file of all the leaf development seed proteins combining data from different sources

import pandas as pd

# File paths
amigo_input_file = "C:\\Users\\BhagyaWijeratne\\Desktop\\Level IV\\Research\\methodolgy\\PPI network\\GO annotations downloading\\Amigo_leaf_development_annotations.txt"
quickGo_input_file = "C:\\Users\\BhagyaWijeratne\\Desktop\\Level IV\\Research\\methodolgy\\PPI network\\GO annotations downloading\\QuickGO-annotations-1733995125434-20241212.tsv"
planteome_input_file = "C:\\Users\\BhagyaWijeratne\\Desktop\\Level IV\\Research\\methodolgy\\PPI network\\GO annotations downloading\\Planteome_annotations.txt"
output_file = "C:\\Users\\BhagyaWijeratne\\Desktop\\Level IV\\Research\\methodolgy\\PPI network\\GO annotations downloading\\seed_proteins_leaf_dev.xlsx"

# Read the data
df_amigo = pd.read_csv(amigo_input_file, sep="\t", header=None)
# Read the QuickGO data
df_quickGo = pd.read_csv(quickGo_input_file, sep="\t", header=0)
# Read the Planteome data
df_planteome = pd.read_csv(planteome_input_file, sep="\t", header = None)

# Assign column names
df_amigo.columns = [
    "Gene/product (bioentity_label)", "GO class (direct) (annotation_class)",
    "GO class (direct) (annotation_class_label)", "Type", "Evidence (evidence_type)",
    "Reference", "Source", "Evidence with (evidence_with)",
    "Gene/product name (bioentity_name)", "Isoform (bioentity_isoform)", "Synonym"
]

# Convert all columns to upper case to avoid mismatches
df_amigo = df_amigo.applymap(lambda x: str(x).upper())
df_quickGo = df_quickGo.applymap(lambda x: str(x).upper())
df_planteome = df_planteome.applymap(lambda x: str(x).upper())

# Write to Excel with separate worksheets
with pd.ExcelWriter(output_file, engine="openpyxl") as writer:
    df_amigo.to_excel(writer, sheet_name="Amigo Data", index=False)
    df_quickGo.to_excel(writer, sheet_name="QuickGO Data", index=False)
    df_planteome.to_excel(writer, sheet_name="Planteome Data", index= False)

print(f"Data has been read from \n,"
      f"{amigo_input_file}\n,"
      f"{quickGo_input_file}\n,"
      f"{planteome_input_file}\n"
      f" and written to {output_file}")

#=================================================================================
