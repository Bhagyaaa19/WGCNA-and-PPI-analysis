# 27/12/2024
# Author - Bhagya WIjeratne

# Input: The excel file containing the planteome IDs with unnecessary suffix "_T001"
# Output: The excel file containing an additional column where the unnecessary suffix is removed
#- such that these IDs can be searched against NCBI batch entrez for the gene Ids

# Removes the unnecessary suffix on the planteome ID

# ==========================================================================================

import pandas as pd

# Load the Excel file and specify the worksheet
file_path = "C:\\Users\\BhagyaWijeratne\\Desktop\\Level IV\\Research\\methodolgy\\PPI network\\GO annotations downloading\\seed_proteins_leaf_dev.xlsx"
sheet_name = 'Planteome Data'

# Read the specified sheet
df = pd.read_excel(file_path, sheet_name=sheet_name)

# Remove the `_T001` part from the first column and create a new column
df['Modified'] = df.iloc[:, 0].str.replace('_T001', '', regex=False)

# Save the updated DataFrame to a new Excel file, including the updated worksheet
output_file_path = "C:\\Users\\BhagyaWijeratne\\Desktop\\Level IV\\Research\\methodolgy\\PPI network\\GO annotations downloading\\seed_proteins_leaf_dev_updated.xlsx"
with pd.ExcelWriter(output_file_path, engine='openpyxl') as writer:
    for sheet in pd.ExcelFile(file_path).sheet_names:
        if sheet == sheet_name:
            df.to_excel(writer, sheet_name=sheet, index=False)
        else:
            pd.read_excel(file_path, sheet_name=sheet).to_excel(writer, sheet_name=sheet, index=False)

print(f"Modified file saved with \"_updated\" suffix")
