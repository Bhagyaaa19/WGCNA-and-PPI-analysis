'''
# 02/01/2025
Author - Bhagya Wijeratne

Input- Raw datasets downloaded from planteome Database
Output- Modified dataset with uniprot ids written to an excel file

Retrieves the uniprot ids for ensemble ids in Planteome dataset
'''

import pandas as pd
import requests, json

taxon_id = 4577

# maize taxon id = 4577

planteome_file = "C:\\Users\\BhagyaWijeratne\\Desktop\\Level IV\\Research\\methodolgy\\PPI network\\GO annotations downloading\\seed_proteins_leaf_dev_updated.xlsx"
sheet_name = "Planteome Data"
# Load the Excel file
data = pd.read_excel(planteome_file, sheet_name=sheet_name)

updated_data = pd.DataFrame()  # Initialize an empty DataFrame to store the updated data

# Define the base URL for the UniProt API
url = "https://rest.uniprot.org/uniparc/search?query="

uniprot_id_list = []

for uni_id in data.iloc[:,0]:
    try:

        request = requests.get(
            f"{url}" + uni_id + "&size=1&&format=json")  # Send a request to the UniProt API to get UniParc cross-references for the given ID
        record = json.loads(request.text)  # Load the response JSON
        results = record["results"]  # Extract the results from the JSON response

        # Print the record to debug
        print(f"Processing ID: {id}")
        print(json.dumps(record, indent=4))

         # Check if "results" key exists and is not empty
        if "results" in record and record["results"]:

            cross_references = results[0]["uniParcCrossReferences",[]]  # Extract cross-references from the first result

            uniprot_id = []
            for ref in cross_references:
                # Check if the reference is from UniProtKB/TrEMBL or UniProtKB/Swiss-Prot and if it matches the taxon ID
                if (ref["database"] == "UniProtKB/TrEMBL" or ref[
                    'database'] == "UniProtKB/Swiss-Prot") and "organism" in ref:
                    if ref["organism"]["taxonId"] == taxon_id:
                        uniprot_id.append(ref["id"])
            print(uniprot_id)

            # Append the list of UniProt IDs to the main list
            uniprot_id_list.append(uniprot_id)

        else:
            print(f"No results for ID: {id}")
            # Append None if no results were found
            uniprot_id_list.append(None)

    except Exception as e:
        print(f"Error processing ID {id}: {e}")
        uniprot_id_list.append(None)

# Add the UniProt IDs as a new column in the DataFrame
data['uniprot_ids'] = uniprot_id_list

# Save the updated DataFrame to a new Excel file, including the updated worksheet
output_file_path = "C:\\Users\\BhagyaWijeratne\\Desktop\\Level IV\\Research\\methodolgy\\PPI network\\GO annotations downloading\\seed_proteins_leaf_dev_updated.xlsx"
with data.ExcelWriter(output_file_path, engine='openpyxl') as writer:
    for sheet in data.ExcelFile(planteome_file).sheet_names:
        if sheet == sheet_name:
            data.to_excel(writer, sheet_name=sheet, index=False)
        else:
            pd.read_excel(output_file_path, sheet_name=sheet).to_excel(writer, sheet_name=sheet, index=False)

print(f"Modified file saved with \"_updated_uniprotIDs\" suffix")
