# postpro
This folder has one script and one subfolder. 
- postprocessing.py
- _tab/_

## postprocessing.py
This script has some functions to post-process the results from the classical and quantum implementation of the N-1 contingency analysis based on the DC-PF. These are:
- postpro_hhl_dc: computes the pre and post-contingency power flows of the quantum approach
- lodf_dc: computes the pre and post-contingency power flows of the classical approach
- branch_list: auxiliary function for _generate_latex_table_pfc_
- get_timestamp_str: auxiliary function for _generate_latex_table_pfc_
- generate_latex_table_pfc: generates a latex table from the results of the N-1 contingency analysis

## Subfolders
The _tab/_ subfolder stores the .tex tables yield by _generate_latex_table_pfc_
