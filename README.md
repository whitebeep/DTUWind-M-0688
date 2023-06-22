# DTUWind-M-0688
This repository contains the code developed during the master's thesis titled _**Quantum Computations for N-1 Secure Power Systems (DTU Wind-M-0688)**_ at the Technical University of Denmark  

Eduard Antolí-Gil  
Friday, June 23rd 2023  

## Code Structure  

The main folder (_N_1_) is structured as follows:  
- N_1.py
- HHL.py
- _LODF/_
  - LODF.py
  - PTDF.py
  - circ_func.py
  - aux_func.py
  - figures_to_export/
  - draw/
  - hist/
  - sv/
- _hhl/_
  - hhl_ac_pf.py
  - hhl_dc_pf.py
  - IBM_credentials.py   
- _other/_
  - qiskit_code_example.py 
- _postpro/_
  - postprocessing.py
  - tab/
- _power/_
  - power_flow.py
  - power_systems.py

Each subfolder has a dedicated README file, except those only storing figures, tables, etc. 

Last of all, the scripts _N_1.py_ and _HHL.py_ are introduced.

## N_1.py
This script gathers and uses the functions in _power_systems.py_, _power_flow.py_, _postprocessing.py_, _aux_func.py_ and _circ_func.py_ to perform a DC-PF and an N-1 contingency analysis via the classical method and the proposed quantum approach. The user can control several options (export results, print variables, etc.)

First, the network data is processed, and the DC-PF+LODF is performed by the classical method. Then the code follows with the quantum implementation. Upon building the final QC, the subsequent code cells encode the input, PTDF, LODF, etc., as a QC. For each of them, a QC is built,  using state vector simulation and post-processed.

## HHL.py
This script gathers and uses the functions in _power_systems.py_, _power_flow.py_, _postprocessing.py_, _hhl_dc_pf.py_, _hhl_ac_pf.py_ and _IBM_credentials.py_ to perform a DC-PF and an AC-FDLF via the classical method and a quantum approach based on the HHL algorithm. The user can control several options (export results, enable DC-PF, enable AC-FDLF, enable HHL algorithm, etc.)

First, the network data is processed. Then,  the DC-PF and the AC-FDLF are solved by the classical method. After that, these are also solved by the quantum approach using the HHL algorithm. Finally, the results are post-processed.