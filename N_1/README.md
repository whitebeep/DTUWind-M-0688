# Main Functions
The main functions are explained next:

## N_1.py
This script gathers and uses the functions in _power_systems.py_, _power_flow.py_, _postprocessing.py_, _aux_func.py_ and _circ_func.py_ to perform a DC-PF and an N-1 contingency analysis via the classical method and the proposed quantum approach. The user can control several options (export results, print variables, etc.)

First, the network data is processed, and the DC-PF+LODF is performed by the classical method. Then the code follows with the quantum implementation. Upon building the final QC, the subsequent code cells encode the input, PTDF, LODF, etc., as a QC. For each of them, a QC is built,  using state vector simulation and post-processed.

## HHL.py
This script gathers and uses the functions in _power_systems.py_, _power_flow.py_, _postprocessing.py_, _hhl_dc_pf.py_, _hhl_ac_pf.py_ and _IBM_credentials.py_ to perform a DC-PF and an AC-FDLF via the classical method and a quantum approach based on the HHL algorithm. The user can control several options (export results, enable DC-PF, enable AC-FDLF, enable HHL algorithm, etc.)

First, the network data is processed. Then,  the DC-PF and the AC-FDLF are solved by the classical method. After that, these are also solved by the quantum approach using the HHL algorithm. Finally, the results are post-processed.
