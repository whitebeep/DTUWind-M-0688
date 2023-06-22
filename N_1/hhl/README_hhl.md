# hhl
This folder has three scripts. 
- hhl_ac_pf.py
- hhl_dc_pf.py
- IBM_credentials.py
These are described next.
  
## IBM_credentials.py
This script implements two functions.
- load_account: given a token, saves and loads its IBMQ account
- get_real_backend: retrieves a device from IBMQ based on the specified filters

## hhl_dc_pf.py
This script has the functions that perform a DC-PF using the HHL algorithm.
- matrix_vector: prepares the A matrix and b vector for the HHL algorithm
- solve_hhl_sv: given A and b, solves the HHL algorithm via state vector simulation

## hhl_ac_pf.py
This script has the functions that perform an AC-FDLF using the HHL algorithm.
- get_delta: calculates the scaling factor to represent exactly lambda_min on nl binary bits
- construct_circuit: given A and b, constructs the quantum circuit for the HHL algorithm
- update_circuit: Update the input circuit (b) and combine it with the constant part (A) of the HHL circuit
- solve_hhl: function to solve the HHL algorithm. It uses the _update_circuit_ function.
- HHL_FDLF: function to solve the AC-FDLF using the HHL algorithm. It uses _construct_circuit_ for the first iteration and then _solve_hhl_
