# power
This folder has two scripts.
- power_systems.py
- power_flow.py

## power_systems.py
This script has all the necessary functions to initialize the network data for an AC power flow (Newton-Raphson, AC-NR), AC-FDLF and DC-PF. It contains main functions, auxiliary functions and network data
- _main_functions_
 - init_network_data: creates the dictionaries to store the network data
 - load_network_data: given some network data, initializes the dictionaries that store the network data
 - load_pf_data: computes the required admittance or susceptance matrices depending on the power flow problem (AC-NR, AC-FDLF or DC-PF)
 - load_contingency_data: computes the PTDF and LODF matrices
- _auxiliary_functions_
 - bus_indices: gets the PQ, PV and Slack indices of the network
 - make_Cft: calculates the line-bus incidence matrix
 - make_Bdc: calculates the susceptance matrix for the DC-PF problem
 - make_B: calculates the susceptance matrices for the AC-FDLF problem
 - make_Ybus: calculates the admittance matrices for the AC-NR problem
 - voltage_power_data: initializes the voltage and apparent power of the network
 - calculate_ptdf: calculates the PTDF matrix
 - calculate_lodf: calculates the LODF matrix
- Network data
 - slides_3bus: system used during the project
 - IEEE_5bus
 - IEEE_5bus_mp
 - IEEE_5bus_mod
 
## power_flow.py
This script has the functions to solve the AC-NR, AC-FDLF and DC-PF using the classical method.
- PowerFlowDC: solves the DC-PF
- PowerFlowFD: solves the AC-FDLF. Uses the bus admittance matrix and the susceptance matrices (B' and B'')
- PowerFlowFD_V2: solves the AC-FDLF. Derives the susceptance matrices from the bus admittance matrix
- PowerFlowNewton: solves the classic power flow problem using the Newton-Raphson method (AC-NR)
- calculate_F: calculates the apparent power mismatch between the calculated and specified values. Auxiliary function for AC-NR
- CheckTolerance: checks if the power mismatch is lower than the defined threshold. Auxiliary function for AC-NR
- generate_Derivatives0: generates the derivatives of the apparent power with respect to the voltage angle and magnitude. Auxiliary function for AC-NR
- generate_Derivatives: generates the derivatives of the apparent power with respect to the voltage angle and magnitude. Auxiliary function for AC-NR
- generate_Jacobian: given the previous derivatives, assembles the Jacobian matrix
- Update_Voltages: updates the solution vector after computing the increment using the Jacobian matrix
- ac_apparent_power: calculates the injected apparent power at the buses, the injected apparent power flowing from a bus, and the injected apparent power to a bus
- DisplayResults: prints the results of the AC-NR power flow. It includes voltage magnitude and angle, apparent power injected, etc.