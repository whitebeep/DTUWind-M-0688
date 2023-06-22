# LODF
This folder has four scripts and four subfolders.
- LODF.py
- PTDF.py
- circ_func.py
- aux_func.py
- figures_to_export/
- draw/
- hist/
- sv/
These are explained below.
  
## LODF.py
This script uses the functions in _aux_func.py_ and _circ_func.py_ to build, simulate and post-process the quantum circuits (QCs) that implement:
- The PTDF matrix
- A row of the LODF matrix
- The PTDF matrix plus a row of the LODF matrix
- The LODF matrix
- The PFC matrix
- The PTDF, LODF and PFC matrix 
  
## PTDF.py
This script uses the functions in _aux_func.py_ and _circ_func.py_ to build, simulate and post-process the QC that implement:
- A loop to construct the PTDF matrix step-by-step to see the evolution of the state vector
- The PTDF matrix

## circ_func.py
This script uses functions in _aux_func.py_ and defines the following functions:
- construct_lcu: performs a LCU to a diagonal matrix
- vector_to_quantum_circuit: to encode an input into a QC via an isometry (amplitude encoding)
- arbitrary_matrix_to_quantum_circuit: to factorize and encode a rectangular matrix into a QC
- diagonal_matrix_to_quantum_circuit: to factorize and encode a diagonal matrix into a QC
- matrix_to_quantum_circuit: to factorize and encode any matrix into a QC
- compose_quantum_circuit: recursive function to compose several QC
- simulate_quantum_circuit: transpiles and runs a QC for a given back-end; and returns its state vector and quasi-probabilities

## aux_func.py
This script comprises several ancillary functions which are used in _LODF.py_, _PTDF.py_ and _circ_func.py_. These are:
- pad_to_square_matrix: pads a rectangular matrix into a square matrix. This function is used in _prepare_lodf_matrix_ and _prepare_pfc_matrix_
- prepare_lodf_matrix: pre-processes the LODF matrix which then will be decomposed into a QC
- prepare_pfc_matrix: pre-processes the PFC matrix which then will be decomposed into a QC
- postpro_statevector: post-processes the state vector and returns its real part
- postpro_counts: post-processes the quasi-probabilities to extract the absolute value of the associated amplitudes
- pfc_from_quantum_circuit: extracts the pre and post-contingency line power flows from a state vector
- pfc_from_quantum_circuit_old: extracts the pre-contingency line power flows and the increments of line power flow due to the contingency. Finally, it computes the post-contingency line power flows
- save_quantum_circuit: exports and draws the QC
- save_quasiprob: exports and plots an histogram of the QC
- save_statevector_as_latex_table: exports a latex table from a state vector
- chop: truncates small values of a complex array
- is_diagonal: checks if a matrix is diagonal
- ishermitian: checks if a matrix is hermitian
- issquare: checks if a matrix is square
- decompose_to_unitary: decomposes a hermitian matrix with norm=1 to four unitary matrices
- compute_control_matrix: adds a control qubit to a unitary matrix (from Qiskit's source code)
- pad_to_n_qubits: pad a matrix to a specified size with ones on the diagonal (default)

## Subfolders
- _figures_to_export/_ has the figures used for documentation
- _draw/_ stores the QC figures generated with _save_quantum_circuit_
- _hist/_ stores the quasi-probability histograms generated with _save_quasiprob_
- _sv/_ stores the .tex generated with _save_statevector_as_latex_table_