This represitory contains the code and data for arxiv:2505.24146, which will be published in PRL (see https://doi.org/10.1103/rwd7-92z9).

The "code" folder contains julia codes for exact diagonalalization calculations and matlab codes for treating structure factor S(q), pair correlation functions G(r), and emergent pseudospin quantum numbers of the SU(2) Hall ferromagnet phase.

To use the julia codes, one need to set correct sublattice potential "AB_potential", twist angle "theta_", dielectric constant "epsilon", corrugation "waa", and finite clusters geometries defined by (V1,V2) in the table S1 of the Supplemental Materials. 

"Helical_Trilayer_ED.jl" calculates the many-body energy spectrum;
"Helical_Trilayer_PairCorrel_Sq.jl" calculates the structure factor and pair correlation functions;
"Helical_Trilayer_Chern.jl" calculates the many-body Chern number;
"Helical_Trilayer_PES.jl" calculates the particle-cut entanglement spectrum for FCI phase;
"Helical_Trilayer_correlation_matrix.jl" calculates the momentum space correlation matrix used for extracting the emergent pseudospins and SU(2) quantum numbers.


The "figs" folder contains matlab codes to plot the figures in main text and Supplemental Materials.
