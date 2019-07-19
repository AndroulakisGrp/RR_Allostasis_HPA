RR_2019_Allostasis_HPA.m is the model file describing the equations for the HPA Axis model used to simulate results published in Rao et al. Allostatic adaptation and personalized physiological trade-offs in the circadian regulation of the HPA axis: A mathematical modeling approach Scientific Reports (accepted July 2019)

RR_2019_Allostasis_Minimal_Results.m simulates a minimal set of results used in the manuscript:
	entrained corticosterone profiles for a randomly chosen point on the nominal 		parametric surface
	the nominal, intermediate stress, and high stress parametric surfaces with the 	9 	representative points chosen from each of the respective surfaces overlaid

RR_Allostasis_minimal_workspace.mat consists of a minimal set of variables used to generate the results
	[variable name] -> [description]
	
	nominal_space -> representative sample of parameters used in the construction of 	the nominal parametric subspace 
	inter_space -> representative sample of parameters used in the construction of the 	intermediate stress parametric subspace
	high_space -> representative sample of parameters used in the construction of 		the high stress parametric subspace
	nominal_grid -> parameter sets on the nominal parametric subspace used in the 		construction of Arnold tongues
	inter_grid -> parameter sets on the intermediate stress parametric subspace used 	in the construction of Arnold tongues
	high_grid -> parameter sets on the high stress parametric subspace used in the 		construction of Arnold tongues
