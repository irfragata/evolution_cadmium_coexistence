Analyses folder
The summary section details general information about the files/folders in the Analyses folder. The section Files explains in detail the different files in this folder.
#### Summary ######

This folder has a series of intermediate files that are needed to run the code or that are produced as intermediate steps to then do figures.

Bootstrap_stats_tests.RData - workspace data saved from running the bootstrap for the statistical tests. This data is used in the Main_analyses script, if the evaluation parameter is set to FALSE.

cxr_normal_Pooled	- folder containing the parameter estimates produced when running the cxr for the pooled data. This data are needed to run the main analyses and for figures 1, S3 to S8.

cxr_normal_Replicates - folder containing the parameter estimates produced when running the cxr for each replicate. This data are needed to run the main analyses and for figures S3 to S8. 

MethodComparison - within this folder is stored the data from running the method comparison. This is to ensure that the script Method_Comparison.R is able to run.

min_distance_per_replicate.csv - Intermediate file produce when estimating the minimal distances to the edges for each replicate. This data file is needed to do figures 3 and S9.

min_distance_pooled.csv	- intermediate file produce when estimating the minimal distances to the edges for the pooled data. This data file is needed to do figures 3 and S9.

Parameter_exploration.RData - intermediate file with the information about the likelihood and estimates when starting the cxr from different initial conditions.

Parameter_exploration_RMSE.RData - intermediate file with the information about the root mean squared error estimates when starting the cxr from different initial conditions

popDyn_pooledData.csv - intermediate file produced when running the predicting section of the main analyses. This corresponds to a data frame with the predicted and observed data from the population dynamics experiment, for the pooled data. This file is needed for figure 4.

popDyn_replicates.csv - intermediate file produced when running the predicting section of the main analyses. This corresponds to a data frame with the predicted and observed data from the population dynamics experiment, for the replicate data. This file is needed for figure 4.

structural_Pooled.csv - intermediate file produced when estimating the structural niche and fitness differences for the different selection regimes, using pooled data. This file is used for figure 2.

structural_Replicates.csv - intermediate file produced when estimating the structural niche and fitness differences for the different selection regimes, using data for each replicate.

####### FILES #######
min_distance_per_replicate.csv
Columns in the file
Tu_Regime: identification of the Tu regime, SR1 - evolved in no cadmium, SR2 - evolved in cadmium
Te_Regime:identification of the Te regime, SR4 - evolved in no cadmium, SR5 - evolved in cadmium
Replicate: Experimental evolution replicate
Environment: Environment in which estimates were obtained
Tu_distance_lower, Tu_distance, Tu_distance_upper: lower, mean and upper distance between the growth rate vector and the edge of the feasibility domain for Tu
Te_distance_lower, Tu_distance, Tu_distance_upper: lower, mean and upper distance between the growth rate vector and the edge of the feasibility domain for Te
minDistance_L, minDistance, minDistance_U - lower and upper boundary for the minimum distance between the vector of growth rates and the edge of the feasibility domain.

min_distance_pooled.csv - 
Columns in the file
Tu_Regime: identification of the Tu regime, SR1 - evolved in no cadmium, SR2 - evolved in cadmium
Te_Regime:identification of the Te regime, SR4 - evolved in no cadmium, SR5 - evolved in cadmium
Environment: Environment in which estimates were obtained
Tu_distance_lower, Tu_distance, Tu_distance_upper: lower, mean and upper distance between the growth rate vector and the edge of the feasibility domain for Tu
Te_distance_lower, Tu_distance, Tu_distance_upper: lower, mean and upper distance between the growth rate vector and the edge of the feasibility domain for Te
minDistance_L, minDistance, minDistance_U - lower and upper boundary for the minimum distance between the vector of growth rates and the edge of the feasibility domain.


structural_Pooled.csv
Columns in the file:
u_Regime: identification of the Tu regime, SR1 - evolved in no cadmium, SR2 - evolved in cadmium
Te_Regime:identification of the Te regime, SR4 - evolved in no cadmium, SR5 - evolved in cadmium
Environment: Environment in which estimates were obtained
Tu_lambda, Te_lambda: estimated growth rate for Tu and Te
Tu_intra, Te_intra: estimated intraspecific competition for Tu and Te
Tu_inter, Te_inter: estimated interspecific competition for Tu and Te
ND: estimated structural niche differences
FD: estimated structural fitness differences
Feasibility: predicted feasibility of the system
ND_L, ND_U: lower and upper boundaries for estimated structural niche differences
FD_L, FD_U: lower and upper boundaries forestimated structural fitness differences
Feasibility_L, Feasibility_U: lower and upper boundaries for predicted feasibility of the system
a21_a11, a22_a12: slopes of the vectors defining the feasibility cone
a21_a11_lower, a22_a12_lower, a21_a11_upper, a22_a12_upper: lower and upper boundaries of the slopes of the vectors defining the feasibility cone
Tu_lambda_lower, Te_lambda_lower, Tu_lambda_upper, Te_lambda_upper: lower and upper boundaries for the estimated growth rate for Tu and Te
distanceTu_lower, distanceTu, distanceTu_upper: lower, mean and upper distance between the growth rate vector and the edge of the feasibility domain for Tu
distanceTe_lower, distanceTe, distanceTe_upper: lower, mean and upper distance between the growth rate vector and the edge of the feasibility domain for Te
min_a21_a11, min_a22_a12, max_a21_a11, max_a22_a12: minimum and maximum slopes of the vectors defining the feasibility cone


structural_Replicates.csv
Columns in the file:
Tu_Regime: identification of the Tu regime, SR1 - evolved in no cadmium, SR2 - evolved in cadmium
Te_Regime:identification of the Te regime, SR4 - evolved in no cadmium, SR5 - evolved in cadmium
Replicate: Experimental evolution replicate
Environment: Environment in which estimates were obtained
Tu_lambda, Te_lambda: estimated growth rate for Tu and Te
Tu_intra, Te_intra: estimated intraspecific competition for Tu and Te
Tu_inter, Te_inter: estimated interspecific competition for Tu and Te
ND: estimated structural niche differences
FD: estimated structural fitness differences
Feasibility: predicted feasibility of the system
ND_L, ND_U: lower and upper boundaries for estimated structural niche differences
FD_L, FD_U: lower and upper boundaries forestimated structural fitness differences
Feasibility_L, Feasibility_U: lower and upper boundaries for predicted feasibility of the system
a21_a11, a22_a12: slopes of the vectors defining the feasibility cone
a21_a11_lower, a22_a12_lower, a21_a11_upper, a22_a12_upper: lower and upper boundaries of the slopes of the vectors defining the feasibility cone
Tu_lambda_lower, Te_lambda_lower, Tu_lambda_upper, Te_lambda_upper: lower and upper boundaries for the estimated growth rate for Tu and Te
distanceTu_lower, distanceTu, distanceTu_upper: lower, mean and upper distance between the growth rate vector and the edge of the feasibility domain for Tu
distanceTe_lower, distanceTe, distanceTe_upper: lower, mean and upper distance between the growth rate vector and the edge of the feasibility domain for Te
min_a21_a11, min_a22_a12, max_a21_a11, max_a22_a12: minimum and maximum slopes of the vectors defining the feasibility cone

popDyn_pooledData.csv
Tu_Regime: identification of the Tu regime, SR1 - evolved in no cadmium, SR2 - evolved in cadmium
Te_Regime:identification of the Te regime, SR4 - evolved in no cadmium, SR5 - evolved in cadmium
Environment: Environment in which estimates were obtained
obs_TeRatio: mean observed estimated ratio of Te in the population after two generations
SE_obs: standard error of the observed ratio
pred_T1, T1_L, T1_U: lower, mean and upper predicted ratio after one generation
meanTe, meanTu: mean number of Te or Tu after two generations
predTu2, predTe2: predicted number of Tu or Te females after two generations

popDyn_replicates.csv - 
Tu_Regime: identification of the Tu regime, SR1 - evolved in no cadmium, SR2 - evolved in cadmium
Te_Regime:identification of the Te regime, SR4 - evolved in no cadmium, SR5 - evolved in cadmium
Environment: Environment in which estimates were obtained
meanRatio, mean_pred: mean observed or predicted Te ratio (Te/Total)
sdRatio, sdPred: mean observed or predicted standard error of the Te ratio (Te/Total)
meanTe, meanTu: average number of females observed

#### FOLDERS #####

MethodComparison - Inside of this folder there are 5 folders, where the estimates of the data following 5 different methods are stored. For each folder the name indicates the type of method used: cxr_lambda_fixed_log, cxr_lambda_fixed_nested, cxr_normal, optim_lambda_fixed_nested, Optim_normal


	cxr_lambda_fixed_log - folder containing the estimates using the cxr package with the lambda fixed.

		parameters_cxr_lambda_fixed_lower.csv, parameters_cxr_lambda_fixed_upper.csv, parameters_cxr_lambda_fixed.csv - Files containing the estimates obtained the lower, upper and mean estimates. The files all have the same structure.

		Tu_Regime: identification of the Tu regime, SR1 - evolved in no cadmium, SR2 - evolved in cadmium
		Te_Regime:identification of the Te regime, SR4 - evolved in no cadmium, SR5 - evolved in cadmium
		Replicate: Experimental evolution replicate
		Environment: Environment in which estimates were obtained
		Tu_lambda, Te_lambda: estimated growth rate for Tu and Te
		Tu_intra, Te_intra: estimated intraspecific competition for Tu and Te
		Tu_inter, Te_inter: estimated interspecific competition for Tu and Te

	cxr_lambda_fixed_nested - folder containing the estimates using the cxr package with the nested approach.

		parameters_cxr_lambda_fixed_lower.csv, parameters_cxr_lambda_fixed_upper.csv, parameters_cxr_lambda_fixed.csv - Files containing the estimates obtained the lower, upper and mean estimates. The files all have the same structure.

		Tu_Regime: identification of the Tu regime, SR1 - evolved in no cadmium, SR2 - evolved in cadmium
		Te_Regime:identification of the Te regime, SR4 - evolved in no cadmium, SR5 - evolved in cadmium
		Replicate: Experimental evolution replicate
		Environment: Environment in which estimates were obtained
		Tu_lambda, Te_lambda: estimated growth rate for Tu and Te
		Tu_intra, Te_intra: estimated intraspecific competition for Tu and Te
		Tu_inter, Te_inter: estimated interspecific competition for Tu and Te

	cxr_normal - folder containing the estimates using the cxr package using the cxr approach.
	
		parameters_cxr_normal_lower.csv, parameters_cxr_normal_upper.csv, parameters_cxr_normal.csv - Files containing the estimates obtained the lower, upper and mean estimates. The files all have the same structure.

		Tu_Regime: identification of the Tu regime, SR1 - evolved in no cadmium, SR2 - evolved in cadmium
		Te_Regime:identification of the Te regime, SR4 - evolved in no cadmium, SR5 - evolved in cadmium
		Replicate: Experimental evolution replicate
		Environment: Environment in which estimates were obtained
		Tu_lambda, Te_lambda: estimated growth rate for Tu and Te
		Tu_intra, Te_intra: estimated intraspecific competition for Tu and Te
		Tu_inter, Te_inter: estimated interspecific competition for Tu and Te

	optim_lambda_fixed_nested

		alpha_estimates_row_is_target_1Cd.csv, alpha_estimates_row_is_target_1N.csv, alpha_estimates_row_is_target_2Cd.csv, alpha_estimates_row_is_target_2N.csv, alpha_estimates_row_is_target_3Cd.csv, alpha_estimates_row_is_target_3N.csv, alpha_estimates_row_is_target_4Cd.csv, alpha_estimates_row_is_target_4N.csv, alpha_estimates_row_is_target_5Cd.csv, alpha_estimates_row_is_target_5N.csv - Files with the estimates obtained using the nested approach with intrinsic growth rate estimated directly from  the data.

		alpha_lower_errors_1Cd.csv, alpha_lower_errors_1N.csv, alpha_lower_errors_2Cd.csv, alpha_lower_errors_2N.csv, alpha_lower_errors_3Cd.csv, alpha_lower_errors_3N.csv, alpha_lower_errors_4Cd.csv, alpha_lower_errors_4N.csv, alpha_lower_errors_5Cd.csv, alpha_lower_errors_5N.csv, alpha_upper_errors_1Cd.csv, alpha_upper_errors_1N.csv, alpha_upper_errors_2Cd.csv, alpha_upper_errors_2N.csv, alpha_upper_errors_3Cd.csv, alpha_upper_errors_3N.csv, alpha_upper_errors_4Cd.csv, alpha_upper_errors_4N.csv, alpha_upper_errors_5Cd.csv, alpha_upper_errors_5N.csv - Files with the lower and upper estimates obtained using the nested approach with intrinsic growth rate estimated directly from  the data.

		All the files have the same structure
		 Intra - column with the estimates of intraspecific competition
		 Inter1 - column with the estimates of interspecific competition with the first selection regime
		 Inter2 - column with the estimates of interspecific competition with the second selection regime

		lambda_estimates_1Cd.csv, lambda_estimates_1N.csv, lambda_estimates_2Cd.csv, lambda_estimates_2N.csv, lambda_estimates_3Cd.csv, lambda_estimates_3N.csv, lambda_estimates_4Cd.csv, lambda_estimates_4N.csv, lambda_estimates_5Cd.csv, lambda_estimates_5N.csv - files with the growth rate estimates obtained directly from the data.

		species - column with species name
		lambda - estimates of intrinsic growth rate
		lower_error - lower boundary of estimates
		upper_error - upper boundary of estimates
		sigma - error associated with parameter estimation
		convergence - likelihood of the model

	Optim_normal

		alpha_estimates_row_is_target_1Cd.csv, alpha_estimates_row_is_target_1N.csv, alpha_estimates_row_is_target_2Cd.csv, alpha_estimates_row_is_target_2N.csv, alpha_estimates_row_is_target_3Cd.csv, alpha_estimates_row_is_target_3N.csv, alpha_estimates_row_is_target_4Cd.csv, alpha_estimates_row_is_target_4N.csv, alpha_estimates_row_is_target_5Cd.csv, alpha_estimates_row_is_target_5N.csv - Files with the estimates obtained using the nested approach.

		alpha_lower_errors_1Cd.csv, alpha_lower_errors_1N.csv, alpha_lower_errors_2Cd.csv, alpha_lower_errors_2N.csv, alpha_lower_errors_3Cd.csv, alpha_lower_errors_3N.csv, alpha_lower_errors_4Cd.csv, alpha_lower_errors_4N.csv, alpha_lower_errors_5Cd.csv, alpha_lower_errors_5N.csv, alpha_upper_errors_1Cd.csv, alpha_upper_errors_1N.csv, alpha_upper_errors_2Cd.csv, alpha_upper_errors_2N.csv, alpha_upper_errors_3Cd.csv, alpha_upper_errors_3N.csv, alpha_upper_errors_4Cd.csv, alpha_upper_errors_4N.csv, alpha_upper_errors_5Cd.csv, alpha_upper_errors_5N.csv - Files with the lower and upper estimates using the nested approach.

		All the files have the same structure
		 Intra - column with the estimates of intraspecific competition
		 Inter1 - column with the estimates of interspecific competition with the first selection regime
		 Inter2 - column with the estimates of interspecific competition with the second selection regime

		lambda_estimates_1Cd.csv, lambda_estimates_1N.csv, lambda_estimates_2Cd.csv, lambda_estimates_2N.csv, lambda_estimates_3Cd.csv, lambda_estimates_3N.csv, lambda_estimates_4Cd.csv, lambda_estimates_4N.csv, lambda_estimates_5Cd.csv, lambda_estimates_5N.csv - files with the growth rate estimates obtained directly from the data.

		species - column with species name
		lambda - estimates of intrinsic growth rate
		lower_error - lower boundary of estimates
		upper_error - upper boundary of estimates
		sigma - error associated with parameter estimation
		likelihood - likelihood of the model

cxr_normal_Pooled - Inside of this folder there 2 Rdata files that store the data after running the Running_cxr_pooledData.R script and 3 csv files were the upper, mean and lower parameter estimates are stored.

	cxr_Cd_allEqual_best.RData and cxr_N_allEqual_best.RData - intermediate files that store the parameters estimates from the cxr

	parameters_cxr_normal_Pooled_lower.csv, parameters_cxr_normal_Pooled_upper.csv, parameters_cxr_normal_Pooled.csv - All files have the same structure and same columns. The estimates presented are the lower, upper or mean, respectively

	Tu_Regime: identification of the Tu regime, SR1 - evolved in no cadmium, SR2 - evolved in cadmium
	Te_Regime:identification of the Te regime, SR4 - evolved in no cadmium, SR5 - evolved in cadmium
	Replicate: Experimental evolution replicate
	Environment: Environment in which estimates were obtained
	Tu_lambda, Te_lambda: estimated growth rate for Tu and Te
	Tu_intra, Te_intra: estimated intraspecific competition for Tu and Te
	Tu_inter, Te_inter: estimated interspecific competition for Tu and Te

cxr_normal_Replicates - Inside of this folder there 17 Rdata files that store the data after running the Running_cxr_Replicates.R script and 9 csv files were the upper, mean and lower parameter estimates are stored for no cadmium and cadmium and environments (N or Cd, respectively) and one with everything.

	cxr_R1_Cd.RData, cxr_R1_N.RData, cxr_R2_SR1_Cd.RData, cxr_R2_SR1_N.RData, cxr_R2_SR4_Cd.RData, cxr_R2_SR4_inter_Cd.RData, cxr_R2_SR4_inter_N.RData, cxr_R2_SR4_N.RData, cxr_R2_SR5_Cd.RData, cxr_R2_SR5_inter_Cd.RData, cxr_R2_SR5_inter_N.RData, cxr_R2_SR5_N.RData,  cxr_R3_Cd.RData, cxr_R3_N.RData, cxr_R4_Cd.RData, cxr_R4_N.RData, cxr_R5_Cd.RData, cxr_R5_N.RData - intermediate files that store the estimates obtained in cxr for each of the replicates and environments.


	parameters_cxr_normal_Cd.csv, parameters_cxr_normal_lower_Cd.csv, parameters_cxr_normal_upper_Cd.csv - Files storing the mean, lower and upper estimates for the cadmium environment
	parameters_cxr_normal_N.csv, parameters_cxr_normal_lower_N.csv, parameters_cxr_normal_upper_N.csv - Files storing the mean, lower and upper estimates for the no cadmium environment
	parameters_cxr_normal.csv, parameters_cxr_normal_lower.csv, parameters_cxr_normal_upper.csv  - Files storing the mean, lower and upper estimates for both environments

	Each file has the same structure and type of columns

	Tu_Regime: identification of the Tu regime, SR1 - evolved in no cadmium, SR2 - evolved in cadmium
	Te_Regime:identification of the Te regime, SR4 - evolved in no cadmium, SR5 - evolved in cadmium
	Replicate: Experimental evolution replicate
	Environment: Environment in which estimates were obtained
	Tu_lambda, Te_lambda: estimated growth rate for Tu and Te
	Tu_intra, Te_intra: estimated intraspecific competition for Tu and Te
	Tu_inter, Te_inter: estimated interspecific competition for Tu and Te







