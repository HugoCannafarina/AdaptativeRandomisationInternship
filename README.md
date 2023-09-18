# AdaptativeRandomisationInternship
The code used for my report called : "Statiscal properties of Minimal Sufficient Balance, Minimization and Stratified Permuted Blocks in a stratification context, a simulation study"

List of R (version 4.2.3) packages required :	“openxlsx”, “dplyr”, “ggplot2”, “MASS”,” corrplot”,” pracma”,” carat”,” VGAM”, “parallel”

List of R scripts :
-	beta.r : compute the required sample sizes for the simulations, as well as the coefficients necessary to simulate the outcome. The coefficients are stored in the file beta.rda.
-	gen_data.r : select the covariate associations to model in the NORTA method, and apply the method to generate the datasets for each scenario. 
-	res_cor_300.r / res_uncor_300.r / res_cor_820.r / res_uncor_820.r: apply the randomization algorithms to each dataset. Each script corresponds to one of the scenarios. Forking is used for parallel computing. Pierre L'Ecuyer's RngStreams is used to set a seed and ensure reproductibility. 
-	final_analysis_cor_300 / final_analysis_uncor_300 / final_analysis_cor_820/ final_analysis_uncor_820 : final analysis of the generated data for each scenario.

Moreover, three script store the functions that are applied in the other scripts :
cov_functions.r : 
-	R : returns the empirical rank correlation between two variables
-	table_cor : returns a rank correlation matrix
-	table_pval_rank : returns the rank_correlation test p-value for every pair of covariates
-	cor_point : returns the obtained correlation between two covariates after applying the NORTA method. fct1 and fct2 are functions (quantile for continuous variables, qbinom for continuous ones). param1 and param2 are list of parameters for these functions (the vector of values and na.rm = TRUE for continuous variables and the probability parameter for discrete ones). The mean rank correlation obtained by generating n_iter NORTA vectors of size n_rep is returned. 
-	dist_cor : returns the absolute difference between a desired correlation and the obtained correlation after specifying a correlation in NORTA.
-	find_cor_sim : applies the “optimize” function of R on dist_cor to find the optimal correlation input for a given desired correlation.
-	g_bern : returns the g function described in the report in the case of a continuous-binary covariate relationship.
-	f_integral : returns the function integrated to compute the previous quantity
-	g_discrete : returns the g function described in the report in the case of a binary-binary covariate relationship.
-	esp_sdra / esp_squared_sdra / esp_immunodep / esp_squared_immunodep / esp_quanti / esp_squared_quanti : when integrated, return the first and second moments of the distribution functions of X for ARDS, immunodeficience and any quantitative variable, as described in the report.
-	objective_f / objective_f_discrete : the objective functions on which to apply Brent’s algorithm, respectively for continuous-discrete relationships and discrete-discrete relationships.
-	find_cor : solves the NORTA problem numerically for all the covariate pairs specified in var_list.
-	noise_cov : adds Gaussian noise to a covariance matrix.
-	simul_data_cor : simulates data using NORTA for a given correlation matrix (use the matrix found with find_cor). Noise is added using noise_cov, and the number of decimals from the original dataset is kept.
-	simul_data : simulates uncorrelated data using bootstrap.
-	n_wilcox_pvals : counts the number of rejected wilcoxon signed rank test between the covariates of two datasets for a given threshold. 
-	n_decimals : counts the number of decimals for a given variable.
-	name_format : replaces variables name in a string vector.
-	is_positive_definite : returns a Boolean that indicates if a matrix is positive definite.

sim_functions.r : 
-	MSB : applies the unstratified MSB algorithm to a dataset. Every vote cast is stored in the dataset
-	MSB_strat : applies the stratified MSB algorithm
-	vote_f : gathers the vote for a covariate
-	statistics_quanti_f : returns the mean, length and variance of a continuous covariate
-	vote_f_quanti : computes the vote for a continuous covariate for given previous treatment allocations
-	vote_f_quali : computes the vote for a qualitative covariate with a low amount of categories for given previous treatment allocations
-	vote_f_center : computes the vote for a qualitative covariate with a large amount of categories for given previous treatment allocations
-	vote_result : applies the MSB’s voting procedure to a vector of votes 
-	rando_blocks : implements’ block randomization
-	rando_pocock_strat : implements minimization
-	sim_outcome : simulates the treatment outcome to a dataset. The randomly generated treatment effect is stored at the first row of the “beta_treatment” variable.
-	inv_logit : returns the inverse of the logit function.

desc_functions.r : 
-	pval_khi2 : returns the p-value of a chi-squared test between two variables.
-	pval_compute_student : computes the p-values of student’s t-tests between the two treatment groups for continuous covariates and the p-values of chi-squared tests between the two treatment groups for categorical covariates for every dataset in a list of datasets. 
-	pval_compute_wilcox : computes the p-values of Wilcoxon signed rank tests between the two treatment groups for continuous covariates and the p-values of chi-squared tests between the two treatment groups for categorical covariates for every dataset in a list of datasets.
-	pval_compute_subgroup_student : computes the p-values of stratified student’s t-tests between the two treatment groups for continuous covariates and the p-values of stratified chi-squared tests between the two treatment groups for categorical covariates for every dataset in a list of datasets.
-	pval_compute_subgroup_wilcox : computes the p-values of stratified Wilcoxon signed rank tests between the two treatment groups for continuous covariates and the p-values of stratified chi-squared tests between the two treatment groups for categorical covariates for every dataset in a list of datasets.
-	MSE : computes the Mean-Squared Error of a list of dataset
-	Power :  computes the power, details of the parameters can be found in the report. 
-	Adjusted power/unadjusted power : return the mean, lower bound and upper bound of the adjusted and unadjusted power to show treatment effect, calculated on a list of datasets.
-	Significant_betas : returns the proportion of datasets from a list of datasets where the treatment effect, the stratified treatment effect and the closed-testing procedures are significative.

