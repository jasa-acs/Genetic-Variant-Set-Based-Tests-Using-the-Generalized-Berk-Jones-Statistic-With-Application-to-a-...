# Adapted from gBJ_power_sim_final.R

# The purpose of this script is to run the simulations which produced Supplementary Tables 3 and 4
# in the submitted manuscript.
# We are interested in knowing the lowest number of causal SNPs that allows the
# different set-based tests (GBJ, GHC, SKAT, minP, OMNI) to reach a power of 80%
# when the effect size is constant at 0.1 or 0.15.

# Here we generate the snp data according to specific parameters, simulate the outcome from 
# a given true model, and apply each of the tests to the simulated data.

# Each job should produce one file with a name like 'power_results_S1_aID1_beta15.txt'
# where the numbers after Snum and aID will vary with each job.
# The parameters for each run are controlled by the aID and Snum variables, which are parameters passed
# in on the command line. 
# Our job submission scripts (*.sbatch and *.bsub) control these parameters
# so that all necessary files are generated.

# After all runs from this script as well as gen_sim_data_supptab3.R have been completed, use 
# make_supptab34.R to reproduce Figures 2 and 3 of the submitted manuscript.


#####################################################################
# Change this location to the directory that holds the necessary R libraries
# (must be installed first).

libs <- '/n/home/user/Rlibrary/3.3.1'

###################################################################
# No manual changes need to be made after this point
###################################################################

library(mvtnorm, lib=libs)
library(SKAT, lib=libs)
library(bindata, lib=libs)
library(data.table, lib=libs)
library(GBJ, lib=libs)
########################################################################################


########################################################################################
# Set parameters for just this specific run.
# If finding the omnibus null distribution, change this value to true.
OMNI_RUN <- FALSE

# Read input arguments (passed by the job submission files on cluster, or if running on a local
# machine, need to be passed in manually).
args <- commandArgs(trailingOnly=TRUE)
aID <- as.numeric(args[1])
Snum <- as.numeric(args[2])
set.seed(03082017 + Snum*2000 + aID + 1000)

# Build a parameter matrix, all possible settings that we could investigate.
rho1_vec <- c(0, 0.3)
rho2_vec <- c(0, 0.3)
rho3_vec <- c(0, 0.3)
max_ncausal <- 10
d_vec <- c(100)
param_mat <- expand.grid(rho1_vec, rho2_vec, rho3_vec, d_vec)
colnames(param_mat) <- c('rho1', 'rho2', 'rho3', 'd')

# Set the parameters for this setting based on input arguments.
runs_per_setting <- 5
if (OMNI_RUN == TRUE) {
	ncausal <- aID
} else {
	ncausal <- ceiling(aID / runs_per_setting)
}
param_row <- Snum
sims_per_setting <- 100					# Runs_per_setting * sims_per_setting = total sims at each setting
rho1 <- param_mat$rho1[param_row]
rho2 <- param_mat$rho2[param_row]
rho3 <- param_mat$rho3[param_row]
d <- param_mat$d[param_row]
print(param_mat[param_row,])
cat('\n')
########################################################################################


########################################################################################
# Parameters that stay constant across all runs.
num_sub <- 2000
omni_cor_boots <- 100
probG <- 0.3
alpha_vec <- c(0)
HAPGEN_FLAG <- FALSE
out_type <- 'C'
link_function <- 'linear'

# Effect sizes
beta_causal <- 0.15
beta_vec <- rep(beta_causal, ncausal)
########################################################################################

	
######################################################################
# Define the expit function
expit <- function(x) {
 	out = exp(x)/(1 + exp(x))
    out[x > 100] = 1
    out
}
########################################################################################



########################################################################################
# Function to generate either binary or continuous outcome data.
gen_data <- function(HAPGEN_FLAG, leg_file=NULL, d, num_sub, probG=NULL, sigma_struct=NULL, 
					ncausal, alpha_vec, beta_vec, out_type)
{
		# If binary, sim 10 times what we need to account for case-control sampling
		# There is a difference between num_sim and num_sub!
		if (out_type == 'D') {
			num_sim <- num_sub * 10
		} else if (out_type == 'C') {
			num_sim <- num_sub
		}
		
		############################################
		# Simulate G
		if (HAPGEN_FLAG)
		{
			G_mat <- -1
			while (class(G_mat) != 'matrix') {
				G_mat <- generate_hapgen(n=num_sim, d=d, aID=aID, Snum=Snum, chr=1, leg_file=leg_file)
			}
			G_mat <- G_mat[,1:d]
		} else 
		{
			G_mat <- rmvbin(n=num_sim, margprob=rep(probG, d), sigma=sigma_struct) + 
				rmvbin(n=num_sim, margprob=rep(probG, d), sigma=sigma_struct)
		}	
						
		# Flip to minor alleles for skat (generally only for HAPGEN SNPs)
		to_flip <- which( apply(G_mat,2,mean)>1 )
		if ( length(to_flip)>0 ) {
			G_mat[,to_flip] <- 2 - G_mat[,to_flip]
		}

		# Do we need other covariates? Look at alpha_vec
		if (length(alpha_vec) > 1)
		{
			X_mat <- cbind( 1, rmvnorm(n=num_sim, mean=rep(0, alpha_vec)) )
		} else {
			X_mat <- as.matrix(rep(1, num_sim))
		}
				
		
		############################################
		# Simulate covariates and outcome. 
		# If binary outcome, then use case-control sampling.
		if (out_type == 'D')
		{
			sim_eta <- X_mat %*% alpha_vec
			if (ncausal > 0) {
				sim_eta <- sim_eta + G_mat[ ,1:ncausal] %*% beta_vec
			}
			sim_mu <- expit(sim_eta)
	
			# Keep making datasets until we have enough cases and controls
			check_CC <- 0
			times_checked <- 0
			while (check_CC < num_sub/2) {
				sim_Y <- rbinom(n=num_sim, size=1, prob=sim_mu)
				check_CC <- min(length(which(sim_Y == 1)), length(which(sim_Y == 0)))
				times_checked <- times_checked + 1
				if (times_checked > 10) {stop('Seems like your average mu is too small')}
			}
			
			# Now use only num_sub/2 of each
			case_ind <- which(sim_Y == 1)[1:(num_sub/2)]
			control_ind <- which(sim_Y == 0)[1:(num_sub/2)]
			X_mat <- rbind(X_mat[case_ind,], X_mat[control_ind,])
			G_mat <- rbind(G_mat[case_ind,], G_mat[control_ind,])
			sim_Y <- c(sim_Y[case_ind], sim_Y[control_ind])
			
			# Fitted values
			null_mod <- glm(sim_Y~X_mat - 1, family=binomial(link='logit'))
			mu_hat_0 <- null_mod$fitted.values
			W_vec <- mu_hat_0 * (1-mu_hat_0)
		} else if (out_type == 'C')
		{
			sim_mu <- X_mat %*% alpha_vec
			if (ncausal > 0) {
				sim_mu <- sim_mu + as.matrix(G_mat[ ,1:ncausal]) %*% beta_vec
			}
			sim_Y <- sim_mu + rnorm(num_sim)			# variance 1
			
			# Fitted values
			null_mod <- glm(sim_Y~X_mat - 1, family=gaussian)
			mu_hat_0 <- null_mod$fitted.values
			W_vec <- rep(summary(null_mod)$sigma^2, num_sim)
		}

		return (list(G_mat=G_mat, 
					X_mat=X_mat, 
					sim_Y=sim_Y,
					mu_hat_0=mu_hat_0,
					null_mod=null_mod))
}
########################################################################################


########################################################################################
# Get correlation between set-based test p-values for use with omnibus test.  
# Needs the correlation between the individual-level test statistics.
# Simulate G matrix using the same correlation structure with bindata package.
# If we can't simulate G for some reason, just return -1.
omnibus_correlations <- function(num_boots, out_type, G_mat, X_mat, fitted_values)
{
	# For each set of test statistics run all the tests
	T_mat <- matrix(data=NA, nrow=num_boots, ncol=4)
	pval_vec <- rep(NA, 4)
	num_sub <- nrow(G_mat)
	bad_boots <- rep(NA, num_boots)
	
	# Bootstrap the test p-values under the null
	for (temp_boot in 1:num_boots)
	{
		# Checkpoint
		cat('omni: ', temp_boot, '\n')

		# Simulate a new outcome according to the null model
		if(out_type == 'D') {
			sim_Y <- rbinom(n=num_sub, size=1, prob=fitted_values)
			link_function <- 'logit'
			# New null model for GOF tests
			null_mod <- glm(sim_Y~X_mat - 1, family=binomial(link='logit'))
		} else if (out_type == 'C') {
			sim_Y <- fitted_values + rnorm(num_sub)
			link_function <- 'linear'
			null_mod <- glm(sim_Y~X_mat - 1, family=gaussian)
		}
		
		# Calculate test stats
		test_stat_output <- calc_score_stats(null_model=null_mod, factor_matrix=G_mat, 
							link_function=link_function)
		test_stats <- test_stat_output$test_stats
		
		
		# Run the 4 tests
		GBJ_omni <- GBJ(test_stats=test_stat_output$test_stats, cor_mat=test_stat_output$cor_mat)			
		GHC_omni <- GHC(test_stats=test_stat_output$test_stats, cor_mat=test_stat_output$cor_mat)	
		minP_omni <- minP(test_stats=test_stat_output$test_stats, cor_mat=test_stat_output$cor_mat)	
		skat_null_obj <- SKAT_Null_Model(sim_Y ~ X_mat - 1, out_type=out_type)
		skat_omni <- SKAT(Z=G_mat, obj=skat_null_obj, weights.beta=c(1,1), method='davies')					 										
		pval_vec[1] <- GBJ_omni$GBJ_pvalue
		pval_vec[2] <- GHC_omni$GHC_pvalue
		pval_vec[3] <- skat_omni$p.value
		pval_vec[4] <- minP_omni$minP_pvalue
		
		# Remember we want the correlation between the transformed T, not the actual p-values.
		# P-values of 1 get transformed to -Inf, so make these 0.999 first
		if (length(which(pval_vec == 1)) > 0)
		{
			pval_vec[which(pval_vec == 1)] <- 0.999
		}
		
		# Sometimes we get NA for minP p-value, skip and remove this later
		if (length(which(is.na(pval_vec))) > 0) {
			bad_boots[temp_boot] <- 1
			next
		} else {
			T_mat[temp_boot, ] <- qnorm(1-pval_vec)
		}
	}
	
	# Remove bad boots
	if (length(which(bad_boots==1)) > 0) {
		T_mat <- T_mat[-which(bad_boots==1),]
	}
	
	# Return the correlation between the transformed group testing statistics
	omni_cor_mat <- cor(T_mat)
	return(omni_cor_mat)
} 
########################################################################################




########################################################################################
# How are we generating the genotypes for this run?
if (HAPGEN_FLAG) {
		leg_file <- fread('hapmap3_r2_b36_chr1.legend')
		sigma_struct <- NULL
} else {
	# Make the simulation structure for new \rho_2
	cor_struct <- matrix(data=NA, nrow=d, ncol=d)
	cor_struct[1:ncausal, 1:ncausal] <- rho1
	cor_struct[1:ncausal, (ncausal+1):d] <- rho2
	cor_struct[(ncausal+1):d, 1:ncausal] <- rho2
	
	# Only 50% of the noncausal SNPs have an exchangeable correlation structure,
	# unless we have rho2 > 0 for positive definiteness
	if (rho2 > 0) {
		cor_struct[(ncausal+1):d, (ncausal+1):d] <- rho3
	} else {
		cor_struct[(ncausal+1):d, (ncausal+1):d] <- 0
	}
	num_noncausal <- d- ncausal
	half_noncausal <- floor(num_noncausal/2)
	cor_struct[(ncausal+half_noncausal+1):d, (ncausal+half_noncausal+1):d] <- rho3
	diag(cor_struct) <- 1
	
	cprob <- bincorr2commonprob( margprob=rep(probG,d), bincorr=cor_struct)
	sigma_struct <- commonprob2sigma(commonprob=cprob)
	leg_file <- NULL
}
########################################################################################



########################################################################################
# Loop through all simulations for this particular run
cat('ncausal: ', ncausal, '\n')
pow_results <- matrix(data=NA, nrow=(sims_per_setting), ncol=18)
pow_results <- data.frame(pow_results)
colnames(pow_results) <- c('ncausal', 'beta_causal',
							'gbj', 'GBJ_p', 'GBJ_err', 'ghc', 'GHC_p', 'GHC_err',
							'minP', 'minP_pvalue', 'skat_p', 'omni_stat',
							'BJ', 'BJ_p', 'rho1', 'rho2', 'rho3', 'd')
for (i in 1:sims_per_setting)
{				
	# Simulate covariates, SNPs, and outcome
	sim_data <- gen_data(HAPGEN_FLAG=HAPGEN_FLAG, leg_file=leg_file,
					d=d, num_sub=num_sub, probG=probG, 
					sigma_struct=sigma_struct, ncausal=ncausal, 
					alpha_vec=alpha_vec, beta_vec=beta_vec, out_type=out_type)
					
	sim_Y <- sim_data$sim_Y
	fitted_Y <- sim_data$null_mod$fitted.values
	G_mat <- sim_data$G_mat
	X_mat <- sim_data$X_mat
	null_mod <- sim_data$null_mod
		
	# Calculate score stats
	score_stats_output <- calc_score_stats(null_model=null_mod, factor_matrix=G_mat, 
			link_function=link_function)
		
	# Run all the tests
	BJ_output <- BJ(test_stats=score_stats_output$test_stats, cor_mat=score_stats_output$cor_mat)
	GBJ_output <- GBJ(test_stats=score_stats_output$test_stats, cor_mat=score_stats_output$cor_mat)
	GHC_output <- GHC(test_stats=score_stats_output$test_stats, cor_mat=score_stats_output$cor_mat)
	minP_output <- minP(test_stats=score_stats_output$test_stats, cor_mat=score_stats_output$cor_mat)
	skat_null_obj <- SKAT_Null_Model(formula=sim_Y ~ X_mat - 1, out_type=out_type)
	skat_output <- SKAT(Z=G_mat, obj= skat_null_obj, weights.beta=c(1,1), method='davies')
								
	#################################################################################
	# Run omnibus test.		
	# Only get the correlations the first time to save time.
	if (OMNI_RUN == TRUE)
	{
		# Get the correlations of the test statistics through bootstrapping
		omni_cor_mat <- omnibus_correlations(num_boots=omni_cor_boots,out_type=out_type,
						G_mat=G_mat, X_mat=X_mat, fitted_values=null_mod$fitted.values)
		
		# Just write it and stop
		colnames(omni_cor_mat) <- c('GBJ', 'GHC', 'SKAT', 'minP')
		omni_cor_name <- paste('d', d, '_rho1_', 10*rho1, '_rho2_', 10*rho2, '_rho3_', 
					10*rho3, '_ncausal_', ncausal, '.txt', sep='')
		write.table(omni_cor_mat, omni_cor_name, append=F, quote=F, row.names=F, col.names=T, sep='\t')
		stop()
	}
					
	# Calculate omnibus test p-value.
	# Doesn't actually matter that we have the correct order since all the bounds are the same.
	allp_vec <- c(GBJ_output$GBJ_pvalue, GHC_output$GHC_pvalue,
						skat_output$p.value, minP_output$minP_pvalue)
	omni_stat <- min(allp_vec)
		
	# Make bounds and calculate
	omni_bound <- tryCatch(qnorm(1-omni_stat), warning=function(w) w, error=function(e) e)
		
	# Record
	pow_results[i,1] <- ncausal
	pow_results[i,2] <- beta_causal
	pow_results[i,3:5] <- unlist(GBJ_output)
	pow_results[i,6:8] <- unlist(GHC_output)
	pow_results[i,9:10] <- unlist(minP_output)
	pow_results[i,11] <- skat_output$p.value
	pow_results[i,12] <- omni_stat
	pow_results[i,13] <- BJ_output$BJ
	pow_results[i,14] <- BJ_output$BJ_pvalue
	pow_results[i,15:18] <- as.numeric(param_mat[param_row,])
	
	# Checkpointing
	cat('i=', i, '\n')
}
########################################################################################


########################################################################################
# Write down the results
out_name <- paste('power_results_S', Snum, '_', aID, 'beta15.txt', sep='')
write.table(pow_results, out_name, append=F, quote=F, col.names=T, row.names=F, sep='\t')



















