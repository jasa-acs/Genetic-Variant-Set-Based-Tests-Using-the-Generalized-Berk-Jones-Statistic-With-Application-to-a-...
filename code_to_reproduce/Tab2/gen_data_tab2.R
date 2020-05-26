# Adapted from cgem_rerun.R

# The purpose of this script is to perform the analysis of the CGEMS dataset that
# is given in Table 2 of the submitted manuscript.
# We are interested in finding gene-level p-values for association with breast cancer.
# We test each gene in the genome with GBJ, GHC, SKAT, minP, and OMNI.

# This code controls the extraction of data for each SNP as well as some data cleaning
# and application of all tests to each gene.

# Each job should produce one file with a name like 'cgems_results_S1_aID1.txt'
# where the numbers after each parameter will index the jobs we are running.
# Our job submission scripts (*.sbatch and *.bsub) control Snum and ID
# so that all necessary files are generated.

# Here the code refers to our attached fake dataset because we are not allowed to release the actual data.
# See the ACC Form for more details.
# If the fake dataset is replaced by the real dataset (and the corresponding files in the real dataset 
# are given the same names that we have given the fake data), then make_tab2.R will 
# reproduce Table 2 of the submitted manuscript (once all runs have finished).

# Obviously running this code on the fake dataset will not reproduce our results.
# In fact, running our code on the fake dataset will produce mostly NA for the results, because the fake
# dataset we have attached only contains genotypes at a small number of markers (in the interest of 
# not sending overly large files), so most genes will have no data.
# To see more numerical results, one can edit the glist-hg18.txt file to define more fake genes at the
# markers we have included (chromosome 10, between 121001183-121876520).


########################################################################################
# Change this location to the directory that holds the necessary R libraries
# (must be installed first).
libs <- '/n/home/user/Rlibrary/3.3.1'

###################################################################
# No manual changes need to be made after this point
###################################################################

library(mvtnorm, lib=libs)
library(data.table, lib=libs)
library(SKAT, lib=libs)
library(GBJ, lib=libs)
########################################################################################


########################################################################################
# Read in the input arguments (given automatically by the job submission scripts)
args <- commandArgs(trailingOnly=TRUE)
aID <- as.numeric(args[1])
Snum <- as.numeric(args[2])

# Two parameters - how many genes to test in every job?
# How much upstream and downstream buffer to put around the gene?
down_buffer <- 20000
up_buffer <- 20000
num_genes <- 100

# Read and order the covariates.
cgems_covar <- read.table('cgems_fake_covar.txt', header=T)
cgems_covar <- cgems_covar[order(cgems_covar[,1]), ]

# Read the gene list.
glist <- read.table('glist-hg18.txt', header=F)
colnames(glist) <- c('chr', 'start', 'end', 'gene')

# Order the gene list by chromosome and then base pair.
# Not necessary, but just for organizational purposes.
glist <- glist[order(glist[,1], glist[,2]),]

# Don't test the X/Y chromosomes - so just use the first 18011 after ordering
start_row <- num_genes * (aID-1) + 1
end_row <- num_genes * aID
if (end_row > 18011) {end_row=18011}
glist <- glist[start_row:end_row,]

# Fit the null model once first.
Y <- cgems_covar$AFFECTED
null_X <- as.matrix( cbind(1, cgems_covar$AGE_CAT0, cgems_covar$AGE_CAT2, cgems_covar$AGE_CAT3,
cgems_covar$AGE_CAT4, cgems_covar$AGE_CAT5, cgems_covar$EV1, cgems_covar$EV2, cgems_covar$EV3) )
null_mod <- glm(Y ~ null_X - 1, family=binomial(link='logit'))
null_mu <- null_mod$fitted.values

# W matrix and P matrix come from the null model and are the same for all genes.
W_vec <- null_mu * (1-null_mu)
W_mat <- diag(W_vec)
P_mat <- W_mat - W_mat %*% null_X %*% solve(t(null_X) %*% W_mat %*% null_X) %*% 
		t(null_X) %*% W_mat
########################################################################################


##################################################################################################
# After we have extracted the relevant set of genotypes using PLINK using 1/2 recoding, 
# read from disk, parse, and then erase from disk.
# Read in and order the .ped file to match covariates.
parse_snpset <- function(plink_fname)
{
	# Read in the ped file, order it to match covariates
	ped_fname <- paste( plink_fname, '.ped', sep='')
	ped_file <- read.table(ped_fname)
	ped_file <- ped_file[order(ped_file[,1]), ]
	
	# 1 or less SNPs in interval
	num_snps <- (ncol(ped_file) - 6)/2
	if (num_snps <= 1 | num_snps > 900) {
		system2(command='rm', args=paste(plink_fname, '.log', sep=''))
		system2(command='rm', args=paste(plink_fname, '.ped', sep=''))
		system2(command='rm', args=paste(plink_fname, '.map', sep=''))
		return(-1)
	}
	
	# Loop over the SNPs to fill the genotype matrix
	G_mat <- matrix(data=NA, nrow=nrow(ped_file), ncol=num_snps)
	for (j in 1:num_snps)
	{
		# First check for missingness, impute if needed
		left <- ped_file[, (5+2*j)]
		right <- ped_file[, (6+2*j)] 
		missing_left <- which(left<0 | left>2)
		missing_right <- which(right<0 | right>2)
		if (length(missing_left) > 0) {left[missing_left] <- NA}
		if (length(missing_right) > 0) {right[missing_right] <- NA}
		
		# Minor/major alleles coded as 1/2.
		temp_G <- left + right - 2
		
		# Impute if any missing.
		to_impute <- which(is.na(temp_G))
		if (length(to_impute) > 0) {
			temp_mean <- mean(temp_G, na.rm=TRUE) / 2
			temp_G[to_impute] <- rbinom(n=length(to_impute), size=2, prob=temp_mean)
		}
		
		# Record the SNP
		G_mat[,j] <- temp_G
	}
	
	# Flip everything to the minor allele
	freqs <- apply(G_mat, 2, mean) / 2
	to_flip <- which(freqs > 0.5)
	if (length(to_flip) > 0) {
		G_mat[,to_flip] <- 2 - G_mat[,to_flip]
	}
	
	# Remove MAF < 5% 
	freqs <- apply(G_mat, 2, mean) / 2
	MAF_remove <- which(freqs < 0.05)
	if (length(MAF_remove) > 0) {
		G_mat <- G_mat[,-MAF_remove]
	}
	
	# Check again if we have enough SNPs to still be a matrix
	num_snps <- num_snps - length(MAF_remove)
	if (num_snps < 2) {
		system2(command='rm', args=paste(plink_fname, '.log', sep=''))
		system2(command='rm', args=paste(plink_fname, '.map', sep=''))
		system2(command='rm', args=paste(plink_fname, '.ped', sep=''))
		return (-1)
	}
	
	# Do we have SNPs with perfect correlation?  If so, remove.
	# This threshold could possibly depend on n.
	duplicate_snps <- which( duplicated(t(G_mat)) )
	if (length(duplicate_snps) > 0) {
		G_mat <- G_mat[,-duplicate_snps]
	}	
	
	# Now check again for only 1 SNP
	num_snps <- num_snps - length(duplicate_snps)
	if (num_snps < 2) {
		system2(command='rm', args=paste(plink_fname, '.log', sep=''))
		system2(command='rm', args=paste(plink_fname, '.map', sep=''))
		system2(command='rm', args=paste(plink_fname, '.ped', sep=''))
		return (-1)
	}
	
	return(G_mat)
}
########################################################################################


########################################################################################
# Run all tests - HC, BJ, GHC, GBJ, SKAT, burden, minP.
# Starting from the input of just the G matrix, sig_vec, and test stats.
# (We could technically start from just the G_mat but it's easier with the test stats and sig_vec
# so we don't have to pass extra matrices in and in case we want to do bootstrapping with just one sig_vec)
run_all_tests <- function(outcome, G_mat, test_stats, sig_vec)
{
	# Sometimes the stats are too big for pchisq, we have to round it down to 8.2
	# 8.19 gives us some breathing room for the optimization procedures, which can go to 8.2
	Z_large_flag <- 0
	if( length(which(test_stats > 8.19))>0 ) {
		Z_large_flag <- 1
		test_stats[which(test_stats > 8.19)] = 8.19
	}
	
	########################
	# Calc  BJ and HC
	BJ_output <- BJ(test_stats=test_stats, pairwise_cors=sig_vec)
	HC_output <- HC(test_stats=test_stats, pairwise_cors=sig_vec)
	
	########################
	# Calculate BB_GBJ and GHC
	GBJ_output <- GBJ(test_stats=test_stats, pairwise_cors=sig_vec)
	GHC_output <- GHC(test_stats=test_stats, pairwise_cors=sig_vec)

	# minP test
	minP_output <- minP(test_stats=test_stats, pairwise_cors=sig_vec)
	
	# Lastly, SKAT
	# Default weights are bad for common variants
	# Fit the null object for SKAT
	skat_null_obj <- SKAT_Null_Model(outcome ~ null_X - 1, out_type="D")
	skat_p <- SKAT(G_mat, skat_null_obj, weights.beta=c(1,1), method='davies')$p.value

	# Return
	return( list(Z_large_flag=Z_large_flag, 
				b=BJ_output$BJ, BJ_p=BJ_output$BJ_pvalue,
				h=HC_output$HC, HC_p=HC_output$HC_pvalue,
				ghc=GHC_output$GHC, GHC_p=GHC_output$GHC_pvalue,
				gbj=GBJ_output$GBJ, gBJ_p=GBJ_output$GBJ_pvalue, 
				gBJ_pvalue_flag=GBJ_output$err_code,
				minP_pvalue=minP_output$minP_pvalue,
				skat_p=skat_p) )
}
########################################################################################


########################################################################################
# Get correlation between set-based test p-values for use with omnibus test.  
# Needs the correlation between the individual-level test statistics.
# Simulate G matrix using the same correlation structure with bindata package.
# If we can't simulate G for some reason, just return -1.
omnibus_correlations <- function(num_boots, null_model, G_mat, sig_vec)
{
	fitted_mu <- null_model$fitted.values

	# For each set of test statistics run all the tests
	T_mat <- matrix(data=NA, nrow=num_boots, ncol=4)
	pval_vec <- rep(NA, 4)
	for (temp_boot in 1:num_boots)
	{
		# Simulate a new outcome according to the null model
		sim_Y <- rbinom(n=length(fitted_mu), size=1, prob=fitted_mu)

		# Bootstrap null model uses new outcomes, careful!
		boot_null_mod <- glm(sim_Y ~ null_X - 1, family=binomial(link='logit'))
		
		# Get new test statistics
		sim_stats <- score_stats_only(null_model=boot_null_mod, factor_matrix=G_mat, 
							link_function='logit', P_mat=P_mat)
				 
		omni_sim_pvalues <- run_all_tests(outcome=sim_Y, G_mat=G_mat,
										test_stats=sim_stats,
										sig_vec=sig_vec)
										
		pval_vec[1] <- omni_sim_pvalues$gBJ_p
		pval_vec[2] <- omni_sim_pvalues$GHC_p
		pval_vec[3] <- omni_sim_pvalues$skat_p
		pval_vec[4] <- omni_sim_pvalues$minP_pvalue
		
		# Remember we want the correlation between the transformed T, not the actual p-values.
		# P-values of 1 get transformed to -Inf, so make these 0.999 first
		if (length(which(pval_vec == 1)) > 0)
		{
			pval_vec[which(pval_vec == 1)] <- 0.99
		}
		T_mat[temp_boot, ] <- qnorm(1-pval_vec)
	}
	
	# Return the correlation between the transformed group testing statistics
	omni_cor_mat <- cor(T_mat)
	omni_cor_vec <- omni_cor_mat[upper.tri(omni_cor_mat)]
	return(omni_cor_vec)
} 
########################################################################################



########################################################################################
# Loop through the gene list and do the testing.
cgems_results <- matrix(data=NA, nrow=num_genes, ncol=16)
for (i in 1:nrow(glist))
{
	start_bp <- as.numeric(glist$start[i]) - down_buffer
	end_bp <- as.numeric(glist$end[i]) + up_buffer
	
	# Some gene names are really weird.
	gene_name <- gsub('/', '.', glist$gene[i])
	
	# Plink to recode as 12, important to as.character(chr)!
	plink_flag <- tryCatch(expr=system2(command='./plink', args=c('--noweb', '--bfile', 'cgems_genotypes', 
									'--recode12', '--chr', as.character(glist$chr[i]), '--from-bp', 
									start_bp, '--to-bp', end_bp, 
									'--out', gene_name), wait=TRUE), 
								warning=function(w) w, error=function(e) e)
	
	# plink_flag is 0 for success and 1 for failure.
	# Avoid removing with '*'
	if (plink_flag != 0) {
		system2(command='rm', args=paste(gene_name, '.log', sep=''))
		system2(command='rm', args=paste(gene_name, '.map', sep=''))
		system2(command='rm', args=paste(gene_name, '.ped', sep=''))
		next
	}
	
	G_mat <- parse_snpset(plink_fname=gene_name)
	
	# If returned -1, not enough SNPs
	if (class(G_mat) != 'matrix') {next}
	# Make sure we have the same number subjects as in covariate table.
	if (nrow(G_mat) != nrow(cgems_covar)) { stop('Mismatch number of subjects') }

	# Loop over the SNPs to get marginal test statistics
	score_stats_output <- calc_score_stats(null_model=null_mod, factor_matrix=G_mat, 
											link_function='logit', P_mat=P_mat)
	Z_stats <- score_stats_output$test_stats
	cor_mat <- score_stats_output$cor_mat
	sig_vec <- cor_mat[upper.tri(cor_mat)]	
	
	# Get all the p-values from HC, BJ, GHC, GBJ, burden, minP, SKAT
	all_pvalues <- run_all_tests(outcome=Y,
								G_mat=G_mat, test_stats=Z_stats, 
								sig_vec=sig_vec)

	#################################################################################
	# Run omnibus test with minP.
	num_boots <- 100
	
	# Get the correlations of the test statistics through bootstrapping
	omni_cor_vec <- omnibus_correlations(num_boots=num_boots, null_model=null_mod, G_mat=G_mat,
										sig_vec=sig_vec)
										
	# It could be that we weren't able to generate SNPs from this correlation structure 
	# for some reason.
	if (length(omni_cor_vec) == 1 | length(which(is.na(omni_cor_vec))) > 0)
	{
		omni_stat <- NA
		omni_p <- NA
	} else 
	{
		# Make the omnibus sigma matrix for use with pmvnorm()
		omni_sig_mat <- matrix(data=1, nrow=4, ncol=4)
		omni_sig_mat[upper.tri(omni_sig_mat)] <- omni_cor_vec
		omni_sig_mat[lower.tri(omni_sig_mat)] <- t(omni_sig_mat)[lower.tri(omni_sig_mat)]
								
		# Calculate omnibus test p-value.
		# Doesn't actually matter that we have the correct order since all the bounds are the same.
		allp_vec <- c(all_pvalues$gBJ_p, all_pvalues$GHC_p,
						all_pvalues$skat_p, all_pvalues$minP_pvalue)
		omni_stat <- min(allp_vec)
		
		# Make bounds and calculate
		omni_Z_bounds <- rep(qnorm(1-omni_stat), length(allp_vec))
		omni_p <- 1 - pmvnorm(lower=rep(-Inf, length(allp_vec)), upper=omni_Z_bounds,
								sigma=omni_sig_mat)[1]
	}
	
	# Record
	cgems_results[i,1] <- gene_name
	cgems_results[i,2:13] <- unlist(all_pvalues)
	cgems_results[i,14] <- omni_stat
	cgems_results[i,15] <- omni_p
	cgems_results[i,16] <- length(Z_stats)

	# Remove the plink files
	system2(command='rm', args=paste(gene_name, '.log', sep=''))
	system2(command='rm', args=paste(gene_name, '.map', sep=''))
	system2(command='rm', args=paste(gene_name, '.ped', sep=''))
	
	# Checkpointing
	cat(i)
}
########################################################################################


########################################################################################
# Write down the results
out_name <- paste('cgems_results_S', Snum, '_', aID, '.txt', sep='')
colnames(cgems_results) <- c('Gene', 'Z_large_flag', 'b', 'BJ_p', 'h', 'HC_p',
							'ghc', 'GHC_p', 'gbj', 'gBJ_p', 'gBJ_pvalue_flag',
							'minP_pvalue', 'skat_p', 'omni_stat',
							'omni_p', 'num_snps')
write.table(cgems_results, out_name, append=F, quote=F, col.names=T, row.names=F)

