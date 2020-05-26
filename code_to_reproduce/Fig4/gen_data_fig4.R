# Adapted from gBJ_power_sim_HAPGEN.R

# The purpose of this script is to run the simulations that are displayed in Figure 4
# of the submitted manuscript.
# We are interested in the power of different set-based tests (GBJ, GHC, SKAT, minP, OMNI)
# on SNP-sets which mimic the structure of genotypes found in real data.

# We generate the data using HAPGEN, simulate outcomes according to the model given in the 
# paper, and apply tests to the simluated data.

# Each job should produce one file with a name like 'chr5_power_S1_aID1.txt'
# where the numbers after Snum and aID will vary with each job.
# The parameters for each run are controlled by the aID and Snum variables, which are parameters passed
# in on the command line. 
# Our job submission scripts (*.sbatch and *.bsub) control these parameters
# so that all necessary files are generated.

# After all of these jobs have output their data, use plot_fig4.R to create a plot similar to Figure 4.
# Since HAPGEN sets a seed for its random number generator according to the time of day, 
# it is not possible to recreate the exact simulations that led to Figure 4.
# However since we are running so many simulations, your output should be very similar.

# Make sure that the HAPGEN2 binary is in the working directory, along with the chromosome 5 reference
# data from the HapMap3 project.


########################################################################################
# Change this location to the directory that holds the necessary R libraries
# (must be installed first).
libs <- '/n/home/user/Rlibrary/3.3.1'

###################################################################
# No manual changes need to be made after this point
###################################################################

library(mvtnorm, lib=libs)
library(SKAT, lib=libs)
library(data.table, lib=libs)
library(GBJ, lib=libs)
########################################################################################


########################################################################################
# Read input arguments (fed automatically by job submission script) to determine simulation parameters 
args <- commandArgs(trailingOnly=TRUE)
aID <- as.numeric(args[1])
Snum <- as.numeric(args[2])
set.seed(05012017 + Snum*1000 + aID)
options(scipen=9)

# Parameters for the run
runs_per_numcausal <- 40
num_causal <- ceiling(aID / runs_per_numcausal)
CHR <- 5
gene_size <- 40
genes_per_run <- 50
omni_cor_boots <- 100
hapgen_per_gene <- 1
runs_per_hapgen <- 10			
sims_per_gene <- hapgen_per_gene*runs_per_hapgen # Number of sims per gene is hapgen_per_gene*runs_per_hapgen
num_sub <- 2000
leg_file <- fread('hapmap3.r2.b36.chr5.legend')

# The effect sizes shrink as number of effects increases
beta_possible <- c(0.16, 0.14, 0.12, 0.12, 0.11, 0.11, 0.10, 0.10)
beta_causal <- beta_possible[num_causal]
########################################################################################



########################################################################################
# While in R, create haps from HAPGEN using system2.
# Pass it the start location, end location, and chromosome.
# Have to pass it the actual legend file (too large to read each time).
# Then remove the bad MAF SNPs, add together 2 haplotypes to make a 'person', transpose to get
# SNPs on the columns.
# If we pass in the SNPs to keep, then use just keep that subset, otherwise cut by MAF=5%
generate_hapgen <- function(n, gene_size, aID, Snum, start_row, end_row, chr, leg_file, snps_to_keep=NULL) {
	
	# HAPGEN output name, make it different for each Snum and each aID so don't overwrite each other
	out_name <- paste('S', Snum, '.ch', aID, '.out', sep='')
	
	# Find start BP and end BP
	start_bp <- leg_file$position[start_row]
	end_bp <- leg_file$position[end_row]
	
	# Disease locus is just the one after the start position
	dl <- leg_file$position[start_row+1]
	
	# Map, legend, and haplotype (CEU only) file name
	map_fname <- paste('genetic_map_chr', chr, '_combined_b36.txt', sep='')
	leg_fname <- paste('hapmap3.r2.b36.chr', chr, '.legend', sep='')
	CEUhap_fname <- paste('CEU.chr', chr, '.hap', sep='')
	
	# Create the haps, return b=0 for success, b=1 for failure (no suitable disease locus)
	b <- 1
	counter <- 2
	while(b==1) {
		# generating n 'controls' and 2 'cases'
		b <- system2("./hapgen2", args=c("-m", map_fname, "-l", leg_fname, "-h", CEUhap_fname, "-o", out_name, "-dl", dl, "1", "1", "1", "-n", n, "2", "-int", start_bp, end_bp, "-no_gens_output"))		
		# move the disease locus
		dl <- leg_file$position[start_row+counter]
		counter <- counter + 1
		if ( (start_row+counter) > end_row) {
			return( list(G_mat=-1, kept_snps=-1) )
		}
	}
		
	# Change the name of the haps we need
	orig_name <- paste(out_name, '.controls.haps', sep='')
	new_name <- paste('S', Snum, '.aID', aID, '.haps', sep='')
	system2("mv", args=c(orig_name, new_name))
	
	# Remove everything else
	remove_name <- paste('S', Snum, '.ch', aID, '.out.*', sep='')
	system2("rm", args=remove_name)
	
	# Read in the haplotype file
	# Every column is one 'strand', add two together for a person.
	# Every row is a different SNP.
	haps <- fread(new_name)
	system2("rm", args=new_name)
	
	# Combine to make 'person'
	left <- seq(from=1, to=(2*n-1), by=2)
	right <- seq(from=2, to=2*n, by=2)
	combined_haps <- haps[,left, with=F] + haps[,right, with=F]
	rm(haps)
	
	# If didn't specify SNPs_to_keep, then remove for MF
	if (is.null(snps_to_keep)) {
		d <- nrow(combined_haps)
		MAF <- 0.05
		freqs <- apply(combined_haps,1,sum) / (2*n)
		bad_snps <- c(which(freqs>(1-MAF)),which(freqs<MAF))
		kept_snps <- 1:d
		if (length(bad_snps) > 0)
		{
			combined_haps <- combined_haps[-bad_snps,]
			kept_snps <- kept_snps[-bad_snps]
		}
		 # Enough snps generated
       		if (nrow(combined_haps) >=  gene_size) {
               		kept_snps <- kept_snps[1:gene_size]
			combined_haps <- combined_haps[1:gene_size, ]	
        	}
	} else {			# Specified the SNPs to keep
		combined_haps <- combined_haps[snps_to_keep,] 
		kept_snps <- snps_to_keep
	}
	
	# Not enough snps generated
	if (nrow(combined_haps) < gene_size) {
		return( list(G_mat=-1, kept_snps=-1) )
	}
	
	# Transpose and return
	combined_haps <- t(combined_haps)
	return( list(G_mat=combined_haps, kept_snps=kept_snps) )
}
########################################################################################



########################################################################################
# Get correlation between set-based test p-values for use with omnibus test.  
# Input the null hypothesis fitted outcomes Y_hat, resimulate data using those fitted values,
# each time record 4 test statistics, repeat ~100 times.
omnibus_correlations <- function(G_mat, num_boots, fitted_values)
{
	# For each set of test statistics run all the tests
	T_mat <- matrix(data=NA, nrow=num_boots, ncol=4)
	
	# Bootstrap the test p-values under the null
	for (temp_boot in 1:num_boots)
	{
		# Checkpoint
		if (temp_boot%%10 == 0) {
			cat('omni: ', temp_boot, '\n')
		}
		
		# New outcome under the null
		sim_Y <- fitted_values + rnorm(length(fitted_values))
		link_function <- 'linear'
		null_mod <- glm(sim_Y~ 1, family=gaussian)
		
		# Calculate test stats
		test_stat_output <- calc_score_stats(null_model=null_mod, factor_matrix=G_mat, 
							link_function=link_function)
		test_stats <- test_stat_output$test_stats
		
		
		# Run the 4 tests
		GBJ_omni <- GBJ(test_stats=test_stat_output$test_stats, cor_mat=test_stat_output$cor_mat)			
		GHC_omni <- GHC(test_stats=test_stat_output$test_stats, cor_mat=test_stat_output$cor_mat)	
		minP_omni <- minP(test_stats=test_stat_output$test_stats, cor_mat=test_stat_output$cor_mat)	
		skat_null_obj <- SKAT_Null_Model(sim_Y ~  1, out_type='C')
		skat_omni <- SKAT(Z=G_mat, obj=skat_null_obj, weights.beta=c(1,1), method='davies')					 						
		# Put them all in a vector	
		pval_vec <- rep(NA, 4)			
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
		
		# Do the inverse normal transform because ultimately we want to apply a Gaussian copula.
		T_mat[temp_boot, ] <- qnorm(1-pval_vec)
	}
	
	# Return the correlation between the transformed group testing statistics
	omni_cor_mat <- cor(T_mat)
	return(omni_cor_mat)
} 


########################################################################################
# Start looping across genes
pow_results <- matrix(data=NA, nrow=genes_per_run*sims_per_gene, ncol=26)
pow_results <- data.frame(pow_results)
colnames(pow_results) <- c('start_row', 'end_row', 'd', 'ncausal', 'beta_causal', 
							'median_rho1', 'median_rho1_abs', 'median_rho3', 'median_rho3_abs', 'median_rho2',
							'median_rho2_abs', 'median_rho', 'median_rho_abs',
							'GBJ', 'GBJ_p', 'GBJ_err', 'GHC', 'GHC_p', 'GHC_err',
							'minP', 'minP_pvalue', 'skat_p', 'omni_stat',
							'omni_p', 'BJ', 'BJ_p')
for (i in 1:genes_per_run) {
	cat('i=', i, '\n')
	
	# Get the new 'gene'
	start_row <- ceiling(runif(n=1, min=100, max=nrow(leg_file)-10*gene_size))
	end_row <- ceiling(start_row + 2.5*gene_size)
	
	# Figure out which SNPs we are going to keep
	hapgen_output <- generate_hapgen(n=num_sub, gene_size=gene_size, aID=aID, Snum=Snum, start_row=start_row,
                                                                        end_row=end_row, chr=CHR, leg_file=leg_file)	
	# Not enough SNPs in the interval?
	while ( class(hapgen_output$G_mat) == 'numeric' ) {
		hapgen_output <- generate_hapgen(n=num_sub, gene_size=gene_size, aID=aID, Snum=Snum, start_row=start_row,
									end_row=end_row, chr=CHR, leg_file=leg_file)
	}
	G_mat <- hapgen_output$G_mat
	snps_to_keep <- hapgen_output$kept_snps
	freqs <- apply(G_mat, 2, mean) / 2
	to_flip <- which(freqs > 0.5) 
	if (length(to_flip) > 0) 
	{
		G_mat[, to_flip] <- 2 - G_mat[, to_flip]
	}
	
	# Record preliminary information
	beta_vec <- rep(beta_causal, num_causal)
	pow_results$d[((i-1)*sims_per_gene+1):(i*sims_per_gene)] <- gene_size
	pow_results$ncausal[((i-1)*sims_per_gene+1):(i*sims_per_gene)] <- num_causal
	pow_results$beta_causal[((i-1)*sims_per_gene+1):(i*sims_per_gene)] <- beta_causal
	pow_results$start_row[((i-1)*sims_per_gene+1):(i*sims_per_gene)] <- start_row
	pow_results$end_row[((i-1)*sims_per_gene+1):(i*sims_per_gene)] <- end_row
	
	# Get the omnibus correlation matrix
	Y <- rnorm(n=num_sub)
	null_mod <- glm(Y ~ 1, family=gaussian)
	omni_cor_mat <- omnibus_correlations(G_mat=G_mat, num_boots=omni_cor_boots, fitted_values=null_mod$fitted.values)

	# We'll have to do a bunch of hapgens for the same gene
	for (j in 1:hapgen_per_gene) {
		
		cat('j=', j, '\n')
		
		# Get the new hapgen data
		hapgen_output <- generate_hapgen(n=(num_sub*runs_per_hapgen), gene_size=gene_size, aID=aID, Snum=Snum, 
									start_row=start_row,end_row=end_row, chr=CHR, leg_file=leg_file,
									snps_to_keep=snps_to_keep)
									
		for (k in 1:runs_per_hapgen) {
						
			# Get the new G
			G_mat <- hapgen_output$G_mat[((k-1)*num_sub+1):(k*num_sub), ]	
			freqs <- apply(G_mat, 2, mean) / 2
			to_flip <- which(freqs > 0.5) 
			if (length(to_flip) > 0) {
				G_mat[, to_flip] <- 2 - G_mat[, to_flip]
			}
			
			# Generate data under the alternative
			causal_pos <- sample(1:gene_size, size=num_causal, replace=FALSE)
			mu <- as.matrix(G_mat[, causal_pos]) %*% beta_vec 
			Y <- mu + rnorm(n=num_sub)
			null_mod <- glm(Y ~ 1, family=gaussian)
			fitted_Y <- null_mod$fitted.values

			# Calculate score stats
			score_stats_output <- calc_score_stats(null_model=null_mod, factor_matrix=G_mat, 
				link_function='linear')
		
			# Run all the tests
			BJ_output <- BJ(test_stats=score_stats_output$test_stats, cor_mat=score_stats_output$cor_mat)
			GBJ_output <- GBJ(test_stats=score_stats_output$test_stats, cor_mat=score_stats_output$cor_mat)
			GHC_output <- GHC(test_stats=score_stats_output$test_stats, cor_mat=score_stats_output$cor_mat)
			minP_output <- minP(test_stats=score_stats_output$test_stats, cor_mat=score_stats_output$cor_mat)
			skat_null_obj <- SKAT_Null_Model(formula=Y ~  1, out_type='C')
			skat_output <- SKAT(Z=G_mat, obj= skat_null_obj, weights.beta=c(1,1), method='davies')
			
			# Calculate omnibus test p-value.
			# Doesn't actually matter that we have the correct order since all the bounds are the same.
			allp_vec <- c(GBJ_output$GBJ_pvalue, GHC_output$GHC_pvalue,
						skat_output$p.value, minP_output$minP_pvalue)
			omni_stat <- min(allp_vec)
		
			# Make bounds and calculate
			omni_bound <- qnorm(1-omni_stat)
			omni_pvalue <- 1 - pmvnorm(lower=rep(-Inf, 4), upper=rep(omni_bound, 4), 
								sigma=omni_cor_mat)

			# Record various statistics and diagnostic information
			record_row <- (i-1)*sims_per_gene + (j-1)*runs_per_hapgen + k
			cor_mat <- score_stats_output$cor_mat
			rho1_mat <- as.matrix(cor_mat[, causal_pos])
			rho1_mat <- rho1_mat[causal_pos, ]
			rho3_mat <- as.matrix(cor_mat[, -causal_pos])
			rho3_mat <- rho3_mat[-causal_pos, ]
			rho2_mat <- as.matrix(cor_mat[, causal_pos])
			rho2_mat <- rho2_mat[-causal_pos, ]
			pow_results$median_rho1[record_row] <- median(rho1_mat)
			pow_results$median_rho1_abs[record_row] <- median(abs(rho1_mat))
			pow_results$median_rho3[record_row] <- median(rho3_mat)
			pow_results$median_rho3_abs[record_row] <- median(abs(rho3_mat))
			pow_results$median_rho2[record_row] <- median(rho2_mat)
			pow_results$median_rho2_abs[record_row] <- median(abs(rho2_mat))
			pow_results$median_rho[record_row] <- median(cor_mat[lower.tri(cor_mat)]) 
			pow_results$median_rho_abs[record_row] <- median(abs(cor_mat[lower.tri(cor_mat)]))
			pow_results$GBJ[record_row] <- GBJ_output$GBJ
			pow_results$GBJ_p[record_row] <- GBJ_output$GBJ_pvalue
			pow_results$GBJ_err[record_row] <- GBJ_output$err_code
			pow_results$GHC[record_row] <- GHC_output$GHC
			pow_results$GHC_p[record_row] <- GHC_output$GHC_pvalue
			pow_results$GHC_err[record_row] <- GHC_output$err_code
			pow_results$minP[record_row] <- minP_output$minP
			pow_results$minP_pvalue[record_row] <- minP_output$minP_pvalue
			pow_results$skat_p[record_row] <- skat_output$p.value
			pow_results$BJ[record_row] <- BJ_output$BJ
			pow_results$BJ_p[record_row] <- BJ_output$BJ_pvalue
			pow_results$omni_stat[record_row] <- omni_stat
			pow_results$omni_p[record_row] <- omni_pvalue
		}
	}
}
########################################################################################


########################################################################################
# Write down the results
out_name <- paste('chr5_power_S', Snum, '_', aID, '.txt', sep='')
write.table(pow_results, out_name, append=F, quote=F, col.names=T, row.names=F, sep='\t')



















