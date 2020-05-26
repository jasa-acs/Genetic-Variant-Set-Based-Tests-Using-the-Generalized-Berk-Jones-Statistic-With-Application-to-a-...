# Adapted from gBJ_TypeIerr_final.R

# The purpose of this script is to generate a portion of the data for the simulations 
# which produced Table 1 in the submitted manuscript.
# Here we need to get the correlation matrix of the OMNI test when running on the FGFR2
# regions.
# This output is necessary to calculate the p-value of the OMNI test in make_tab1.R

# We generate predefined sets of SNPs in FGFR2 using HAPGEN2, 
# generate an outcome from under the null, apply GBJ/GHC/SKAT/minP to the simulated data,
# and get the correlation of the transformed p-values from these tests.  
# See the omnibus section of the manuscript for more details.

# This script will produce the two files: low_cor_mat.txt and high_cor_mat.txt.

# After this script as well as all runs from gen_data_tab1_part1.R and gen_data_part2.R have been completed, use 
# make_tab1.R to make Table 1 of the submitted manuscript.
# Since HAPGEN sets a seed for its random number generator according to the time of day, 
# it is not possible to recreate the exact simulations that led to Table 1.
# However since we are running so many simulations, your output should be very similar.

# Make sure that the HAPGEN2 binary is in the working directory, along with the chromosome 5 and
# chromosome 10 reference data from the HapMap3 project.

########################################################################################
# Change this location to the directory that holds the necessary R libraries
# (must be installed first).
libs <- '/n/home/user/Rlibrary/3.3.1'

###################################################################
# No manual changes need to be made after this point
###################################################################

library(data.table, lib=libs)
library(SKAT, lib=libs)
library(GBJ, lib=libs)

###########################################################################
# Function to generate FGFR2 gene of n people.
# Return a low LD subset and a high LD subset
gen_FGFR2 <- function(num_to_generate, aID, Snum)
{
	# Predetermined
	low_LD_snps <- c('rs3135814', 'rs2433759', 'rs3135772', 'rs7090018', 'rs3135761',
					'rs2912787', 'rs2981428', 'rs7899765')
	high_LD_snps <- c('rs3750817', 'rs11200014', 'rs17542768', 'rs1219648',
					'rs17102287', 'rs2420946', 'rs3135715', 'rs1047111')


	# Use hapgen2 to make the haplotypes.
	gen_success <- system2("./hapgen2", args=c("-m", "genetic_map_chr10_combined_b36.txt", "-l", 
						"hapmap3.r2.b36.chr10.legend", "-h", "CEU.chr10.hap", "-o", 
						paste("FGFR2_S", Snum, "_aID", aID, sep=""), "-dl", 123228042, 1, 1, 1, "-n",
						as.character(num_to_generate), 2, "-int", 123228042, 123347962), wait=TRUE)
						
	if (gen_success != 0) {stop('Problem generating FGFR2')}
	
	# Read the map file.
	FGFR2_map <- read.table(paste("FGFR2_S", Snum, "_aID", aID, ".legend", sep=''), header=T)
	low_rows <- which(FGFR2_map[,1] %in% low_LD_snps)
	high_rows <- which(FGFR2_map[,1] %in% high_LD_snps)
	
	# Read the haps file
	haps_file <- fread(paste("FGFR2_S", Snum, "_aID", aID, ".controls.haps", sep=''), header=F)
	
	# Subset into high and low LD, transpose.
	left <- seq(from=1, to=(ncol(haps_file)-1), by=2)
	right <- seq(from=2, to=ncol(haps_file), by=2)
	low_LD_haps <- t( haps_file[low_rows, left, with=F] + haps_file[low_rows, right, with=F] )
	high_LD_haps <- t( haps_file[high_rows, left, with=F] + haps_file[high_rows, right, with=F] )
	
	# Remove everything we wrote to disk
	system2(command="rm", args=paste("FGFR2_S", Snum, "_aID", aID, ".*", sep=""))
	
	# Now we have n people by 8 SNPs for each subset
	return( list(low_LD_haps=low_LD_haps, high_LD_haps=high_LD_haps) )
}
###########################################################################


########################################################################################
# Perform many bootstrap samples to get the necessary correlation matrix.
# Can't pass scientific notation to HAPGEN2
options(scipen=9)
num_subjects <- 2000
num_boots <- 100

# Just generate one set of G, and simulate Y num_boots times
set.seed(100)
all_data <- gen_FGFR2(num_to_generate=num_subjects, aID=1, Snum=1)

set.seed(1000)
typeIerr_results <- matrix(data=NA, nrow=num_boots, ncol=18)
typeIerr_results <- data.frame(typeIerr_results)
colnames(typeIerr_results) <- c('low_GBJ', 'low_GBJ_p', 'low_GBJ_err', 'low_GHC',
								'low_GHC_p', 'low_GHC_err', 'low_minP', 'low_minP_p', 'low_SKAT_p', 
								'high_GBJ', 'high_GBJ_p', 'high_GBJ_err', 'high_GHC', 
								'high_GHC_p', 'high_GHC_err', 'high_minP', 'high_minP_p', 'high_SKAT_p')
for (jjj in 1:num_boots) {
	if (jjj%%10 == 0) {cat(jjj)}
	
	# Simulate outcome.
	Y <- rnorm(n=num_subjects)			# Mean 0, Variance 1
	null_mod <- glm(Y~1, family=gaussian)
	skat_null_obj <- SKAT_Null_Model(Y ~ 1, out_type='C')
	
	# First run with low LD
	G_mat <- all_data$low_LD_haps
			
	# Calculate score statistics
	low_score_stats <- calc_score_stats(null_model=null_mod, 
					factor_matrix=G_mat, link_function='linear')
	test_stats <- low_score_stats$test_stats
	cor_mat <- low_score_stats$cor_mat
	# Run tests
	low_GBJ_list <- GBJ(test_stats=test_stats, cor_mat=cor_mat)
	low_GHC_list <- GHC(test_stats=test_stats, cor_mat=cor_mat)
	low_minP_list <- minP(test_stats=test_stats, cor_mat=cor_mat)
	low_skat_list <- SKAT(Z=G_mat, obj=skat_null_obj, weights.beta=c(1,1), method='davies')			
			
	# Record low LD results
	record_row <- jjj
	typeIerr_results$low_GBJ[record_row] <- low_GBJ_list$GBJ
	typeIerr_results$low_GBJ_p[record_row] <- low_GBJ_list$GBJ_pvalue
	typeIerr_results$low_GBJ_err[record_row] <- low_GBJ_list$err_code
	typeIerr_results$low_GHC[record_row] <- low_GHC_list$GHC
	typeIerr_results$low_GHC_p[record_row] <- low_GHC_list$GHC_pvalue
	typeIerr_results$low_GHC_err[record_row] <- low_GHC_list$err_code
	typeIerr_results$low_minP[record_row] <- low_minP_list$minP
	typeIerr_results$low_minP_p[record_row] <- low_minP_list$minP_pvalue
	typeIerr_results$low_SKAT_p[record_row] <- low_skat_list$p.value
			
	# Now do the high LD
	G_mat <- all_data$high_LD_haps
			
	# Calculate score statistics
	high_score_stats <- calc_score_stats(null_model=null_mod, 
						factor_matrix=G_mat, link_function='linear')
	test_stats <- high_score_stats$test_stats
	cor_mat <- high_score_stats$cor_mat
	# Run tests
	high_GBJ_list <- GBJ(test_stats=test_stats, cor_mat=cor_mat)
	high_GHC_list <- GHC(test_stats=test_stats, cor_mat=cor_mat)
	high_minP_list <- minP(test_stats=test_stats, cor_mat=cor_mat)
	high_skat_list <- SKAT(Z=G_mat, obj=skat_null_obj, weights.beta=c(1,1), method='davies')
			
	# Record high LD results
	typeIerr_results$high_GBJ[record_row] <- high_GBJ_list$GBJ
	typeIerr_results$high_GBJ_p[record_row] <- high_GBJ_list$GBJ_pvalue
	typeIerr_results$high_GBJ_err[record_row] <- high_GBJ_list$err_code
	typeIerr_results$high_GHC[record_row] <- high_GHC_list$GHC
	typeIerr_results$high_GHC_p[record_row] <- high_GHC_list$GHC_pvalue
	typeIerr_results$high_GHC_err[record_row] <- high_GHC_list$err_code
	typeIerr_results$high_minP[record_row] <- high_minP_list$minP
	typeIerr_results$high_minP_p[record_row] <- high_minP_list$minP_pvalue
	typeIerr_results$high_SKAT_p[record_row] <- high_skat_list$p.value
}
########################################################################################



########################################################################################
# Transform the data slightly
# Transform the 1 to 0.999 for qnorm()
p_columns <- c(2, 5, 8, 9, 11, 14, 17, 18)
p_mat <- typeIerr_results[, p_columns]
for (i in 1:ncol(p_mat)) {
	temp_col <- p_mat[,i]
	one_elements <- which(p_mat[, i] == 1) 
	if (length(one_elements) > 0) {
		p_mat[one_elements,i] <- 0.999
	}
}

low_T_mat <- p_mat[,1:4]
for (i in 1:ncol(low_T_mat)) {
	low_T_mat[,i] <- qnorm(1-low_T_mat[,i])
}
low_cor_mat <- cor(low_T_mat)
high_T_mat <- p_mat[,5:8]
for (i in 1:ncol(high_T_mat)) {
	high_T_mat[,i] <- qnorm(1-high_T_mat[,i])
}
high_cor_mat <- cor(high_T_mat)
########################################################################################


########################################################################################
# Write the correlation matrices
write.table(low_cor_mat, file='low_cor_mat.txt', append=F, quote=F,
			row.names=F, col.names=F, sep='\t')
write.table(high_cor_mat, file='high_cor_mat.txt', append=F, quote=F,
			row.names=F, col.names=F, sep='\t')