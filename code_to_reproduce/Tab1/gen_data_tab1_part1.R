# Adapted from gBJ_TypeIerr_final.R

# The purpose of this script is to generate a portion of the data for the simulations 
# which produced Table 1 in the submitted manuscript.
# We are interested in knowing the Type I error of GBJ and OMNI
# across high and low LD regions of FGFR2.

# We generate FGFR2 SNPs using HAPGEN2, generate an outcome from under the null, and apply
# GBJ and OMNI to the simulated data.

# Each job should produce one file with a name like 'FGFR2_typeIerr_S1_aID1.txt'
# where the numbers after each parameter will index the jobs we are running.
# Our job submission scripts (*.sbatch and *.bsub) control Snum and ID
# so that all necessary files are generated.

# After all runs from this script as well as gen_data_tab1_part2.R and gen_data_part3.R have been completed, use 
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

library(mvtnorm, lib=libs)
library(data.table, lib=libs)
library(SKAT, lib=libs)
library(GBJ, lib=libs)
########################################################################################


########################################################################################
# Read input arguments (automatically provided by job submission scripts)
FGFR2_flag <- TRUE
args <- commandArgs(trailingOnly=TRUE)
aID <- as.numeric(args[1])
Snum <- as.numeric(args[2])

# Can't pass scientific notation to HAPGEN2
options(scipen=9)

# How many subjects, how many genes per hapgen activation, how many genes_total
num_subjects <- 2000
num_genes <- 50000
genes_per_hapgen <- 10
num_hapgen <- num_genes / genes_per_hapgen 			# Make sure this is an integer!
########################################################################################


########################################################################################
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
########################################################################################


########################################################################################
# Function to generate n people with random gene.
# Pass in the gene list and legend_file so we don't have to load every time.
gen_random_gene <- function(num_to_generate, aID, Snum, chr, glist, leg_file)
{
	# Map, legend, and haplotype (CEU only) file name, for input into hapgen
	map_fname <- paste('genetic_map_chr', chr, '_combined_b36.txt', sep='')
	leg_fname <- paste('hapmap3.r2.b36.chr', chr, '.legend', sep='')
	CEUhap_fname <- paste('CEU.chr', chr, '.hap', sep='')
	
	# HAPGEN output name, make it different for each Snum and each aID so no overwrites
	out_name <- paste("S", Snum, "_aID", aID, sep="")
	
	# Keep going until we get a suitable gene with >=5 SNPs
	keep_searching <- TRUE
	while (keep_searching)
	{
		random_gene <- sample(x=1:nrow(glist), size=1)
		start_pos <- glist[random_gene,2]
		end_pos <- glist[random_gene,3]
		gene_name <- as.character(glist[random_gene,4])
		rows_in_gene <- which(leg_file$position>=start_pos & leg_file$position<=end_pos)
		num_snps <- length(rows_in_gene)
		
		# Needs to be at least 5 SNPs in the gene
		if (num_snps < 5) {next}
		
		# Sometimes even if we have 3 SNPs, we will get 'no suitable disease locus'
		keep_generating <- TRUE
		counter <- 0
		while (keep_generating)
		{
			dl <- leg_file$position[rows_in_gene[1]+counter]
			b <- system2("./hapgen2", args=c("-m", map_fname, "-l", leg_fname, 
									"-h", CEUhap_fname, "-o", out_name, "-dl", dl, "1", "1", "1",
									"-n", as.character(num_to_generate), "2", "-int", start_pos, end_pos,
									 "-no_gens_output"), wait=TRUE)
									
			# If good generation, then break out of both loops
			if (b==0) {
				keep_generating <- FALSE
				keep_searching <- FALSE
			}	
			# Move the disease locus.
			counter <- counter + 1
			if (counter+1 > num_snps) {keep_generating <- FALSE}	
		}
	}

	# Read the haps file
	haps_file <- fread(paste(out_name, ".controls.haps", sep=''), header=F)
	
	# Subset into high and low LD, transpose.
	left <- seq(from=1, to=(ncol(haps_file)-1), by=2)
	right <- seq(from=2, to=ncol(haps_file), by=2)
	random_g_haps <- t( haps_file[, left, with=F] + haps_file[, right, with=F] )
	
	# Remove everything we wrote to disk
	system2(command="rm", args=paste("S", Snum, "_aID", aID, "*", sep=""))
	
	# Check for MAF and sd, if we fail too many, then return -1
	freqs <- apply(random_g_haps, 2, mean)/2
	std_devs <- apply(random_g_haps, 2, sd)
	bad_snps <- unique(c(which(freqs<0.05), which(freqs>0.95), which(std_devs==0)))
	if (length(bad_snps)>0)
	{
		if ( (num_snps-length(bad_snps)) >= 3)			# Still 3 snps left? OK	
		{
			random_g_haps <- random_g_haps[,-bad_snps]
		} else {
			return (-1)
		}
	}
	
	# Now we have n people by some # of SNPs for each subset
	return( list(random_g_haps=random_g_haps, gene_name=gene_name) )
}
########################################################################################


########################################################################################
# Run the simulation

# Check if doing FGFR or random genes.
if (FGFR2_flag) {
	typeIerr_results <- matrix(data=NA, nrow=num_genes, ncol=18)
	typeIerr_results <- data.frame(typeIerr_results)
	colnames(typeIerr_results) <- c('low_GBJ', 'low_GBJ_p', 'low_GBJ_err', 'low_GHC',
									'low_GHC_p', 'low_GHC_err', 'low_minP', 'low_minP_p', 'low_SKAT_p', 
									'high_GBJ', 'high_GBJ_p', 'high_GBJ_err', 'high_GHC', 
									'high_GHC_p', 'high_GHC_err', 'high_minP', 'high_minP_p', 'high_SKAT_p')
} else {
	glist <- read.table('glist-hg18.txt', header=F)
	glist <- glist[which(glist[,1]==10),]
	leg_file <- read.table('hapmap3.r2.b36.chr10.legend', header=T)
	typeIerr_results <- matrix(data=NA, nrow=num_genes, ncol=8)
	colnames(typeIerr_results) <- c('Zlargeflag', 'GBJpflag', 'gbj', 'GBJp', 'ghc', 'GHCp', 
									'burdenP','gene', 'num_snps')
}

# Start looping through genes
for (iii in 1:num_hapgen)
{
	# Get the genotypes, we have num_gens sets of them
	all_data <- -1
	while (class(all_data) == 'numeric')
	{
		if (FGFR2_flag) {
			all_data <- gen_FGFR2(num_to_generate=(num_subjects*genes_per_hapgen), aID=aID, Snum=Snum)
		} else {
			all_data <- gen_random_gene(num_to_generate=(num_subjects*genes_per_hapgen), aID=aID, Snum=Snum, chr=10, 
						glist=glist, leg_file=leg_file)
		}
	}
	
	for (jjj in 1:genes_per_hapgen)
	{
		# Simulate outcome.
		Y <- rnorm(n=num_subjects)			# Mean 0, Variance 1
		null_mod <- glm(Y~1, family=gaussian)
		skat_null_obj <- SKAT_Null_Model(Y ~ 1, out_type='C')
		
		record_row <- genes_per_hapgen * (iii-1) + jjj
	
		# First run with low LD
		G_mat <- all_data$low_LD_haps[((jjj-1)*num_subjects+1):(jjj*num_subjects), ]
			
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
		G_mat <- all_data$high_LD_haps[((jjj-1)*num_subjects+1):(jjj*num_subjects), ]
			
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
}
########################################################################################

########################################################################################
# Write down the results
out_name <- paste('FGFR2_typeIerr_S', Snum, '_', aID, '.txt', sep='')
write.table(typeIerr_results, out_name, append=F, quote=F, col.names=T, row.names=F, sep='\t')

