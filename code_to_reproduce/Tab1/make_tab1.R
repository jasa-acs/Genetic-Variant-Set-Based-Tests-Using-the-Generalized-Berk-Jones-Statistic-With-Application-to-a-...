# The purpose of this script is to recreate Table 1 of the submitted manuscript.

# We read in the output produced by gen_data_tab1_part1.R-gen_data_tab1_part3.R, 
# calculate the omnibus statistic p-values for certain runs, and format the data 
# into tables.
# Simply type in the appropriate variable names (labeled below) to see certain columns of the table. 

# Since HAPGEN sets a seed for its random number generator according to the time of day, 
# it is not possible to recreate the exact simulations that led to Table 1.
# However since we are running so many simulations, your output should be very similar.

# Make sure all the jobs from gen_data_tab1_part1.R-gen_data_tab1_part3.R have finished running,
# otherwise the necessary data will not exist and there will be errors.



########################################################################################
# Change this location to the directory that holds the necessary R libraries
# (must be installed first).
libs <- '/n/home/user/Rlibrary/3.3.1'
library(mvtnorm, lib=libs)

# Change data_dir your working directory directory which holds all the necessary data files.
data_dir <- '/users/user/desktop'

###################################################################
# No manual changes need to be made after this point
###################################################################

##################################################################
# Make columns 2 and 3 of table 1
FGFR2_typeIerr <- read.table('FGFR2_typeIerr_S1_1.txt', header=T)
for(i in 2:400)
{
	temp_fname <- paste('FGFR2_typeIerr_S1_', i, '.txt', sep='')
	temp_file <- read.table(temp_fname, header=T)
	FGFR2_typeIerr <- rbind(FGFR2_typeIerr, temp_file)
	if (i%%100 == 0) {cat(i)}
}

# Make the first two columns of table 1
alpha_vec <- c(0.01, 0.001, 0.0001, 0.00001)
results_tab1 <- matrix(data=NA, nrow=length(alpha_vec), ncol=2)
for (i in 1:length(alpha_vec)) {
	alpha <- alpha_vec[i]
	results_tab1[i, ] <- c(length(which(FGFR2_typeIerr$low_GBJ_p<alpha)) / nrow(FGFR2_typeIerr),
							length(which(FGFR2_typeIerr$high_GBJ_p<alpha)) / nrow(FGFR2_typeIerr))
}
########################################################################################



########################################################################################
# Make columns 5 and 6 of table 1

# We presimulated the omnibus test statistic correlation matrix
# for FGFR2.
# Ensure that these files are also in your data directory.
low_cor_mat <- read.table('low_cor_mat.txt', header=F)
low_cor_mat <- unname(as.matrix(low_cor_mat))
high_cor_mat <- read.table('high_cor_mat.txt', header=F)
high_cor_mat <- unname(as.matrix(high_cor_mat))


# Calculate omnibus statistic
all_low_pvalues <- cbind(FGFR2_typeIerr$low_GBJ_p, FGFR2_typeIerr$low_GHC_p,
						FGFR2_typeIerr$low_minP_p, FGFR2_typeIerr$low_SKAT_p)
all_high_pvalues <- cbind(FGFR2_typeIerr$high_GBJ_p, FGFR2_typeIerr$high_GHC_p,
						FGFR2_typeIerr$high_minP_p, FGFR2_typeIerr$high_SKAT_p)
omni_low_stat <- apply(all_low_pvalues, 1, min)	
omni_low_stat <- qnorm(1-omni_low_stat)	
omni_high_stat <- apply(all_high_pvalues, 1, min)			
omni_high_stat <- qnorm(1-omni_high_stat)

# Omnibus p-values
omni_low_p <- rep(NA, length(omni_low_stat))
omni_high_p <- rep(NA, length(omni_high_stat))
for (i in 1:length(omni_low_p)) {
	omni_low_p[i] <- 1-pmvnorm(lower=rep(-Inf, 4), upper=rep(omni_low_stat[i], 4), sigma=low_cor_mat)
	omni_high_p[i] <- 1-pmvnorm(lower=rep(-Inf, 4), upper=rep(omni_high_stat[i], 4), sigma=high_cor_mat)
}

# Make columns 5 and 6 of table 1
alpha_vec <- c(0.01, 0.001, 0.0001, 0.00001)
results_tab2 <- matrix(data=NA, nrow=length(alpha_vec), ncol=2)
for (i in 1:length(alpha_vec)) {
	alpha <- alpha_vec[i]
	results_tab2[i, ] <- c(length(which(omni_low_p <alpha)) / length(omni_low_p),
							length(which(omni_high_p<alpha)) / length(omni_high_p))
}
########################################################################################


########################################################################################
# Make columns 4 and 7 of table 1

# Have to read in Snum=1, 2, and 3
chr5_typeIerr <- read.table('chr5_typeIerr_S1_1.txt', header=T)
for(i in 2:2000)
{
	temp_fname <- paste('chr5_typeIerr_S1_', i, '.txt', sep='')
	temp_file <- read.table(temp_fname, header=T)
	chr5_typeIerr <- rbind(chr5_typeIerr, temp_file)
	if (i%%100 == 0) {cat(i)}
}

# Cut to only the first 20 million rows
chr5_typeIerr <- chr5_typeIerr[1:20000000, ]

# Make the first two columns of table 1
alpha_vec <- c(0.01, 0.001, 0.0001, 0.00001)
results_tab3 <- matrix(data=NA, nrow=length(alpha_vec), ncol=2)
for (i in 1:length(alpha_vec)) {
	alpha <- alpha_vec[i]
	results_tab3[i, ] <- c(length(which(chr5_typeIerr$GBJ_p<alpha)) / nrow(chr5_typeIerr),
							length(which(chr5_typeIerr$omni_p<alpha)) / nrow(chr5_typeIerr))
}





