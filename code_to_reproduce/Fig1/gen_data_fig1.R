# Adapted from GOF_bound_finder.R

# This script is used to reproduce Figure 1 from our submitted manuscript.
# For a given SNP-set structure, the code below will find the rejection region
# of the Higher Criticism, Berk-Jones, Generalized Higher Criticism, and Generalized
# Berk-Jones tests at level alpha.

# The SNP-set structure and alpha will change for each run, and they are controlled
# by parameters (aID and Snum) passed in on the command line.
# Our job submission scripts (*.sbatch and *.bsub) will control these parameters so that
# the correct values are used.
# Output is a text file containing the bounds; plot_fig1.R will read these files to make
# the plots in Figure 1.

# You will need to have the compiled binary "ebb_crossprob_cor" in the
# same directory as this file before running!


###################################################################
# These parameters are specific to each computing environment, 
# will need to be changed depending on where your libraries are!
###################################################################
libs <- '/n/home/user/R/library/3.3.1'
library(GBJ, lib=libs)

###################################################################
# No manual changes need to be made after this point
###################################################################


########################################################################################
# Takes in 2 arguments at runtime.  Our job submission script controls these values.
args <- commandArgs(trailingOnly=TRUE)
aID <- as.numeric(args[1])
Snum <- as.numeric(args[2])

# Parameter settings we want to test
param_row <- aID
d_vec <- c(20, 50, 100, 200)
alpha_vec <- c(0.01, 0.05)
rho_vec <- c(0.3)
percent_cor <- c(0.25, 0.5, 0.75, 1)
param_mat <- expand.grid(alpha_vec, rho_vec, percent_cor, d_vec)
param_mat <- data.frame(param_mat)
colnames(param_mat) <- c('alpha', 'rho', 'percent', 'd')

# Parameters for each run are set here
d <- param_mat$d[param_row]
alpha <- param_mat$alpha[param_row]
rho <- param_mat$rho[param_row]
percent_cor <- param_mat$percent[param_row]
sig_vec <- rep(0, (d*(d-1)/2))
num_cor <- percent_cor * d
sig_vec[1:(num_cor*(num_cor-1)/2)] <- rho

# Have to write the correlation vector to disk so ebb_crossprob_cor can see it later.
sig_vec_name <- paste('sigvec_S', Snum, '_aID', aID, '.txt', sep='')
write.table(sig_vec, sig_vec_name, append=F, quote=F, row.names=F, col.names=F)
########################################################################################


########################################################################################
# The following functions are defined to 'invert' the observed statistic value and give us
# rejection region in terms of the order statistics of the absolute values of the marginal
# test statistics.  For GHC and BJ we need to use uniroot(), for HC the inversion is analytic.

########################################################################################
# Inverse of the HC objective.
HC_inv <- function(h, d) {
	i_vec <- 1:d
	((2*i_vec+h^2)/d - sqrt((2*i_vec/d+h^2/d)^2 - 4*i_vec^2/d^2 - 4*i_vec^2*h^2/d^3))/(2*(1+h^2/d))
}
########################################################################################

########################################################################################
# The GHC objective function minus the observed value, put it into uniroot to find the bounds.		
GHC_func <- function(x, kkk, d, ghc, sig_vec) {
	# we get errors in qnorm once x is too small
	if(x<10^(-16)) {
		temp_Z <- qnorm(1-10^(-15)/2)
	} else {
		temp_Z <- qnorm(1-x/2)
	}
	(kkk - d*x) / sqrt(calc_var_mu_nonzero(d=d, t=temp_Z, mu=0, sig_vec=sig_vec)) - ghc
}
########################################################################################


########################################################################################
# The BJ objective function minus the observed value, put it into uniroot to find the bounds.
BJ_func <- function(x, k, d, b) {
	k*log(k/(d*x)) + (d-k)*log((1-k/d)/(1-x)) - b
}
########################################################################################



########################################################################################
# (1) Start by finding the HC bounds

# Use root-finding on this function to find the value of h that gets us to the level \alpha.
HC_pvalue_root <- function(h) {
	HC_p_bounds <- HC_inv(h, d)
	HC_z_bounds <- qnorm(1-HC_p_bounds/2)
	HC_z_bounds <- sort(HC_z_bounds, decreasing=F)
	HC(test_stats=HC_z_bounds, pairwise_cors=sig_vec)$HC_pvalue - alpha
}

# Find the necessary h
h <- tryCatch(uniroot(HC_pvalue_root, interval=c(1, 20))$root, warning=function(w) w,
				error=function(e) e)
if (length(class(h)) > 1) {
	h <- uniroot(HC_pvalue_root, interval=c(1, 100))$root
}

# Define the rejection region at the given h
HC_p_bounds <- HC_inv(h, d)
HC_z_bounds <- qnorm(1-HC_p_bounds/2)
HC_z_bounds <- sort(HC_z_bounds, decreasing=F)

# Check our answer with the C++ binary p-value calculation
hcbounds_name <- paste('HCbounds_S', Snum, '_aID', aID, '.txt', sep='')
write.table(HC_z_bounds, file=hcbounds_name, append=F, quote=F, row.names=F, col.name=F, sep=', ')
pvalue <- as.numeric(system2(command="./ebb_crossprob_cor", args=c(d, hcbounds_name, sig_vec_name), stdout=TRUE))
if (abs(pvalue - alpha) > 10^(-4)) {
	stop('Bad calculation for HC')
}
		
# Remove bounds so we don't get mixed up for the next iteration
system2(command="rm", args=hcbounds_name, wait=TRUE)
########################################################################################


########################################################################################
# (2) Now do the BJ bounds

# Use root-finding on this function to find the value of b that gets us to the level \alpha.
BJ_pvalue_root <- function(b) {
	BJ_p_bounds <- rep(NA, d)
	
	# Nested uniroot to 'invert' the test statistic and find the bounds for pvalue calculation
	for ( jjj in 1:(floor(d/2)) ) {
		BJ_p_bounds[jjj] <- uniroot(BJ_func, k=jjj, d=d, b=b, lower=0, upper=jjj/d, tol=(10^(-12)))$root
	}
	# The last half of the order statistic bounds are the same
	BJ_p_bounds[(floor(d/2)+1):d] <- BJ_p_bounds[floor(d/2)]
	
	# now put the bounds in terms of the Z statistics
	BJ_z_bounds <- qnorm(1 - BJ_p_bounds/2)
	BJ_z_bounds <- sort(BJ_z_bounds, decreasing=F)

	# qnorm can't handle more precision than ~1*10^-15
	if(length(which(BJ_z_bounds==Inf))>0) {
		BJ_z_bounds[which(BJ_z_bounds==Inf)]= 8.209536
	}
	
	BJ(test_stats=BJ_z_bounds, pairwise_cors=sig_vec)$BJ_pvalue - alpha
}

# Use uniroot to get the observed BJ statistic which would result in our desired p-value.
b <- uniroot(BJ_pvalue_root, interval=c(1, 50))$root
if (length(class(b)) > 1) {
	b <- uniroot(BJ_pvalue_root, interval=c(1, 200))$root
}
BJ_p_bounds <- rep(NA, d)

# Then use uniroot to 'invert' the observed test statistic and find the bounds for p-value calculation.
for ( jjj in 1:(floor(d/2)) ) {
	BJ_p_bounds[jjj] <- uniroot(BJ_func, k=jjj, d=d, b=b, lower=0, upper=jjj/d, tol=(10^(-12)))$root
}
# The last half of the order statistic bounds are the same
BJ_p_bounds[(floor(d/2)+1):d] <- BJ_p_bounds[floor(d/2)]
	
# Now put the bounds in terms of the Z statistics
BJ_z_bounds <- qnorm(1 - BJ_p_bounds/2)
BJ_z_bounds <- sort(BJ_z_bounds, decreasing=F)

# qnorm can't handle more precision than ~1*10^-15
if(length(which(BJ_z_bounds==Inf))>0) {
	BJ_z_bounds[which(BJ_z_bounds==Inf)]= 8.209536
}

# Check our answer with the p-value binary
bjbounds_name <- paste('BJbounds_S', Snum, '_aID', aID, '.txt', sep='')
write.table(BJ_z_bounds, file=bjbounds_name, append=F, quote=F, row.names=F, col.name=F, sep=', ')
pvalue <- as.numeric(system2(command="./ebb_crossprob_cor", args=c(d, bjbounds_name, sig_vec_name), stdout=TRUE))
if (abs(pvalue - alpha) > 10^(-4)) {
	stop('Bad calculation for BJ')
}
		
# Remove so we don't get mixed up for the next iteration
system2(command="rm", args=bjbounds_name, wait=TRUE)
########################################################################################

	
########################################################################################
# (3) Next find the GHC bounds

# Use root-finding on this function to find the value of b that gets us to the level \alpha.
my_tol <- (-12)
GHC_pvalue_root <- function(ghc) {
	GHC_p_bounds <- rep(NA, d)
	
	# Use uniroot to 'invert' the test statistic and find bounds for p-value calculation.
	GHC_lowerbound <- 10^(-20)
	for(temp_it in 1:d) {
		temp_ghc <- uniroot(GHC_objective, k=temp_it, d=d, offset=ghc, 
				pairwise_cors=sig_vec, interval=c(GHC_lowerbound,(1-10^(-12))), 
				tol=(10^(my_tol)))	
			
		# if it doesn't work, just return NA
		if(length(class(temp_ghc))>1) {
			return (NA)
		} 
		GHC_p_bounds[temp_it] <- temp_ghc$root
			
		# small security measure to ensure that GHC bounds are increasing
		GHC_lowerbound <- GHC_p_bounds[temp_it]
	}
		
	# Now put the bounds in terms of the Z statistics
	GHC_z_bounds <- qnorm(1-GHC_p_bounds/2)
	GHC_z_bounds <- sort(GHC_z_bounds, decreasing=F)
	
	# qnorm can't handle more precision than ~1*10^-16
	# Also crossprob_cor can only handle Z up to 8.2
	GHC_z_bounds[which(GHC_z_bounds > 8.2)]= 8.2

	GHC(test_stats=GHC_z_bounds, pairwise_cors=sig_vec)$GHC_pvalue - alpha
}

# Use uniroot to get the observed BJ statistic which would result in our desired p-value.
ghc <- uniroot(GHC_pvalue_root, interval=c(1, 20))$root
if (length(class(ghc)) > 1) {
	ghc <- uniroot(GHC_pvalue_root, interval=c(1, 200))$root
}
		
# Use uniroot to 'invert' the observed ghc statistic and get the bounds
GHC_p_bounds <- rep(NA, d)
GHC_lowerbound <- 10^(-20)
for(temp_it in 1:d) {
	
	# Shouldn't have precision errors
	temp_ghc <- uniroot(GHC_objective, k=temp_it, d=d, offset=ghc, 
				pairwise_cors=sig_vec, interval=c(GHC_lowerbound,(1-10^(-12))), 
				tol=(10^(my_tol)))			
	
	GHC_p_bounds[temp_it] <- temp_ghc$root
			
	# small security measure to ensure that GHC bounds are increasing
	GHC_lowerbound <- GHC_p_bounds[temp_it]
}
		
# now put the bounds in terms of the Z statistics
GHC_z_bounds <- qnorm(1-GHC_p_bounds/2)
GHC_z_bounds <- sort(GHC_z_bounds, decreasing=F)
	
# qnorm can't handle more precision than ~1*10^-16
# Also crossprob_cor can only handle Z up to 8.2
GHC_z_bounds[which(GHC_z_bounds > 8.2)]= 8.2
	
# Check work on these bounds
ghcbounds_name <-  paste('GHCbounds_S', Snum, '_aID', aID, '.txt', sep='')
write.table(GHC_z_bounds, file=ghcbounds_name, append=F, quote=F, row.names=F, col.name=F, sep=', ')
pvalue <- as.numeric(system2(command="./ebb_crossprob_cor", args=c(d, ghcbounds_name, sig_vec_name), stdout=TRUE))
if (abs(pvalue - alpha) > 10^(-4)) {
	stop('Bad calculation for GHC')
}

# Remove so we don't get mixed up for the next iteration
system2(command="rm", args=ghcbounds_name, wait=TRUE)
########################################################################################

		
########################################################################################
# (4) Finally find the GBJ bounds

# Function to pass to uniroot, find the value of gbj that gives us a p-value of \alpha
GBJ_pvalue_root <- function(gbj) {
	gBJ_z_bounds <- rep(NA,d)
	prev_bound <- 8.2
	
	# Use uniroot to 'invert' the observed value and get bounds for p-value calculation
	for ( kkk in 1:(floor(d/2)) ) 
	{
		temp_gbj <- tryCatch(uniroot(GBJ_objective, interval = c(0, 
                  prev_bound), d = d, k_vec = kkk, pairwise_cors = sig_vec, 
                  offset = gbj), error = function(e) e, 
                  warning = function(w) w)
        if (length(class(temp_gbj)) > 1) {
       		break
        }
        else {
        	gBJ_z_bounds[kkk] <- temp_gbj$root
        }
        prev_bound <- gBJ_z_bounds[kkk]
    }
			
	# The last half of the GBJ bounds are the same.
	gBJ_z_bounds[(floor(d/2)+1):d] <- gBJ_z_bounds[floor(d/2)]
	gBJ_z_bounds <- sort(gBJ_z_bounds, decreasing=FALSE)

	# qnorm can't handle more precision than 10^-16
	# Also crossprob_cor can only handle Z up to 8.2
	gBJ_z_bounds[which(gBJ_z_bounds > 8.2)]= 8.2
	
	GBJ(test_stats = gBJ_z_bounds, pairwise_cors=sig_vec)$GBJ_pvalue - alpha
}
	
# Find the observed gbj value which will give us the desired p-value.
gbj <- uniroot(GBJ_pvalue_root, interval=c(1, 10))$root
if (length(class(gbj)) > 1) {
	gbj <- uniroot(GBJ_pvalue_root, interval=c(1, 200))$root
}
		
# Use uniroot to 'invert' the observed statistic and make the rejection region
gBJ_z_bounds <- rep(NA,d)
prev_bound <- 8.2
for ( kkk in 1:(floor(d/2)) ) { 
		temp_gbj <- tryCatch(uniroot(GBJ_objective, interval = c(0, 
                  prev_bound), d = d, k_vec = kkk, pairwise_cors = sig_vec, 
                  offset = gbj), error = function(e) e, 
                  warning = function(w) w)
        if (length(class(temp_gbj)) > 1) {
       		break
        }
        else {
        	gBJ_z_bounds[kkk] <- temp_gbj$root
        }
        prev_bound <- gBJ_z_bounds[kkk]
}
		
# Only need the first half 
# Make sure to sort in increasing order for crossprob_cor	
gBJ_z_bounds[(floor(d/2)+1):d] <- gBJ_z_bounds[floor(d/2)]
gBJ_z_bounds <- sort(gBJ_z_bounds, decreasing=FALSE)

# qnorm can't handle more precision than ~1*10^-16
# Also crossprob_cor can only handle Z up to 8.2
gBJ_z_bounds[which(gBJ_z_bounds > 8.2)]= 8.2

# Check work with c++ binary
gBJ_bounds_name <- paste('GBJbounds_S', Snum, '_aID', aID, '.txt', sep='')
write.table(gBJ_z_bounds, file=gBJ_bounds_name, append=F, quote=F, row.names=F, col.name=F, sep=', ')
pvalue <- as.numeric(system2(command="./ebb_crossprob_cor", args=c(d, gBJ_bounds_name, sig_vec_name), stdout=TRUE))
if (abs(pvalue - alpha) > 10^(-4)) {
	stop('Bad calculation for GBJ')
}

# Remove so we don't get mixed up for next iteration
system2(command="rm", args=gBJ_bounds_name, wait=TRUE)

# Remove sigvec
system2(command="rm", args=sig_vec_name, wait=TRUE)
########################################################################################


########################################################################################
# Write to disk the table with all of the bounds
bounds_data <- cbind(1:d, BJ_z_bounds, gBJ_z_bounds, HC_z_bounds, GHC_z_bounds)
bounds_data <- data.frame(bounds_data)
colnames(bounds_data) <- c('Index', 'BJ', 'GBJ', 'HC', 'GHC')
out_name <- paste('d', d, '_alpha', alpha*100, '_rho', rho*10, '_pct', percent_cor*100, '_bounds.txt', sep='')
write.table(bounds_data, out_name, append=F, quote=F, row.names=F, col.names=T, sep='\t')

