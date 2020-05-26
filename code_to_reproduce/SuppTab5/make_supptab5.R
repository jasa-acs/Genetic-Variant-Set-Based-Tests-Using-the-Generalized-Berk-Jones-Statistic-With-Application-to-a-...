# Adapted from cgem_rerun.R

# The purpose of this script is to check the approximation provided for the correlation
# structure of summary statistics in Sec 2.2.
# We compare the estimated correlation structure using (1) all the individual-level data from
# CGEMS and (2) data from the 1000 Genomes Project reference panel (CEU population).

# This code controls the generation of summary statistics as well as the correlation estimation
# using both methods. A number of intermediate files will be produced.

# Here the code refers to our attached fake dataset because we are not allowed to release the actual data.
# See the ACC Form for more details.
# If the fake dataset is replaced by the real dataset (and the corresponding files in the real dataset 
# are given the same names that we have given the fake data), then make_supptab5.R will 
# reproduce supplementary table 5 of the submitted manuscript (once all runs have finished).

# Make sure you have the necessary 1000 Genomes files (genotypes and eigenvectors, provided with this code)
# in your working directory.


########################################################################################
# Change this location to the directory that holds the necessary R libraries
# (must be installed first).
libs <- '/n/home/user/Rlibrary/3.3.1'

###################################################################
# No manual changes need to be made after this point
###################################################################

library(dplyr, lib=libs)
library(data.table, lib=libs)
library(magrittr, lib=libs)
library(GBJ, lib=libs)
########################################################################################


########################################################################################
# Two parameters - how many genes to test in every job?
# How much upstream and downstream buffer to put around the gene?
down_buffer <- 20000
up_buffer <- 20000
num_genes <- 5

# Read and order the covariates.
cgem_covar <- read.table('cgems_fake_covar.txt', header=T)
cgem_covar <- cgem_covar[order(cgem_covar[,1]), ]

# Read the gene list.
glist <- read.table('glist-hg18.txt', header=F)
colnames(glist) <- c('chr', 'start', 'end', 'gene')

# Order the gene list by chromosome and then base pair.
# Not necessary, but just for organizational purposes.
glist <- glist[order(glist[,1], glist[,2]),]
glist <- glist[which(glist$gene %in% c('FGFR2', 'CNGA3', 'PTCD3', 'POLR1A', 'ZNF263', 'VWA3B', 'TBK1', 'ABCA1', 'MMRN1', 'TIGD7')), ]

# Fit the null model once first.
Y <- cgem_covar$AFFECTED
null_X <- as.matrix( cbind(1, cgem_covar[,3:7], cgem_covar[,10:12]) )
null_mod <- glm(Y ~ null_X - 1, family=binomial(link='logit'))
null_mu <- null_mod$fitted.values

# W matrix and P matrix come from the null model and are the same for all genes.
W_vec <- null_mu * (1-null_mu)
W_mat <- diag(W_vec)
P_mat <- W_mat - W_mat %*% null_X %*% solve(t(null_X) %*% W_mat %*% null_X) %*% 
		t(null_X) %*% W_mat


##################################################################################################
# After we have extracted the relevant set of genotypes using PLINK using 1/2 recoding, 
# read from disk, parse, and then erase from disk.
# Don't need to do much error checking since we know we are only working with the 10 top genes
# that have already passed QC.
parse_snpset <- function(plink_fname)
{
	# Read in the ped file, order it to match covariates
	ped_fname <- paste( plink_fname, '.ped', sep='')
	ped_file <- read.table(ped_fname)
	ped_file <- ped_file[order(ped_file[,1]), ]
	
	# Loop over the SNPs to fill the genotype matrix
	num_snps <- (ncol(ped_file) - 6)/2
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

	# Read in the map file
    map_fname <- paste( plink_fname, '.map', sep='')
    map_file <- read.table(map_fname)
	colnames(map_file) <- c('Chr', 'RS', 'o1', 'BP')

	# Sort the map file and G_mat
	map_file_order <- order(map_file$BP)
	map_file <- map_file %>% arrange(BP)
	G_mat <- G_mat[, map_file_order]
	
	# Read in the bim file to get the allele information
	bim_file <- fread('cgems_genotypes.bim', header=F) %>%
		set_colnames(c('Chr', 'RS', 'o1', 'BP', 'A0', 'A1')) %>%
		select('RS', 'A0', 'A1') %>%
		filter(RS %in% map_file$RS)

	# Merge bim and map file, sort again
	map_file <- merge(map_file, bim_file, by='RS') %>% arrange(BP)

	# Flip everything to the minor allele
	freqs <- apply(G_mat, 2, mean) / 2
	to_flip <- which(freqs > 0.5)
	if (length(to_flip) > 0) {
		G_mat[,to_flip] <- 2 - G_mat[,to_flip]

		# Flip the map file
		temp_A0 <- map_file$A0
		map_file$A0[to_flip] <- map_file$A1[to_flip]
		map_file$A1[to_flip] <- temp_A0[to_flip]
	}

	# Remove MAF < 5% 
	freqs <- apply(G_mat, 2, mean) / 2
	MAF_remove <- which(freqs < 0.05)
	if (length(MAF_remove) > 0) {
		G_mat <- G_mat[,-MAF_remove]
		map_file <- map_file[-MAF_remove, ]
	}
	
	# Add MAF to map_file
	map_file <- map_file %>% mutate(MAF = apply(G_mat,2 , mean)/2)
	
	# Do we have SNPs with perfect correlation?  If so, remove.
	# This threshold could possibly depend on n.
	duplicate_snps <- which( duplicated(t(G_mat)) )
	if (length(duplicate_snps) > 0) {
		G_mat <- G_mat[,-duplicate_snps]
		map_file <- map_file[-duplicate_snps, ]
	}	
	
	return(list(G_mat=G_mat, map_file=map_file))
}


#################################################################################
# Loop through the gene list, save the summary statistics and correlation structures.

for (gene_it in 1:nrow(glist))
{
    # Start/end points for plink
	start_bp <- as.numeric(glist$start[gene_it]) - down_buffer
	end_bp <- as.numeric(glist$end[gene_it]) + up_buffer
	gene_name <- paste0('cgems', as.character(glist$gene[gene_it]))
	
	# Plink to recode as 12, important to as.character(chr)!
	plink_flag <- tryCatch(expr=system2(command='./plink', args=c('--bfile', 'cgems_genotypes', 
									'--recode12', '--chr', as.character(glist$chr[gene_it]), '--from-bp', 
									start_bp, '--to-bp', end_bp, '--noweb',
									'--out', gene_name), wait=TRUE), 
								warning=function(w) w, error=function(e) e)

	# Load the SNPs and map file.
	snp_output <- parse_snpset(plink_fname=gene_name)
	map_file <- snp_output$map_file
	G_mat <- snp_output$G_mat
	
	# Calculate the summary statistics and correlation matrix.
	score_stats_output <- calc_score_stats(null_model=null_mod, factor_matrix=G_mat, 
											link_function='logit', P_mat=P_mat)
	Z_stats <- score_stats_output$test_stats
	ss_file <- cbind(map_file, Z_stats)
	colnames(ss_file) <- c('RS', 'Chr', 'o1', 'BP', 'A0', 'A1', 'MAF', 'Z')	
	cor_mat <- score_stats_output$cor_mat
	sig_vec <- cor_mat[upper.tri(cor_mat)]	

	# Remove the plink files
	system2(command='rm', args=paste(gene_name, '.log', sep=''))
	system2(command='rm', args=paste(gene_name, '.map', sep=''))
	system2(command='rm', args=paste(gene_name, '.ped', sep=''))

	# Save the test statistics and correlation matrix	
	write.table(ss_file, paste0(gene_name, '_stats.txt'), append=F, quote=F, row.names=F, col.names=T)
	write.table(cor_mat, paste0(gene_name, '_cormat.txt'), append=F, quote=F, row.names=F, col.names=F)
	
	# Checkpointing
	cat(gene_it)
}



#################################################################################
# Next we "pretend" that we do not have the individual-level data and attempt to 
# estimate the correlation structure of the summary statistics using publicly available 
# data from the 1000 Genomes Project.

# Read in the precalculated eigenvectors
evecs_tab <- read.table('PCs_1000G_Euro.txt', header=T)
evec_cols_to_use <- paste("PC", 1:3, sep = "")
X_mat <- as.matrix(subset(evecs_tab, select = evec_cols_to_use))
X_mat <- cbind(1, X_mat)
W_mat <- diag(x = 1, nrow = nrow(evecs_tab), ncol = nrow(evecs_tab))
P_mat <- tryCatch(W_mat - X_mat %*% solve(t(X_mat) %*% X_mat) %*% 
                      t(X_mat), warning = function(w) w, error = function(e) e)


# Loop through the genes and estimate the correlation using reference data
for (gene_it in 1:nrow(glist)) {
    gene_name <- glist$gene[gene_it]
    start_bp <- glist$start[gene_it]
    end_bp <- glist$end[gene_it] 
    CHR <- glist$chr[gene_it]
                                                                 
    # Read 1000 Genomes ped file and order it to match eigenvectors 
    fname_root <- gene_name 
    ped_fname <- paste(fname_root, ".ped", sep = "")
    ped_file <- tryCatch(read.table(ped_fname, header = F), warning = function(w) w, 
                         error = function(e) e)
    ped_file <- ped_file[order(ped_file[, 1], decreasing = FALSE),  ]
    # Read the map file
    map_fname <- paste(fname_root, ".map", sep = "")
    map_file <- tryCatch(read.table(map_fname, header = F), warning = function(w) w, 
                         error = function(e) e)
    colnames(map_file) <- c("Chr", "RS", "Other", "BP")
    
    # Read previously calculated summary statistics for this gene
    ss_file <- read.table(paste0('cgems', gene_name, '_stats.txt'), header=T)
    
    # A few of the CGEMS SNPs may not be in the 1000 Genomes dataset, will have to remove 
    # these for comparison
    keep_map_idx <- which(map_file$RS %in% ss_file$RS)
    map_file <- map_file[keep_map_idx, ]
    keep_ped_cols <- sort(c(keep_map_idx*2 + 5, keep_map_idx*2 + 6))
    ped_file <- ped_file[, keep_ped_cols]
    
    # Important to make sure we are using the correct reference alleles in this
    # analysis as the original analysis, or else our correlations will be opposite signs.
    temp_Gmat <- matrix(data=NA, nrow=nrow(ped_file), ncol=nrow(map_file))
    flip_df <- data.frame(Orig=c('A', 'C', 'G', 'T'), Flip=c('T', 'G', 'C', 'A'))	
    bad_snps <- rep(NA, nrow(map_file))	
    for (snp_it in 1:nrow(map_file)) {	
        # Find the allele information from the original CGEMS analysis.
        temp_ss_row <- which(as.character(ss_file$RS) == as.character(map_file$RS[snp_it]))
        temp_A0 <- ss_file$A0[temp_ss_row]	
        temp_A1 <- ss_file$A1[temp_ss_row]
        temp_MAF <- ss_file$MAF[temp_ss_row]
        
        # Find the alleles used in the 1000 Genomes data.
        ped_column_1 <- snp_it*2 - 1
        ped_column_2 <- snp_it*2
        all_alleles_1000G <- c(as.character(ped_file[, ped_column_1]), as.character(ped_file[, ped_column_2]))
        unique_alleles_1000G <- unique(all_alleles_1000G)
        freq1 <- length(which(all_alleles_1000G == unique_alleles_1000G[1])) / length(all_alleles_1000G)
        freq2 <-  length(which(all_alleles_1000G == unique_alleles_1000G[2])) / length(all_alleles_1000G)
        # First try to set the effect alelle to the minor allele
        if (freq1 < freq2) {
            A0_1000G <- unique_alleles_1000G[2]; A1_1000G <- unique_alleles_1000G[1]
        } else {
            A0_1000G <- unique_alleles_1000G[1]; A1_1000G <- unique_alleles_1000G[2]		   
        }
        # If the minor alleles don't match between 1000 Genomes and CGEMS
        if (A0_1000G != temp_A0) {
            flipped_A0_1000G <- as.character(as.matrix(flip_df %>% filter(Orig == A0_1000G) %>% select(Flip)))
            flipped_A1_1000G <- as.character(as.matrix(flip_df %>% filter(Orig == A1_1000G) %>% select(Flip)))
            
            # If one or both of the 1000 Genomes alleles are not also CGEMS alleles
            if (length(which(c(as.character(temp_A0), as.character(temp_A1)) %in% c(A0_1000G, A1_1000G))) < 2) {
                # If not a simple strand flip
                if (flipped_A0_1000G != as.character(temp_A0) | flipped_A1_1000G != as.character(temp_A1)) {
                    # If a strand flip and minor allele flip
                    if (flipped_A0_1000G == as.character(temp_A1) & flipped_A1_1000G == as.character(temp_A0) & temp_MAF > 0.4) {	
                        temp_holder <- A1_1000G
                        A1_1000G <- A0_1000G
                        A0_1000G <- temp_holder
                    } else {bad_snps[snp_it] <- 1}   # Not a strand and/or minor allele flip, need to remove this one.
                }
            } else {		# If the alleles are correct, just swapped
                # If it not's possibly an ambigous flip and the MAFs are quite similar.
                if (flipped_A0_1000G != temp_A1 & temp_MAF > 0.4)  {
                    temp_holder <- A1_1000G
                    A1_1000G <- A0_1000G
                    A0_1000G <- temp_holder
                } else {	# Ambigous flip, can't tell, need to remove this one
                    bad_snps[snp_it] <- 1
                }
            }
        }
        
        # Now that the alleles are known, fill in the effect allele that matches the CGEMS effect allele 
        temp_snp1 <- rep(0, nrow(ped_file))
        temp_snp2 <- rep(0, nrow(ped_file))
        side1 <- which(ped_file[, ped_column_1] == A1_1000G)
        side2 <- which(ped_file[, ped_column_2] == A1_1000G)
        temp_snp1[side1] <- 1
        temp_snp2[side2] <- 1
        temp_Gmat[, snp_it] <- temp_snp1 + temp_snp2
    }
    
    # Remove the alleles where we can't be sure of the strand orientation.
    total_bad_snps <- length(which(bad_snps == 1))
    if (total_bad_snps > 0) {
        temp_Gmat <- temp_Gmat[, -which(bad_snps == 1)]
        map_file <- map_file[-which(bad_snps == 1), ]
    }	
    
    # Order the map file and temp_Gmat by BP, same as CGEMS analysis
    map_file_order <- order(map_file$BP)
    temp_Gmat <- temp_Gmat[, map_file_order]
    map_file <- map_file %>% arrange(BP)
    
    # Estimate the correlation matrix using the reference data
    cor_mat <- matrix(data = NA, nrow = ncol(temp_Gmat), ncol = ncol(temp_Gmat))
    denominators <- rep(NA, ncol(temp_Gmat))
    for (i in 1:ncol(temp_Gmat)) {
        temp_G <- temp_Gmat[, i]
        denominators[i] <- sqrt(t(temp_G) %*% P_mat %*% temp_G)
    }
    for (i in 2:ncol(temp_Gmat)) {
        for (j in 1:(i - 1)) {
            temp_G1 <- temp_Gmat[, i]
            temp_G2 <- temp_Gmat[, j]
            numerator <- t(temp_G1) %*% P_mat %*% temp_G2
            sig_hat <- numerator/(denominators[i] * denominators[j])
            cor_mat[i, j] <- sig_hat
            cor_mat[j, i] <- sig_hat
        }
    }
    
    # Write it
    write.table(cor_mat, paste0(gene_name, '_1000G_cor.txt'), append=F, quote=F, row.names=F, col.names=F)
    write.table(map_file, paste0(gene_name, '_1000G_map.txt'), append=F, quote=F, row.names=F, col.names=F)
    
  	cat(gene_it)
}


#################################################################################
# Next compare the two sets of estimated correlation matrices.

results_df <- data.frame(Gene=glist$gene, GBJ_1000G=NA, 
                         Frob=NA, Mat_l1=NA, med_abs=NA, med=NA)	
for (gene_it in 1:nrow(results_df)) {
    
    # Load 1000G data
    temp_gene <- results_df$Gene[gene_it]
    cor_1000G <- read.table(paste0(temp_gene, '_1000G_cor.txt'), header=F) %>%
        as.matrix(.)
    map_1000G <- read.table(paste0(temp_gene, '_1000G_map.txt'), header=F) %>%
        set_colnames(c('Chr', 'RS', 'o1', 'BP'))	
    
    # Load summary stats from CGEMS	
    cor_SS <- read.table(paste0('cgems', temp_gene, '_cormat.txt'), header=F) %>%
        as.matrix(.)
    map_SS <- read.table(paste0('cgems', temp_gene, '_stats.txt'), header=T)
    
    # Keep only CGEMS summary statistics that exist in 1000G
    keep_SS <- which(map_SS$RS %in% map_1000G$RS)
    cut_cor_SS <- cor_SS[keep_SS, keep_SS]
    cut_map_SS <- map_SS[keep_SS, ]
    
    # Run GBJ using the correlations estimated from 1000 Genomes.
    GBJ_1000G <- GBJ(test_stats=cut_map_SS$Z, pairwise_cors=cor_1000G[lower.tri(cor_1000G)])
    results_df$GBJ_1000G[gene_it] <- GBJ_1000G$GBJ_pvalue
   
    # Difference in correlation matrices norms
    diff_cor <- cor_1000G - cut_cor_SS
    diag(diff_cor) <- 0
    diff_cor_vec <- as.vector(diff_cor)
    frob_norm <- sqrt(sum(diff_cor_vec^2))
    matrix_l1_norm <- max(apply(abs(diff_cor), 2, sum))
    
    # Record
    results_df$Frob[gene_it] <- frob_norm
    results_df$Mat_l1[gene_it] <- matrix_l1_norm
    results_df$med_abs[gene_it] <- median(abs(diff_cor))
    results_df$med[gene_it] <- median(diff_cor)
    
    # Checkpoint
    cat(gene_it)
}

# Print
results_df