# Adapted from gBJ_power_sim_HAPGEN.R

# The purpose of this script is to run the simulations that are displayed in Supplementary Figures 1 and 2
# of the submitted manuscript.
# We want to find the accuracy of the GBJ p-value calculation for genes on chromosome 5.
# The code is very similar to that of gen_data_fig4.R, except here we focus on four random genes
# and do a large number of simulations for just those genes.
# Also we only test with GBJ, not the other methods, and we only test under the null.

# We generate the data using HAPGEN, simulate outcomes according to the model given in the 
# paper, and apply tests to the simluated data.

# Each job should produce one file with a name like 'chr5_pval_acc_S1_aID1.txt'
# where the numbers after Snum and aID will vary with each job.
# The parameters for each run are controlled by the aID and Snum variables, which are parameters passed
# in on the command line. 
# Our job submission scripts (*.sbatch and *.bsub) control these parameters
# so that all necessary files are generated.

# After all of these jobs have output their data, use plot_suppfig1.R to create plots similar to Supplementary Figures 1 and 2.
# Since HAPGEN sets a seed for its random number generator according to the time of day, 
# it is not possible to recreate the exact simulations that led to Supplementary Figures 1 and 2.
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

library(dplyr, lib=libs)
library(data.table, lib=libs)
library(GBJ, lib=libs)
########################################################################################


########################################################################################
# Read input arguments (fed automatically by job submission script) to determine simulation parameters 
args <- commandArgs(trailingOnly=TRUE)
aID <- as.numeric(args[1])
Snum <- as.numeric(args[2])
set.seed(05012017 + Snum*1000 + aID)

# load the legend file
leg_file <- fread('hapmap3.r2.b36.chr5.legend', header=T)

# figure out which gene we are going to use
if (Snum == 1) {                    # just recording the correlation structure
    runs_per_gene <- 1
    runs_per_sim <- 500
} else if (Snum == 2) {             # performing the simulation
    runs_per_gene <- 1000
    runs_per_sim <- 1000
}

# Parameters for the run
num_subjects <- 2000
param_row <- ceiling(aID / runs_per_gene)
param_mat <- data.frame(start_row=c(1, 1, 75878,83510,22969,65400), 
                        end_row=c(1, 1, 75897,83520,22987,65407))
gene_name_list <- c('FGFR2_low', 'FGFR2_high', 'RNF145', 'BNIP1', 'SEPP1', 'HARS2')
keep_snps_list <- list(c(1), c(1), 2:20, c(1,3,4,5,6,7,8,10,11),  c(3,4,8,10,12,14,19), c(2,3,4,8))
start_row <- param_mat$start_row[param_row]
end_row <- param_mat$end_row[param_row]
keep_snps <- keep_snps_list[param_row][[1]]


###########################################################################
# Function to generate FGFR2 gene of n people.
# Return a low LD subset and a high LD subset.
# Must have chr10 legend/map/hap files along with hapgen2 in directory.
gen_FGFR2 <- function(num_to_generate, aID, Snum, low=TRUE)
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
    if (low) {
        return(low_LD_haps)
    } else {
        return(high_LD_haps)
    }
}
###########################################################################


########################################################################################
# Function to generate genotypes of Chr5 genes for n people.
# Must have chr5 legend/map/hap files along with hapgen2 in directory.
# Pass it the start row, end row, and leg_file.
# Remove the bad MAF SNPs, add together 2 haplotypes to make a 'person', transpose to get
# SNPs on the columns.
generate_hapgen <- function(n, aID, Snum, leg_file, start_row, end_row, keep_snps) {
    
    # HAPGEN output name, make it different for each Snum and each aID so don't overwrite each other
    out_name <- paste('S', Snum, '.ch', aID, '.out', sep='')
    
    # Find start BP and end BP
    start_bp <- leg_file$position[start_row]
    end_bp <- leg_file$position[end_row]
    
    # Disease locus is just the one after the start position
    dl <- leg_file$position[start_row]
    
    # Map, legend, and haplotype (CEU only) file name
    map_fname <- paste('genetic_map_chr5_combined_b36.txt', sep='')
    leg_fname <- paste('hapmap3.r2.b36.chr5.legend', sep='')
    CEUhap_fname <- paste('CEU.chr5.hap', sep='')
    
    # Create the haps, return b=0 for success, b=1 for failure (no suitable disease locus)
    b <- 1
    counter <- 0
    while(b==1) {
        if ( (start_row+counter) > end_row) {
            remove_name <- paste('S', Snum, '.ch', aID, '.out.*', sep='')
            system2("rm", args=remove_name)
            return( list(G_mat=-1, kept_snps=-1) )
        }
        
        # generating n 'controls' and 2 'cases'
        b <- system2("./hapgen2", args=c("-m", map_fname, "-l", leg_fname, "-h", CEUhap_fname, "-o", out_name, "-dl", dl, "1", "1", "1", "-n", n, "2", "-int", start_bp, end_bp, "-no_gens_output"))
        # move the disease locus
        counter <- counter + 1
        dl <- leg_file$position[start_row+counter]
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
    
    # Keep same SNPs every time to ensure consistency across chr 5 random gene simulations
    combined_haps <- combined_haps[keep_snps, ]
    
    # Transpose and return
    combined_haps <- t(combined_haps)
    return( G_mat=combined_haps )
}
########################################################################################



########################################################################################
# Start running the simulation
gbj_results <- data.frame(gbj=rep(NA, runs_per_sim), gbj_p=NA, gbj_err_code=NA,	
                          start_bp=NA, end_bp=NA,  size=NA, Gene=NA)
for (sim_it in 1:runs_per_sim) {
    # Generate SNPs with HAPGEN
    if (param_row == 1) {
        gen_Gmat <- gen_FGFR2(num_to_generate=num_subjects, aID=aID, Snum=Snum, low=TRUE)
    } else if (param_row == 2) {
        gen_Gmat <- gen_FGFR2(num_to_generate=num_subjects, aID=aID, Snum=Snum, low=FALSE)
    } else {
        gen_Gmat <- generate_hapgen(n=num_subjects, aID=aID, Snum=Snum,
                                    leg_file=leg_file, start_row=start_row,
                                    end_row=end_row, keep_snps=keep_snps)
    }
   
    
    # Generate outcome and score statistics
    Y <- rnorm(n=num_subjects)
    null_mod <- glm(Y ~ 1)
    score_stats_output <- calc_score_stats(null_model=null_mod, factor_matrix=gen_Gmat,
                                           link_function='linear')
    
    # save the correlation matrix or gbj p-value
    if (Snum == 1 ) {
        if (sim_it == 1) {
            cor_mat <- score_stats_output$cor_mat / runs_per_sim
        } else {
            cor_mat <- cor_mat + score_stats_output$cor_mat / runs_per_sim
        }
    } else if (Snum == 2) {
        GBJ_output <- GBJ(test_stats=score_stats_output$test_stats, cor_mat=score_stats_output$cor_mat)
        gbj_results$gbj[sim_it] <- GBJ_output$GBJ
        gbj_results$gbj_p[sim_it] <- GBJ_output$GBJ_pvalue
        gbj_results$gbj_err_code[sim_it] <- GBJ_output$err_code	
        gbj_results$start_bp[sim_it] <- leg_file$position[start_row]
        gbj_results$end_bp[sim_it] <- leg_file$position[end_row]
        gbj_results$size[sim_it] <- length(score_stats_output$test_stats)
        gbj_results$Gene[sim_it] <- gene_name_list[param_row]
    }
    
    cat(sim_it, '-', nrow(score_stats_output$cor_mat), '\n')
}

# Record results
if (Snum == 1 ) {
    out_name <- paste0('cor_mat_S', Snum, '_aID', aID, '.txt')
    write.table(cor_mat, out_name, append=F, quote=F, row.name=F, col.names=F)
} else if (Snum == 2) {
    out_name <- paste0('chr5_pval_acc_S', Snum, '_aID', aID, '.txt')                                                                                                      
    write.table(gbj_results, out_name, append=F, quote=F, row.name=F, col.names=T)  
}



















