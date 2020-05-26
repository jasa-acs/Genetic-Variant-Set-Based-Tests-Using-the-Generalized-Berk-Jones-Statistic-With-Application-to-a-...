# The purpose of this script is to recreate Table 2 of the submitted manuscript.

# We read in the output produced by gen_data_tab2.R and format the data 
# into tables.
# Simply view the variable 'concise_results' to see the entire table.

# Here we have attached a fake dataset because we are not allowed to release the actual data.
# See the ACC Form for more details.
# If the fake dataset is replaced by the real dataset (and the corresponding files in the real dataset 
# are given the same names that we have given the fake data) then make_tab2.R will 
# reproduce Table 2 of the submitted manuscript (once all runs have finished).

# Obviously running this code on the fake dataset will not reproduce our results.
# In fact, running our code on the fake dataset will produce mostly NA for the results, because the fake
# dataset we have attached only contains genotypes at a small number of markers (in the interest of 
# not sending overly large files), so most genes will have no data.
# To see more numerical results, one can edit the glist-hg18.txt file to define more fake genes at the
# markers we have included (chromosome 10, between 121001183-121876520).

########################################################################################
# Change data_dir your working directory directory which holds all the necessary data files.
data_dir <- '/users/user/desktop'

###################################################################
# No manual changes need to be made after this point
###################################################################

##################################################################
##################################################################
# Read in data and make Table 2
cgems_results <- read.table('cgems_results_S1_1.txt', header=T)
for(i in 2:181)
{
	temp_fname <- paste('cgems_results_S1_', i, '.txt', sep='')
	temp_file <- read.table(temp_fname, header=T)
	cgems_results <- rbind(cgems_results, temp_file)
}

# Output version you see in the main manuscript
concise_results <- cgems_results[, c(1,8,10,12,13,15,16)]
smallest_pvalue <- apply(concise_results[, 2:5], 1, min)
concise_results <- concise_results[order(smallest_pvalue), ]
concise_results

