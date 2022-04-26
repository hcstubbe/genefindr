# Search for occureneces in sex biased genes

# Source functions
source('scripts/genefindr.R')

# Select abstract files
abstractFileNames = c("corona.txt")

names(abstractFileNames) = abstractFileNames

# Run genefindr for counting occurrences
runGeneFindr(abstractFileNames, 
             n.cores = 1, # to run on widonws
             chunks = 10, # number of chunks needs to be smaller or equal to number of investigated genes
             GEOsearchFile = "files/X12864_2013_5693_MOESM1_ESM.csv", # investigate only sex biased genes
             subset_data = TRUE
             )

# Count numbers of abstracts
lapply(abstractFileNames, count.abstracts)