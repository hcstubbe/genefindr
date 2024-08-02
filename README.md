[![DOI](https://zenodo.org/badge/485934720.svg)](https://zenodo.org/doi/10.5281/zenodo.13167120)

# genefindr

*Genefindr is currently converted into a Shiny app.*

## Example
 
``` r
# Select abstract files
abstractFileNames = c("hiv 110422.txt",
                      "bcells 28793.txt",
                      "tcells 84002.txt")
 
names(abstractFileNames) = abstractFileNames
 
# Run genefindr for counting occurrences
runGeneFindr(abstractFileNames, 
             n.cores = 1, # to run on widonws
             chunks = 10, # number of chunks needs to be smaller or equal to number of investigated genes
             GEOsearchFile = "files/result.strict.brisbane.basel.durham.csv" # investigate only sex biased genes
             )

# Count numbers of abstracts
lapply(abstractFileNames, count.abstracts)
```
