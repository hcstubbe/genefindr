# Version 20190928

# Read deepsearch/genefindr results
# handling of genefindr results
hiv <- read.csv("results/result_hiv 110422.csv", row.names = 1)
bcells <- read.csv("results/result_bcells 26902.csv", row.names = 1)
tcells <- read.csv("results/result_tcells 84002.csv", row.names = 1)

bool = c()
for(i in list(hiv, bcells, tcells)){
  id=identical(rownames(hiv), rownames(i))
  bool=c(bool, id)
}
bool = all(bool)
bool
if(!bool){stop("Rownames not matching")}

# make data frame with proper colnames and rownames. subset data fram on intersection
genefindrResult <- data.frame(bcells=bcells$occurences,
                  hiv=hiv$occurences,
                  tcells=tcells$occurences,
                  row.names = hiv$symbol)
# merge with additional symbols
genefindrResult = cbind(genefindrResult, hiv)

# Read GEO results
geoResult <- read.csv("files/result.strict.brisbane.basel.durham.csv")[,-1]

# Merge results
result <- merge(geoResult, genefindrResult, by.x="Gene.symbol", by.y = "symbol", all.x = TRUE)
write.csv(result, "reports/deepserachGeoCombinedSignificantStrict.csv")

result$mean.pval = apply(result[,grep("P.Value", colnames(result))], 1, mean)
result$mean.logFC = apply(result[,grep("logFC", colnames(result))], 1, mean)


# Subset result
resultFiltered <- subset(result, (bcells >=1 & tcells >=1) & hiv >=1)
resultFiltered <- resultFiltered[,c("Gene.symbol" , "mean.pval", "mean.logFC", "tcells", "bcells", "hiv")]
write.csv(resultFiltered, "reports/deepserachGeoCombinedFiltered.csv")