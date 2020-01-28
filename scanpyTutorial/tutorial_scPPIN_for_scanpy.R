# example script that demonstrates the usage of scPPIN

# 1) Load the library
#library(scPPIN)


### A) Small Example
# 1) Load data
# the biogrid protein--protein interaction network
ppin <- loadPPIN()
# some example p-values
pValuesRaw <- read.csv("./data/PBMCpValuesCluster0.csv")
# rewrite them as a named vector of doubles
pValues <- pValuesRaw$pVal
names(pValues) <- pValuesRaw$gene
# if the gene symbols are not in upper case, use the `toupper` function
#names(pValues) <- toupper(pValuesRaw$gene)

# you should make sure that the p-values are in the semi-open interval (0,1], so this checksum should be zero
checksum = sum(pValues<0) + sum(pValues>1)
if (checksum>0){
  warning('p-values are not in (0,1]!')
}

# 2) Run our algorithms
FDR <- 10^{-2} # choose the false discovery rate
functionalModule <-detectFunctionalModule(ppin,pValues,FDR)


# 3) Plot the results
plotFunctionalModule(functionalModule,FDR)



