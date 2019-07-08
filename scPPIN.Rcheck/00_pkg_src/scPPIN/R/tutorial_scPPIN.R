# example script that demonstrates the usage of scPPIN

# 1) Load the library
#library(scPPIN)


### A) Small Example
# 1) Load data
# the biogrid protein--protein interaction network
ppin <- loadPPIN()
# some example p-values
pValuesRaw <- read.csv("./inst/extdata/examplePvalues.csv")
# rewrite them as a named vector of doubles
pValues <- pValuesRaw$pVal
names(pValues) <- pValuesRaw$gene

# 2) Run our algorithms
FDR <- 10^{-2} # choose the false discovery rate
#functionalModule <-detectFunctionalModule(ppin,pValues,FDR)


# 3) Plot the results
#plotFunctionalModule(functionalModule,FDR)

### B) Larger, real-world example
# 1) Load  data
# the biogrid protein--protein interaction network (not necessary if we loaded it before)
ppin <- loadPPIN()
# load p-values
pValuesH1vsH3 <- readRDS("./inst/extdata/pValuesClusterH1vsClusterH3.rds")

FDR <- 10^{-2}# choose the false discovery rate

# 2) Computation of the maximum-weight spanning tree
functionalModuleH1H3  <- detectFunctionalModule(ppin,pValuesH1vsH3,FDR)
# 3) Plot the result
# a) using the normal potting function
plotFunctionalModule(functionalModuleH1H3,FDR)

# b) using qgraph for a nicer layout of the nodes
e <- get.edgelist(functionalModuleH1H3,names=FALSE)
computedNodePositions <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(functionalModuleH1H3))

plotFunctionalModule(functionalModuleH1H3,fdr=FDR,layoutGraph=computedNodePositions)

#4) We can also keep proteins without expression information
functionalModuleH1H3_missingData  <- detectFunctionalModule(ppin,pValuesH1vsH3,FDR,missingDataScore=TRUE)
e2 <- get.edgelist(functionalModuleH1H3_missingData,names=FALSE)
computedNodePositions2 <- qgraph.layout.fruchtermanreingold(e2,vcount=vcount(functionalModuleH1H3_missingData))
plotFunctionalModule(functionalModuleH1H3_missingData,fdr=FDR,layoutGraph=computedNodePositions2)
# 5) If you want to check the fit of the BUM model, you can use the following command with the plot option set to TRUE
fitBUM(pValuesH1vsH3,plot=TRUE)






