# scPPIN Tutorial 2
# this tutorial shows how to use the underlying functions in a more step-by-step manner

# 1) load library
library(scPPIN)


# 2) Load  data
# the biogrid protein--protein interaction network (not necessary if we loaded it before)
ppin <- loadPPIN()
# load p-values
pValues <- readRDS("./inst/extdata/pValuesClusterH1vsClusterH3.rds")

FDR <- 10^{-2}# choose the false discovery rate


geneSubset <- intersect(V(ppin)$name,names(pValues) )
# 2) Construct node-weighted network
networkFiltered <- induced_subgraph(ppin,  geneSubset)
indexTemp = match(names(pValues),geneSubset)
networkNodeWeighted <- set_vertex_attr(networkFiltered, "pVal",indexTemp[is.finite(indexTemp)], pValues[is.finite(indexTemp )] )
# just to check that there is no vertex without a weight
networkNodeWeighted <- delete_vertices(networkNodeWeighted, which(is.na(V(networkNodeWeighted)$pVal)))

# 3) Run the analysis to find the maximal node-weighted spanning tree

pValues <- V(networkNodeWeighted)$pVal
names(pValues) <- V(networkNodeWeighted)$name

# 4) Calculate the scores from the p-values (we exclude those that are 2)

fitModel <- fitBUM(pValues[which(is.finite(pValues)) ], plot=FALSE)
nodeScores<- computeNodeScores(networkNodeWeighted,pValues,fitModel,fdr=FDR,missingDataScore=TRUE)

# 5) Shift the scores to formulate the MWST to a PCST
minimumScore <- min(nodeScores)
shiftedNodeScores <- nodeScores - minimumScore

# 6) Write a graph
graphForSTCP <- networkNodeWeighted %>%
  set_edge_attr("edgeWeight", value = -minimumScore) %>%
  set_vertex_attr("nodeWeight", value = shiftedNodeScores)

activeModule <- calculatePrizeCollectingSteinerTree(graphForSTCP)


activeModules <- activeModule
activeModulesScores<- nodeScores[activeModule]

# 4) create a subgraph and return it
activeSubgraph <- induced_subgraph(networkNodeWeighted, activeModules)
activeSubgraph <- set_vertex_attr(activeSubgraph,"nodeScore", value = activeModulesScores)

# Plotting
plotFunctionalModule(activeSubgraph,FDR)
