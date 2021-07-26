# scPPIN library for iGraph version for functional module detection
# If used please cite Klimm et al.
# "Functional module detection through integration of single-cell RNA sequencing data with protein–protein interaction networks"


# load necessary other libraries
library(igraph) # for handlign graphs/networks
library(qgraph) # for nicer plotting of the graphs
library(RColorBrewer) # for colourmaps for plotting
library(MASS) # for the fitting of the density function

loadPPIN <-function(){
  # loads biogrid protein--protein interaction network and returns it as a graphml
  biogridNetwork <- read_graph('./inst/extdata/biogridHomoSapiens3.5.166.graphml', format='graphml')
  biogridNetworkSimple <- simplify(biogridNetwork) # removing parallel edges and self-loops
  return(biogridNetwork)
}



plotFunctionalModule <-function(functionalModule,fdr,plotColorbar,layoutGraph,nodeScale){
  # function that plots the functional modules
  # input : module -- iGraph object
  #         fdr -- false discovery rate
  #         plotColorbar -- Bool (whether to plot the colorbar)
  if(missing(fdr)) {
    fdr = 0
  }
  # default is to plot the color bar
  if(missing(plotColorbar)) {
    plotColorbar = TRUE
  }
  # if no layout of the graph is given we compute a simple one
  if(missing(layoutGraph)) {
    layoutGraph <- layout_nicely(functionalModule, dim = 2)
  }

  # we might ahve to scale the size of the nodes and the labels
  if(missing(nodeScale)) {
    nodeScale <- 1
  }


  pValues <- log10(V(functionalModule)$pVal)
  pValues[is.nan(pValues)] <- 0  # for plotting purposes we set pValues that are uknown to 1

  # shape of proteins above FDR is a square
  if(fdr!=0){
    shapeVec <- ifelse(pValues < log10(fdr) , "circle" , "square") # above FDR is a square
    labelColorVec <- ifelse(pValues < log10(fdr) , "black" , "red") # above FDR has red font
    boldFontVec <- ifelse(pValues < log10(fdr) , 1 , 2) # above FDR has bold font
  } else {
    shapeVec <-"circle"
    labelColorVec <- "black"
    boldFontVec<- 1
  }

  # frame only for proteins that had no expression information
  fracmeVec <- ifelse(pValues == 0 , "red" , "NA")


  # colour of nodes is p-value
  pMax<- max(pValues) #max(log10(vertex_attr(functionalModule, "pVal")))
  pMin<- min(pValues) # min(log10(vertex_attr(functionalModule, "pVal")))
  pValBorders <- seq(from = pMin-1, to = pMax+0.0001,  length.out = 9)
  V(functionalModule)$pValDiscrete <-cut(log10(vertex_attr(functionalModule, "pVal")),pValBorders)


  pal <-  brewer.pal(length(pValBorders), "Purples")
  vColor <- pal[vertex_attr(functionalModule, "pValDiscrete")]

  # prepare the plotting
  if (plotColorbar){
    par(fig=c(0,0.8,0,0.9)) # if we have a colorbar we need space for it
  }
  else {
    par(fig=c(0,0.99,0,0.99)) # if not, we can use the whole area to plot the graph
  }


  # plot the graph
  plot(functionalModule,
       layout=layoutGraph,
      # vertex properties
      vertex.shape=shapeVec,
      vertex.color = vColor,vertex.label.color=labelColorVec,
      vertex.frame.color=fracmeVec,
      vertex.label.font=boldFontVec	,
      vertex.size	= nodeScale*12,
      vertex.label.cex=nodeScale*0.3,

      # edge properties
      edge.width=2,
      edge.curved=0.3,
      edge.color="black"
       )

  if (plotColorbar){
    # add a colorbar
    par(fig=c(0.8,0.99,0.10,0.70), new=TRUE) # location of bar
    colorbar(colorRampPalette(pal)(100)) # the actual plot
    # add labels to the colorbar
    mtext(toString(round(pMin, digits=0)), side=1) # bottom
    mtext(toString(round(pMax, digits=0)), side=3) # top
  }


}
#
colorbar <- function(colors) {
  count <- length(colors) # number of colours
  m <- matrix(1:count, 1,count) # create matrix with colour numbers
  image(m, col=colors, ylab="", axes=FALSE)
  mtext("log10(p-values)", side=2) # left
}


detectFunctionalModule <-function(network,pValues,FDR,missingDataScore=FALSE){
  # Input:
  #   network -- iGraph network
  #   pValues -- names list with pValues
  #   FDR -- False Discovery Rate

  if(missingDataScore==FALSE){
  # 1) intersection between the gene names in the network and the pValues
  geneSubset <- intersect(V(network)$name,names(pValues) )
  # 2) Construct node-weighted network
  networkFiltered <- induced_subgraph(network,  geneSubset)
  indexTemp = match(names(pValues),geneSubset)
  networkNodeWeighted <- set_vertex_attr(networkFiltered, "pVal",indexTemp[is.finite(indexTemp)], pValues[is.finite(indexTemp )] )
  # just to check that there is no vertex without a weight
  networkNodeWeighted <- delete_vertices(networkNodeWeighted, which(is.na(V(networkNodeWeighted)$pVal)))

  # 3) Run the analysis to find the maximal node-weighted spanning tree
  outputList <- FunctionalModuleComputation(networkNodeWeighted,FDR)
  activeModules <- outputList[[1]]
  activeModulesScores<- outputList[[2]]

  # 4) create a subgraph and return it
  activeSubgraph <- induced_subgraph(networkNodeWeighted, activeModules)
  activeSubgraph <- set_vertex_attr(activeSubgraph,"nodeScore", value = activeModulesScores)

  return(activeSubgraph)
  }
  else{
    # 1) intersection between the gene names in the network and the pValues
    geneSubset <- names(pValues)[names(pValues) %in% V(network)$name]
    # 2) Construct node-weighted network
    indexTemp = match(V(network)$name,names(pValues))
    V(network)$pVal <- NaN # this is a placeholder, these will be overwritten later for those with p-values
    #networkNodeWeighted <- set_vertex_attr(networkFiltered, "pVal",indexTemp[which(is.finite(indexTemp))], pValues[which(is.finite(indexTemp))])
    networkNodeWeighted <- set_vertex_attr(network, "pVal",which(is.finite(indexTemp)) ,  pValues[indexTemp[is.finite(indexTemp)]]     )

    # 3) Run the analysis to find the maximal node-weighted spanning tree
    outputList <- FunctionalModuleComputationMD(networkNodeWeighted,FDR)
    activeModules <- outputList[[1]]
    activeModulesScores<- outputList[[2]]

    # 4) create a subgraph and return it
    activeSubgraph <- induced_subgraph(networkNodeWeighted, activeModules)
    activeSubgraph <- set_vertex_attr(activeSubgraph,"nodeScore", value = activeModulesScores)

  }
}


FunctionalModuleComputationMD<- function(graph,falseDiscoveryRate,suboptimalCost){
  # detecting a functional module with missing expression data


  if(missing(suboptimalCost)) {
    suboptimalCost = Inf
  }

  N <- length(V(graph)) # number of nodes

  pValues <- V(graph)$pVal
  names(pValues) <- V(graph)$name


  # 1) Calculate the scores from the p-values (we exclude those that are 2)
  fitModel <- fitBUM(pValues[which(is.finite(pValues)) ], plot=FALSE)
  nodeScores<- computeNodeScores(graph,pValues,fitModel,fdr=falseDiscoveryRate,missingDataScore=TRUE)

  # 2) Shift the scores to formulate the MWST to a PCST
  minimumScore <- min(nodeScores)
  shiftedNodeScores <- nodeScores - minimumScore

  if (is.infinite(suboptimalCost)){
    # this is for the normal case of no suboptimal solutions
    # create a new graph object with the node scores and edge weights
    graphForSTCP <- graph %>%
      set_edge_attr("edgeWeight", value = -minimumScore) %>%
      set_vertex_attr("nodeWeight", value = shiftedNodeScores)

    activeModule <- calculatePrizeCollectingSteinerTree(graphForSTCP)
  }
  else {
    # if we want suboptimal solutions, we create a `root` node that connects with all nodes
    # first we create the normal graph
    graphForSTCP <- graph %>%
      set_edge_attr("edgeWeight", value = -minimumScore) %>%
      set_vertex_attr("nodeWeight", value = shiftedNodeScores) %>%
      add_vertices(N+1,name="rootnode",nodeWeight=0)  %>%
      add_edges(c(rbind(1:N, N+1)),edgeWeight=suboptimalCost)  # edges to the root node

    activeModule <- calculatePrizeCollectingSteinerTree(graphForSTCP)
    activeModule<- setdiff(activeModule,N+1) # remove the root node
  }


  returnList = list(activeModule,nodeScores[activeModule])

  return(returnList)
}


FunctionalModuleComputation<- function(graph,falseDiscoveryRate,suboptimalCost){
  # detecting a functional module in a monolayer network


  if(missing(suboptimalCost)) {
    suboptimalCost = Inf
  }

  N <- length(V(graph)) # number of nodes

  # 1) Calculate the scores from the p-values
  pValues <- V(graph)$pVal
  names(pValues) <- V(graph)$name

  fitModel <- fitBUM(pValues[which(is.finite(pValues)) ], plot=FALSE)
  nodeScores<- computeNodeScores(graph,pValues,fitModel,fdr=falseDiscoveryRate,missingDataScore=FALSE)

  # 2) Shift the scores to formulate the MWST to a PCST
  minimumScore <- min(nodeScores)
  shiftedNodeScores <- nodeScores - minimumScore

  if (is.infinite(suboptimalCost)){
    # this is for the normal case of no suboptimal solutions
  # create a new graph object with the node scores and edge weights
  graphForSTCP <- graph %>%
    set_edge_attr("edgeWeight", value = -minimumScore) %>%
    set_vertex_attr("nodeWeight", value = shiftedNodeScores)

  activeModule <- calculatePrizeCollectingSteinerTree(graphForSTCP)
  }
  else {
    # if we want suboptimal solutions, we create a `root` node that connects with all nodes
    # first we create the normal graph
    graphForSTCP <- graph %>%
      set_edge_attr("edgeWeight", value = -minimumScore) %>%
      set_vertex_attr("nodeWeight", value = shiftedNodeScores) %>%
      add_vertices(N+1,name="rootnode",nodeWeight=0)  %>%
      add_edges(c(rbind(1:N, N+1)),edgeWeight=suboptimalCost)  # edges to the root node

  activeModule <- calculatePrizeCollectingSteinerTree(graphForSTCP)
  activeModule<- setdiff(activeModule,N+1) # remove the root node
  }


  returnList = list(activeModule,nodeScores[activeModule])

  return(returnList)
}

calculatePrizeCollectingSteinerTree <- function(graph){
  # takes a graph in iGraph format
  # calls dapcstp to calculate a prize-collecting steiner tree
  # input: graph -  iGraph object with edge-weights and node-weigths

  fileName <- 'DAPCSTP_temp.stp'

  # 1) Writes an input file
  fileConn<-file(fileName,"w") # creates a file in write mode
  writeHeaderDAPCSTP(fileConn,"\"Prize-Collecting Steiner Problem in Graphs\"")
  writeGraphDAPCSTP(fileConn,graph)
  close(fileConn)
  # 2) Calls DAPCSTP
  # if it is a mac we have to call a different executable
  if (get_os()=="osx"){
    system("./dapcstpMAC DAPCSTP_temp.stp --type pcstp -o DAPCSTP_temp.sol")
    }
  else{
    system("./dapcstp DAPCSTP_temp.stp --type pcstp -o DAPCSTP_temp.sol")
  }
  # 3) reads DAPCSTP output file
  vertexOptimalSolution <- readSTPoutfileMWCS("DAPCSTP_temp.sol")

  return(vertexOptimalSolution)
}


writeHeaderDAPCSTP <- function(fileHandle,problemName){
  # writing the header for a STP file
  if(missing(problemName)) {
    problemName = "\"Maximum Node Weight Connected Subgraph\""
  }
  # the actual writing of the header lines
  writeLines(c("33D32945 STP File, STP Format Version 1.0",""), fileHandle)
  writeLines(c("SECTION Comments","Name \"Testfile\"","Creator \"Florian Klimm\"",paste("Problem ",problemName),"END",""), fileHandle)

  return(0)
}



writeGraphDAPCSTP <- function(fileHandle,graph){
  # write the graph section for a STP file
  write("SECTION Graph", fileHandle)
  # get the size of the graph
  N = length(V(graph)) # number of nodes
  E = gsize(graph) # number of edges
  # and write it
  write(paste("Nodes",N), fileHandle)
  write(paste("Edges",E), fileHandle)

  # writing the edges
  edgeList <- cbind('E', as_edgelist(graph, names = FALSE), E(graph)$edgeWeight) # edge list with edge weights
  write.table(edgeList,fileHandle, row.names=FALSE, col.names=FALSE, quote=FALSE)
  write("END", fileHandle)
  write("", fileHandle)

  # write the nodes
  write("SECTION Terminals", fileHandle)
  write(paste("Terminals",N), fileHandle)
  nodeIDs <- cbind('TP', as.numeric(V(graph)),V(graph)$nodeWeight)
  write.table(nodeIDs,fileHandle, row.names=FALSE, col.names=FALSE, quote=FALSE)
  write("END", fileHandle)
  write("", fileHandle)
  write("EOF", fileHandle)

  return(0)
}


readSTPoutfileMWCS <- function(fileName){
  # function that reads the output from dapcstp for solving a MWCS problem
  loadedText <- readLines(fileName)
  numNertexIndexSolution <- scan(fileName, skip = 11, nlines = 1, what = list("","") ) # read how many vertexes the optimal solution has
  vertexIndexSolution <- scan(fileName, skip = 12, nlines = numNertexIndexSolution[[2]], what = list("","") ) # read the vertex indexes
  return(strtoi(vertexIndexSolution[[2]]))
}


# fitting BUM

# the BUM density as a function
BUMdensity <- function(x,lambda,alpha){ lambda + (1-lambda)*alpha*x^(alpha-1) }

BUMdensityCummulative <- function(x,lambda,alpha){ lambda*x + (1-lambda)*x^(alpha) }


# fitting a a beta-uniform mixture to the observed, uncorrected p-values
fitBUM <- function(pValues, plot = FALSE){
  # fitting the distribution
  mlFit<-fitdistr(pValues, BUMdensity,start=list(lambda=0.50,alpha=0.50))

  # if wanted, we can also show a distirbution of the p-values and a plot
  if(plot){
    plot.new()
    par(mfrow=c(1,1))
    # plot a histogram of the p-values and add the estimations
    histogram <- hist(pValues, freq = FALSE,xlab="p-values",ylab="distribution",main="Histogram of p-values",breaks=50)
    maxDensity <- max(histogram$density)
    lines(density(pValues,from=0, to=1), col= "red",lwd=3) # density estimator
    # fitted BUM model
    curve(BUMdensity(x,coef(mlFit)['lambda'],coef(mlFit)['alpha']),add=TRUE, col="blue",lty=2,xlim=c(histogram$mids[1]/2,1),lwd=3)

    # value for x=1
    yInf <- coef(mlFit)['lambda'] +(1-coef(mlFit)['lambda'])*coef(mlFit)['alpha']

    segments(x0=0,y0=yInf,x1=1,y1=yInf, col="orange",lty=3,lwd=3) # plot this value
    legend(0.5, 1.9, legend=c("density estimation", "BUM model estimation"),col=c("red", "blue"), lty=1:2, cex=0.8,lwd=c(3,3))
    text(1.05,yInf,expression(tilde(f)),col="orange")

  }
  return(mlFit)
}

# give each node a weight
computeNodeScores <- function(network,pValues,bumFit,fdr,missingDataScore=FALSE){
# compute the node score for each protein
  # input:
  #       network -- an iGraph object with nodes named as proteins
  #       pValues -- named vector with all pValues
  #       bumFit -- the fit of the Beta-uniform mixture
  #       fdr -- false discovery rate
  #       missingDataScore -- if TRUE it sets the node score of proteins without p-value to zero

  # compute tau from the fit
  alpha <- coef(bumFit)['alpha']
  lambda <- coef(bumFit)['lambda']
  pihat <- lambda + (1 - lambda) * alpha
  tau <- ((pihat - (fdr * lambda))/(fdr * (1 - lambda)))^(1/(alpha - 1))

  # compute the node scores for all pValues
  nodeScores <- nodeScoreFunction(pValues,alpha,tau)

  # if we give missing data a zero scores for these proteins
  if(missingDataScore==FALSE){
    return(nodeScores)
  }
  else{
    # set node scores that don't have a p-value to a small negative number, e.g. -1, or the minimal value
    minimalScore <- median(nodeScores[is.finite(nodeScores)]  )
    nodeScores[which(is.nan(nodeScores))] <- -1
    return(nodeScores)
  }
}


nodeScoreFunction <- function(pValues,alpha,tau){
# compute the node score from the fit parameters and tau
# returns a named vector with the scores for every protein
  return((alpha - 1) * (log(pValues) - log(tau)))
}



#### Some further helper functions
get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}

