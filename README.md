# scPPIN
Functional module detection through integration of single-cell RNA sequencing data with protein–protein interaction networks

This R library allows the computation of active modules in protein-protein interaction networks. The method is outlined in our manuscript

> Functional module detection through integration of single-cell RNA sequencing data with protein–protein interaction networks.
>
> Florian Klimm, Enrique M. Toledo, Thomas Monfeuga, Fang Zhang, Charlotte M. Deane, and Gesine Reinert

The article is available HERE (link to be included)

## Dependencies
* R
* some standard R libraries:
    * igraph
    * qgraph
    * RColorBrewer
    * MASS
    * jsonlite
* dapcstp [(available on GitHub)](https://github.com/mluipersbeck/dapcstp)

**The solver dapcstp can be installed with the linked source code. We also provide pre-compiled binaries for Unix (Fedora 30) and Mac. It is likely, however, that you have to compile it for your system.**

## Usage

The pipeline is as follows
1. Preprocessing of single-cell RNA-seq data
2. Detection of cell clusters e.g., `FindClusters` function in [SEURAT](https://satijalab.org/seurat/))
3. Computation of differentially expressed genes p-values with an approach of your choice (e.g., `FindAllMarkers` function in [SEURAT](https://satijalab.org/seurat/)) **Use the option `return.thresh=1` to obtain all p-values**
4. Load a protein--protein interaction network (we here provide a PPIN for *Homo sapiens* that was constructed from [BioGRID](https://thebiogrid.org/) and can be loaded with the `loadPPIN()' function)
5. Use the function `detectFunctionalModule(ppin,pValues,FDR)` to compute the functional module
6. Illustrate the detected modules with the function `plotFunctionalModule(functionalModule,FDR)`

In Step 4 all computations are executed:
* Fitting of a beta-uniform model to the observed p-values,
* Construction of a node-weighted graph, 
* Rewriting of the maximum-weight subgraph problem as a prize-collecting Steiner tree problem,
* Writing input files for `dapcstp`,
* Solving of the prize-collecting Steiner tree problem by calling the `dapcstp` solver, and
* Reading the solution file into R.

All networks and modules (which are subnetworks) are igraph objects.

## Tutorial

The usage is demonstrated for two examples in *tutorial_scPPIN.R*. In these tutorials the steps 1 to 3 (preprocessing, cluster identification, and computation of differentially expressed genes) replaced by a loading of pre-computed p-values. 

The tutorial also demonstrates the usage of `qgraph` for a nicer plotting of the functional modules and some helper functions (e.g., `fitBUM`)


## License
This project is licensed under the AGPL - see the [LICENSE](https://github.com/floklimm/scPPIN/blob/master/LICENSE) file for details.
