# scPPIN

**There exists also a [webtool](https://floklimm.shinyapps.io/scPPIN-online/).**

This R library allows the computation of active modules in protein-protein interaction networks by integration with single-cell RNA sequencing data. The method is outlined in our manuscript

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

**The solver dapcstp can be installed with the linked source code. We also provide pre-compiled binaries for Unix (Fedora 30) and Mac. It is likely, however, that you have to compile it for your system. If you name the executable differently, you have to change its call in the R function `calculatePrizeCollectingSteinerTree()`**

## Usage

The pipeline is as follows
1. Preprocessing of single-cell RNA-seq data
2. Detection of cell clusters e.g., `FindClusters` function in [SEURAT](https://satijalab.org/seurat/))
3. Computation of differentially expressed genes p-values with an approach of your choice (e.g., `FindAllMarkers` function in [SEURAT](https://satijalab.org/seurat/)) **Use the option `return.thresh=1` to obtain all p-values**
4. Load a protein--protein interaction network (we here provide a PPIN for *Homo sapiens* that was constructed from [BioGRID](https://thebiogrid.org/) and can be loaded with the `loadPPIN()` function)
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

The usage is demonstrated for two examples in *tutorial_scPPIN.R*. In these tutorials the steps 1 to 3 (preprocessing, cell cluster identification, and computation of differentially expressed genes) are replaced by a loading of pre-computed p-values.

For the first, small example the obtained functional module consists of three nodes (APP, SCD, and ALDOB). For the second example, the functional module consists of sixteen nodes.

The tutorial also demonstrates the usage of `qgraph` for a nicer plotting of the functional modules and some helper functions (e.g., `fitBUM`)

The function `functionalModuleH1H3_missingData` allows the computation of functional modules while keeping proteins without gene-expression information.

![alt text][ppinModule]

[ppinModule]: https://github.com/floklimm/scPPIN/blob/master/images/activeModuleExampleMissingData.png "Example functional module with missing gene-expression information"

## FAQ

1. When loading the PPIN with the `load_ppin()` function I get an error *Can not open GraphML file*
> Most likely this is happening because R does not find the file. Make sure that you are in the correct working directory `./scPPIN-master/R`

2. I receive an error when executing the function `detectFunctionalModule`.
> This occurs often when `dapcstp` is not properly installed. Make sure that it works by calling it directly from the terminal.

3. I receive a segmentation fault (Core Dumped) when executing `dapcstp`.
> This should not happen. Most likely the input file given to `dapcstp` is not correctly formatted. Inspect the file and make sure that the file is correctly written (including End-of-File).

## License
This project is licensed under the AGPL - see the [LICENSE](https://github.com/floklimm/scPPIN/blob/master/LICENSE) file for details.
