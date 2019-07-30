# scPPIN

**There exists also a [webtool](https://floklimm.shinyapps.io/scPPIN-online/).**

This R library allows the computation of active modules in protein-protein interaction networks by integration with single-cell RNA sequencing data. The method is outlined in our manuscript

> Functional module detection through integration of single-cell RNA sequencing data with protein–protein interaction networks.  
> Florian Klimm, Enrique M. Toledo, Thomas Monfeuga, Fang Zhang, Charlotte M. Deane, and Gesine Reinert  
> bioRxiv 698647; doi: https://doi.org/10.1101/698647

The preprint is available on [bioRxiv](https://www.biorxiv.org/content/10.1101/698647v1).

## Dependencies
* R
* some standard R libraries:
    * igraph
    * qgraph
    * RColorBrewer
    * MASS
* dapcstp [(available on GitHub)](https://github.com/mluipersbeck/dapcstp) as discussed in
> A Dual Ascent-Based Branch-and-Bound Framework for the Prize-Collecting Steiner Tree and Related Problems  
> M. Leitner, I. Ljubic, M. Luipersbeck, M. Sinnl  
> INFORMS Journal on Computing 30(2):402-420, 2018

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
### Tutorial 1: Basic usage

The usage is demonstrated for two examples in *tutorial_scPPIN.R*. In these tutorials the steps 1 to 3 (preprocessing, cell cluster identification, and computation of differentially expressed genes) are replaced by a loading of pre-computed p-values.

For the first, small example the obtained functional module consists of three nodes (APP, SCD, and ALDOB). For the second example, the functional module consists of sixteen nodes.

The tutorial also demonstrates the usage of `qgraph` for a nicer plotting of the functional modules and some helper functions (e.g., `fitBUM`)

The function `detectFunctionalModule(ppin,pValues,FDR)` has also an optional argument `missingDataScore=TRUE`, which allows the computation of functional modules while keeping proteins without gene-expression information (shown as red boxes in the image below).

![alt text][ppinModule]

[ppinModule]: https://github.com/floklimm/scPPIN/blob/master/images/activeModuleExampleMissingData.png "Example functional module with missing gene-expression information"

### Tutorial 2: Step-by-step
In the script *tutorial_scPPIN-stepwise.R* the functionality is shown step-by-step. This might be helpful if user want to adapt some of the steps with their own function (e.g., a different choice of node scores). The result should be the same as in the real-world example in the first tutorial.

## FAQs

1. When loading the PPIN with the `load_ppin()` function I get an error *Can not open GraphML file*
> Most likely this is happening because R does not find the file. Make sure that you are in the correct working directory `./scPPIN-master/R`

2. I receive an error when executing the function `detectFunctionalModule`.
> This occurs often when `dapcstp` is not properly installed. Make sure that it works by calling it directly from the terminal.

3. I receive a segmentation fault (Core Dumped) when executing `dapcstp`.
> This should not happen. Most likely the input file given to `dapcstp` is not correctly formatted. Inspect the file and make sure that the file is correctly written (including End-of-File).

4. I want to run it on a different organims than *Homo sapiens*.
> In the folder `R/inst/extdata/morePPINs` you can find graphML files for all 68 organisms for which BioGRID data is available. Load them with the `biogridNetwork <- read_graph('./inst/extdata/morePPINs/ biogridSaccharomyces_cerevisiae_S288c3.5.169.tab2.txt.graphml', format='graphml')` command.

5. I don't like BioGRID and would rather use my own PPIN.
> You can construct your own network in the igraph format and use the provided functions. But it is important that the gene symbols are the same as the names of the nodes in the PPIN.

6. When using the `fitBUM` function I receive an error.
> This often occurs when the p-values are not in the half-open interval (0,1]. This means that p-values of zero are not allowed. (a first workaround would be to replace all zero p-values with the smallest of all non-zero p-values.)

7. When compyling `dapcstp`, I receive an error.
> Please report such errors to the [(creators)](https://github.com/mluipersbeck/dapcstp) of `dapcstp`.

## License
This project is licensed under the AGPL - see the [LICENSE](https://github.com/floklimm/scPPIN/blob/master/LICENSE) file for details.
