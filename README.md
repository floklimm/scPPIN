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
1. Computation of differentially expressed genes p-values with an approach of your choice (e.g., `FindMarkers` function in [SEURAT](https://satijalab.org/seurat/))
2. Load a protein--protein interaction network (we here provide a PPIN for *Homo sapiens* that was constructed from [BioGRID](https://thebiogrid.org/)
3. Use the function `detectFunctionalModule(ppin,pValues,FDR)` to compute the functional module
4. Illustrate the detected modules with the function `plotFunctionalModule(functionalModule,FDR)`

In Step 3 all computations are executed:
A. Fitting of a beta-uniform model to the observed p-values,
B. Construction of a node-weighted graph, 
C. Rewriting of the maximum-weight subgraph problem as a prize-collecting Steiner tree problem, and
D. Writing input files for `dapcstp',
E. Solving of the prize-collecting Steiner tree problem by calling the dapcstp solver, and
F. Reading the solution file into R.

## Tutorial

The usage is demonstrated for two examples in the *tutorial_scPPIN.R*. In these tutorials the step 1 (computation of differentially expressed genes) replaced by a loading of pre-computed p-values. 

The tutorial also demonstrates the usage of `qgraph' for a nicer plotting of the functional modules and some helper functions (e.g., `fitBUM')


## License
This project is licensed under the AGPL - see the [LICENSE](https://github.com/floklimm/scPPIN/blob/master/LICENSE) file for details.
