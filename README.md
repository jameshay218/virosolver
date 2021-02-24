# A method to infer epidemic dynamics using viral load data: virosolver
------------
Documentation for this package is a work in progress, but all code is working correctly. This particular version accompanies the git repository at https://github.com/jameshay218/virosolver_paper, where use cases and instructions are provided.
 
## Setup
------------
This package uses the `lazymcmc` R package, which is used for the MCMC procedure. This is easy to do with `devtools::install_github("jameshay218/lazymcmc")`. *However*, for many of the analyses, a separate branch implementing parallel tempering is needed. This can be installed using `devtools::install_github("jameshay218/lazymcmc@parallel_tempering")`.
  
A number of generic R packages are also used throughout:
```r
c("tidyverse","ggthemes","ggpubr","ggsci","data.table","patchwork",
"fitdistrplus","deSolve","lazymcmc","odin","doParallel","coda")
```
