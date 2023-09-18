# KIG_v2
Updated version of the Known Ionomic Gene list using [OrthoFinder v2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1832-y). Visit the [KIG v1](https://github.com/baxterlab/KIG/) for the original publication using [InParanoid](https://www.sciencedirect.com/science/article/pii/S0022283600951970?via%3Dihub) to call orthologs.

# Known Ionomics Gene List pipeline (v2.0) 

This Rscript takes an input list of previously known ionomic genes primary genes and finds orthologs from other species and appends them to the list as inferred genes. Each species will have its own csv output file that contains all of its primary genes and genes that have been inferred as orthologs to other species' primary genes, along with each gene's corresponding orthologs in every other species. Current primary species, include Arabidopsis thaliana, Oryza sativa, Zea mays, Triticum aestivum, and Medicago truncatula. Current inferred species include all the organisms listed in the primary list, as well as Sorghum bicolor, Glycine max, Setaria viridis, Setaria italica, Populus trichocarpa, Vitis vinifera.

R version 4.1.0 (2021-05-18)

#### R Packages: The Ionomics Known Gene List script runs with these additional packages:
    1.1   plyr v1.8.8 and dplyr v1.1.2 (Wickham 2011 & Wickham et al. 2023)
    1.2   data.table v1.14.8 (Dowle & Srinivasan 2023)
    1.3   doParallel v1.0.17 (Microsoft & Weston 2017)
    1.4   reader v1.0.6 (Cooper 2017)
    1.5   stringr v1.5.0 (Wickham 2022)
