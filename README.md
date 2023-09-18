# KIG_v2
Updated version of the Known Ionomic Gene list using [OrthoFinder v2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1832-y). Visit the [KIG v1](https://github.com/baxterlab/KIG/) for the original publication using [InParanoid](https://www.sciencedirect.com/science/article/pii/S0022283600951970?via%3Dihub) to call orthologs.

# Known Ionomics Gene List pipeline (v2.0) 

This Rscript takes an input list of previously known ionomic genes primary genes and finds orthologs from other species and appends them to the list as inferred genes. Each species will have its own csv output file that contains all of its primary genes and genes that have been inferred as orthologs to other species' primary genes, along with each gene's corresponding orthologs in every other species. Current primary species, include Arabidopsis thaliana, Oryza sativa, Zea mays, Triticum aestivum, and Medicago truncatula. Note: T. aestivum was removed from analysis due to its genome annotation quality which impacts ortholog calling. Current inferred species include all the organisms listed in the primary list (except for T. aestivum), as well as Sorghum bicolor, Glycine max, Setaria viridis, Setaria italica, Populus trichocarpa, Vitis vinifera.

R version 4.1.0 (2021-05-18)

#### R Packages: The Ionomics Known Gene List script runs with these additional packages:
    1.1   plyr v1.8.8 and dplyr v1.1.2 (Wickham 2011 & Wickham et al. 2023)
    1.2   data.table v1.14.8 (Dowle & Srinivasan 2023)
    1.3   doParallel v1.0.17 (Microsoft & Weston 2017)
    1.4   reader v1.0.6 (Cooper 2017)
    1.5   stringr v1.5.0 (Wickham 2022)

#### Input list
The pipeline takes in the known ionome list of submissions from collaborators as csv input to start the pipeline. The file should have one gene per line, with no duplicates, and include information about the gene ID, the species, elements the gene was associated with, the tissue type of the analysis, and a citation of the study characterizing this gene/element interaction. Gene name, closest ortholog species, and additional comments are optional but, if included, will be copied along with the gene to its respective species table in the output.

### Ortholog table
Ortholog table was generated using OrthoFinder v2.0 (Emms and Kelly 2019) with protein fasta files obtained from [Phytozome v13](https://phytozome-next.jgi.doe.gov/phytomine/begin.do) (Goodstein et al. 2012) for all the species in the list (except T. aestivum), plus Liriodendron tulipifera as an appropriate outgroup. Best practices for using OrthoFinder to build your orthogroups have been stated by the authors [here](https://davidemms.github.io/orthofinder_tutorials/orthofinder-best-practices.html).

#### Function
The beginning of the pipeline reads in the primary list and sorts it by species and gene ID. White spaces are removed from the elements file to prepare for string manipulations later in the pipeline. The OrthoFinder ortholog table is read in and converted to a list of strings, where each row is pasted together into one big string - this is to make searching faster. The function get_orthogroups searches for each geneID and saves the entire row, which contains all the orthologs for each species. This function also removes the protein suffixes from each ID, which may need to be updated with new species additions, depending on their naming conventions. After the ortholog search has been executed for each gene in the list, the orthologs of all the genes are sorted into data frames by species. A list of all the orthologs for a gene is condensed into lists by species and appended to their species column in the list. For example, if you wanted to see all the G. max orthologs for AT1G01580, you would go to the A.thaliana_knownIonomicsGenesWOrthologs table, look at the row for AT1G01580 and the column "G. max orthologs". All of the species in the pipeline will have their primary and inferred genes combined into the same data table and saved to a csv with the suffix "_knownIonomicsGenesWOrthologs".

#### Output files
Each output file (species.name_knownIonomicsGenesWOrthologs.csv) will retain the same format as the primary list, except it will only include genes for its respective species, and it will have additional columns from the pipeline. After the original columns, you will find the generated columns "Primary/Inferred", "Inferred from", "Inferred.Elements" and ortholog columns for each species.
