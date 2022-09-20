
<!-- README.md is generated from README.Rmd. Please edit that file -->

# FerrenaSCRNAseq

<!-- badges: start -->
<!-- badges: end -->

The goal of FerrenaSCRNAseq is to aid with QC, preprocessing, and analysis
of scRNAseq data.

## Installation

You can install from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("FerrenaAlexander/FerrenaSCRNAseq")
```

## Usage

This package works with Seurat objects.

The core function is called `FerrenaSCRNAseq::autofilter()`. This function selects poor quality cells that are outliers relative to the rest of the data.

* baseline cutoffs of minimum number of UMIs (default is 1000); minimum number of genes (default is 200); maximum percent mitochondrial content (max 25%); max percent hemoglobin content (max 25%)

* complexity analysis, default as number of genes / number of UMIs. The pipeline detects abnormally low complexity genes, ie those with lower number of unique genes expected given the number of UMIs captured. Usually, this is either poor quality cells or extreme cells like RBCs.

* data driven mitochondrial cutoffs: based on the distribtuion of mitochondrial gene content per cell, learn an upper cutoff. It uses median + (median absolute deviation * 2.5) by default.

* data driven nUMI cutoff: based on the distribtuion of UMIs per cell, learn an lower cutoff. It uses median - (median absolute deviation * 2.5) by default.

See `?FerrenaSCRNAseq::automatedfiltering()`


    ### quick run ###
    # sobj is a seurat object
    
    # identify outliers
    af <- autofilter(sobj)

    # get high quality cells
    goodcells <- af$cellstatus[af$cellstatus$filteredout==F,"barcodes"]

    # filter the seurat object for the high quality cells
    sobj <- sobj[,goodcells]
    
    
I recommend doing this as a first step, before other filtering.
After removal of poor quality outlier cells, you should then do (or re-do) any processing: HVG detection, PCA, graph building, clustering etc.
    
    
## View details of filtering

    ### calculate % mito content and % hemoglobin
    # good to do it beforehand, but the pipeline does it internally if you don't
    
    #percent.mito
    mito.features <- grep(pattern = "^mt-", x = rownames(x = sobj), value = TRUE, ignore.case = T)
    sobj[["percent.mito"]] <- Seurat::PercentageFeatureSet(sobj, features = mito.features)
    
    #percent.hemoglobin, by default checks all mouse/human hemoglobin genes
    sobj[["percent.hemoglobin"]] <- FerrenaSCRNAseq::calculate_percent.hemoglobin(sobj)
    
    
    #detect the outliers
    af <- autofilter(sobj)
    
    #see the classification as outlier or not for each cell, 
    # it also shows the reason each cell may have been called an outlier
    head( af$cellstatus )
    
    # check the complexity filtering step
    af$globalfilter.complexity
    
    #check lower lib size cutoff step
    # with the baseline UMI cutoff, this may not remove many or even remove none at all
    af$globalfilter.libsize
    
    # check the percent mito step
    af$globalfilter.mito
    
    #filter out the outliers
    goodcells <- af$cellstatus[af$cellstatus$filteredout==F,"barcodes"]
    sobj <- sobj[,goodcells]



## DoubletFinder Wrapper for automated doublet calling

I use DoubletFinder a lot, so I added a wrapper of DoubletFinder

-   DoubletFinder paper:
    <https://www.sciencedirect.com/science/article/pii/S2405471219300730>

-   DoubletFinder github
    <https://github.com/chris-mcginnis-ucsf/DoubletFinder>

It assumes processing with SCT. It also uses an estimated doublet rate
from 10X genomics,

see `?FerrenaSCRNAseq::dratedf` and
(<https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled->)

This table is accurate as of 2022 Feb 09.


    ## preprocess with SCT and cluster ##
    suppressWarnings(sobj <- Seurat::SCTransform(sobj, verbose = T, method="glmGamPoi"))

    sobj <- Seurat::RunPCA(object = sobj, verbose = F)
    sobj <- Seurat::FindNeighbors(object = sobj, dims = 1:30, verbose = F)
    sobj <- Seurat::FindClusters(object = sobj, resolution = 0.1, verbose = F, algorithm = 1)


    # get the doubletfinder dataframe
    # num.cores, parallelization is now supported by DF
    dfdf <- doubletfinderwrapper(sobj, #autofilterres = af,
                                        num.cores = 5)
                                        
                                        
    # alternatively, you can also pass the autofilter result, 
    # it will then return an updated autofilter result list including the doubletfinder results
    af <- doubletfinderwrapper(sobj, autofilterres = af,
                           num.cores = 5)
                           

    # filter out doublets
    # from doubletfinder data.frame, if autofilter res was not passed above
    singlets <- dfdf[dfdf$DoubletFinderClassification=='Singlet','cells']
    sobj <- sobj[,singlets]
    
    # or, filter using updated autofilter result object
    goodcells <- af$cellstatus[af$cellstatus$filteredout==F,"barcodes"]
    sobj <- sobj[,goodcells]
    
    

After removal of doublets, I recommend re-processing the data (HVG detection, PCA, graph, clustering).
