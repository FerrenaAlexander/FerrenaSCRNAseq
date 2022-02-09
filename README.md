
<!-- README.md is generated from README.Rmd. Please edit that file -->

# FerrenaSCRNAseq

<!-- badges: start -->
<!-- badges: end -->

The goal of FerrenaSCRNAseq is to perform QC, processing, and analysis
of scRNAseq data.

## Installation

You can install from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("FerrenaAlexander/FerrenaSCRNAseq")
```

## Usage

This package works with Seurat objects.

The idea is to automatically filter out poor quality cells, defined as
cells with high mito content, low UMI, low “complexity” (lower genes
than expected given nUMI)

I suggest you process with SCT pipeline and run clustering on
un-filtered data. This allows cell-type specific QC filtering.
Sometimes, mito conent, number of UMI or number of genes is actually a
function of cell type. So globbal cutoffs suffer from lack of both
specificty and sensitivity, by throwing out cells that are actually
good, and keeping cells that are actually bad. So adjusting for cell
type is pretty important.

It works by multivariable linear models which allow for “outlier
diagnostics”.

After filtering, you should re-process the data.

    ### read in and pre-process
    #read in object
    rawh5 = '/path/to/h5/file/from/cellranger'
    samp = 'sample_label'
    sobj <- CreateSeuratObject(   Read10X_h5(rawh5), min.cells= 3, project = samp)

    #mito content, add to metadata
    mito.features <- grep(pattern = "^mt-", x = rownames(x = sobj), value = TRUE, ignore.case = T)
    sobj[["percent.mito"]] <- Seurat::PercentageFeatureSet(sobj, features = mito.features)


    #normalize and cluster
    suppressWarnings(sobj <- Seurat::SCTransform(sobj, verbose = T))

    sobj <- Seurat::RunPCA(object = sobj, verbose = F)

    sobj <- Seurat::FindNeighbors(object = sobj, dims = 1:20, verbose = F)
    sobj <- Seurat::FindClusters(object = sobj, resolution = 0.1, verbose = F, algorithm = 4)

    sobj <- RunUMAP(sobj, dims = 1:20)


    ### run auto filter ###

    # sobj is a suerat object
    # clusters refers to a column in the seurat@meta.data - here, we use the clsutering computed above.
    # iterative mito filter (cell-wise) is not as good as identifying and removing the mito cluster.
    # see ?FerrenaSCRNAseq::automatedfiltering

    reportlist <- FerrenaSCRNAseq::automatedfiltering(sobj, clusters = 'SCT_snn_res.0.1',
    iterativefilter.mito = F)



    #add autofilter results to metadata
    autofilterres <- reportlist[[1]]
    sobj$filteredout <- autofilterres$filteredout
    sobj$filterreason <- autofilterres$filterreason


    #make an output dir for the report
    outputdir = '.'
    dir.create( paste0(outputdir) )
    dir.create( paste0(outputdir, '/qc') )

    #plot autofilter results
    pdf(  paste0(outputdir, '/qc/autofilter.pdf'), 7,7)
    print( DimPlot(sobj, label = T) )


    print( FeaturePlot(sobj, c('nCount_RNA', 'nFeature_RNA', 'percent.mito'), order = T) + 
    DimPlot(sobj, group.by = 'filteredout')
    )

    print(reportlist)
    dev.off()

    #save autofilter output
    saveRDS( reportlist, paste0(outputdir, '/qc/reportlist-autofilter.rds') )

    #filter
    goodcells <- autofilterres[autofilterres$filteredout == 'No', 'barcodes']
    sobj <- sobj[,goodcells]

    #reprocess
    #normalize and cluster
    suppressWarnings(sobj <- Seurat::SCTransform(sobj, verbose = T))

    sobj <- Seurat::RunPCA(object = sobj, verbose = F)

    sobj <- Seurat::FindNeighbors(object = sobj, dims = 1:20, verbose = F)
    sobj <- Seurat::FindClusters(object = sobj, resolution = 0.1, verbose = F, algorithm = 4)

    sobj <- RunUMAP(sobj, dims = 1:20)

## DoubletFinder Wrapper

I use DouletFinder a lot, so I added a wrapper of DouletFinder

It assumes processing with SCT. It also uses an estimated doublet rate
from 10X genomics, see FerrenaSCRNAseq::dratedf This table is accurate
as of 2022 Feb 09.

    #use doublet filtering
    dfdf <- FerrenaSCRNAseq::doubletfinderwrapper(sobj, clusters = 'SCT_snn_res.0.1')
