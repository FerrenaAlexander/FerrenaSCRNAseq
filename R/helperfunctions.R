# https://r-pkgs.org/whole-game.html

#' DoubletFinder Wrapper
#'
#' Models homotypic doublets using the following table:
#' https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
#'
#' @param seuratobject A Seurat object.
#' @param clusters string, should match a column name of seuratobject metadata with clusters or other cell group annotations. default = 'seurat_clusters'
#'
#' @return a data.frame with barcodes and DoubletFinder classifications
#' @export
#'
#' @examples
doubletfinderwrapper <- function(seuratobject, clusters){

  if( missing( clusters )) {clusters <- "seurat_clusters"}
  if( !(clusters %in% colnames(seuratobject@meta.data)) ) {stop('No clusters detected in Seurat object')}


  message('Running DoubletFinder paramsweep (may take a while)')

  #param sweep
  sweepres <- DoubletFinder::paramSweep_v3(seu = seuratobject, PCs = 1:30, sct = T)

  message('DF Parameter Sweep Completed')

  sweepstats <- DoubletFinder::summarizeSweep(sweepres)

  pdf(NULL) #prevent plotted output
  bcmvn <- DoubletFinder::find.pK(sweepstats)
  dev.off()

  maxscorepk <- bcmvn[bcmvn$BCmetric == max(bcmvn$BCmetric),2]
  maxscorepk <- as.numeric( levels(maxscorepk)[maxscorepk] )





  #homotypic doublet modelling

  ### using 10x table, use linear regression --> important for predicting homotypic / total doublet number
  #tamlabscpipeline::dratedf
  #dratedf <- read.delim('/Users/ferrenaa/Documents/tam/scripts/doublets/doubletrate.txt', header = T)
  #
  #   dratedf[,1] <- as.numeric(sub("%","",dratedf[,1]))/100
  #
  #   names(dratedf) <- c('MultipletRate', 'CellsLoaded_100%Viability', 'CellsRecovered')

  dbmodel <- lm(MultipletRate ~ CellsRecovered, data = dratedf)

  predicteddoubletrate <- as.numeric((dbmodel$coefficients['CellsRecovered'] * ncol(seuratobject)) + dbmodel$coefficients[1])

  #choose annotations to model homotypic doublets
  homotypicprop <- DoubletFinder::modelHomotypic(    as.vector(  seuratobject@meta.data[,clusters] )   )

  nexppoi <- round(predicteddoubletrate * length(rownames(seuratobject@meta.data)))
  nexppoiadj <- round(nexppoi * (1 - homotypicprop))

  #classify doublets
  message('Initiating third-pass clustering for doublet estimation:\n')

  seuratobject <- suppressWarnings(DoubletFinder::doubletFinder_v3(seu = seuratobject, PCs = 1:30, pN = 0.25, pK = maxscorepk, nExp = nexppoi, sct = T))

  ddf <- data.frame(cells = rownames(seuratobject@meta.data),
                    DoubletFinderClassification = seuratobject@meta.data[,ncol(seuratobject@meta.data)])

  ddf

}





#' Create an alluvial plot from long categorical data
#'
#' Wrapper around ggalluvium package for ggplot based alluvial plot, for quickly making alluvial plot from "long", "raw" categorical data (such as Seurat object meta.data), rather than two-way counts of categories.
#'
#'
#'
#'
#' @param labelsdf data.frame with two columns of raw categorical label: for example, each row is a cell (or other observation), and each column is metadata column 1 and metadata column 2
#' @param ggfittext T/F - whether to use ggfittext, to try to squeeze or remove tiny stratum labels
#'
#' @return a ggplot object
#' @export
#'
#' @examples
#' #' labelsdf can look like this, row.names of cells not needed, just two columns of categorical data as a data.frame:
#'                      From      To
#' AACCCAAGCATGCGA-1    Malignant  2
#' AAACCCAAGTAGGTTA-1    Malignant  2
#' AAACCCACAAAGCACG-1   Neutrophil  0
#' AAACCCACAGCAGTAG-1    Malignant  2
#' AAACCCACATACCGTA-1    Malignant  2
alluvialplot <- function(labelsdf, ggfittext){


  if( missing(ggfittext) ){ggfittext = F}

  #if levels not set, get them by ordering hi > lo

  labelsdf2 <- lapply(labelsdf, function(i){
    if( !is.factor(i) ){
      factor(i, levels = names(sort(table(i),decreasing = T)))
    } else{
      i
    }
  })

  labelsdf <- data.frame(labelsdf2, row.names = rownames(labelsdf))

  require(ggalluvial)

  #for ease, we'll set colnames to from and to
  # colnames(labelsdf)[1:2] <- c('From', 'To')
  # do this with .data trick now to keep colname!


  # turn it into a matrix
  mat <- table(labelsdf[,1], labelsdf[,2])


  #make the table long format
  longfreqs <- reshape2::melt(mat)
  colnames(longfreqs) <- c(colnames(labelsdf)[1:2], 'Freq')


  #factorize, using input levels or existing levels
  # inputting levels is mostly about order.


  longfreqs[,1] <- factor(longfreqs[,1], levels = levels(labelsdf[,1]) )
  longfreqs[,2] <- factor(longfreqs[,2], levels = levels(labelsdf[,2]))

  #if both are numerics, it seems to cause an issue, so convert to char vector...
  # if(is.numeric(longfreqs)[1] & is.numeric(longfreqs)[2])
  # can't reporduce that problem...


  if(ggfittext == T){

    require(ggfittext)
    ap <- ggplot(longfreqs, aes(y = Freq, axis1=.data[[colnames(longfreqs[1])]], axis2= .data[[colnames(longfreqs[2])]] ))+
      geom_alluvium(aes(fill= .data[[colnames(longfreqs[1])]] )) +
      geom_stratum()+
      #geom_label(stat = "stratum", aes(label = after_stat(stratum)))+
      ggfittext::geom_fit_text(stat = "stratum",aes(label = after_stat(stratum)), width = 1/4, min.size = 3) +
      theme_void()


  } else{


    ap <- ggplot(longfreqs, aes(y = Freq, axis1=.data[[colnames(longfreqs[1])]], axis2= .data[[colnames(longfreqs[2])]] ) )+
      geom_alluvium(aes(fill= .data[[colnames(longfreqs[1])]] )) +
      geom_stratum()+
      geom_label(stat = "stratum", aes(label = after_stat(stratum)))+
      # ggfittext::geom_fit_text(stat = "stratum",aes(label = after_stat(stratum)), width = 1/4, min.size = 3) +
      theme_void()


  }


  ap

}










#' Pseudobulk Seurat objects at whole sample or celltype level
#'
#' Pseudobulking is performed adding up gene expression values for each cell, for either cell types or whole sample.
#'
#' @param obj - seurat object, or a matrix
#' @param grouping_colname_in_md - optional, a string, the column name of `sobj@meta.data` (if obj is a seurat object) or metadata (if using matrix and metadata input) to use as "cell types" or any other grouping for pseudobulking at group level. if not provided, pseudobulk the entire matrix. default, not used.
#' @param metadata - data.frame with cell metadata, similar to `seuratobject@meta.data`. pass this only if
#' @param rawh5_path - optional, a string, the path to a rawH5 file, if provided will use the raw matrix subsetted by cells in sobj; if not will just pseudobulk the seurat object. useful if some filtering was applied to seurat object but you want to pseudobulk the whole matrix without that filtering, but with only using cells in seurat object. Default is not to use this.
#' @param assay - optional, a string, the name of the Seurat object assay to pull matrix from if rawh5_path is not provided. Default is "RNA" assay
#' @param slot - optional, a string, the name of the Seurat object slot within the designated object assay to pull matrix from if rawh5_path is not provided. Default is "counts" slot
#'
#' @return a data.frame. if grouping_colname_in_md is provided, each celltype will have a pseudobulked column, if not the data.frame will just be one column for the whole sample matrix.
#' @export
#'
#' @examples
pseudobulk <- function(obj, grouping_colname_in_md, metadata, rawh5_path, assay, slot){

  if(missing(assay)){assay = 'RNA'}
  if(missing(slot)){slot = 'counts'}

  require(Seurat)
  require(Matrix)

  #if rawh5_path is given, read in from raw data for all genes
  # if not, just use the seurat object as is

  #if grouping_colname_in_md is given, pseudobulk at celltype level
  # if not, pseudobulk whole object


  #get matrix and md


  if( any(grepl('Seurat', is(obj), ignore.case = T))  ){

    message('Seurat object detected')
    #md from seurat obj
    sobj <- obj
    md <- sobj@meta.data

    #mat: read in H5, or use
    if( !missing(rawh5_path) ){

      message('Reading raw matrix from:\n', rawh5_path)

      #read in raw mat
      mat  <- Read10X_h5(rawh5)

      mat <- mat[, match(colnames(sobj), colnames(mat)) ]

    } else{

      message('Using matrix from Seurat object:',
              '\n - Assay = ', assay,
              '\n - Slot = ', slot)

      mat <- Seurat::GetAssayData(sobj, assay=assay, slot=slot)
    }
  } else{
    message('Assuming input is matrix-like')

    mat <- obj

    if(!missing(grouping_colname_in_md)){
      if(missing(metadata)){stop('Please pass metadata dataframe to "metadata" argument')} else{md = metadata}
    }

  }




  #pseudobulk (at whole or celltype level)

  if( !missing(grouping_colname_in_md) ){
    message('For grouping, using metadata column: "', grouping_colname_in_md, '"')

    #get celltypes by order of number hi-->lo
    if( is.factor(md[,grouping_colname_in_md]) ){cts <- levels(md[,grouping_colname_in_md])} else{
      cts <- names( sort(table(as.vector(md[,grouping_colname_in_md])), decreasing = T) )
    }

    #for each cell type, pseudobulk
    dflist <- lapply(cts, function(ct){

      #subset md
      md_ct <- md[md[,grouping_colname_in_md]==ct,]

      #subset mat
      mat_ct <- mat[,match(rownames(md_ct), colnames(mat))]

      df <- data.frame(Matrix::rowSums(mat_ct))
      colnames(df) <- ct
      df
    })

    df <- dplyr::bind_cols(dflist)

  } else{
    message('Groupings not provided, will pseudobulk whole dataset')


    df <- data.frame(Matrix::rowSums(mat))
    colnames(df) <- NULL

  }


  df


}

