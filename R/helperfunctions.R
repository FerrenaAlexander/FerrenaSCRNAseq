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
#' labelsdf can look like this:
#'                      From      To
#' AACCCAAGCATGCGA-1    Malignant  2
#' AAACCCAAGTAGGTTA-1    Malignant  2
#' AAACCCACAAAGCACG-1   Neutrophil  0
#' AAACCCACACACTTAG-1 FILTERED_OUT  3
#' AAACCCACAGCAGTAG-1    Malignant  2
#' AAACCCACATACCGTA-1    Malignant  2
#'
#'
#' @param labelsdf data.frame with two columns of raw categorical label: for example, each row is a cell (or other observation), and each column is metadata column 1 and metadata column 2
#' @param fromlevels character vector. levels of labelsdf[,1] - can be used to set a desired order from first (top) to last (bottom)
#' @param tolevels character vector. levels of labelsdf[,1] - can be used to set a desired order from first (top) to last (bottom)
#' @param ggfittext T/F - whether to use ggfittext, to try to squeeze or remove tiny stratum labels
#'
#' @return a ggplot object
#' @export
#'
#' @examples
alluvialplot <- function(labelsdf, fromlevels, tolevels, ggfittext){

  if( missing(fromlevels) ) { fromlevels = levels(labelsdf[,1])}
  if( missing(tolevels) ) { tolevels = levels(labelsdf[,2])}
  if( missing(ggfittext) ){ggfittext = F}

  #if levels not set, get them by ordering hi > lo
  if(is.null(fromlevels)){fromlevels <- names(sort(table(labelsdf[,1]),decreasing = T)) }
  if(is.null(tolevels)){tolevels <- names(sort(table(labelsdf[,2]),decreasing = T)) }


  require(ggalluvial)

  #for ease, we'll set colnames to from and to
  colnames(labelsdf) <- c('From', 'To')


  # turn it into a matrix
  mat <- table(labelsdf[,1], labelsdf[,2])


  #make the table long format
  longfreqs <- reshape2::melt(mat)
  colnames(longfreqs) <- c('From', 'To', 'Freq')


  #factorize, using input levels or existing levels
  # inputting levels is mostly about order.


  longfreqs[,1] <- factor(longfreqs[,1], levels = fromlevels)
  longfreqs[,2] <- factor(longfreqs[,2], levels = tolevels)

  #if both are numerics, it seems to cause an issue, so convert to char vector...
  # if(is.numeric(longfreqs)[1] & is.numeric(longfreqs)[2])
  # can't reporduce that problem...


  if(ggfittext == T){

  require(ggfittext)
    ggplot(longfreqs, aes(y = Freq, axis1=From, axis2=To))+
      geom_alluvium(aes(fill=From))+
      geom_stratum()+
      #geom_label(stat = "stratum", aes(label = after_stat(stratum)))+
      ggfittext::geom_fit_text(stat = "stratum",aes(label = after_stat(stratum)), width = 1/4, min.size = 3) +
      theme_void()


  } else{


    ggplot(longfreqs, aes(y = Freq, axis1=From, axis2=To))+
      geom_alluvium(aes(fill=From))+
      geom_stratum()+
      geom_label(stat = "stratum", aes(label = after_stat(stratum)))+
      # ggfittext::geom_fit_text(stat = "stratum",aes(label = after_stat(stratum)), width = 1/4, min.size = 3) +
      theme_void()

  }


}



