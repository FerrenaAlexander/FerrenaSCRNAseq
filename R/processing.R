# https://r-pkgs.org/whole-game.html


#' Perform data-driven filtering of scRNAseq data
#'
#' Apply simple cutoffs and discover data-driven thresholds for poor quality cells in scRNAseq.
#'
#' Simple cutoffs include minimum number of UMIs, minimum number of unique genes detected, maximum percent mito, and maximum percent hemoglobin.
#' More complex cutoffs are learnt for lower than expected complexity (defined for each cell as num unique genes / num UMIs). Additionally, median absolute deviation is used to exclude remaining cells with high mito content or low UMI content.
#'
#' @param sobj seurat object
#' @param min_num_UMI numeric, default is 1000, if no filter is desired set to -Inf
#' @param min_num_Feature numeric, default is 2, if no filter is desired set to -Inf
#' @param max_perc_mito numeric, default is 25, if no filter is desired set to Inf
#' @param max_perc_hemoglobin numeric, default is 25, if no filter is desired set to Inf
#' @param mad.score.threshold numeric, default is 2.5, threshold for median abs deviation thresholding, ie cutoffs set to `median +/- mad * threshold`
#' @param globalfilter.complexity T/F, default T, whether to filter cells with lower than expected number of genes given number of UMIs
#' @param globalfilter.mito T/F, default T, whether to filter cells with higher than normal mito content
#' @param globalfilter.libsize T/F, default T, whether to filter cells with lower than normal UMI content
#'
#' @return
#' @export
#'
#' @examples
autofilter <- function(
    sobj,

    min_num_UMI,
    min_num_Feature,
    max_perc_mito,
    max_perc_hemoglobin,

    mad.score.threshold,

    globalfilter.complexity,
    globalfilter.mito,
    globalfilter.libsize

){

  #check if seurat object is there
  if( missing( sobj )) {stop('please input a Seurat object')}

  #set the basic cutoffs
  if( missing(min_num_UMI ) ){min_num_UMI <- 1000}
  if( missing(min_num_Feature ) ){min_num_Feature <- 200}
  if( missing(max_perc_mito ) ){max_perc_mito <- 25}
  if( missing(max_perc_hemoglobin ) ){max_perc_hemoglobin <- 25}


  #set maximum distance of deviations from median tolerated before outlier classification
  if( missing( mad.score.threshold )) {mad.score.threshold <- 2.5}

  # baseline (global) filtration
  if( missing( globalfilter.complexity )) {globalfilter.complexity <- T}
  if( missing( globalfilter.mito )) {globalfilter.mito <- T}
  if( missing( globalfilter.libsize )) {globalfilter.libsize <- T}



  #output list object; add to this as needed for all filtering.
  reportlist <- list()

  #start keeping a cell status DF
  cellstatus <- data.frame(barcodes = colnames(sobj), filteredout = F, filterreason = 'Unfiltered')

  reportlist[['cellstatus']] <- cellstatus





  #report commands used
  reportlist[['allcommands']] <- data.frame(
    Command = c("mad.score.threshold",
                'min_num_UMI', 'min_num_Feature', 'max_perc_mito', 'max_perc_hemoglobin',
                "globalfilter.complexity", "globalfilter.libsize", "globalfilter.mito"),

    Option = c(mad.score.threshold,
               min_num_UMI, min_num_Feature, max_perc_mito, max_perc_hemoglobin,
               globalfilter.complexity, globalfilter.libsize, globalfilter.mito)

  )


  if( !('percent.mito' %in% colnames(sobj@meta.data)) ){

    message('Calculating percent.mito')

    #mito content, add to metadata
    mito.features <- grep(pattern = "^mt-", x = rownames(x = sobj), value = TRUE, ignore.case = T)
    sobj[["percent.mito"]] <- Seurat::PercentageFeatureSet(sobj, features = mito.features)

  }


  if( !('percent.hemoglobin' %in% colnames(sobj@meta.data)) ){


    sobj[["percent.hemoglobin"]] <- FerrenaSCRNAseq::calculate_percent.hemoglobin(sobj)

  }

  #report baseline stuff
  md <- sobj@meta.data
  reportlist[['baseline_summary']] <- data.frame(summary_nCount_RNA = as.vector(summary(md$nCount_RNA)),
                                                 summary_nFeature_RNA = as.vector(summary(md$nFeature_RNA)),
                                                 summary_perc.mito = as.vector(summary(md$percent.mito)),
                                                 summary_perc.hemoglobin = as.vector(summary(md$percent.hemoglobin)),
                                                 row.names = names(summary(md$nCount_RNA)))


  ### apply cutoffs based on min umi, min gene, max mito, max hemoglobin

  md <- sobj@meta.data

  bad <- md[md$nCount_RNA < min_num_UMI | md$nFeature_RNA < min_num_Feature | md$percent.mito > max_perc_mito | md$percent.hemoglobin > max_perc_hemoglobin,]

  cellstatus[cellstatus$barcodes %in% rownames(bad),"filteredout"] <- T
  cellstatus[cellstatus$barcodes %in% rownames(bad),"filterreason"] <- 'BasicFilter'


  cellstatus$BasicFilter <- 'No'

  #add reasons for filtering
  #min UMI:
  explanation_string <- paste0('_nUMI-below-', min_num_UMI)
  bad_reason <- rownames( bad[bad$nCount_RNA < min_num_UMI,] )
  cellstatus[ cellstatus$barcodes %in% bad_reason , 'BasicFilter'] <- paste0(cellstatus[ cellstatus$barcodes %in% bad_reason , 'filterreason'],
                                                                             explanation_string)

  #min feature:
  explanation_string <- paste0('_nFeature-below-', min_num_Feature)
  bad_reason <- rownames( bad[bad$nFeature_RNA < min_num_Feature,] )
  cellstatus[ cellstatus$barcodes %in% bad_reason , 'BasicFilter'] <- paste0(cellstatus[ cellstatus$barcodes %in% bad_reason , 'filterreason'],
                                                                             explanation_string)

  #max mito:
  explanation_string <- paste0('_percent.mito-above', max_perc_mito)
  bad_reason <- rownames( bad[bad$percent.mito > max_perc_mito,] )
  cellstatus[ cellstatus$barcodes %in% bad_reason , 'BasicFilter'] <- paste0(cellstatus[ cellstatus$barcodes %in% bad_reason , 'filterreason'],
                                                                             explanation_string)

  #max hemoglobin:
  explanation_string <- paste0('_percent.hemoglobin-above', max_perc_hemoglobin)
  bad_reason <- rownames( bad[bad$percent.hemoglobin > max_perc_hemoglobin,] )
  cellstatus[ cellstatus$barcodes %in% bad_reason , 'BasicFilter'] <- paste0(cellstatus[ cellstatus$barcodes %in% bad_reason , 'filterreason'],
                                                                             explanation_string)



  #retain good cells
  good <- cellstatus[cellstatus$filteredout==F,"barcodes"]
  sobj <- sobj[,good]


  if(globalfilter.complexity == T){

    message(' - Initiating baseline complexity filtration')

    # complexity = num Genes given num UMIs.
    # nFeature_RNA / nCount_RNA

    # use a linear model to try to ID outliers with low-complexity


    #set up md
    md <- sobj@meta.data

    #calculate and order by complexity
    md$complexity <-   md$nFeature_RNA / md$nCount_RNA
    md <- md[order(md$complexity),]

    #two-step model:
    # loess: must have low negative residuals < -1
    # lm: must have cooks distance > 4 / n

    #lm
    m <- stats::lm(data = md, formula = log(nFeature_RNA) ~ log(nCount_RNA))

    #calculate cooks distance for each point
    cd <- cooks.distance(m)

    #get scaled residuals for each point
    resid_lm <- scale( residuals(m))

    # log transforms are used to flatten variance, make things more normal-looking, bring to similar scale
    m_loess <- stats::loess(data = md, formula = log(nFeature_RNA) ~ log(nCount_RNA))

    #get scaled residuals for each point from LOESS
    resid_loess <- scale( residuals(m_loess))

    out_decision <- data.frame(cd_lm = cd,
                               resid_lm  =resid_lm,
                               resid_loess = resid_loess)

    #to be an outlier, must have scaled loess residuals < -1 and lm cooks distance > 4/n
    outs <- rownames( out_decision[out_decision$cd_lm > (4 / nrow(out_decision)) & out_decision$resid_loess < -2,] )

    #plot outliers
    md$outlier <- 'nonoutlier'
    md[rownames(md) %in% outs, "outlier"] <- 'outlier'


    #plot relationship with outlier calls
    complexityplot.outliers <- ggplot2::ggplot(md, ggplot2::aes(log(nCount_RNA), log(nFeature_RNA), col = outlier)) +
      ggplot2::geom_point(alpha = 0.7, size = 0.7)+
      ggplot2::scale_color_brewer(palette = 'Set1', direction = -1)+
      ggplot2::labs( caption = 'Outliers = Cooks Distance > (4/n)\n& Loess residual < -2' )


    #plot the relationship; we must use log transforms to see things more clearly.
    complexityplot <- ggplot2::ggplot(md, ggplot2::aes(log(nCount_RNA), log(nFeature_RNA), col = complexity)) +
      ggplot2::geom_point(alpha = 0.7, size = 0.7)+
      ggplot2::labs(caption = 'Complexity = nFeature / nCount')


    #plot loess with residuals and lm with cd
    md <- cbind(md, out_decision)
    md$CooksDistance_lm <- md$cd_lm

    cp_lm_cd <- ggplot2::ggplot(md, ggplot2::aes(log(nCount_RNA), log(nFeature_RNA), col = CooksDistance_lm)) +
      ggplot2::geom_point(alpha = 0.7, size = 0.7)+
      ggplot2::geom_smooth(method = 'lm')+
      ggplot2::scale_color_gradient(high = 'red')+
      ggplot2::labs(caption = paste0('Cutoff = cooks > 4 / N, here = ', round(4/nrow(md), 5) ) )

    cp_ls_rd <- ggplot2::ggplot(md, ggplot2::aes(log(nCount_RNA), log(nFeature_RNA), col = resid_loess)) +
      ggplot2::geom_point(alpha = 0.7, size = 0.7)+
      ggplot2::geom_smooth(method = 'loess')+
      ggplot2::scale_color_gradient2(mid = 'grey', low = 'red', high = 'blue')+
      ggplot2::labs(caption = 'Cutoff = loess residuals < -2' )


    finalplot <- patchwork::wrap_plots(complexityplot,complexityplot.outliers, cp_lm_cd, cp_ls_rd)+
      patchwork::plot_annotation(title = 'Complexity Filtering',
                                 subtitle = paste0('Total Cells: ', ncol(sobj),
                                                   '; Num Outlers: ', length(outs)) )


    #subset and save
    bad <- outs
    filteredcells <- rownames( md[md$outlier == 'nonoutlier',] )

    sobj <- sobj[,filteredcells]

    #report
    cellstatus[cellstatus$barcodes %in% bad,"filteredout"] <- T
    cellstatus[cellstatus$barcodes %in% bad,"filterreason"] <- 'globalfilter.complexity'


    reportlist[["globalfilter.complexity"]] <- list(plot = finalplot,
                                                    table = table(md$outlier) )





  }



  #find libsize cutoff low.
  if(globalfilter.libsize == T){

    message(' - Initiating global baseline lib size low filtration')


    #the use of median absolute deviation in outlier detection is inspired by the following:
    # https://www.sciencedirect.com/science/article/abs/pii/S0022103113000668
    # https://eurekastatistics.com/using-the-median-absolute-deviation-to-find-outliers/

    # don't use long-tailed hack.


    # get variable, in this case libsize
    # log transform to deal with extreme tails and negative values
    x <- log( sobj$nCount_RNA )


    # make a mad cutoff; similar to mean +/- sd*2.5, we use median +/- mad*2.5
    # get lib size cutoffs; use a lower cutoff only
    globalfilter.libsize.cutoff.lo <- median(x) - (mad.score.threshold  * mad(x) )

    #try multimode analysis...
    # mm <- multimode::locmodes(log(sobj$nCount_RNA), 2)
    #globalfilter.libsize.cutoff.lo <-mm$locations[1] - ( mm$locations[2] - mm$locations[1] )


    #get cells outside threshold
    bad <- x[x < globalfilter.libsize.cutoff.lo]
    #bad <- x[x < globalfilter.libsize.cutoff.lo | x > globalfilter.libsize.cutoff.hi]

    #exponentiate to get out of log space
    nonexp <- globalfilter.libsize.cutoff.lo
    globalfilter.libsize.cutoff.lo <- exp(globalfilter.libsize.cutoff.lo)


    #get cells within threshold
    filteredcells <- names(x[!names(x) %in% names(bad)] )


    #plot cutoffs for library size / UMIs
    g_lib_hist <- ggplot2::ggplot(sobj@meta.data) +
      ggplot2::geom_histogram(ggplot2::aes(nCount_RNA),
                              color="black", fill = "steelblue",
                              # binwidth = 0.05
      )+
      #ggplot2::geom_vline(xintercept = globalfilter.libsize.cutoff.hi, linetype = "dashed", colour = "red")+
      ggplot2::geom_vline(xintercept =  globalfilter.libsize.cutoff.lo, linetype = "dashed", colour = "red")+
      ggplot2::geom_vline(xintercept = median(sobj@meta.data$nCount_RNA),
                          linetype = "dotted", colour = "red", size = 1.2)+
      ggplot2::scale_x_log10(labels = scales::comma)+
      ggplot2::labs(x = "", y = "Number of cells" ,
                    caption = paste0('Num cells presubset = ', ncol(sobj),
                                     '\nNum cells remaining = ', length(filteredcells))  )+
      ggplot2::theme_linedraw()+
      ggplot2::coord_flip()

    g_lib_vln <- ggplot2::ggplot(sobj@meta.data, ggplot2::aes(x = 0, y = nCount_RNA))+
      ggplot2::geom_violin(fill='steelblue')+
      ggplot2::geom_jitter(height = 0, width = 0.25, size = 0.1)+
      ggplot2::geom_hline(yintercept =  globalfilter.libsize.cutoff.lo, col = 'red', linetype = 'dashed')+
      ggplot2::geom_hline(yintercept = median(sobj@meta.data$nCount_RNA),
                          linetype = "dotted", colour = "red", size = 1.2)+
      ggplot2::theme_linedraw()+
      ggplot2::scale_y_log10(labels = scales::comma)+
      ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                     axis.ticks.x=ggplot2::element_blank())+
      ggplot2::labs(y = "nUMI, log10 scale",
                    caption = paste0("Low UMI cutoff: ", as.character(round(globalfilter.libsize.cutoff.lo, digits = 3)),
                                     "\nMedian lib size = ", median(sobj@meta.data$nCount_RNA)) )+
      ggplot2::xlab('Cells')



    #report
    cellstatus[cellstatus$barcodes %in% names(bad),"filteredout"] <- T
    cellstatus[cellstatus$barcodes %in% names(bad),"filterreason"] <- 'globalfilter.libsize'


    #set up reportables
    finalplot <- patchwork::wrap_plots(g_lib_vln, g_lib_hist)+
      patchwork::plot_annotation(title = 'Num UMI Filtering',
                                 subtitle = paste0('Total Cells: ', ncol(sobj),
                                                   '; Num Outlers: ', length(bad)) )

    finaltable <- data.frame(var = x, var_cond = x < nonexp)
    finaltable$outlier <- 'nonoutlier'
    finaltable[finaltable$var_cond==T,'outlier'] <- 'outlier'
    finaltable <- table(finaltable$outlier)

    reportlist[['globalfilter.libsize']] <- list(plot = finalplot, table = finaltable)


    sobj <- sobj[,filteredcells]

  } #close baseline libsize loop



  #find mito cutoff high
  if(globalfilter.mito == T){

    message(' - Initiating global baseline mito hi filtration')



    # get variable, in this case mito
    x <- sobj$percent.mito

    # make a mad cutoff; similar to mean +/- sd*2.5, we use median +/- mad*2.5
    # get lib size cutoffs; use a high cutoff only
    globalfilter.mito.cutoff <- median(x) + (mad.score.threshold  * mad(x) )


    #get cells above the threshold
    bad <- x[x > globalfilter.mito.cutoff]

    #get cells within threshold
    filteredcells <- names(x[!names(x) %in% names(bad)] )



    g_mito_hist <- ggplot2::ggplot(sobj@meta.data) +
      ggplot2::geom_histogram(ggplot2::aes(percent.mito),
                              color="black", fill = "steelblue",
                              binwidth = 0.5)+
      ggplot2::geom_vline(xintercept =  globalfilter.mito.cutoff, linetype = "dashed", colour = "red")+
      ggplot2::geom_vline(xintercept = median(sobj@meta.data$percent.mito),
                          linetype = "dotted", colour = "red", size = 1.2)+
      ggplot2::labs(y = "Number of cells" , x = '',
                    caption = paste0('Num cells presubset = ', ncol(sobj),
                                     '\nNum cells remaining = ', length(filteredcells))  )+
      ggplot2::theme_linedraw()+
      ggplot2::coord_flip()

    g_mito_vln <- ggplot2::ggplot(sobj@meta.data, ggplot2::aes(x = 0, y = percent.mito))+
      ggplot2::geom_violin(fill='steelblue')+
      ggplot2::geom_jitter(height = 0, width = 0.25, size = 0.1)+
      ggplot2::geom_hline(yintercept = globalfilter.mito.cutoff, col = 'red', linetype = 'dashed')+
      ggplot2::geom_hline(yintercept = median(sobj@meta.data$percent.mito),
                          linetype = "dotted", colour = "red", size = 1.2)+
      ggplot2::scale_y_continuous(breaks = c(0, 25, 50, 75, 100, round(globalfilter.mito.cutoff, digits = 2)))+
      ggplot2::theme_linedraw()+
      ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                     axis.ticks.x=ggplot2::element_blank())+
      ggplot2::labs(
        caption = paste0('Percent.mito cutoff = ', round(globalfilter.mito.cutoff, 2) ,
                         "\nMedian percent.mito = ", round( median(sobj@meta.data$percent.mito), 2))
      )+
      ggplot2::xlab('Cells')




    #report
    cellstatus[cellstatus$barcodes %in% names(bad),"filteredout"] <- T
    cellstatus[cellstatus$barcodes %in% names(bad),"filterreason"] <- 'globalfilter.mito'



    #set up reportables
    finalplot <- patchwork::wrap_plots(g_mito_vln, g_mito_hist)+
      patchwork::plot_annotation(title = 'Percent Mito Filtering',
                                 subtitle = paste0('Total Cells: ', ncol(sobj),
                                                   '; Num Outlers: ', length(bad)) )

    finaltable <- data.frame(var = x, var_cond = x > globalfilter.mito.cutoff)
    finaltable$outlier <- 'nonoutlier'
    finaltable[finaltable$var_cond==T,'outlier'] <- 'outlier'
    finaltable <- table(finaltable$outlier)

    reportlist[['globalfilter.mito']] <- list(plot = finalplot, table = finaltable)


    sobj <- sobj[,filteredcells]


  } #end mito baseline block



  #update reportlist with finalized cellstatus:
  reportlist[["cellstatus"]] <- cellstatus

  #return whole result list.
  reportlist


}




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








### testing below

# hemoglobin.features <- c('HBA1', 'HBA2', 'HBB', 'HBD', 'HBE1', 'HBG1', 'HBG2', 'HBM', 'HBQ1', 'HBZ',
#                          'Hba', 'Hba-a1', 'Hba-a2', 'Hba-ps3', 'Hba-ps4', 'Hba-x', 'Hbb', 'Hbb-ar', 'Hbb-b1', 'Hbb-b2', 'Hbb-bh0', 'Hbb-bh1', 'Hbb-bh2', 'Hbb-bh3', 'Hbb-bs', 'Hbb-bt', 'Hbb-y')
#
#
# min_num_UMI <- 1000
# min_num_Feature <- 200
# max_perc_mito <- 25
# max_perc_hemoglobin <- 25
#
#
# mad.score.threshold = 2.5
# globalfilter.complexity <- T
# globalfilter.mito <- T
# globalfilter.libsize <- T
#
# rawh5 <- '~/Dropbox/data/bangdata/scrnaseq-TKO-DKOAA-DKO/rawdata/Sample-06_TL494/filtered_feature_bc_matrix.h5'
# rawh5 <- '~/Dropbox/data/bangdata/scrnaseq-TKO-DKOAA-DKO/rawdata/Sample-04_DJ582M11/filtered_feature_bc_matrix.h5'
#
# sobj <- CreateSeuratObject(   Read10X_h5(rawh5), min.cells= 3)
#
# af <- autofilter(sobj)
#
#
# good <- af$cellstatus[af$cellstatus$filteredout==F,"barcodes"]







