# https://r-pkgs.org/whole-game.html





#' Automated dataset-specific quality control.
#'
#' Given a Seurat object, perform QC filtering.
#'
#' This function utilizes "complexity filtration" in the form of  multivariable outlier identification based on modelling number of unique genes by total UMIs and percent mitochondrial content. Cooks distance is used to call outliers in this setting.
#' Also, this function uses one-sided univariate filtering of library size and percent mitochdonrial content. This makes use of median absolute deviation outlier thresholding, with adjustments for long-tailed distributions.
#' An important parameter for univariate median absolute deviaiton cutoff is mad.score.threshold, by default this is set to 2.5; this is how many deviations from the median are tolerated, beyond which cells are classified as outliers.
#'
#' Additionally, this function will attempt to perform partition (or cluster)-specific outlier filtering. This is done by stratifying the data into detected or provided clusters and then applying multivariable complexity detection or univariate outlier detection.
#' The reasoning behind this is that the number of UMIs and Unique Genes may be cell-type specific. Taking this into account during threshold can reduce Type-1 and Type-2 error in outlier calling.
#'
#' For reference, see:
#' univariate outlier calling with Median Absolute Deviation: Leys et al ScienceDirect 2013 (https://www.sciencedirect.com/science/article/abs/pii/S0022103113000668)
#' adjustment of median absolute deviation method for long-tailed distributions: Peter Rosenmai 2013 (https://eurekastatistics.com/using-the-median-absolute-deviation-to-find-outliers/)
#'
#' @param seuratobject A Seurat object. should have per-cell total UMI and Feature columns in metadata titled "nCount_RNA" and "nFeature_RNA"
#' @param mad.score.threshold Numeric. Maximum number of deviations from median tolerated before outlier classification. Default = 2.5. lower is more strict, higher is more lenient.
#' @param baselinefilter T/F; whether to perform "global" QC, ie without pre-clustering. Default = True.
#' @param baselinefilter.complexity T/F; whether to perform global filtering of cells based on complexity. Uses linear model for each cell: Num Unique Genes ~ Num UMI + percent mito. outliers are called based on Cook's distance < 1/4n. Default = T
#' @param baselinefilter.libsize T/F; whether to perform global lower lib size filtration using median absolute deviation threshold; will only work if baselinefilter is set to True. Default = T.
#' @param baselinefilter.mito T/F; whether to perform global upper mitochondrial content filtration using median absolute deviation threshold; will only work if baselinefilter is set to True. Default = T.
#' @param clusters string. Name of seuratobject metadata to be used as clusters. Optimally corresponds to cell type. Or, set to "automate" to call clusters using Seurat pipeline (SCT->30PCs->Leiden with res=0.1). Default = 'automate'
#' @param removemitohiclust T/F ; whether to identify and remove abnormally high mitochondrial content clusters; if no baseline mito filtration is used, there will almost certainly be a mito-driven cluster. Identification is via outliers::grubbs.test(). Default = T.
#' @param iterativefilter T/F ; whether to perform cluster-specific filtering. Default = T.
#' @param iterativefilter.complexity T/F; whether to perform cell type-specific filtering of cells based on complexity. Uses linear model for each cell: Num Unique Genes ~ Num UMI + percent mito. outliers are called based on Cook's distance < 1/4n. Default = T
#' @param iterativefilter.libsize T/F ; whether to perform cluster-specific filtering of lower library size cutoff. iterativefilter must be true. Default = T.
#' @param iterativefilter.mito T/F; whether to perform cluster-specific filtering of upper mito content cutoff. iterativefilter must be true. Default = T.
#'
#' @return a list object. The first element has the barcodes and outlier classification (ie whether the cell should be filtered out or not), along with the step of filtering if called as an outlier. other elements include details of the filtering.
#' @export
#'
#' @examples
#' \dontrun{
#' #Run automated filtering
#' reportlist <- automatedfiltering(seuratobject)
#'
#' #Actually filter based on the calls
#' barcodefilter <- reportlist[[1]]
#'
automatedfiltering <- function(
  seuratobject,

  mad.score.threshold,

  baselinefilter,
  baselinefilter.complexity,
  baselinefilter.mito,
  baselinefilter.libsize,

  clusters,
  removemitohiclust,
  iterativefilter,
  iterativefilter.complexity,
  iterativefilter.libsize,
  iterativefilter.mito){

  #check if seurat object is there
  if( missing( seuratobject )) {stop('please input a Seurat object')}
  if( is(seuratobject) != 'Seurat') {stop('please input a Seurat object')}

  #set maximum distance of deviations from median tolerated before outlier classification
  if( missing( mad.score.threshold )) {mad.score.threshold <- 2.5}

  # baseline (global) filtration
  if( missing( baselinefilter )) {baselinefilter <- T}
  if( missing( baselinefilter.complexity )) {baselinefilter.complexity <- T}
  if( missing( baselinefilter.mito )) {baselinefilter.mito <- T}
  if( missing( baselinefilter.libsize )) {baselinefilter.libsize <- T}

  # cluster based filtration
  if( missing( clusters )) {clusters <- "automate"}
  if( missing( removemitohiclust )) {removemitohiclust <- T}
  if( missing( iterativefilter )) {iterativefilter <- T}
  if( missing( iterativefilter.complexity )) {iterativefilter.complexity <- T}
  if( missing( iterativefilter.libsize )) {iterativefilter.libsize <- T}
  if( missing( iterativefilter.mito )) {iterativefilter.mito <- T}






  message('num genes = ', nrow(seuratobject) , '\n',
          'num cells = ', ncol(seuratobject)  , '\n')


  #output list object; add to this as needed for all filtering.
  reportlist <- list()

  #start keeping a cell status DF
  cellstatus <- data.frame(barcodes = colnames(seuratobject), filteredout = 'No', filterreason = 'Unfiltered_NA')

  reportlist[['cellstatus']] <- cellstatus


  #rport commands used
  reportlist[['allcommands']] <- data.frame(
    Command = c("mad.score.threshold",
                "baselinefilter", "baselinefilter.complexity", "baselinefilter.libsize", "baselinefilter.mito",
                "clusters", "removemitohiclust",
                "iterativefilter", "iterativefilter.complexity", "iterativefilter.libsize", "iterativefilter.mito"),



    Option = c(mad.score.threshold,
               baselinefilter, baselinefilter.complexity, baselinefilter.libsize, baselinefilter.mito,
               clusters, removemitohiclust,
               iterativefilter, iterativefilter.complexity, iterativefilter.libsize, iterativefilter.mito)

  )





  #mito content, add to metadata
  mito.features <- grep(pattern = "^mt-", x = rownames(x = seuratobject), value = TRUE, ignore.case = T)
  seuratobject[["percent.mito"]] <- Seurat::PercentageFeatureSet(seuratobject, features = mito.features)

  #whether to perform baseline filtering at all.
  if(baselinefilter == T) {

    message('Running baseline filtering:')


    if(baselinefilter.complexity == T){

      message(' - Initiating baseline complexity filtration')

      #there is usually a linear relationship between nUMIs and nGenes.
      # high counts + low genes = weird cells...
      # these cells may also have high mito.

      # use a linear model to try to ID outliers with low-complexity


      #set up df; order by mito, for plotting
      df <- seuratobject@meta.data[order(seuratobject$percent.mito),]

      #plot the relationship; we must use log transforms to see things more clearly.
      complexityplot <- ggplot2::ggplot(df, ggplot2::aes(log(nCount_RNA), log(nFeature_RNA), col = percent.mito)) +
        ggplot2::geom_point(alpha = 0.7, size = 0.7)

      #set up the model
      # log transforms are used to flatten variance, make things more normal-looking, bring to similar scale
      # pseudocounts prevent log(0)
      m <- stats::lm(data = df, formula = log(nFeature_RNA+1) ~ log(nCount_RNA+1))

      #calculate cooks distance for each point
      cd <- cooks.distance(m)

      #define outiers using cutoff of cooks distance > 4 / n
      outs <- cd[cd > 4/length(cd)]
      #plot(cd) ; abline(h=4/length(cd))


      #plot outliers
      df$outlier <- 'nonoutlier'
      df[rownames(df) %in% names(outs), "outlier"] <- 'outlier'

      #plot relationship with outlier calls
      complexityplot.outliers <- ggplot2::ggplot(df, ggplot2::aes(log(nCount_RNA), log(nFeature_RNA), col = outlier)) +
        ggplot2::geom_point(alpha = 0.7, size = 0.7)+
        ggplot2::scale_color_brewer(palette = 'Set1', direction = -1)+
        ggplot2::labs(caption = paste0('Total Cells: ', ncol(seuratobject),
                                       '\nNum Outlers: ', length(outs) ) )


      #just to visualize, make the non-log plots
      complexityplot.nolog <- ggplot2::ggplot(df, ggplot2::aes(nCount_RNA, nFeature_RNA, col = percent.mito)) +
        ggplot2::geom_point(alpha = 0.7, size = 0.7)

      complexityplot.outliers.nolog <- ggplot2::ggplot(df, ggplot2::aes(nCount_RNA, nFeature_RNA, col = outlier)) +
        ggplot2::geom_point(alpha = 0.7, size = 0.7)+
        ggplot2::scale_color_brewer(palette = 'Set1', direction = -1)

      #       #testing purposes: maybe use MCD?
      #https://www.sciencedirect.com/science/article/abs/pii/S0022103117302123
      #       df2 <- df[,2:4]
      #       df2 <- apply(df2, MARGIN = 2, FUN = function(x){log(x+1)})
      #       outs <- Routliers::outliers_mcd(df2, alpha = 0.001)
      #       Routliers::plot_outliers_mcd(df2, res = outs)

      #subset and save
      bad <- outs
      filteredcells <- names(cd[!names(cd) %in% names(bad)] )

      seuratobject <- subset(x = seuratobject,
                             cells = filteredcells)

      #report
      cellstatus[cellstatus$barcodes %in% names(bad),"filteredout"] <- 'Yes'
      cellstatus[cellstatus$barcodes %in% names(bad),"filterreason"] <- 'baselinefilter.complexity'


      reportlist[["baselinefilter.complexity"]] <- patchwork::wrap_plots(complexityplot.nolog,complexityplot.outliers.nolog,
                                                                         complexityplot,complexityplot.outliers)



    }



    #find libsize cutoff low.
    if(baselinefilter.libsize == T){

      message(' - Initiating global baseline lib size low filtration')

      #the use of median absolute deviation in outlier detection is inspired by the following:
      # https://www.sciencedirect.com/science/article/abs/pii/S0022103113000668
      # https://eurekastatistics.com/using-the-median-absolute-deviation-to-find-outliers/

      #UPDATE 2021.07.20 - WE WILL NOT USE THE LONG-TAILED MAD HACK.


      # get variable, in this case libsize
      # log transform to deal with extreme tails
      x <- log(seuratobject$nCount_RNA)

      #instead of using normal mad (ie mad(x)),
      # use an adjustment for long-tailed distributions by Peter Rosenmai
      # this is the median of absolute deviations for points above the data median
      left.mad  <- median(   abs(x - median(x))[x < median(x)]    )
      #right.mad <- median(   abs(x - median(x))[x > median(x)]    )
      #only do this for lower tail.

      # get lib size cutoffs; use a right/upper tail only
      # right tail of Z-score, solve for which value matches the
      baselinefilter.libsize.cutoff.lo <- median(x) - (mad.score.threshold  * mad(x) )
      #baselinefilter.libsize.cutoff.hi <- median(x) + (mad.score.threshold  * mad(x) )
      #only do this for lower tail.


      #get cells outside threshold
      bad <- x[x < baselinefilter.libsize.cutoff.lo]
      #bad <- x[x < baselinefilter.libsize.cutoff.lo | x > baselinefilter.libsize.cutoff.hi]

      #exponentiate to get out of log space
      #baselinefilter.libsize.cutoff.hi <- exp(baselinefilter.libsize.cutoff.hi)
      baselinefilter.libsize.cutoff.lo <- exp(baselinefilter.libsize.cutoff.lo)


      #get cells within threshold
      filteredcells <- names(x[!names(x) %in% names(bad)] )


      #plot cutoffs for library size / UMIs
      g_lib_hist <- ggplot2::ggplot(seuratobject@meta.data) +
        ggplot2::geom_histogram(ggplot2::aes(nCount_RNA),
                                color="black", fill = "steelblue",
                                binwidth = 0.05)+
        #ggplot2::geom_vline(xintercept = baselinefilter.libsize.cutoff.hi, linetype = "dashed", colour = "red")+
        ggplot2::geom_vline(xintercept =  baselinefilter.libsize.cutoff.lo, linetype = "dashed", colour = "red")+
        ggplot2::geom_vline(xintercept = median(seuratobject@meta.data$nCount_RNA),
                            linetype = "dotted", colour = "red", size = 1.2)+
        ggplot2::scale_x_log10()+
        ggplot2::labs(x = "nUMI, log10 scale",
                      y = "Number of cells" ,
                      subtitle = paste0("Median lib size = ", median(seuratobject@meta.data$nCount_RNA)),
                      caption = paste0('nUMI cutoff = ', round(baselinefilter.libsize.cutoff.lo, digits = 3),
                                       '\nNum cells presubset = ', ncol(seuratobject),
                                       '\nNum cells remaining = ', length(filteredcells))  )+
        ggplot2::theme_linedraw()

      g_lib_vln <- ggplot2::ggplot(seuratobject@meta.data, ggplot2::aes(x = 0, y = nCount_RNA))+
        ggplot2::geom_violin(fill='steelblue')+
        ggplot2::geom_jitter(height = 0, width = 0.25, size = 0.1)+
        ggplot2::geom_hline(yintercept =  baselinefilter.libsize.cutoff.lo, col = 'red', linetype = 'dashed')+
        ggplot2::theme_linedraw()+
        ggplot2::scale_y_log10()+
        ggplot2::theme(axis.text.x=ggplot2::element_blank(),
              axis.ticks.x=ggplot2::element_blank())+
        ggplot2::labs(title = paste0("nUMI cutoff"),
                      y = "nUMI, log10 scale",
                      subtitle = paste0("low UMI cutoff: ", as.character(round(baselinefilter.libsize.cutoff.lo, digits = 3))) )+
        ggplot2::xlab('Cells')



      #report
      cellstatus[cellstatus$barcodes %in% names(bad),"filteredout"] <- 'Yes'
      cellstatus[cellstatus$barcodes %in% names(bad),"filterreason"] <- 'baselinefilter.libsize'


      reportlist[['baselinefilter.libsize']] <- patchwork::wrap_plots(g_lib_vln, g_lib_hist)


      seuratobject <- subset(x = seuratobject,
                             cells = filteredcells)


    } #close baseline libsize loop



    #find mito cutoff high
    if(baselinefilter.mito == T){

      message(' - Initiating global baseline mito hi filtration')

      #the use of median absolute deviation in outlier detection is inspired by the following:
      # https://www.sciencedirect.com/science/article/abs/pii/S0022103113000668

      #the adjustment for non-symmetric / long-tailed data is inspired by this:
      # https://eurekastatistics.com/using-the-median-absolute-deviation-to-find-outliers/


      # get variable, in this case mito
      x <- seuratobject$percent.mito

      #instead of using normal mad (ie mad(x)),
      # use an adjustment for long-tailed distributions by Peter Rosenmai
      # this is the median of absolute deviations for points above the data median
      right.mad <- median(   abs(x - median(x))[x > median(x)]    )

      # get percent mito cutoff; use a right/upper tail only
      # right tail of Z-score, solve for which value matches the
      baselinefilter.mito.cutoff <- median(x) + (mad.score.threshold  * mad(x) )


      #get cells above the threshold
      bad <- x[x > baselinefilter.mito.cutoff]

      #get cells within threshold
      filteredcells <- names(x[!names(x) %in% names(bad)] )



      gm <- ggplot2::ggplot(seuratobject@meta.data, ggplot2::aes(x = 0, y = percent.mito))+
        ggplot2::geom_violin(fill='steelblue')+
        ggplot2::geom_jitter(height = 0, width = 0.25, size = 0.1)+
        ggplot2::geom_hline(yintercept = baselinefilter.mito.cutoff, col = 'red', linetype = 'dashed')+
        ggplot2::scale_y_continuous(breaks = c(0, 25, 50, 75, 100, round(baselinefilter.mito.cutoff, digits = 2)))+
        ggplot2::theme_linedraw()+
        ggplot2::theme(axis.text.x=ggplot2::element_blank(),
              axis.ticks.x=ggplot2::element_blank())+
        ggplot2::labs(title = paste0("percent mito cutoff"),
                      subtitle = paste0("Hi percent.mitop cutoff: ", as.character(round(baselinefilter.mito.cutoff, digits = 2))),
                      caption = paste0('Percent.mito cutoff = ', round(baselinefilter.mito.cutoff, 2) ,
                                       '\nNum cells presubset = ', ncol(seuratobject),
                                       '\nNum cells remaining = ', length(filteredcells)))+
        ggplot2::xlab('Cells')


      #report
      cellstatus[cellstatus$barcodes %in% names(bad),"filteredout"] <- 'Yes'
      cellstatus[cellstatus$barcodes %in% names(bad),"filterreason"] <- 'baselinefilter.mito'

      reportlist[['baselinefilter.mito']] <- gm

      seuratobject <- subset(x = seuratobject,
                             cells = filteredcells)


    } #end mito baseline block


  } # close baseline filtering loop






  # for cluster-filtering, either cluster the data or get the inputted clusters

  if(clusters == 'automate'){

    message('\nClusters set to automate. Initiating clustering\n')

    #normalize and cluster, first pass, for lib size filtration
    suppressWarnings(seuratobject <- Seurat::SCTransform(seuratobject, verbose = T))

    seuratobject <- Seurat::RunPCA(object = seuratobject, verbose = F)

    seuratobject <- Seurat::FindNeighbors(object = seuratobject, dims = 1:30, verbose = F)
    seuratobject <- Seurat::FindClusters(object = seuratobject, resolution = 0.1, verbose = F, algorithm = 4)

  } else{
    seuratobject$seurat_clusters <- seuratobject@meta.data[,clusters]
  }



  ### remove the mito cluster and recluster ###
  if(removemitohiclust == T){
    message('\nSearching for Mito High clusters')

    #get average mito for each cluster
    md <- seuratobject@meta.data
    avgs <- aggregate(data = md, cbind(nCount_RNA , nFeature_RNA , percent.mito) ~ seurat_clusters, FUN=mean)
    rm(md)
    avgsseuratobject <- avgs
    gt <- outliers::grubbs.test(avgsseuratobject$percent.mito)

    #grubbs test while loop for right-tail outliers, w/ bonferroni correction
    iter=1
    while(
      stringr::str_split(gt$alternative, ' ')[[1]][1] == 'highest' &
      gt$p.value < (0.05 / nrow(avgsseuratobject))
    ) {
      message('\tgrubbs outlier test iteration: ',iter)

      avgsseuratobject <- avgsseuratobject[-which.max(avgsseuratobject$percent.mito),]

      gt <- outliers::grubbs.test(avgsseuratobject$percent.mito)

      iter=iter+1
    }

    mitohiclusts <- as.vector(avgs$seurat_clusters[!(avgs$seurat_clusters %in% avgsseuratobject$seurat_clusters)])

    #sometimes no mito clusts are detected, so this if-esle avoids crashes in that case
    if( length(mitohiclusts) != 0 ) {

      #plot violin of mito clusts
      vln1 <- Seurat::VlnPlot(seuratobject, c('percent.mito'), pt.size = 0.1)+Seurat::NoLegend()

      #get cells, report
      mitoclustcells <- Seurat::WhichCells(seuratobject, idents = mitohiclusts)
      bad <- mitoclustcells
      filteredcells <- filteredcells <- names(x[!names(x) %in% names(bad)] )


      seuratobject <- seuratobject[,!(colnames(seuratobject) %in% mitoclustcells)]

      vln2 <- Seurat::VlnPlot(seuratobject, c('percent.mito'), pt.size = 0.1) + ggplot2::labs(caption=paste0('Mito hi clusts: ', mitohiclusts))+Seurat::NoLegend()


      #report
      cellstatus[cellstatus$barcodes %in% names(bad),"filteredout"] <- 'Yes'
      cellstatus[cellstatus$barcodes %in% names(bad),"filterreason"] <- 'removemitohiclust'

      reportlist[['removemitohiclust']] <-  patchwork::wrap_plots(vln1,vln2, ncol=1)


      seuratobject <- subset(x = seuratobject,
                             cells = filteredcells)



    } else{
      message('No Mito High clusters identified by Grubbs test\n')

      reportlist[['removemitohiclust']] <- paste0("No Mito High clusters identified by Grubbs test")

    } #end if-else loop
  } # end mitohiclust loop


  ########################### ITERATIVE FILTRATION LOOP ####################################



  if(iterativefilter == T) {

    message('\nInitiate cluster-specific filtration')

    #set default ident
    seuratobject <- Seurat::SetIdent(seuratobject, value = seuratobject$seurat_clusters)

    #make list for cluster filtering. will add as one element to reportlist at end.
    iterativefilterresults <- list()


    #loop through clusters and get cells with proper library size, based on MAD

    for(clust in stringr::str_sort(unique(seuratobject$seurat_clusters), numeric = T) ) {        #print(clust) }

      message('Filtering cluster ', clust)
      seuratobjectclust <- subset(seuratobject, idents = clust)


      if(ncol(seuratobjectclust) < 30) {
        warning(paste0('Cluster ', clust ,' has under 30 cells; treat results with caution'))
      }

      #make list for this cluster, will add at end
      clustlist <- list()

      ########### ITERATIVE COMPLEXITY FILTER  ###########
      if(iterativefilter.complexity == T){
        #there is usually a linear relationship between nUMIs and nGenes.
        # high counts + low genes = weird cells...
        # these cells may also have high mito.

        # use a linear model to try to ID outliers with low-complexity


        #set up df; order by mito, for plotting
        df <- seuratobjectclust@meta.data[order(seuratobjectclust$percent.mito),]

        #plot the relationship; we must use log transforms to see things more clearly.
        complexityplot <- ggplot2::ggplot(df, ggplot2::aes(log(nCount_RNA), log(nFeature_RNA), col = percent.mito)) +
          ggplot2::geom_point(alpha = 0.7, size = 0.7)

        #set up the model
        # log transforms are used to flatten variance, make things more normal-looking, bring to similar scale
        # pseudocounts prevent log(0)
        m <- stats::lm(data = df, formula = log(nFeature_RNA+1) ~ log(nCount_RNA+1))

        #calculate cooks distance for each point
        cd <- cooks.distance(m)

        #define outiers using cutoff of cooks distance > 4 / n
        outs <- cd[cd > 4/length(cd)]
        #plot(cd) ; abline(h=4/length(cd))


        #plot outliers
        df$outlier <- 'nonoutlier'
        df[rownames(df) %in% names(outs), "outlier"] <- 'outlier'

        #plot relationship with outlier calls
        complexityplot.outliers <- ggplot2::ggplot(df, ggplot2::aes(log(nCount_RNA), log(nFeature_RNA), col = outlier)) +
          ggplot2::geom_point(alpha = 0.7, size = 0.7)+
          ggplot2::scale_color_brewer(palette = 'Set1', direction = -1)+
          ggplot2::labs(caption = paste0('Total Cells: ', ncol(seuratobjectclust),
                                         '\nNum Outlers: ', length(outs) ) )


        #just to visualize, make the non-log plots
        complexityplot.nolog <- ggplot2::ggplot(df, ggplot2::aes(nCount_RNA, nFeature_RNA, col = percent.mito)) +
          ggplot2::geom_point(alpha = 0.7, size = 0.7)

        complexityplot.outliers.nolog <- ggplot2::ggplot(df, ggplot2::aes(nCount_RNA, nFeature_RNA, col = outlier)) +
          ggplot2::geom_point(alpha = 0.7, size = 0.7)+
          ggplot2::scale_color_brewer(palette = 'Set1', direction = -1)

        #       #testing purposes: maybe use MCD?
        #https://www.sciencedirect.com/science/article/abs/pii/S0022103117302123
        #       df2 <- df[,2:4]
        #       df2 <- apply(df2, MARGIN = 2, FUN = function(x){log(x+1)})
        #       outs <- Routliers::outliers_mcd(df2, alpha = 0.001)
        #       Routliers::plot_outliers_mcd(df2, res = outs)

        #subset and save
        bad <- outs
        filteredcells <- names(cd[!names(cd) %in% names(bad)] )

        seuratobjectclust <- subset(x = seuratobjectclust,
                                    cells = filteredcells)

        #report
        cellstatus[cellstatus$barcodes %in% names(bad),"filteredout"] <- 'Yes'
        cellstatus[cellstatus$barcodes %in% names(bad),"filterreason"] <- 'iterativefilter.complexity'


        clustlist[["iterativefilter.complexity"]] <- patchwork::wrap_plots(complexityplot.nolog+ggplot2::ggtitle(paste0('Cluster ', clust)),complexityplot.outliers.nolog,
                                                                           complexityplot,complexityplot.outliers
        )



      }


      ########### ITERATIVE LIB SIZE FILTER  ###########
      if(iterativefilter.libsize == T){


        #the use of median absolute deviation in outlier detection is inspired by the following:
        # https://www.sciencedirect.com/science/article/abs/pii/S0022103113000668
        # https://eurekastatistics.com/using-the-median-absolute-deviation-to-find-outliers/


        # get variable, in this case libsize
        # log transform to deal with extreme tails
        x <- log(seuratobjectclust$nCount_RNA)

        #instead of using normal mad (ie mad(x)),
        # use an adjustment for long-tailed distributions by Peter Rosenmai
        # this is the median of absolute deviations for points above the data median
        left.mad  <- median(   abs(x - median(x))[x < median(x)]    )
        #right.mad <- median(   abs(x - median(x))[x > median(x)]    )
        #only do this for lower tail.

        # get lib size cutoffs; use a right/upper tail only
        # right tail of Z-score, solve for which value matches the
        iterative.libsize.cutoff.lo <- median(x) - (mad.score.threshold  * mad(x) )
        #iterative.libsize.cutoff.hi <- median(x) + (mad.score.threshold  * mad(x) )
        #only do this for lower tail.


        #get cells outside threshold
        bad <- x[x < iterative.libsize.cutoff.lo]
        #bad <- x[x < iterative.libsize.cutoff.lo | x > iterative.libsize.cutoff.hi]

        #exponentiate to get out of log space
        #iterative.libsize.cutoff.hi <- exp(iterative.libsize.cutoff.hi)
        iterative.libsize.cutoff.lo <- exp(iterative.libsize.cutoff.lo)


        #get cells within threshold
        filteredcells <- names(x[!names(x) %in% names(bad)] )


        #plot cutoffs for library size / UMIs
        g_lib_hist <- ggplot2::ggplot(seuratobjectclust@meta.data) +
          ggplot2::geom_histogram(ggplot2::aes(nCount_RNA),
                                  color="black", fill = "steelblue",
                                  binwidth = 0.05)+
          #ggplot2::geom_vline(xintercept = iterative.libsize.cutoff.hi, linetype = "dashed", colour = "red")+
          ggplot2::geom_vline(xintercept =  iterative.libsize.cutoff.lo, linetype = "dashed", colour = "red")+
          ggplot2::geom_vline(xintercept = median(seuratobjectclust@meta.data$nCount_RNA),
                              linetype = "dotted", colour = "red", size = 1.2)+
          ggplot2::scale_x_log10()+
          ggplot2::labs(x = "nUMI, log10 scale",
                        y = "Number of cells" ,
                        subtitle = paste0("Median lib size = ", median(seuratobject@meta.data$nCount_RNA)),
                        caption = paste0('nUMI cutoff = ', round(baselinefilter.libsize.cutoff.lo, digits = 3),
                                         '\nNum cells presubset = ', ncol(seuratobject),
                                         '\nNum cells remaining = ', length(filteredcells))  )+
          ggplot2::theme_linedraw()

        g_lib_vln <- ggplot2::ggplot(seuratobjectclust@meta.data, ggplot2::aes(x = 0, y = nCount_RNA))+
          ggplot2::geom_violin(fill='steelblue')+
          ggplot2::geom_jitter(height = 0, width = 0.25, size = 0.1)+
          ggplot2::geom_hline(yintercept =  iterative.libsize.cutoff.lo, col = 'red', linetype = 'dashed')+
          ggplot2::theme_linedraw()+
          ggplot2::scale_y_log10()+
          ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                axis.ticks.x=ggplot2::element_blank())+
          ggplot2::labs(title = paste0("nUMI cutoff"),
                        y = "nUMI, log10 scale",
                        subtitle = paste0("low UMI cutoff: ", as.character(round(baselinefilter.libsize.cutoff.lo, digits = 3))) )+
          ggplot2::xlab('Cells')



        #report
        cellstatus[cellstatus$barcodes %in% names(bad),"filteredout"] <- 'Yes'
        cellstatus[cellstatus$barcodes %in% names(bad),"filterreason"] <- 'iterativefilter.libsize'


        clustlist[['iterativefilter.libsize']] <- patchwork::wrap_plots(g_lib_vln+ggplot2::ggtitle(paste0('Cluster ', clust)), g_lib_hist)


        seuratobjectclust <- subset(x = seuratobjectclust,
                                    cells = filteredcells)


      } #close iterative libsize loop


      ########### ITERATIVE MITO FILTER  ###########
      if(iterativefilter.mito == T){


        #the use of median absolute deviation in outlier detection is inspired by the following:
        # https://www.sciencedirect.com/science/article/abs/pii/S0022103113000668

        #the adjustment for non-symmetric / long-tailed data is inspired by this:
        # https://eurekastatistics.com/using-the-median-absolute-deviation-to-find-outliers/


        # get variable, in this case mito
        x <- seuratobjectclust$percent.mito

        #instead of using normal mad (ie mad(x)),
        # use an adjustment for long-tailed distributions by Peter Rosenmai
        # this is the median of absolute deviations for points above the data median
        right.mad <- median(   abs(x - median(x))[x > median(x)]    )

        # get percent mito cutoff; use a right/upper tail only
        # right tail of Z-score, solve for which value matches the
        iterativefilter.mito.cutoff <- median(x) + (mad.score.threshold  * mad(x) )


        #get cells above the threshold
        bad <- x[x > iterativefilter.mito.cutoff]

        #get cells within threshold
        filteredcells <- names(x[!names(x) %in% names(bad)] )



        gm <- ggplot2::ggplot(seuratobjectclust@meta.data, ggplot2::aes(x = 0, y = percent.mito))+
          ggplot2::geom_violin(fill='steelblue')+
          ggplot2::geom_jitter(height = 0, width = 0.25, size = 0.1)+
          ggplot2::geom_hline(yintercept = iterativefilter.mito.cutoff, col = 'red', linetype = 'dashed')+
          ggplot2::theme_linedraw()+
          ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                axis.ticks.x=ggplot2::element_blank())+
          ggplot2::labs(title = paste0("percent mito cutoff"),
                         subtitle = paste0("Hi percent.mitop cutoff: ", as.character(round(baselinefilter.mito.cutoff, digits = 2))),
                         caption = paste0('Percent.mito cutoff = ', round(baselinefilter.mito.cutoff, 2) ,
                                          '\nNum cells presubset = ', ncol(seuratobject),
                                          '\nNum cells remaining = ', length(filteredcells)))+
          ggplot2::xlab('Cells')


        #report
        cellstatus[cellstatus$barcodes %in% names(bad),"filteredout"] <- 'Yes'
        cellstatus[cellstatus$barcodes %in% names(bad),"filterreason"] <- 'iterativefilter.mito'

        clustlist[['iterativefilter.mito']] <- gm

        seuratobjectclust <- subset(x = seuratobjectclust,
                                    cells = filteredcells)


      } #end mito iterative block


      iterativefilterresults[[clust]] <- clustlist

    } #close cluster loop


    reportlist[['iterativefilter']] <- iterativefilterresults


  } # close iterative filter conditional loop





  #update reportlist with finalized cellstatus:
  reportlist[["cellstatus"]] <- cellstatus

  #return whole result list.
  reportlist

}

