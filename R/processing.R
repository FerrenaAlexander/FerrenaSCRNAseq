# https://r-pkgs.org/whole-game.html




# Automated filtering for Seurat object ---------------------------
#'Basic SC Pipeline for QC and clustering
#'
#' Given a Seurat object, perform QC filtering.
#' This is best to do on an object that has not undergone other QC/filtering.
#'
#' @param seuratobject A Seurat object
#' @param baselinefilter.mad T/F; whether to perform "global" QC, ie without pre-clustering. Default = False.
#' @param baseline.mito.filter T/F; whether to perform global maximum mitochondrial content filtration using median absolute deviation threshold; will only work if baselinefilter.mad is set to True. Default is True.
#' @param madmax.dist.percentmito.baseline a numeric, or a string reading 'predict'. If numeric is provided, will use this as median absolute deviation threshold for global mito cutoff. If set to 'predict', will attempt to learn cutoff from data. Default = 'predict'
#' @param baseline.libsize.filter T/F; whether to perform global minimal lib size filtration using median absolute deviation threshold; will only work if baselinefilter.mad is set to True. Default is True.
#' @param madmax.dist.nCount_RNA.baseline a numeric, or a string reading 'predict'. If numeric is provided, will use this as median absolute deviation threshold for global libsize cutoff. If set to 'predict', will attempt to learn cutoff from data. Default = 'predict'
#' @param removemitomaxclust T/F ; whether to identify and remove abnormally high mitochondrial content clusters after first-pass clustering; if no baseline mito filtration is used, there will almost certainly be a mito-driven cluster. Identification is via Grubbs' test for outliers, based on Lukasz Komsta's implementation in the outliers package. Default is True.
#' @param iterativefilter T/F ; whether to perform iterative filtering on first-pass filtering. Default is True.
#' @param iterativefilter.libsize either one of two strings ('twosided' or 'lefttail') or False. Twosided will attempt to learn median abs. dev. cutoffs for both min and max to catch debris and, ostensibly, doublets. Lefttail will attempt to learn median abs. dev. threshold for min to catch debris; doublets may not accurately be captured by max cutoffs as this is more an artifact of sequencing than cell suspension. False will skip. Default is 'lefttail'.
#' @param iterativefilter.mito T/F; whether to learn right-tail median abs. dev. thresholds and filter maximal mitochondrial content from each cluster. May incorrectly remove mito-okay cells while missing true mito-hi cells. Default = F.
#' @return will return a bunch of plots related to QC and an output in the form of a Seurat object to the standard out.
#' @examples
#' \dontrun{
#' pdf('qcplots.pdf')
#' sobj <- seuratpipeline('datafilepath.h5', format=h5)
#' dev.off()
#' }
automatedfiltering <- function(
  seuratobject,
  baselinefilter.mad,
  baseline.mito.filter,
  madmax.dist.percentmito.baseline,
  baseline.libsize.filter,
  madmax.dist.nCount_RNA.baseline,

  removemitomaxclust,
  iterativefilter,
  iterativefilter.libsize,
  iterativefilter.libsize.twosided,
  iterativefilter.libsize.lefttail,
  iterativefilter.mito){

  # baseline (global) filtration
  if( missing( baselinefilter.mad )) {baselinefilter.mad <- F}
  if( missing( baseline.mito.filter )) {baseline.mito.filter <- T}
  if( missing( madmax.dist.percentmito.baseline )) {madmax.dist.percentmito.baseline <- 'predict'}

  if( missing( baseline.libsize.filter )) {baseline.libsize.filter <- T}
  if( missing( madmax.dist.nCount_RNA.baseline )) {madmax.dist.nCount_RNA.baseline <- 'predict'}

  # cluster based filtration
  if( missing( iterativefilter )) {iterativefilter <- T}
  if( missing( removemitomaxclust )) {removemitomaxclust <- T}
  if( missing( iterativefilter.libsize )) {iterativefilter.libsize <- 'lefttail'}
  if( missing( iterativefilter.mito )) {iterativefilter.mito <- F}

  if( missing( cellcycleregression )) {cellcycleregression <- F}





  message('\tnum genes = ', nrow(tmp) , '\n',
          '\tnum cells = ', ncol(tmp)  , '\n')


  tmp@commands$ferrenascpipeline <- list()

  tmp@commands$ferrenascpipeline[['allcommands']] <-
    data.frame(Command = c('baselinefilter.mad',
                           'baseline.mito.filter', 'baseline.libsize.filter',

                           'removemitomaxclust',

                           'iterativefilter',
                           'iterativefilter.libsize', 'iterativefilter.mito',

                           'cellcycleregression'

    ),



    Option = c(baselinefilter.mad,
               baseline.mito.filter, baseline.libsize.filter,

               removemitomaxclust,

               iterativefilter,
               iterativefilter.libsize, iterativefilter.mito,

               cellcycleregression

    )
    )



  #mito and ribo content, add to metadata
  #mito
  mito.features <- grep(pattern = "^mt-", x = rownames(x = tmp), value = TRUE, ignore.case = T)

  tmp[["percent.mito"]] <- PercentageFeatureSet(tmp, features = mito.features)

  #ribo --> remove sincethese are ribo protein genes...
  # ribogenes <- readxl::read_excel('~/Documents/tam/scripts/ribogenes.xlsx')
  # goodribo <- ribogenes$genesymbol[ribogenes$genesymbol %in% rownames(tmp)]
  #
  # tmp[["percent.ribo"]] <- PercentageFeatureSet(tmp, features = goodribo)


  #checkpoint_prebaseline <- tmp



  #whether to perform baseline filtering at all.
  if(baselinefilter.mad == T) {

    #find mito cutoff high
    if(baseline.mito.filter == T){

      message('Initiating global baseline mito hi filtration')

      x <- tmp$percent.mito
      m         <- median(x)
      abs.dev   <- abs(x - m)

      right.mad <- median(abs.dev[x>=m])

      x.mad <- rep(right.mad, length(x))

      mad.distance <- abs(x - m) / x.mad
      mad.distance[x==m] <- 0

      #select bad cells (outside of MAD cutoff)
      if(madmax.dist.percentmito.baseline == 'predict'){

        message('Predicting baseline mito MAD threshold...')

        #predict mad cutoff. for right tail, x > m
        y=sort(mad.distance[x>m])
        ecpout <- ecp::e.divisive(diff(as.matrix(y)), k = 2, min.size = 2)
        #this step is non-trivial... may have to play around with k value and which estimate to use.


        ecpp = ecpout$estimates[2]
        madmax.dist.percentmito.baseline = as.numeric(y[ecpout$estimates[2]])
        ecpplot <- data.frame(x = 1:length(y), y = y) #should have called this ecpdf...

        ecpplot1 <- ggplot(data=ecpplot, aes(x = x, y = y)) +
          geom_point(shape=1, size = 0.1) +
          geom_hline(yintercept = madmax.dist.percentmito.baseline, col = 'red', linetype = 'dashed') +
          geom_vline(xintercept = ecpp, col = 'red', linetype = 'dashed')+
          geom_point(mapping = aes(x = ecpp, y = madmax.dist.percentmito.baseline),col = 'red') +
          theme_linedraw()+
          xlab(label = 'Cells') + ylab(label = 'Med. Abs. Dev.')


        ecpplot2 <- ggplot(data=ecpplot, aes(x = 0, y = y)) +
          geom_violin(fill = 'steelblue')+
          geom_jitter(height = 0, width = 0.25, size = 0.1)+
          geom_hline(yintercept = madmax.dist.percentmito.baseline, col = 'red', linetype = 'dashed')+
          theme_linedraw()+
          xlab(label = 'Cells') + ylab(label = 'Med. Abs. Dev.')+
          scale_x_discrete(limits=0)
        #scaling x this way so the line is at the same level in the cowplot

        title <- ggdraw() +
          draw_label(
            paste0('ECP Prediction = ', ecpp),
            fontface = 'bold',
            x = 0,
            hjust = 0
          ) +
          theme(
            # add margin on the left of the drawing canvas,
            # so title is aligned with left edge of first plot
            plot.margin = margin(0, 0, 0, 7)
          )

        plotrow <- cowplot::plot_grid(ecpplot1, ecpplot2)

        print(
          plot_grid(
            title, plotrow,
            ncol = 1,
            # rel_heights values control vertical title margins
            rel_heights = c(0.1, 1)
          )
        )

        rm(y, ecpout, ecpplot, ecpplot1, ecpplot2, ecpp, plotrow, title)
      } #end predict block

      bad <- x[mad.distance > madmax.dist.percentmito.baseline]

      #for plotting / message output: get mito upper bound
      mitohi <- ifelse( test = length(bad[bad > median(x)]) == 0,
                        yes = Inf,
                        no = min(bad[bad > median(x)])
      )



      gm <-   ggplot(tmp@meta.data, aes(x = 0, y = percent.mito))+
        geom_violin(fill='steelblue')+
        geom_jitter(height = 0, width = 0.25, size = 0.1)+
        geom_hline(yintercept = mitohi, col = 'red', linetype = 'dashed')+
        theme_linedraw()+
        theme(axis.text.x=element_blank(),
              axis.ticks.x=element_blank())+
        labs(title = paste0("percent mito cutoff"),
             subtitle = paste0(length(rownames(tmp@meta.data)), " cells total",
                               "; cutoff for percent mito is ", as.character(round(mitohi, digits = 3))),
             caption = paste0('Mito max MAD = ', round(madmax.dist.percentmito.baseline, 3) ,
                              '\npercent.mito cutoff = ', round(mitohi, digits = 3),
                              '\nNum cells presubset = ', ncol(tmp),
                              '\nNum cells remaining = ', length(tmp$percent.mito[tmp@meta.data$percent.mito < mitohi]), '\n'))+
        xlab('Cells')



      #save to log and print
      tmp@commands$ferrenascpipeline[['baseline_mito_filter']] <- gm
      print(gm)

      filteredcells <- names(x[!names(x) %in% names(bad)] )


      message('\nMito max MAD = ', round(madmax.dist.percentmito.baseline, 3) ,
              '\npercent.mito cutoff = ', round(mitohi, digits = 3),
              '\nNum cells presubset = ', ncol(tmp),
              '\nNum cells remaining = ', length(filteredcells), '\n')


      tmp <- subset(x = tmp,
                    cells = filteredcells)

      rm(abs.dev, bad, right.mad, m, x, x.mad, mad.distance,filteredcells,
         madmax.dist.percentmito.baseline, mitohi)


    } #end mito baseline block



    #find libsize cutoff low.
    if(baseline.libsize.filter == T){

      message('Initiating global baseline lib size low filtration')

      x <- tmp$nCount_RNA
      m         <- median(x)
      abs.dev   <- abs(x - m)

      left.mad <- median(abs.dev[x<=m])

      x.mad <- rep(left.mad, length(x))

      mad.distance <- abs(x - m) / x.mad
      mad.distance[x==m] <- 0

      #select bad cells (outside of MAD cutoff)
      if(madmax.dist.nCount_RNA.baseline == 'predict'){

        message('Predicting baseline lib size threshold...')

        #predict mad cutoff. for left tail, x < m...
        y=sort(mad.distance[x<m])
        ecpout <- ecp::e.divisive(diff(as.matrix(y)), k = 1, min.size = 2)
        #this step is non-trivial... may have to play around with k value and which estimate to use.


        ecpp = ecpout$estimates[2]
        madmax.dist.nCount_RNA.baseline = as.numeric(y[ecpout$estimates[2]])
        ecpplot <- data.frame(x = 1:length(y), y = y) #should have called this ecpdf...

        ecpplot1 <- ggplot(data=ecpplot, aes(x = x, y = y)) +
          geom_point(shape=1, size = 0.1) +
          geom_hline(yintercept = madmax.dist.nCount_RNA.baseline, col = 'red', linetype = 'dashed') +
          geom_vline(xintercept = ecpp, col = 'red', linetype = 'dashed')+
          geom_point(mapping = aes(x = ecpp, y = madmax.dist.nCount_RNA.baseline),col = 'red') +
          theme_linedraw()+
          xlab(label = 'Cells') + ylab(label = 'Med. Abs. Dev.')


        ecpplot2 <- ggplot(data=ecpplot, aes(x = 0, y = y)) +
          geom_violin(fill = 'steelblue')+
          geom_jitter(height = 0, width = 0.25, size = 0.1)+
          geom_hline(yintercept = madmax.dist.nCount_RNA.baseline, col = 'red', linetype = 'dashed')+
          theme_linedraw()+
          xlab(label = 'Cells') + ylab(label = 'Med. Abs. Dev.')+
          scale_x_discrete(limits=0)
        #scaling x this way so the line is at the same level in the cowplot

        title <- ggdraw() +
          draw_label(
            paste0('ECP Prediction = ', ecpp),
            fontface = 'bold',
            x = 0,
            hjust = 0
          ) +
          theme(
            # add margin on the left of the drawing canvas,
            # so title is aligned with left edge of first plot
            plot.margin = margin(0, 0, 0, 7)
          )

        plotrow <- cowplot::plot_grid(ecpplot1, ecpplot2)

        print(
          plot_grid(
            title, plotrow,
            ncol = 1,
            # rel_heights values control vertical title margins
            rel_heights = c(0.1, 1)
          )
        )

        rm(y, ecpout, ecpplot, ecpplot1, ecpplot2, ecpp, plotrow, title)
      }

      bad <- x[mad.distance > madmax.dist.nCount_RNA.baseline]


      ncountslo <- ifelse( test = length(bad[bad < median(x)]) == 0,
                           yes = 0,
                           no = max(bad[bad < median(x)])
      )

      #plot cutoffs for library size / UMIs
      g_lib_hist <- ggplot(tmp@meta.data) +
        geom_histogram(aes(nCount_RNA),
                       color="black", fill = "steelblue",
                       binwidth = 0.05)+
        #geom_vline(xintercept = ncountshi, linetype = "dashed", colour = "red")+
        geom_vline(xintercept = ncountslo, linetype = "dashed", colour = "red")+
        geom_vline(xintercept = median(tmp@meta.data$nCount_RNA),
                   linetype = "dotted", colour = "red", size = 1.2)+
        scale_x_log10()+
        labs(title = paste0("Library Size per Cell"),
             x = "Library size (UMI, aka 'nCount'), log10 scale",
             y = "Number of cells" ,
             subtitle = paste0(length(rownames(tmp@meta.data)), " cells; median lib size = ",
                               median(tmp@meta.data$nCount_RNA)),
             caption = paste0("cells remaining = ",
                              length(x) - length(bad)
             )
        )+
        theme_linedraw()

      g_lib_vln <- ggplot(tmp@meta.data, aes(x = 0, y = nCount_RNA))+
        geom_violin(fill='steelblue')+
        geom_jitter(height = 0, width = 0.25, size = 0.1)+
        geom_hline(yintercept = ncountslo, col = 'red', linetype = 'dashed')+
        theme_linedraw()+
        scale_y_log10()+
        theme(axis.text.x=element_blank(),
              axis.ticks.x=element_blank())+
        labs(title = paste0("LibSize Cutoff"),
             subtitle = paste0(length(rownames(tmp@meta.data)), " cells total",
                               "; cutoff for libsize low is ", as.character(round(ncountslo, digits = 3))),
             caption = paste0('Mito max MAD = ', round(madmax.dist.nCount_RNA.baseline, 3) ,
                              '\npercent.mito cutoff = ', round(ncountslo, digits = 3),
                              '\nNum cells presubset = ', ncol(tmp),
                              '\nNum cells remaining = ', length(tmp$percent.mito[tmp@meta.data$nCount_RNA > ncountslo]), '\n'))+
        xlab('Cells')

      print(g_lib_hist)
      print(g_lib_vln)

      tmp@commands$ferrenascpipeline[['baseline_libsize_filer_hist']] <- g_lib_hist
      tmp@commands$ferrenascpipeline[['baseline_libsize_filer_vln']] <- g_lib_vln



      filteredcells <- names(x[!names(x) %in% names(bad)] )


      message('Libsize global max MAD = ', round(madmax.dist.nCount_RNA.baseline, 3) ,
              '\nLibSize Low Cutoff = ', ncountslo,
              '\nNum cells presubset = ', ncol(tmp),
              '\nNum cells remaining = ', length(filteredcells), '\n')

      #save cutoff params...


      tmp <- subset(x = tmp,
                    cells = filteredcells)

      rm(abs.dev, bad, left.mad, m, x, x.mad, mad.distance, ncountslo, madmax.dist.nCount_RNA.baseline, filteredcells)


    }

  }













  message('Initiating first-pass clustering\n')

  #normalize and cluster, first pass, for lib size filtration
  suppressWarnings(seuratobject <- SCTransform(seuratobject, verbose = T))

  seuratobject <- RunPCA(object = seuratobject, verbose = F)

  seuratobject <- FindNeighbors(object = seuratobject, dims = 1:30, verbose = F)
  seuratobject <- FindClusters(object = seuratobject, resolution = 0.1, verbose = F)


  #record params for all clusters...
  md <- seuratobject@meta.data
  avgs <- aggregate(data = md, cbind(nCount_RNA , nFeature_RNA , percent.mito) ~ seurat_clusters, FUN=mean)
  rm(md)

  seuratobject@commands$ferrenascpipeline[['pre_iterative-filter_vln']] <- VlnPlot(seuratobject, c('percent.mito', 'nCount_RNA', 'nFeature_RNA'), ncol=1, pt.size = 0.1)


  #checkpoint_postbaseline <- seuratobject



  ### remove the mito cluster and recluster ###
  if(removemitomaxclust == T){
    message('\nSearching for Mito High clusters')

    avgsseuratobject <- avgs
    gt <- grubbs(avgsseuratobject$percent.mito)

    #grubbs test while loop for right-tail outliers, w/ bonferroni correction
    iter=1
    while(
      str_split(gt$alternative, ' ')[[1]][1] == 'highest' &
      gt$p.value < (0.05 / nrow(avgsseuratobject))
    ) {
      message('\tgrubbs outlier test iteration: ',iter)

      avgsseuratobject <- avgsseuratobject[-which.max(avgsseuratobject$percent.mito),]

      gt <- grubbs(avgsseuratobject$percent.mito)

      iter=iter+1
    }

    mitohiclusts <- as.vector(avgs$seurat_clusters[!(avgs$seurat_clusters %in% avgsseuratobject$seurat_clusters)])

    #sometimes no mito clusts are detected, so this if-esle avoids crashes in that case
    if( length(mitohiclusts) != 0 ) {
      mitoclustcells <- WhichCells(seuratobject, idents = mitohiclusts)

      seuratobject <- seuratobject[,!(colnames(seuratobject) %in% mitoclustcells)]

      message('Reclustering without mito clusters\n')
      suppressWarnings(seuratobject <- SCTransform(seuratobject, verbose = T))

      seuratobject <- RunPCA(object = seuratobject, verbose = F)

      seuratobject <- FindNeighbors(object = seuratobject, dims = 1:30, verbose = F)
      seuratobject <- FindClusters(object = seuratobject, resolution = 1.0, verbose = T)

      seuratobject@commands$ferrenascpipeline[['mitohiclust-result']] <- VlnPlot(seuratobject, c('percent.mito', 'nCount_RNA', 'nFeature_RNA'), ncol=1, pt.size = 0.1)


    } else{
      message('\nNo Mito High clusters identified by Grubbs test\n')

      seuratobject@commands$ferrenascpipeline[['mitohiclust-result']] <- paste0("No Mito High clusters identified by Grubbs test")

    }
  }







  message('Initiating lib-size QC of first-pass clustering\n')


  ########################### filtration loop ####################################



  if(iterativefilter == T) {





    #start recording
    params <- list()




    #loop through clusters and get cells with proper library size, based on MAD
    filteredcells <- c()

    for(clust in levels(seuratobject$seurat_clusters) ) {        #print(clust) }

      message('Filtering cluster ', clust)
      seuratobjectclust <- subset(seuratobject, idents = clust)


      if(ncol(seuratobjectclust) < 100) {
        message('  -Under 100 cells; will not subset')
      }


      ### LIBSIZE FILTER ###

      #two sided, left sided, or False;
      #if false, else statement will just add all cells to filtered cells obj


      if(iterativefilter.libsize ==  'twosided'){

        #calculate median absolute deviation (MAD) of library size (nCount_RNA)
        x <- seuratobjectclust$nCount_RNA
        m <- median(x)

        #get absolute deviation from the median
        abs.dev   <- abs(x - m)

        #get left and right side of absolute deviation
        left.mad  <- median(abs.dev[x<=m])
        right.mad <- median(abs.dev[x>=m])
        two.sided.mad <- c(left.mad, right.mad)

        #for each cell, assign to tail of median
        x.mad <- rep(two.sided.mad[1], length(x))
        x.mad[x > m] <- two.sided.mad[2]

        #for each cell, get distance from the MAD for the cell's assigned tail
        mad.distance <- abs(x - m) / x.mad
        mad.distance[x==m] <- 0

        #predict mad cutoff. for left tail, x < m...
        y=sort(mad.distance[x<m])
        ecpout <- ecp::e.divisive(diff(as.matrix(y)), k = 1, min.size = 2)
        #this step is non-trivial... may have to play around with k value and which estimate to use.


        ecpp = ecpout$estimates[2]
        madmax.dist.nCount_RNA = as.numeric(y[ecpout$estimates[2]])
        ecpplot <- data.frame(x = 1:length(y), y = y) #should have called this ecpdf...

        ecpplot1 <- ggplot(data=ecpplot, aes(x = x, y = y)) +
          geom_point(shape=1, size = 0.1) +
          geom_hline(yintercept = madmax.dist.nCount_RNA, col = 'red', linetype = 'dashed') +
          geom_vline(xintercept = ecpp, col = 'red', linetype = 'dashed')+
          geom_point(mapping = aes(x = ecpp, y = madmax.dist.nCount_RNA),col = 'red') +
          theme_linedraw()+
          xlab(label = 'Cells') + ylab(label = 'Med. Abs. Dev.')


        ecpplot2 <- ggplot(data=ecpplot, aes(x = 0, y = y)) +
          geom_violin(fill = 'steelblue')+
          geom_jitter(height = 0, width = 0.25, size = 0.1)+
          geom_hline(yintercept = madmax.dist.nCount_RNA, col = 'red', linetype = 'dashed')+
          theme_linedraw()+
          xlab(label = 'Cells') + ylab(label = 'Med. Abs. Dev.')+
          scale_x_discrete(limits=0)
        #scaling x this way so the line is at the same level in the cowplot

        title <- ggdraw() +
          draw_label(
            paste0('Cluster: ', clust, '\nLibsize MAD threshold prediction = ', ecpp),
            fontface = 'bold',
            x = 0,
            hjust = 0
          ) +
          theme(
            # add margin on the left of the drawing canvas,
            # so title is aligned with left edge of first plot
            plot.margin = margin(0, 0, 0, 7)
          )

        plotrow <- cowplot::plot_grid(ecpplot1, ecpplot2)

        print(
          plot_grid(
            title, plotrow,
            ncol = 1,
            # rel_heights values control vertical title margins
            rel_heights = c(0.1, 1)
          )
        )

        rm(y, ecpout, ecpplot, ecpplot1, ecpplot2, ecpp, plotrow, title)


        #select bad cells (outside of MAD cutoff)
        bad <- x[mad.distance > madmax.dist.nCount_RNA]


        #for plotting: get hi and lo bounds for histogram based on hi and lo cells
        ncountshi <- ifelse( test = length(bad[bad > median(x)]) == 0,
                             yes = Inf,
                             no = min(bad[bad > median(x)])
        )

        ncountslo <- ifelse( test = length(bad[bad < median(x)]) == 0,
                             yes = 0,
                             no = max(bad[bad < median(x)])
        )

        #plot cutoffs for library size / UMIs
        g1 <- ggplot(seuratobjectclust@meta.data) +
          geom_histogram(aes(nCount_RNA),
                         color="black", fill = "steelblue",
                         binwidth = 0.05)+
          geom_vline(xintercept = ncountshi, linetype = "dashed", colour = "red")+
          geom_vline(xintercept = ncountslo, linetype = "dashed", colour = "red")+
          geom_vline(xintercept = median(seuratobjectclust@meta.data$nCount_RNA),
                     linetype = "dotted", colour = "red", size = 1.2)+
          scale_x_log10()+
          labs(title = paste0("Library Size per Cell for Cluster", clust),
               x = "Library size (UMI, aka 'nCount'), log10 scale",
               y = "Number of cells" ,
               subtitle = paste0(length(rownames(seuratobjectclust@meta.data)), " cells; median lib size = ",
                                 median(seuratobjectclust@meta.data$nCount_RNA)),
               caption = paste0("cells remaining = ",
                                length(x) - length(bad)
               )
          )+
          theme_linedraw()




        g2 <-  ggplot(seuratobjectclust@meta.data, aes(x = 0, y = nCount_RNA))+
          geom_violin(fill='steelblue')+
          geom_jitter(height = 0, width = 0.25, size = 0.1)+
          geom_hline(yintercept = ncountslo, col = 'red', linetype = 'dashed')+
          geom_hline(yintercept = ncountshi, col = 'red', linetype = 'dashed')+
          geom_hline(yintercept = median(seuratobjectclust@meta.data$nCount_RNA),
                     linetype = "dotted", colour = "red", size = 1.2)+
          theme_linedraw()+
          scale_y_log10()+
          theme(axis.text.x=element_blank(),
                axis.ticks.x=element_blank())+
          labs(title = paste0("Library Size per Cell for Cluster", clust),
               subtitle = paste0(length(rownames(seuratobjectclust@meta.data)), " cells total",
                                 "; cutoff for percent mito is ", as.character(round(ncountslo, digits = 3))),
               caption = paste0('libsize max MAD = ', round(madmax.dist.nCount_RNA, 3) ,
                                '\nNum cells presubset = ', ncol(seuratobjectclust),
                                '\nNum cells remaining = ',
                                length(x) - length(bad), '\n'))+
          xlab('Cells') + ylab('LibSize, Log10 scale')



        title <- ggdraw() +
          draw_label(
            paste0('Cluster: ', clust, '\nLibsize cutoffs'),
            fontface = 'bold',
            x = 0,
            hjust = 0
          ) +
          theme(
            # add margin on the left of the drawing canvas,
            # so title is aligned with left edge of first plot
            plot.margin = margin(0, 0, 0, 7)
          )

        plotrow <- cowplot::plot_grid(g1, g2)

        print(
          plot_grid(
            title, plotrow,
            ncol = 1,
            # rel_heights values control vertical title margins
            rel_heights = c(0.1, 1)
          )
        )

        rm(g1,g2,title, plotrow)

        #add good cells to iterator list
        filteredclust <- names(x[!names(x) %in% names(bad)] )

        #save params
        #params[[clust]] <- data.frame(libsizemad = madmax.dist.nCount_RNA,
        #                              libsizelo = ncountslo,
        #                              libsizehi = ncountshi)

        #remove objects
        rm(abs.dev, bad, left.mad, right.mad, m, x, x.mad,
           two.sided.mad, mad.distance,
           ncountshi, ncountslo,
           madmax.dist.nCount_RNA)


        cells_after_libsizefilt <- length(filteredclust)

      } #end two sided iterative lib size filt





      if(iterativefilter.libsize == 'lefttail'){

        #calculate median absolute deviation (MAD) of library size (nCount_RNA)
        x <- seuratobjectclust$nCount_RNA
        m <- median(x)

        #get absolute deviation from the median
        abs.dev   <- abs(x - m)

        left.mad  <- median(abs.dev[x<=m])

        x.mad <- rep(left.mad, length(x))

        #for each cell, get distance from the MAD for the cell's assigned tail
        mad.distance <- abs(x - m) / x.mad
        mad.distance[x==m] <- 0

        #predict mad cutoff. for left tail, x < m...
        y=sort(mad.distance[x<m])
        ecpout <- ecp::e.divisive(diff(as.matrix(y)), k = 1, min.size = 2)
        #this step is non-trivial... may have to play around with k value and which estimate to use.


        ecpp = ecpout$estimates[2]
        madmax.dist.nCount_RNA = as.numeric(y[ecpout$estimates[2]])
        ecpplot <- data.frame(x = 1:length(y), y = y) #should have called this ecpdf...

        ecpplot1 <- ggplot(data=ecpplot, aes(x = x, y = y)) +
          geom_point(shape=1, size = 0.1) +
          geom_hline(yintercept = madmax.dist.nCount_RNA, col = 'red', linetype = 'dashed') +
          geom_vline(xintercept = ecpp, col = 'red', linetype = 'dashed')+
          geom_point(mapping = aes(x = ecpp, y = madmax.dist.nCount_RNA),col = 'red') +
          theme_linedraw()+
          xlab(label = 'Cells') + ylab(label = 'Med. Abs. Dev.')


        ecpplot2 <- ggplot(data=ecpplot, aes(x = 0, y = y)) +
          geom_violin(fill = 'steelblue')+
          geom_jitter(height = 0, width = 0.25, size = 0.1)+
          geom_hline(yintercept = madmax.dist.nCount_RNA, col = 'red', linetype = 'dashed')+
          theme_linedraw()+
          xlab(label = 'Cells') + ylab(label = 'Med. Abs. Dev.')+
          scale_x_discrete(limits=0)
        #scaling x this way so the line is at the same level in the cowplot

        title <- ggdraw() +
          draw_label(
            paste0('Cluster: ', clust, '\nLibsize MAD threshold prediction = ', ecpp),
            fontface = 'bold',
            x = 0,
            hjust = 0
          ) +
          theme(
            # add margin on the left of the drawing canvas,
            # so title is aligned with left edge of first plot
            plot.margin = margin(0, 0, 0, 7)
          )

        plotrow <- cowplot::plot_grid(ecpplot1, ecpplot2)

        print(
          plot_grid(
            title, plotrow,
            ncol = 1,
            # rel_heights values control vertical title margins
            rel_heights = c(0.1, 1)
          )
        )

        rm(y, ecpout, ecpplot, ecpplot1, ecpplot2, ecpp, plotrow, title)


        #select bad cells (outside of MAD cutoff)
        bad <- x[mad.distance > madmax.dist.nCount_RNA]


        #for plotting: get lo bounds for histogram based on lo cells

        ncountslo <- ifelse( test = length(bad[bad < median(x)]) == 0,
                             yes = 0,
                             no = max(bad[bad < median(x)])
        )

        #plot cutoffs for library size / UMIs
        g1 <- ggplot(seuratobjectclust@meta.data) +
          geom_histogram(aes(nCount_RNA),
                         color="black", fill = "steelblue",
                         binwidth = 0.05)+
          geom_vline(xintercept = ncountslo, linetype = "dashed", colour = "red")+
          geom_vline(xintercept = median(seuratobjectclust@meta.data$nCount_RNA),
                     linetype = "dotted", colour = "red", size = 1.2)+
          scale_x_log10()+
          labs(title = paste0("Library Size per Cell for Cluster", clust),
               x = "Library size (UMI, aka 'nCount'), log10 scale",
               y = "Number of cells" ,
               subtitle = paste0(length(rownames(seuratobjectclust@meta.data)), " cells; median lib size = ",
                                 median(seuratobjectclust@meta.data$nCount_RNA)),
               caption = paste0("cells remaining = ",
                                length(x) - length(bad)
               )
          )+
          theme_linedraw()




        g2 <-  ggplot(seuratobjectclust@meta.data, aes(x = 0, y = nCount_RNA))+
          geom_violin(fill='steelblue')+
          geom_jitter(height = 0, width = 0.25, size = 0.1)+
          geom_hline(yintercept = ncountslo, col = 'red', linetype = 'dashed')+
          geom_hline(yintercept = median(seuratobjectclust@meta.data$nCount_RNA),
                     linetype = "dotted", colour = "red", size = 1.2)+
          theme_linedraw()+
          scale_y_log10()+
          theme(axis.text.x=element_blank(),
                axis.ticks.x=element_blank())+
          labs(title = paste0("Library Size per Cell for Cluster", clust),
               subtitle = paste0(length(rownames(seuratobjectclust@meta.data)), " cells total"),
               caption = paste0('libsize max MAD = ', round(madmax.dist.nCount_RNA, 3) ,
                                '\nNum cells presubset = ', ncol(seuratobjectclust),
                                '\nNum cells remaining = ',
                                length(x) - length(bad), '\n'))+
          xlab('Cells') + ylab('LibSize, Log10 scale')



        title <- ggdraw() +
          draw_label(
            paste0('Cluster: ', clust, '\nLibsize cutoffs'),
            fontface = 'bold',
            x = 0,
            hjust = 0
          ) +
          theme(
            # add margin on the left of the drawing canvas,
            # so title is aligned with left edge of first plot
            plot.margin = margin(0, 0, 0, 7)
          )

        plotrow <- cowplot::plot_grid(g1, g2)

        print(
          plot_grid(
            title, plotrow,
            ncol = 1,
            # rel_heights values control vertical title margins
            rel_heights = c(0.1, 1)
          )
        )

        rm(g1,g2,title, plotrow)

        #add good cells to iterator list
        filteredclust <- names(x[!names(x) %in% names(bad)] )

        cells_after_libsizefilt <- length(filteredclust)

        #save params
        # params[[clust]] <- data.frame(libsizemad = madmax.dist.nCount_RNA,
        #                                libsizelo = ncountslo)

        #remove objects
        rm(abs.dev, bad, left.mad, m, x, x.mad,
           mad.distance,  ncountslo,
           madmax.dist.nCount_RNA)

      } #end left iterative lib size LEFT ONLY filt

      if(iterativefilter.libsize ==  F){
        filteredclust <- colnames(seuratobject)
      }

      #end iterative lib size filt








      ### MITO FILTER ###

      if(iterativefilter.mito == T){

        ### find mito cutoff high ###
        x <- seuratobjectclust$percent.mito
        m         <- median(x)
        abs.dev   <- abs(x - m)

        right.mad <- median(abs.dev[x>=m])

        x.mad <- rep(right.mad, length(x))

        mad.distance <- abs(x - m) / x.mad
        mad.distance[x==m] <- 0

        #predict mad cutoff. for right tail, x > m
        y=sort(mad.distance[x>m])
        ecpout <- ecp::e.divisive(diff(as.matrix(y)), k = 3, min.size = 2)
        #this step is non-trivial... may have to play around with k value and which estimate to use.


        ecpp = ecpout$estimates[2]
        madmax.dist.percentmito = as.numeric(y[ecpout$estimates[2]])
        ecpplot <- data.frame(x = 1:length(y), y = y) #should have called this ecpdf...

        ecpplot1 <- ggplot(data=ecpplot, aes(x = x, y = y)) +
          geom_point(shape=1, size = 0.1) +
          geom_hline(yintercept = madmax.dist.percentmito, col = 'red', linetype = 'dashed') +
          geom_vline(xintercept = ecpp, col = 'red', linetype = 'dashed')+
          geom_point(mapping = aes(x = ecpp, y = madmax.dist.percentmito),col = 'red') +
          theme_linedraw()+
          xlab(label = 'Cells') + ylab(label = 'Med. Abs. Dev.')


        ecpplot2 <- ggplot(data=ecpplot, aes(x = 0, y = y)) +
          geom_violin(fill = 'steelblue')+
          geom_jitter(height = 0, width = 0.25, size = 0.1)+
          geom_hline(yintercept = madmax.dist.percentmito, col = 'red', linetype = 'dashed')+
          theme_linedraw()+
          xlab(label = 'Cells') + ylab(label = 'Med. Abs. Dev.')+
          scale_x_discrete(limits=0)
        #scaling x this way so the line is at the same level in the cowplot

        title <- ggdraw() +
          draw_label(
            paste0('Cluster: ', clust, '\nMito MAD threshold prediction = ', ecpp),
            fontface = 'bold',
            x = 0,
            hjust = 0
          ) +
          theme(
            # add margin on the left of the drawing canvas,
            # so title is aligned with left edge of first plot
            plot.margin = margin(0, 0, 0, 7)
          )

        plotrow <- cowplot::plot_grid(ecpplot1, ecpplot2)

        print(
          plot_grid(
            title, plotrow,
            ncol = 1,
            # rel_heights values control vertical title margins
            rel_heights = c(0.1, 1)
          )
        )

        rm(y, ecpout, ecpplot, ecpplot1, ecpplot2, ecpp, plotrow, title)


        #select bad cells (outside of MAD cutoff)
        bad <- x[mad.distance > madmax.dist.percentmito]



        #for plotting / message output: get mito upper bound
        mitohi <- ifelse( test = length(bad[bad > median(x)]) == 0,
                          yes = Inf,
                          no = min(bad[bad > median(x)])
        )



        m1 <- ggplot(seuratobjectclust@meta.data, aes(x = 0, y = percent.mito))+
          geom_violin(fill='steelblue')+
          geom_jitter(height = 0, width = 0.25, size = 0.1)+
          geom_hline(yintercept = mitohi, col = 'red', linetype = 'dashed')+
          theme_linedraw()+
          theme(axis.text.x=element_blank(),
                axis.ticks.x=element_blank())+
          labs(title = paste0("percent mito cutoff for cluster ", clust),
               subtitle = paste0(length(rownames(seuratobjectclust@meta.data)), " cells total",
                                 "; cutoff for percent mito is ", as.character(round(mitohi, digits = 3))),
               caption = paste0('Mito max MAD = ', round(madmax.dist.percentmito, 3) ,
                                '\npercent.mito cutoff = ', round(mitohi, digits = 3),
                                '\nNum cells presubset = ', ncol(seuratobjectclust),
                                '\nNum cells remaining = ', length(seuratobjectclust$percent.mito[seuratobjectclust@meta.data$percent.mito < mitohi]), '\n'))+
          xlab('Cells')

        title <- ggdraw() +
          draw_label(
            paste0('Cluster: ', clust, '\nMito Hi Cutoff'),
            fontface = 'bold',
            x = 0,
            hjust = 0
          ) +
          theme(
            # add margin on the left of the drawing canvas,
            # so title is aligned with left edge of first plot
            plot.margin = margin(0, 0, 0, 7)
          )

        plotrow <- cowplot::plot_grid(m1)

        print(
          plot_grid(
            title, plotrow,
            ncol = 1,
            # rel_heights values control vertical title margins
            rel_heights = c(0.1, 1)
          )
        )

        rm(m1, title, plotrow)

        #save mito filtered cells
        mitofilteredcells <- names(x[!names(x) %in% names(bad)] )


        #remove objects
        rm(abs.dev, bad, right.mad, m, x, x.mad, mad.distance,mitohi, madmax.dist.percentmito)


        cells_after_mitofilt <- length(mitofilteredcells)

      } else{mitofilteredcells <- colnames(seuratobject)} #end iterative mito filt



      #get the good cell barcodes for this cluster
      filteredclust <- intersect(filteredclust, mitofilteredcells)

      cells_after_bothfilt <- length(filteredclust)

      paramdf <- data.frame(Cluster = clust,
                            Cells_prefilt = ncol(seuratobjectclust)
      )

      #save cell numbers pre and post
      if(iterativefilter.libsize != F){
        paramdf$Cells_after_libsizefilt <- cells_after_libsizefilt
      }
      if(iterativefilter.mito != F) {
        paramdf$Cells_after_mitofilt <- cells_after_mitofilt
      }

      paramdf$Cells_after_bothfilt_theoretical <- cells_after_bothfilt





      #if under 100, do not actually subset.
      if(ncol(seuratobjectclust) > 100) {

        #barcodes of this cluster to global good barcode list
        filteredcells <- c(filteredcells, filteredclust)

        paramdf$Cells_after_bothfilt_true <- length(filteredclust)

      } else{

        filteredcells <- c(filteredcells, colnames(seuratobjectclust))

        paramdf$Cells_after_bothfilt_true <- ncol(seuratobjectclust)
      }

      rm(filteredclust, mitofilteredcells, seuratobjectclust)


      params[[as.character(clust)]] <- paramdf

    } #close for clust loop






    finalparamdf <- bind_rows(params)

    reportingdf <- data.frame(Cells_prefilt = sum(finalparamdf$Cells_prefilt),
                              cells_postfilt = sum(finalparamdf$Cells_after_bothfilt_true))

    seuratobject@commands$ferrenascpipeline[['iterativefilter_results']] <- reportingdf

    seuratobject@commands$ferrenascpipeline[['iterativefilter_details']] <- finalparamdf

    seuratobject <- seuratobject[,filteredcells]

    seuratobject@commands$ferrenascpipeline[['iterativefilter_post-vln']] <- VlnPlot(seuratobject, c('percent.mito', 'nCount_RNA', 'nFeature_RNA'), ncol=1, pt.size = 0.1)





  } #close iterative filter conditional loop


}
