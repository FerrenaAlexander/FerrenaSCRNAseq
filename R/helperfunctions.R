# https://r-pkgs.org/whole-game.html

doubletfinderwrapper <- function(seuratobject){


  message('DF Parameter Sweep Completed')

  #param sweep
  sweepres <- DoubletFinder::paramSweep_v3(seu = tmp, PCs = 1:30, sct = T)

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

  dratedf[,1] <- as.numeric(sub("%","",dratedf[,1]))/100

  names(dratedf) <- c('MultipletRate', 'CellsLoaded_100%Viability', 'CellsRecovered')

  dbmodel <- lm(MultipletRate ~ CellsRecovered, data = dratedf)

  predicteddoubletrate <- as.numeric((dbmodel$coefficients['CellsRecovered'] * ncol(tmp)) + dbmodel$coefficients[1])

  homotypicprop <- modelHomotypic(tmp$seurat_clusters)
  nexppoi <- round(predicteddoubletrate * length(rownames(tmp@meta.data)))
  nexppoiadj <- round(nexppoi * (1 - homotypicprop))

  #classify doublets
  message('Initiating third-pass clustering for doublet estimation:\n')

  tmp <- suppressWarnings(doubletFinder_v3(seu = tmp, PCs = 1:30, pN = 0.25, pK = maxscorepk, nExp = nexppoi, sct = T))

  ddf <- data.frame(cells = rownames(tmp@meta.data),
                    orig.ident = tmp$orig.ident,
                    classification = tmp@meta.data[,ncol(tmp@meta.data)])

  #print( table(ddf$classification ))

  ddf <- ddf[ddf$classification == 'Singlet',]

  doublets <- table(ddf$classification )['Doublet']
  singlets <- table(ddf$classification )['Singlet']


  dfoutput <- data.frame(cells_pre_df = cells_pre_df,
                         doublets = doublets,
                         singlets = singlets,
                         row.names = NULL)


  tmp@commands$tamlabscpipeline[['DoubletFinder_results']] <- dfoutput

  tmp <- tmp[,ddf$cells]



}
