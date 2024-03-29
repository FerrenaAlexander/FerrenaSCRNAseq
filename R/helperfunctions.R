# https://r-pkgs.org/whole-game.html


#' Create an alluvial plot from long categorical data
#'
#' Wrapper around ggalluvium package for ggplot based alluvial plot, for quickly making alluvial plot from "long", "raw" categorical data (such as Seurat object meta.data), rather than two-way counts of categories.
#'
#'
#'
#'
#'
#' @param labelsdf data.frame with two columns of raw categorical label: for example, each row is a cell (or other observation), and each column is metadata column 1 and metadata column 2
#' @param repel T/F, whether to repel the labels, default = T
#' @param nudge_x numeric, default nudge to the left and right of the repel labels, default 0.2
#' @param ggfittext T/F - whether to use ggfittext, to try to squeeze or remove tiny stratum labels
#' @param ...
#'
#' @return a ggplot object
#' @export
#'
#' @examples
#' #labelsdf can look like this, row.names of cells not needed,
#' # just two columns of categorical data as a data.frame:
#'                      From      To
#' AACCCAAGCATGCGA-1    Malignant  2
#' AAACCCAAGTAGGTTA-1    Malignant  2
#' AAACCCACAAAGCACG-1   Neutrophil  0
#' AAACCCACAGCAGTAG-1    Malignant  2
#' AAACCCACATACCGTA-1    Malignant  2
alluvialplot <- function(labelsdf, repel, nudge_x, ggfittext, ...){


  if( missing(repel) ){ repel = T}
  if( missing(nudge_x) ){nudge_x = 0.3}

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


  } else if(repel==T){

    require(ggrepel)

    ap <- ggplot(longfreqs, aes(y = Freq, axis1=.data[[colnames(longfreqs[1])]], axis2= .data[[colnames(longfreqs[2])]] ) )+
      scale_x_discrete(expand = c(.4, 0))+
      geom_alluvium( aes(fill= .data[[colnames(longfreqs[1])]] ), width = 1/4 ) +
      geom_stratum(width = 1/4) +
      scale_linetype_manual(values = c("blank", "solid")) +

      ggrepel::geom_label_repel(
        aes(label = .data[[colnames(longfreqs[1])]] ),
        stat = "stratum", nudge_x = nudge_x * -1, ...) +

      ggrepel::geom_label_repel(
        aes(label = .data[[colnames(longfreqs[2])]]),
        stat = "stratum", nudge_x = nudge_x, ...) +
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
      mat_ct <- mat[,match(rownames(md_ct), colnames(mat)), drop=F]

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









#' Calculate hemoglobin features in Seurat object
#'
#' @param sobj Seurat object
#' @param hemoglobin.features character vector, hemoglobin gene symbols to search in dataset, default will search all mouse and human hemoglobin gene symbols
#'
#' @return data.frame with percent.hemoglobin in all cells
#' @export
#'
#' @examples
#' sobj[['percent.hemoglobin']] <- Seurat::PercentageFeatureSet(sobj, features = hemoglobin.features)
calculate_percent.hemoglobin <- function(sobj, hemoglobin.features){


  # is missing, check all human / mouse hemoglobin genes
  if( missing(hemoglobin.features) ){
    hemoglobin.features <- c('HBA1', 'HBA2', 'HBB', 'HBD', 'HBE1', 'HBG1', 'HBG2', 'HBM', 'HBQ1', 'HBZ',
                             'Hba', 'Hba-a1', 'Hba-a2', 'Hba-ps3', 'Hba-ps4', 'Hba-x', 'Hbb', 'Hbb-ar', 'Hbb-b1', 'Hbb-b2', 'Hbb-bh0', 'Hbb-bh1', 'Hbb-bh2', 'Hbb-bh3', 'Hbb-bs', 'Hbb-bt', 'Hbb-y')
  }


  message('Calculating percent.hemoglobin')
  #hemoglobin content, add to metadata

  #hemoglobin features, all mouse and human hemoglobin genes are searched for, or hemoglobins are user-provided
  hemoglobin.features <- hemoglobin.features[hemoglobin.features %in% rownames(sobj)]

  hb <- Seurat::PercentageFeatureSet(sobj, features = hemoglobin.features)

  message('Returning percent.hemoglobin as data.frame column\nMake sure to add to sobj$percent.hemoglobin')

  return(hb)

}





#' Make a table that can print to a pdf page
#'
#' @param tabledf data.frame to put as table on pdf
#' @param title title for table
#' @param titlesize title font size, default is 15
#' @param padding whitespace between title and table, default=1
#'
#' @return
#' @export
#'
#' @examples
pdftable <- function(tabledf, title, titlesize, padding){

  if(missing(titlesize)){titlesize <- 15}
  if(missing(padding)){padding <- 1}


  table <- gridExtra::tableGrob(tabledf)


  if(missing(title)){

    grid::grid.newpage()

    return(grid::grid.draw(table))

  } else{


    #set up title and table
    title <- grid::textGrob(label = title,
                            gp=gpar(fontsize=titlesize) )


    #add padding and table
    # https://stackoverflow.com/a/33738678
    padding <- unit(padding,"line")

    table <- gtable::gtable_add_rows(
      table, heights = grid::grobHeight(title) + padding, pos = 0
    )
    table <- gtable::gtable_add_grob(
      table, list(title),
      t = 1, l = 1, r = ncol(table)
    )

    grid::grid.newpage()

    return(grid::grid.draw(table))

  }
}






#' Plot marker heatmap with labelled genes
#'
#' Plot a heatmap where both cell clusters (columns) AND genes (rows) are labelled. The row labels are taken from the output of `Seurat::FindAllMarkers()`.
#'
#' @param sobj seurat object
#' @param markers data.frame in formal of `Seurat::FindAllMarkers()`
#' @param grouping.var column name in `sobj@meta.data` that contains cell grouping info such as cluster or celltype; it should match the cluster. default is 'seurat_clusters'
#' @param assay
#' @param slot
#' @param numgenes
#' @param maxval
#' @param minval
#' @param cellgroup_color_palette
#' @param genegroup_color_palette
#' @param show_heatmap_legend
#' @param show_cellgroup_legend
#' @param show_genegroup_legend
#' @param heatmaptitle
#' @param legendname
#' @param row_names_gp
#' @param row_title_gp
#' @param row_gap
#' @param column_gap
#' @param row_title_rot
#' @param column_title_rot
#' @param column_title_gp
#' @param cluster_columns
#' @param cluster_rows
#' @param use_raster
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
markerheatmap <- function(sobj,
                          markers,
                          grouping.var,
                          assay,
                          slot,
                          numgenes,
                          maxval,
                          minval,
                          cellgroup_color_palette,
                          genegroup_color_palette,
                          show_heatmap_legend,
                          show_cellgroup_legend,
                          show_genegroup_legend,
                          # heatmaptitle,
                          legendname,
                          row_names_gp,
                          row_title_gp,
                          row_gap,
                          column_gap,
                          row_title_rot,
                          column_title_rot,
                          column_title_gp,
                          cluster_columns,
                          cluster_rows,
                          use_raster,
                          ...

){


  if(missing(sobj)){stop('Provide Seurat object')}
  if(missing(markers)){stop("Provide markers in format of Seurat::FindAllMarkers() output")}
  if(missing(grouping.var)){grouping.var <- 'seurat_clusters'}
  if(missing(assay)){assay <- DefaultAssay(sobj)}
  if(missing(slot)){slot <- 'scale.data'}
  if(missing(numgenes)){numgenes <- 5}

  if(missing(maxval)){maxval <- 5}
  if(missing(minval)){minval <- -5}

  if(missing(show_heatmap_legend)){show_heatmap_legend = T}
  if(missing(show_cellgroup_legend)){show_cellgroup_legend = T}
  if(missing(show_genegroup_legend)){show_genegroup_legend = F}

  # if(missing(heatmaptitle)){heatmaptitle = ''}
  if(missing(legendname)){
    if(slot == 'scale.data'){legendname = 'Scaled\nExpression'}
    if(slot == 'data'){legendname = 'Normalized\nExpression'}
    if(slot == 'counta'){legendname = 'Expression'}
  }

  if(missing(row_names_gp)) {row_names_gp = grid::gpar(fontsize = 5)}
  if(missing(row_title_gp)) {row_title_gp = grid::gpar(fontsize = 5)}
  if(missing(row_gap)){ row_gap = unit(0.8, "mm")}
  if(missing(column_gap)){column_gap = unit(0.8, "mm")}
  if(missing(row_title_rot)){row_title_rot = 0}
  if(missing(column_title_rot)){column_title_rot = 45}
  if(missing(column_title_gp)){column_title_gp = grid::gpar(fontsize = 7)}
  if(missing(cluster_columns)){cluster_columns = F}
  if(missing(cluster_rows)){cluster_rows = F}
  if(missing(use_raster)){use_raster = F}



  if(assay == 'integrated'){warning('Using assay "integrated" may cause issues')}

  #set defualt assay as the one given
  DefaultAssay(sobj) <- assay

  require(magrittr)
  require(dplyr)
  require(ComplexHeatmap)

  #get top genes
  n <- numgenes
  top <- markers %>% group_by(cluster) %>% top_n(n = n, wt = avg_log2FC)
  genes <- top$gene

  #heatmap of cluster markers


  #make sure all genes are in if using scale.data
  if(slot == 'scale.data'){

    #if using SCT, use "getresidual" to add the gnenes to scale.data
    if(assay == 'SCT'){

      if( any( !(genes %in% rownames(sobj@assays$SCT@scale.data)) ) ){

        #get missing genes
        missinggenes <- genes[!(genes %in% rownames(sobj@assays$SCT@scale.data))]

        warning('The following genes are missing from assay SCT slot scale.data, will try to recover using GetResidual() function:',
                '\n',
                paste(missinggenes,collapse = ', ') )

        #try getresidual...
        sobj <- GetResidual(sobj, missinggenes, na.rm = F, replace.value = T)

        #it can be complicated doing this after integration, some genes are NAs...
        scgem <- sobj@assays$SCT@scale.data

        if( any( !complete.cases(scgem) ) ){
          scgem <- scgem[complete.cases(scgem),]
          top <- top[top$gene %in% rownames(scgem),]
          sobj@assays$SCT@scale.data <- scgem
        }
        rm(scgem)

      }

    }

    if(assay == 'RNA'){

      if( any( !(genes %in% rownames(sobj@assays$RNA@scale.data)) ) ){

        #try get missing genes
        missinggenes <- genes[!(genes %in% rownames(sobj@assays$RNA@scale.data))]

        warning('The following genes are missing from assay RNA slot scale.data, will try to recover using ScaleData() function:',
                '\n',
                paste(missinggenes,collapse = ', ') )

        # try to recover
        sobj <- ScaleData(sobj, features = missinggenes)

        #it can be complicated doing this after integration, some genes are NAs...
        scgem <- sobj@assays$RNA@scale.data

        if( any( !complete.cases(scgem) ) ){
          scgem <- scgem[complete.cases(scgem),]
          top <- top[top$gene %in% rownames(scgem),]
          sobj@assays$RNA@scale.data <- scgem
        }
        rm(scgem)

      }

    }


  }



  #prep markers
  top <- top[top$gene %in% rownames(sobj),]

  # get gem
  gem <- GetAssayData(sobj, assay = assay, slot = slot)

  #subset gem for just the markers
  gem <- gem[match(top$gene, rownames(gem)),]

  #annot for clusters
  #set up grouping var tmp column in md
  md <- sobj@meta.data
  md$AF_groupingvar <- md[,grouping.var]

  #factorize if not, using order from marker data
  if(!is.factor(md$AF_groupingvar)){
    md$AF_groupingvar <- factor(md$AF_groupingvar, levels = unique(markers$cluster) )
  }


  #order md by grouping var
  md <- md[order(md$AF_groupingvar),]

  #order gem using ordered md
  gem <- gem[,match(rownames(md), colnames(gem))]


  ### heatmap column and row annotatons for cells and genes
  ## column annots
  #get named vector: values are grouping,s names are barcodes
  clust_bc <- setNames(md$AF_groupingvar,
                       nm = rownames(md)
  )

  #set up column color
  if(missing(cellgroup_color_palette)){
    COLORPAL <- scales::hue_pal()(length(levels(md$AF_groupingvar)))
  } else{
    COLORPAL <- cellgroup_color_palette
  }


  col_clust <- setNames(COLORPAL,
                        nm = levels(md$AF_groupingvar))

  #set up column cluster color

  ha_clust <- ComplexHeatmap::HeatmapAnnotation(Group = clust_bc, col = list(Group = col_clust),
                                                name = grouping.var,
                                                show_legend = show_cellgroup_legend)



  ##annot for rows (markers)

  #facotrize if not
  if(!is.factor(top$cluster)){
    top$cluster <- factor(top$cluster, levels=unique(top$cluster))
  }

  #make sure gene order matches features
  gem <- gem[match(top$gene, rownames(gem)),]

  #set up named vector: calues are clusters, names are genes
  ct_gene <- setNames(top$cluster,
                      nm=top$gene)

  #set up row colors;
  # if not provided, we use same as column colors,
  # but make sure to remove missing cluster marker colors if cluster had no markers
  if(missing(genegroup_color_palette)){
    col_gene <- col_clust
    col_gene <- col_gene[names(col_gene) %in% top$cluster]
  } else{

    COLORPAL <- genegroup_color_palette
    col_gene <- setNames(COLORPAL,
                         nm = levels(top$cluster))

  }



  #get row annot
  ha_genes <- ComplexHeatmap::rowAnnotation(Group = ct_gene, col = list(Group = col_gene),
                                            show_annotation_name=F,
                                            show_legend = show_genegroup_legend)




  #restrict range
  gem[gem>maxval] <- maxval
  gem[gem<minval] <- minval



  #actual heatmap
  hm <- ComplexHeatmap::Heatmap(gem,

                                column_labels = rep('', ncol(gem)),
                                column_split = md$AF_groupingvar,
                                row_split = top$cluster,

                                top_annotation = ha_clust,
                                left_annotation = ha_genes,


                                column_title = '',
                                name = legendname,

                                row_names_gp = row_names_gp,
                                row_title_gp = row_title_gp,
                                row_gap = row_gap,
                                column_gap = column_gap,
                                row_title_rot = row_title_rot,
                                column_title_rot = column_title_rot,
                                column_title_gp = column_title_gp,
                                cluster_columns = cluster_columns,
                                cluster_rows = cluster_rows,
                                use_raster = use_raster,
                                ...
  )



  pdf(NULL) #prevent plotted output
  hm <- draw(hm, merge_legend=T,
             show_heatmap_legend = show_heatmap_legend)
  dev.off()

  return(hm)

}





### try to collapse clusters in seurat object with some cutoff of correlation?
# what's the reference for this?
# "extend" the collapsing?
# A may be cor with B, and B with C, but A not with C...
# collapse A, B and C?

# #collapse clusters: make cor matrix
# avgs <- AverageExpression(sobj, assays = 'SCT', return.seurat = F, slot = 'data')$SCT
#
# #pairwise correlation
# cormat <- cor(as.matrix(avgs))
#
# #collapse them:
# mask <- cormat>correlation_threshold_collapse_hires_clusters
#
# #for each cluster, record which ones have cor > threshold
# clusts <- colnames(mask)
# clusts_to_collapse <- lapply(clusts, function(clust){
#
#   #get correlations of this cluster
#   maskcol <- mask[,clust]
#
#
#   #get any clusts above thres
#   if( any(maskcol) == T ){
#     return( names(maskcol[maskcol==T]) )
#   } else{
#     return()
#   }
#
# })
# names(clusts_to_collapse) <- clusts
#
#
# #check if list elements are null; if so they have no clusters above thres, so drop them
# clusts_to_collapse <- clusts_to_collapse[ sapply(clusts_to_collapse, function(x){length(x)>1}) ]
#
#
# # if there are clusters to collapse, we should get them...
# identdf <- data.frame(orig = names(clusts_to_collapse),
#                       collapse = sapply(clusts_to_collapse,function(x){paste(x,collapse = '_')}) )
#
# identdf
