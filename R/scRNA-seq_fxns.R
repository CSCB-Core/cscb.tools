library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(SoupX)
library(sva)

#' Load Cell Ranger Multi Results into SoupX
#'
#' Takes in 10x Cell Ranger-generated `multi` results and produces a modified table of counts, with background contamination removed.
#' Learn more here \url{https://cran.r-project.org/web/packages/SoupX/SoupX.pdf} or here \url{https://github.com/constantAmateur/SoupX}
#'
#' @param raw file path of raw multi results; found in outs/multi/count/raw_feature_bc_matrix/
#' @param filtered file path of filtered multi results; found in count/sample_filtered_feature_bc_matrix/
#' @param clusters file path of clusters.csv; found in count/analysis/clustering/gene_expression_graphclust/clusters.csv
#'
#' @return SoupX object adjusted for ambient RNA contamination
#' @export
#'
#' @examples
#' gex_multi_soup_raw <-
#'    R.utils::getAbsolutePath(Sys.glob(
#'   file.path(
#'     "../../SCTC-VP-16/Data_Analysis/vinodh-vdj-gex/code/02-multi/runs/*_VDJ_GEX/outs/multi/count/raw_feature_bc_matrix/"
#'   )
#' ))
#' gex_multi_soup_filtered <-
#'   R.utils::getAbsolutePath(Sys.glob(
#'     file.path(
#'       "../../SCTC-VP-16/Data_Analysis/vinodh-vdj-gex/code/02-multi/runs/*_VDJ_GEX/outs/per_sample_outs/*_VDJ_GEX/count/sample_filtered_feature_bc_matrix/"
#'     )
#'   ))
#' gex_multi_soup_clusters <-
#'   R.utils::getAbsolutePath(Sys.glob(
#'     file.path(
#'       "../../SCTC-VP-16/Data_Analysis/vinodh-vdj-gex/code/02-multi/runs/*_VDJ_GEX/outs/per_sample_outs/*_VDJ_GEX/count/analysis/clustering/gene_expression_graphclust/clusters.csv"
#'     )
#'   ))
#'
#' names(gex_multi_soup_raw) <- NULL
#' names(gex_multi_soup_filtered) <- NULL
#' names(gex_multi_soup_clusters) <- NULL
#'
#' gex_obj <- list()
#'
#' for (i in 1:length(gex_multi_soup_raw)) {
#'   gex_obj[[i]] <-
#'     soupx.load.multi(gex_multi_soup_raw[i],
#'                      gex_multi_soup_filtered[i],
#'                      gex_multi_soup_clusters[i])
#' }
soupx.load.multi <- function(raw, filtered, clusters) {
  drops <- Read10X(raw)
  dat.filtered <-
    Read10X(filtered)
  
  sc = SoupChannel(drops, dat.filtered)
  
  cr.clusters <-
    read.csv(clusters)
  clust <- cr.clusters$Cluster
  names(clust) <- cr.clusters$Barcode
  sc = setClusters(sc, clust)
  
  sc = autoEstCont(sc, soupQuantile = 0.5)
  out = adjustCounts(sc)
  return(out)
}

#' Adds layers to Seurat object, percent mitochondrial and ribosomal genes
#'
#' Matches string prefixes to gene symbols ("MT" and "RPS") to find percent mitochondrial and ribosomal genes per cell
#'
#' @param seuratObj A Seurat object with RNA assay
#'
#' @return A Seurat object with two new layers
#' @export
#'
#'
#' @examples
#' for (i in 1:length(gex_obj)) {
#'  gex_obj[[i]] <- get.percent(gex_obj[[i]])
#' }
get.percent <- function(seuratObj) {
  seuratObj[["percent.mt"]] <-
    PercentageFeatureSet(object = seuratObj,
                         pattern = "^MT-",
                         assay = "RNA")
  
  seuratObj[["percent.rps"]] <-
    PercentageFeatureSet(object = seuratObj,
                         pattern = "^RPS*",
                         assay = "RNA")
  
  return(seuratObj)
}


#' Plots mitochondrial and ribosomal genes with ggplot2
#'
#' @param seuratObj A Seurat object with percent mitochondrial and ribosomal genes calculated with get.percent()
#' @param mito.cutoff Cutoff of percent mitochondrial genes, default 10
#' @param rps.cutoff Cutoff of percent ribosomal genes, default 10
#' @param title Title of ggplot, default NULL
#'
#' @return A dual-faceted ggplot printed to active device
#' @export
#'
#' @examples
#'
#' plot.mito(gex_obj[[1]], title = gex_names[[1]])
#' plot.mito(gex_obj[[2]], title = gex_names[[2]])
#' plot.mito(gex_obj[[3]], title = gex_names[[3]])
#' plot.mito(gex_obj[[4]], title = gex_names[[4]])
#' plot.mito(gex_obj[[6]], title = gex_names[[6]])
#' plot.mito(gex_obj[[5]], title = gex_names[[5]])
plot.mito <-
  function(seuratObj,
           mito.cutoff = 10,
           rps.cutoff = 10,
           title = NULL) {
    count_mt <-
      FeatureScatter(seuratObj, feature1 = "nCount_RNA", feature2 = "percent.mt") +
      geom_hline(yintercept = 0.05,
                 linetype = "dashed",
                 color = "black") +
      theme(legend.position = "none") +
      ggtitle(title) + geom_hline(yintercept = mito.cutoff) +
      geom_point(aes(colour = cut(percent.mt, c(0, mito.cutoff))))
    
    count_rps <-
      FeatureScatter(seuratObj, feature1 = "nCount_RNA", feature2 = "percent.rps") +
      geom_hline(yintercept = 0.05,
                 linetype = "dashed",
                 color = "black") +
      theme(legend.position = "none") +
      ggtitle(NULL) + geom_hline(yintercept = rps.cutoff) +
      geom_point(aes(colour = cut(percent.rps, c(0, rps.cutoff))))
    
    
    return((count_mt + count_rps))
  }

#' Generate feature plots of transcript counts vs gene counts
#'
#' @param seuratObj A Seurat object
#' @param count.cutoff A cutoff indicating max counts per cell; default 10k
#' @param downsample A number representing number of cells to subset
#'
#' @return A FeatureScatter printed to device
#' @export
#'
#' @examples
#'
#' seurat.counts(gex_obj[[1]], count.cutoff = 22000) + ggtitle(gex_names[[1]])
#' seurat.counts(gex_obj[[2]], count.cutoff = 52000) + ggtitle(gex_names[[2]])
#' seurat.counts(gex_obj[[3]], count.cutoff = 10000) + ggtitle(gex_names[[3]])
#' seurat.counts(gex_obj[[4]], count.cutoff = 70000) + ggtitle(gex_names[[4]])
#' seurat.counts(gex_obj[[5]], count.cutoff = 65000) + ggtitle(gex_names[[5]])
#' seurat.counts(gex_obj[[6]], count.cutoff = 60000) + ggtitle(gex_names[[6]])
#'
seurat.counts <-
  function(seuratObj,
           count.cutoff = 10000,
           downsample = 0) {
    if (downsample > 0) {
      seuratCounts <- FeatureScatter(
        subset(seuratObj, downsample = downsample),
        feature1 = "nCount_RNA",
        feature2 = "nFeature_RNA"
      ) +
        theme(legend.position = "none") +
        ggtitle(NULL) +
        geom_vline(xintercept = count.cutoff) +
        geom_point(aes(colour = cut(nCount_RNA, c(
          -Inf, count.cutoff
        ))))
      return(seuratCounts)
    }
    else {
      seuratCounts <-
        FeatureScatter(seuratObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
        theme(legend.position = "none") +
        ggtitle(NULL) +
        geom_vline(xintercept = count.cutoff) +
        geom_point(aes(colour = cut(nCount_RNA, c(
          -Inf, count.cutoff
        ))))
      return(seuratCounts)
    }
  }

#' Preprocess Seurat object
#'
#' @param counts Seurat object
#' @param counts_name Name of sample corresponding to seurat object
#' @param count.cutoff Cutoff of counts used in seurat.counts, default 10k
#' @param mito.cutoff Cutoff of percent mitochondrial genes per cell, default 10
#' @param rps.cutoff Cutoff of percent ribosomal genes per cell, default 10
#'
#' @return Seurat object that has been QC'd, SCTransformed, PCA'd, Clustered and UMAP.
#' @export
#'
#' @examples
#' for (i in 1:length(gex_obj)) {
#' gex_obj[[i]] <-
#'   seurat.process(gex_obj[[i]], counts_name = gex_names[[i]])
#' }
seurat.process <-
  function(counts,
           counts_name,
           count.cutoff = 10000,
           mito.cutoff = 10,
           rps.cutoff = 10) {
    if (class(counts)[1] != "Seurat") {
      seuratObj <- CreateSeuratObject(
        counts = counts,
        project = counts_name,
        min.cells = 10,
        min.features = 200
      )
    } else {
      seuratObj <- counts
      seuratObj@active.assay <- "RNA"
    }
    
    # Absolutely imperative that these thresholds are set PER SAMPLE
    seuratObj <-
      subset(
        seuratObj,
        subset = nFeature_RNA > 200 &
          nCount_RNA < count.cutoff &
          percent.mt < mito.cutoff &
          percent.rps < rps.cutoff
      )
    
    # normalize
    
    # TODO expanding limit to prevent some arbitrary error msg
    # options(future.globals.maxSize = 20971520000)
    
    seuratObj <- SCTransform(seuratObj, ncells = 2000)
    seuratObj <- RunPCA(seuratObj)
    
    eb <- ElbowPlot(seuratObj, reduction = "pca")
    
    seuratObj <-
      FindNeighbors(seuratObj, dims = 1:10) # we pick this based on elbow plot
    seuratObj <- FindClusters(seuratObj, resolution = 0.5)
    seuratObj <- RunUMAP(seuratObj, dims = 1:10)
    return(seuratObj)
  }



#' Running paramSweep and summarizeSweep from DoubletFinder
#'
#' @param sample Seurat object
#' @param cores Number of cores, default 1
#'
#' @return Summary states from summarizeSweep
#' @export
#'
#' @examples
#' for (i in 1:length(gex_obj)) {
#'   message(crayon::red$underline$bold(
#'   paste0("Finding Parameter Sweep with DoubletFinder for ", gex_names[[i]])
#' ))
#' gex_pk[[i]] <-
#'   sum.sweep(gex_obj[[i]])
#'   message(crayon::red$underline$bold(paste0("Saving ", gex_names[[i]])))
#'   save(gex_pk, file = here("gex_obj_postDoublet.RData"))
#' }
sum.sweep <- function(sample, cores = 1) {
  message("Parameter sweep...")
  
  if ("future" %in% .packages(TRUE))
    message("Package {future} loaded, setting cores to 1...")
  
  system.time({
    sweep.res.sample <-
      paramSweep_v3(
        sample,
        PCs = 1:10,
        sct = TRUE,
        num.cores = if ("future" %in% .packages(TRUE))
          1
        else
          cores
      )
  })
  
  sweep.stats.sample <-
    summarizeSweep(sweep.res.sample, GT = FALSE)
  
  return(sweep.stats.sample)
}


#' Plots result from find.pk
#'
#' @param sweep.stats.sample Summary stats generated by sum.sweep
#' @param title Plot title, default NA
#'
#' @return R plot printed to active device
#' @export
#'
#' @examples
#'
#' pdf(here("processed/01-5_3-comp/ParamSweepPeaks.pdf"))
#'
#' plot.pk(gex_pk[[1]], title = gex_names[[1]])
#' plot.pk(gex_pk[[2]], title = gex_names[[2]])
#' plot.pk(gex_pk[[3]], title = gex_names[[3]])
#' plot.pk(gex_pk[[4]], title = gex_names[[4]])
#' plot.pk(gex_pk[[5]], title = gex_names[[5]])
#' plot.pk(gex_pk[[6]], title = gex_names[[6]])
#'
#' dev.off()
plot.pk <- function(sweep.stats.sample,  title = NA) {
  bcmvn_sample <- find.pK(sweep.stats.sample)
  x <-
    plot(
      x = bcmvn_sample$ParamID,
      y = bcmvn_sample$MeanAUC,
      pch = 18,
      col = "black",
      cex = 0.75,
      xlab = NA,
      ylab = NA,
      main = title
    )
  x <-
    lines(
      x = bcmvn_sample$ParamID,
      y = bcmvn_sample$MeanAUC,
      col = "black",
      lty = 2
    )
  par(new = TRUE)
  x <-
    plot(
      x = bcmvn_sample$ParamID,
      y = bcmvn_sample$BCmetric,
      pch = 16,
      col = "#41b6c4",
      cex = 0.75
    )
  axis(side = 4)
  x <- lines(x = bcmvn_sample$ParamID,
             y = bcmvn_sample$BCmetric,
             col = "#41b6c4")
  
  return(x)
  
}

#' Predict and remove doublets from Seurat object
#'
#' @param sample Seurat object
#' @param sweep.stats.sample Summary stats from sum.sweep()
#' @param rm.doubs logical; default FALSE; removes doublets from `sample` if TRUE
#'
#' @return A Seurat object with identified doublets
#' @export
#'
#' @examples
#'
#' for (i in 1:length(gex_obj)) {
#'   message(crayon::red$underline$bold(paste0(
#'     "Finding Doublets with DoubletFinder for ", gex_names[[i]]
#'   )))
#'   gex_doubs[[i]] <-
#'     doubFinder(gex_obj[[i]], sweep.stats.sample = gex_pk[[i]])
#'   message(crayon::red$underline$bold(paste0("Saving ", gex_names[[i]])))
#'   save(gex_doubs, file = here("gex_obj_postDoublet_actual.RData"))
#' }
doubFinder <-
  function(sample, sweep.stats.sample, rm.doubs = FALSE) {
    bcmvn_sample <- find.pK(sweep.stats.sample) # save this plot
    
    message("Selecting pK Value...")
    
    peaks <-
      bcmvn_sample[bcmvn_sample$BCmetric > (max(bcmvn_sample$BCmetric) / 2),]
    
    # consider adding trycatch
    # https://stackoverflow.com/questions/12193779/how-to-write-trycatch-in-r
    
    # pay special attention to the different types of boolean operators in R
    # https://medium.com/biosyntax/single-or-double-and-operator-and-or-operator-in-r-442f00332d5b
    for (v in 1:nrow(peaks)) {
      if (nrow(peaks) == 1) {
        pK_val <- as.numeric(as.vector(peaks[v, ]$pK))
        break
      }
      else if (peaks[v, ]$BCmetric > peaks[v + 1, ]$BCmetric |
               is.na(peaks[v + 1, ]$BCmetric)) {
        pK_val <- as.numeric(as.vector(peaks[v, ]$pK))
        break
      }
      else {
        next
      }
    }
    
    message("Homotypic Doublet Proportion Estimation")
    homotypic.prop <-
      modelHomotypic(sample@meta.data$SCT_snn_res.0.5) ## ex: annotations <- sample@meta.data$ClusteringResults
    nExp_poi <-
      round(0.075 * nrow(sample@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
    nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
    
    ## Run DoubletFinder with varying classification stringencies
    ## Remember to update with appropriate params
    # options(future.globals.maxSize = 2097152000)
    
    message("Running doubletFinder...")
    sample <-
      doubletFinder_v3(
        sample,
        PCs = 1:10,
        pN = 0.25,
        pK = pK_val,
        nExp = nExp_poi,
        reuse.pANN = FALSE,
        sct = TRUE
      )
    
    # grep for the column name that has pANN at the start
    pANN_col <-
      colnames(sample@meta.data)[grep("pANN*", colnames(sample@meta.data))]
    
    sample <-
      doubletFinder_v3(
        sample,
        PCs = 1:10,
        pN = 0.25,
        pK = pK_val,
        nExp = nExp_poi.adj,
        reuse.pANN = pANN_col,
        sct = TRUE
      )
    
    class_col <-
      colnames(sample@meta.data)[grep("*classification*", colnames(sample@meta.data))]
    
    class_col <- class_col[length(class_col)]
    
    colnames(sample@meta.data)[ncol(sample@meta.data)] <- "doublet"
    
    if (rm.doubs) {
      sample <- subset(sample, subset = doublet == "Singlet")
    }
    
    return(sample)
    
  }

#' Visualize doublets on UMAP embedding w/ RNA counts for reference
#'
#' @param seuratObj Seurat object with predicted doublets
#' @param title Sample name
#'
#' @return Printed ggplot2 to active device
#' @export
#'
#' @examples
#' pdf(
#'   here("processed/01-5_3-comp/GEX_doubletClusters.pdf"),
#'     height = 8,
#'     width = 14
#'   )
#'
#' viz.doubs(gex_doubs[[1]], title = gex_names[[1]])
#' viz.doubs(gex_doubs[[2]], title = gex_names[[2]])
#' viz.doubs(gex_doubs[[3]], title = gex_names[[3]])
#' viz.doubs(gex_doubs[[4]], title = gex_names[[4]])
#' viz.doubs(gex_doubs[[5]], title = gex_names[[5]])
#' viz.doubs(gex_doubs[[6]], title = gex_names[[6]])
#'
#' dev.off()

viz.doubs <- function(seuratObj, title = "Sample") {
  print(
    DimPlot(seuratObj, reduction = 'umap', group.by = 'doublet') +
      ggtitle(paste0("Doublet Clusters for ", title)) + FeaturePlot(seuratObj,
                                                                    reduction = "umap",
                                                                    features = "nCount_RNA")
  )
}


#' ComBat Batch-Correct SCT Counts in Seurat Object
#'
#' @param seuratObj A Seurat object with an SCT assay
#' @param seuratBatch A Vector of batch labels for each
#'
#' @return Seurat object with combatBatch Assay
#' @export
#'
#' @examples
#' ComBat_merge <-
#'   merge(
#'     rna_seq[[1]],
#'     y = c(
#'       rna_seq[[2]],
#'       rna_seq[[3]],
#'       rna_seq[[4]],
#'       rna_seq[[5]],
#'       rna_seq[[6]],
#'       gex_doubs[[1]],
#'       gex_doubs[[2]],
#'       gex_doubs[[3]],
#'       gex_doubs[[4]],
#'       gex_doubs[[5]],
#'       gex_doubs[[6]]
#'     )
#'   )
#'
#' combat_proccess <-
#'   process.combat(ComBat_merge, ComBat_merge$comp.ident)
process.combat <- function(seuratObj, seuratBatch) {
  pheno <- seuratObj@meta.data
  
  batch <- seuratBatch
  
  edata <-
    as.matrix(GetAssayData(seuratObj, assay = "SCT", slot = "data"))
  
  mod0 <- model.matrix( ~ 1, data = pheno)
  
  combat_edata = ComBat(
    dat = edata,
    batch = batch,
    mod = mod0,
    par.prior = TRUE,
    prior.plots = FALSE
  )
  
  seuratObj[["combatBatch"]] <-
    CreateAssayObject(data = combat_edata)
  
  DefaultAssay(seuratObj) <- "combatBatch"
  
  seuratObj <-
    ScaleData(seuratObj, verbose = TRUE, assay = "combatBatch")
  
  return(seuratObj)
  
}
