% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scRNA-seq_fxns.R
\name{soupx_load_multi}
\alias{soupx_load_multi}
\title{Load Cell Ranger Multi Results into SoupX}
\usage{
soupx_load_multi(raw, filtered, clusters)
}
\arguments{
\item{raw}{file path of raw multi results; found in outs/multi/count/raw_feature_bc_matrix/}

\item{filtered}{file path of filtered multi results; found in count/sample_filtered_feature_bc_matrix/}

\item{clusters}{file path of clusters.csv; found in count/analysis/clustering/gene_expression_graphclust/clusters.csv}
}
\value{
SoupX object adjusted for ambient RNA contamination
}
\description{
Takes in 10x Cell Ranger-generated \code{multi} results and produces a modified table of counts, with background contamination removed.
Learn more here \url{https://cran.r-project.org/web/packages/SoupX/SoupX.pdf} or here \url{https://github.com/constantAmateur/SoupX}
}
\examples{
gex_multi_soup_raw <-
   R.utils::getAbsolutePath(Sys.glob(
  file.path(
    "../../SCTC-VP-16/Data_Analysis/vinodh-vdj-gex/code/02-multi/runs/*_VDJ_GEX/outs/multi/count/raw_feature_bc_matrix/"
  )
))
gex_multi_soup_filtered <-
  R.utils::getAbsolutePath(Sys.glob(
    file.path(
      "../../SCTC-VP-16/Data_Analysis/vinodh-vdj-gex/code/02-multi/runs/*_VDJ_GEX/outs/per_sample_outs/*_VDJ_GEX/count/sample_filtered_feature_bc_matrix/"
    )
  ))
gex_multi_soup_clusters <-
  R.utils::getAbsolutePath(Sys.glob(
    file.path(
      "../../SCTC-VP-16/Data_Analysis/vinodh-vdj-gex/code/02-multi/runs/*_VDJ_GEX/outs/per_sample_outs/*_VDJ_GEX/count/analysis/clustering/gene_expression_graphclust/clusters.csv"
    )
  ))

names(gex_multi_soup_raw) <- NULL
names(gex_multi_soup_filtered) <- NULL
names(gex_multi_soup_clusters) <- NULL

gex_obj <- list()

for (i in 1:length(gex_multi_soup_raw)) {
  gex_obj[[i]] <-
    soupx_load_multi(gex_multi_soup_raw[i],
                     gex_multi_soup_filtered[i],
                     gex_multi_soup_clusters[i])
}
}
