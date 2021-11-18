#' Jaccard Analysis Function
#'
#' This function conducts Jaccard analysis to calculate the Jaccard similarity between cells.
#' @param celltag.obj A CellTag object with the counts filtered based on metrics
#' @param plot.corr Would you like to plot the correlation matrix?
#' @return A CellTag object with attribute (jaccard.mtx) filled
#' @keywords single-cell RNA-seq data, CellTagging
#' @export
#' @examples
#' JaccardAnalysis(bam.test.obj)
#'
JaccardAnalysis <- function(celltag.obj, plot.corr = TRUE, fast = FALSE) {
  filtered.whitelised.data <- GetCellTagCurrentVersionWorkingMatrix(celltag.obj, "metric.filtered.count")
  # Calculating the Jaccard matrix
  if (fast) {
    Jac <- proxyC::simil(filtered.whitelised.data, method = "jaccard")
  } else {
    Jac <- proxy::simil(as.matrix(filtered.whitelised.data), method = "Jaccard")
    Jac <- as(Jac, "dsTMatrix")
  }
  
  if ((!fast) & plot.corr) {
    diag(Jac) <- 1
    corrplot(Jac, method="color", order="hclust", hclust.method ="ward.D2", cl.lim=c(0,1), tl.cex=0.1)
  }
  
  celltag.obj@jaccard.mtx <- Jac
  return(celltag.obj)
}

#' Clone Calling Function
#'
#' This function conducts clone calling based on the Jaccard results.
#' @param celltag.obj A CellTag object with the jaccard matrix generated
#' @param correlation.cutoff Correlation cutoff for clone membership
#' @return A CellTag object with attributes (clone.composition & clone.size.info) filled.
#' @keywords single-cell RNA-seq data, CellTagging
#' @export
#' @examples
#' CloneCalling(bam.test.obj, 0.7)
#'
CloneCalling <- function(celltag.obj, correlation.cutoff) {
  Jaccard.Matrix <- celltag.obj@jaccard.mtx
  
  # Using the igraph package to facilitate the identification of membership to each clone
  jac.summ <- Matrix::summary(Jaccard.Matrix)
  jac.lower.i <- jac.summ$j
  jac.summ$j <- jac.summ$i
  jac.summ$i <- jac.lower.i
  lower.tri.summ <- subset(jac.summ, i>j) # Exclude diagnol
  
  test <- sparseMatrix(i = lower.tri.summ$i,
                       j = lower.tri.summ$j,
                       x = lower.tri.summ$x,
                       dims = dim(Jaccard.Matrix),
                       dimnames = dimnames(Jaccard.Matrix))
  
  test.df <- as.data.frame(Matrix::summary(test))
  test.df.sub <- test.df[which(test.df$x > correlation.cutoff), ]
  
  check.corelation <- test.df.sub[,c(1,2)]
  colnames(check.corelation) <- c("row", "col")
  check.corelation <- as.matrix(check.corelation)

  graph.cor <- graph.data.frame(check.corelation, directed = FALSE)
  groups.cor <- split(unique(as.vector(check.corelation)), clusters(graph.cor)$membership)
  conv.groups.cor <- lapply(groups.cor,
                            function(list.cor){
                              rownames(test)[list.cor]})
  
  # Put clones into tables
  l <- seq(1, length(groups.cor))
  df.conv <- apply(as.matrix(l), 1, 
                   function(x) {
                     data.table(clone.id = x, 
                                cell.barcode = conv.groups.cor[[x]])
                   }
  )
  
  df.comb <- rbindlist(df.conv)
  
  # Calculate the size of each clone
  counts <- table(df.comb$clone.id)
  counts <- as.data.frame(counts)
  colnames(counts) <- c("Clone.ID", "Frequency")
  
  celltag.obj@clone.composition[[celltag.obj@curr.version]] <- df.comb
  celltag.obj@clone.size.info[[celltag.obj@curr.version]] <- counts
  return(celltag.obj)
}

