#' Metric-Base Filtering Function
#'
#' This function applies further filtering on scRNA-seq data with CellTags based on cutoff values identified from the metric plots.
#' @param celltag.obj A CellTag Object with count matrix generated
#' @param cutoff The cutoff decided from the metric plots
#' @param comparison Would you like to maintain the part less than/greater than the cutoff? Default to less. Choices can be greater or less.
#' @return A CellTag Object with attribute (metric.filtered.count) filled
#' @keywords single-cell RNA-seq data, CellTagging
#' @export
#' @examples
#' MetricBasedFiltering(bam.test.object, 20, "less")
#'
MetricBasedFiltering <- function(celltag.obj, cutoff, comparison = "less") {
  whitelisted.ct.data <- GetCellTagCurrentVersionWorkingMatrix(celltag.obj, "whitelisted.count")
  metric.filter.ct.data <- GetCellTagCurrentVersionWorkingMatrix(celltag.obj, "metric.filtered.count")
  if (ncol(metric.filter.ct.data) <= 0) {
    whitelisted.celltag.data <- as.matrix(whitelisted.ct.data)
  } else {
    whitelisted.celltag.data <- as.matrix(metric.filter.ct.data)
  }
  # Set up the filtering data frame
  CellTags.per.cell.whitelisted.pf <- as.data.frame(Matrix::rowSums(whitelisted.celltag.data))
  
  # Set up the filtered celltag dataset object
  if (comparison == "less") {
    cell.filter <- subset(CellTags.per.cell.whitelisted.pf, CellTags.per.cell.whitelisted.pf < (cutoff))
  } else {
    cell.filter <- subset(CellTags.per.cell.whitelisted.pf, CellTags.per.cell.whitelisted.pf > (cutoff))
  }
  cell.bc.filter <- row.names(cell.filter)
  # Filter celltag dataset
  celltags.whitelisted.new <- whitelisted.celltag.data[cell.bc.filter, ]

  new.obj <- SetCellTagCurrentVersionWorkingMatrix(celltag.obj, "metric.filtered.count", as(celltags.whitelisted.new, "dgCMatrix"))

  return(new.obj)
}

#' CellTag Metric Plotting Function
#'
#' This function provides some metric plots for further downstream celltag filtering in the scRNA-seq dataset.
#' @param celltag.obj A CellTag Object
#' @keywords single-cell RNA-seq data, CellTagging
#' @export
#' @examples
#' MetricPlots(bam.test.obj)
#'
MetricPlots <- function(celltag.obj) {
  
  obj.metric.filtered.count <- GetCellTagCurrentVersionWorkingMatrix(celltag.obj, "metric.filtered.count")
  obj.whitelisted.count <- GetCellTagCurrentVersionWorkingMatrix(celltag.obj, "whitelisted.count")
  
  if (ncol(obj.metric.filtered.count) <= 0) {
    if (ncol(obj.whitelisted.count) <= 0) {
      celltag.data <- GetCellTagCurrentVersionWorkingMatrix(celltag.obj, "binary.mtx")
    } else {
      celltag.data <- obj.whitelisted.count
    }
  } else {
    celltag.data <- obj.metric.filtered.count
  }
  
  CellTags.per.cell.whitelisted.pf <- Matrix::rowSums(celltag.data)
  CellTags.per.cell.avg <- mean(CellTags.per.cell.whitelisted.pf)
  CellTags.frequency.whitelisted.pf <- Matrix::colSums(celltag.data)
  CellTags.freq.avg <- mean(CellTags.frequency.whitelisted.pf)
  par(mfrow=c(2,2))
  plot(CellTags.per.cell.whitelisted.pf, main = "CellTag Counts of Individual Cells", xlab = "Cell Index", ylab = "CellTag Counts")
  plot(CellTags.frequency.whitelisted.pf, main = "CellTag Occurrence Frequency Across All Cells", xlab = "Cell Index", ylab = "CellTag Frequency")
  hist(CellTags.per.cell.whitelisted.pf, main = "Histogram of CellTag Counts of Individual Cells", xlab = "CellTag Counts", ylab = "Count")
  hist(CellTags.frequency.whitelisted.pf, main = "Histogram of CellTag Occurrence Frequency Across All Cells", xlab = "CellTag Occurrence Frequency", ylab = "Count")
  cat("Average: ", CellTags.per.cell.avg, "\n")
  cat("Frequency: ", CellTags.freq.avg, "\n")
}

