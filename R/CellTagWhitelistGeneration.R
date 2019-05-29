#' CellTag Whitelist Filtering Function
#'
#' This function conducts whitelist filtering such that only CellTags with count number over their certain percentile would be considered for clone calling
#' @param celltag.obj A CellTag Object with CellTag frequency table counted and sorted
#' @param percentile A fraction cutoff percentile for filtering the CellTags e.g. 0.9 for 90th percentile
#' @param output.dir Which directory would you like to store these files? If NULL, save to the same directory as the fastq/bam file
#' @return A CellTag Object with attribute (whitelist) filled.
#' @keywords single-cell RNA-seq data, CellTagging
#' @export
#' @examples
#' CellTagWhitelistFiltering(bam.test.obj, 0.9)
#' 
CellTagWhitelistFiltering <- function(celltag.obj, percentile, output.dir = NULL) {
  # Load table and calculate cutoff
  count.sorted.table <- celltag.obj@celltag.freq.stats[[celltag.obj@curr.version]]
  count.cutoff <- quantile(count.sorted.table$Count, probs = percentile)
  count.true.cut <- floor(count.cutoff/10)
  
  # Plot
  plot(count.sorted.table$Count, main="CellTag Whitelist",xlab="CellTag",ylab="Reads")
  abline(v=sum(count.sorted.table$Count >= count.true.cut), col="red", lty=2)
  cat(paste0("Abline Threshold: ", sum(count.sorted.table$Count >= count.true.cut)), "\n")
  
  # Subset the ones pass filtering
  whitelist <- subset(count.sorted.table, Count>=count.true.cut)
  
  if (is.null(output.dir)) output.dir <- paste0(dirname(celltag.obj@fastq.bam.dir), "/", celltag.obj@curr.version, "_whitelist.csv")
  write.csv(whitelist, output.dir, quote = F, row.names = F)
  
  cat("File is saved: ", output.dir, "\n")
  
  celltag.obj@whitelist[[celltag.obj@curr.version]] <- whitelist
  return(celltag.obj)
}

#' CellTag Frequency Sort Table
#'
#' This function counts and sorts the identified CellTags from Fastq file
#' @param celltag.obj A CellTag Object with CellTags extracted
#' @return A CellTag Object with attribute (celltag.freq.stats) filled.
#' @keywords single-cell RNA-seq data, CellTagging
#' @export
#' @examples
#' CellTagWhitelistFiltering(bam.test.obj)
#' 
AddCellTagFreqSort <- function(celltag.obj) {
  # Count the occurrence of each CellTag
  cell.tag.count <- as.data.table(table(celltag.obj@fastq.only.celltag[[celltag.obj@curr.version]]), stringsAsFactors = F)
  # Sort the CellTags in descending order of occurrence
  cell.tag.count.sort <- cell.tag.count[order(-cell.tag.count$N), ]
  colnames(cell.tag.count.sort) <- c("CellTag", "Count")
  # Add to the slot in celltag object
  celltag.obj@celltag.freq.stats[[celltag.obj@curr.version]] <- cell.tag.count.sort
  return(celltag.obj)
}
