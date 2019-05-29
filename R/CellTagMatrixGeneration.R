#' CellTag Matrix Generation Function
#'
#' This function uses the extract information from data processed before and generate a Cell Barcode x CellTag matrix
#' @param celltag.obj A CellTag object with bam file result filled
#' @param barcodes.file A .tsv output file from 10x CellRanger pipeline. It contains a list of all cell barcodes identified in the filtered dataset.
#' @return A CellTag object with the attribute (raw.count) filled
#' @keywords single-cell RNA-seq data, CellTagging
#' @export
#' @examples
#' CellTagMatrixCount(bam.test.obj, "barcodes.tsv")
#' 
CellTagMatrixCount <- function(celltag.obj, barcodes.file) {
  # Read in the cell barcodes identified during alignment
  barcodeList <- fread(barcodes.file, header = FALSE)[[1]]
  celltagData <- celltag.obj@bam.parse.rslt[[celltag.obj@curr.version]]
  # Filter based on filtered barcodes
  celltagData <- celltagData[which(celltagData$Cell.BC %in% barcodeList), ]
  
  #With the parsed CellTag reads loaded we can then easily filter the data and generate UMI Counts for each Cell Barcode/Cell Tag combination.
  #-Groups the data.table by Cell Barcode/Cell Tag combination and creates a new column "UMI.Count" which has the number of unique UMI associated with each Cell Barcode/Cell.Tag combination. uniqueN is equivalent to length(unique(UMI))
  celltagCounts <- celltagData[, .(UMI.Count = uniqueN(UMI)), .(Cell.BC, Cell.Tag)]          
  #The data is now in a long format and needs to be reshaped. We will cast the long data into a wide format resembling a matrix.
  celltagCountsWide <- dcast(data = celltagCounts, formula = Cell.BC ~ Cell.Tag, value.var = "UMI.Count", fill = 0 )
  
  #Now we have the data we want in the correct format. Next we can add Cells from the barcode list that were not in the celltagData.
  missingCells <- barcodeList[!(barcodeList %in% celltagCountsWide$Cell.BC)]
  #Lets make a data.table with one column Cell.BC which will contain a list of the missing cells. This can then be merged with the UMI Count data table.
  missingCells <- data.table(Cell.BC = missingCells)
  #Bind the missing cells to the data.table containing the Cell Tag UMI Counts.
  alltagCounts <- rbind(celltagCountsWide, missingCells, fill = TRUE)
  #We have added the missing cells, whose values now need to be changed from NA to 0.
  alltagCounts[is.na(alltagCounts)] <- 0
  #Now we can filter out cells which are not in our barcode list.
  alltagCounts <- alltagCounts[Cell.BC %in% barcodeList, ]
  #Generate dgCMatrix
  cols <- colnames(alltagCounts)
  cols.tags <- cols[2:length(cols)]
  dgc.mtx <- as.matrix(alltagCounts[, ..cols.tags])
  rownames(dgc.mtx) <- alltagCounts$Cell.BC
  
  #Lets also filter Cell Tags in which no UMIs are counted.
  celltagExpr <- colSums(alltagCounts[, -1])
  tagsRemove <- names(celltagExpr)[celltagExpr == 0]
  alltagCounts[, (tagsRemove):= NULL]
  
  #We now have a final matrix. Next lets generate some stats about the Cell Tags. 
  celltagExpr <- summary(colSums(alltagCounts[, -1]))
  cellsPerTag <- summary(colSums(alltagCounts[, -1] > 0))
  cellExpr <- summary(rowSums(alltagCounts[, -1]))
  
  tagsPerCell <- rowSums(alltagCounts[, -1] > 0)
  tagsPerCellSum <- summary(tagsPerCell)
  
  stats.df <- rbind(celltagExpr, cellsPerTag, cellExpr, tagsPerCellSum)
  rownames(stats.df) <- c("CellTag.UMI.Counts", "Cells.per.CellTag", "Cell.UMI.Counts", "CellTags.per.Cell")
  stats.df <- as.data.frame(stats.df)
  
  celltag.obj@celltag.stats[[celltag.obj@curr.version]] <- stats.df
  
  dgc.mtx.filter <- dgc.mtx[, which(colSums(dgc.mtx) > 0)]
  new.obj <- SetCellTagCurrentVersionWorkingMatrix(celltag.obj, "raw.count", as(dgc.mtx.filter, "dgCMatrix"))
  return(new.obj)
}
