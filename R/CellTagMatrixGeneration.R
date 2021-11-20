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
CellTagMatrixCount <- function(celltag.obj, barcodes.file, replace.option = FALSE) {
  # Read in the cell barcodes identified during alignment
  barcodeList <- fread(barcodes.file, header = FALSE)[[1]]
  celltagData <- celltag.obj@bam.parse.rslt[[celltag.obj@curr.version]]
  # Filter based on filtered barcodes
  celltagData <- celltagData[which(celltagData$Cell.BC %in% barcodeList), ]

  #With the parsed CellTag reads loaded we can then easily filter the data and generate UMI Counts for each Cell Barcode/Cell Tag combination.
  #-Groups the data.table by Cell Barcode/Cell Tag combination and creates a new column "UMI.Count" which has the number of unique UMI associated with each Cell Barcode/Cell.Tag combination. uniqueN is equivalent to length(unique(UMI))
  celltagCounts <- celltagData[, .(UMI.Count = uniqueN(UMI)), .(Cell.BC, Cell.Tag)]          
  # The data is now in a long format and needs to be reshaped. We will cast the long data into a wide format resembling a matrix.
  # celltagCountsWide <- dcast(data = celltagCounts, formula = Cell.BC ~ Cell.Tag, value.var = "UMI.Count", fill = 0 )

  #Now we have the data we want in the correct format. Next we can add Cells from the barcode list that were not in the celltagData.
  missingCells <- barcodeList[!(barcodeList %in% celltagCounts$Cell.BC)]
  #Lets make a data.table with one column Cell.BC which will contain a list of the missing cells. This can then be merged with the UMI Count data table.
  missingCells <- setDT(expand.grid(Cell.BC = missingCells, Cell.Tag = unique(celltagCounts$Cell.Tag)))
  missingCells$UMI.Count <- 0
  #Bind the missing cells to the data.table containing the Cell Tag UMI Counts.
  alltagCounts <- rbind(celltagCounts, missingCells, fill = TRUE)
  #Now we can filter out cells which are not in our barcode list.
  alltagCounts <- alltagCounts[Cell.BC %in% barcodeList, ]

  #Generate dgCMatrix
  ### Reference for code
  ## https://datawookie.netlify.app/blog/2016/01/casting-a-wide-and-sparse-matrix-in-r/

  df <- as.data.frame(alltagCounts)
  df <- transform(df, Cell.BC = factor(Cell.BC), Cell.Tag = factor(Cell.Tag))

  celltag.count.sparse <- sparseMatrix(as.integer(df$Cell.BC), as.integer(df$Cell.Tag), x = df$UMI.Count)
  colnames(celltag.count.sparse) <- levels(df$Cell.Tag)
  rownames(celltag.count.sparse) <- levels(df$Cell.BC)

  #Lets also filter Cell Tags in which no UMIs are counted.
  celltagExpr <- Matrix::colSums(celltag.count.sparse)
  tagsRemove <- names(celltagExpr)[celltagExpr == 0]
  alltagCounts[, (tagsRemove):= NULL]

  ## Let's make the dgc matrix again with the tags removed
  df <- as.data.frame(alltagCounts)
  df <- transform(df, Cell.BC = factor(Cell.BC), Cell.Tag = factor(Cell.Tag))

  celltag.count.sparse <- sparseMatrix(as.integer(df$Cell.BC), as.integer(df$Cell.Tag), x = df$UMI.Count)
  colnames(celltag.count.sparse) <- levels(df$Cell.Tag)
  rownames(celltag.count.sparse) <- levels(df$Cell.BC)

  #We now have a final matrix. Next lets generate some stats about the Cell Tags. 
  celltagExpr <- summary(Matrix::colSums(celltag.count.sparse))
  cellsPerTag <- summary(Matrix::colSums(celltag.count.sparse > 0))
  cellExpr <- summary(Matrix::rowSums(celltag.count.sparse))
  
  tagsPerCell <- Matrix::rowSums(celltag.count.sparse > 0)
  tagsPerCellSum <- summary(tagsPerCell)
  
  stats.df <- rbind(celltagExpr, cellsPerTag, cellExpr, tagsPerCellSum)
  rownames(stats.df) <- c("CellTag.UMI.Counts", "Cells.per.CellTag", "Cell.UMI.Counts", "CellTags.per.Cell")
  stats.df <- as.data.frame(stats.df)
  

  dgc.mtx.filter <- celltag.count.sparse
  new.obj <- SetCellTagCurrentVersionWorkingMatrix(celltag.obj, "raw.count", as(dgc.mtx.filter, "dgCMatrix"), replace = replace.option)

  return(new.obj)
}
