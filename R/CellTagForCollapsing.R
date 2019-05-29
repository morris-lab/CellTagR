#' CellTag Starcode Prior Collapsing
#'
#' This function generate the .txt file that will be fed into starcode - https://github.com/gui11aume/starcode - to collapse similar CellTags.
#' @param celltag.obj A CellTag object with the raw count matrix filled.
#' @param output.file The filepath and name to save the table for collapsing (usually a .txt file)
#' @return A CellTag object with collapsing mapping table stored in pre.starcode slot
#' @keywords single-cell RNA-seq data, CellTagging
#' @export
#' @examples
#' CellTagDataForCollapsing(bam.test.obj, "./collapsing.txt")
#' 
CellTagDataForCollapsing <- function(celltag.obj, output.file) {
  # Get the data out from the CellTag object
  umi.matrix <- GetCellTagCurrentVersionWorkingMatrix(celltag.obj, "raw.count")
  for.collapse <- t(as.matrix(umi.matrix))
  # Melt the matrix
  for.collapse <- melt(for.collapse)
  # Subset the matrix to only contain tags with positive UMI numbers
  for.collapse <- subset(for.collapse, value > 0)
  for.collapse$X1 <- as.character(for.collapse$X1)
  for.collapse$X2 <- as.character(for.collapse$X2)
  # Create the contatenation column
  for.collapse$concat <- paste0(for.collapse$X1, unlist(lapply(strsplit(for.collapse$X2, "-"), function(x) x[1])))
  for.collapse.sub <- for.collapse[, c("concat", "value")]
  write.table(for.collapse.sub, output.file, sep = "\t", row.names = F, quote = F, col.names = F)
  # Set CellTag object
  celltag.obj@pre.starcode[[celltag.obj@curr.version]] <- for.collapse
  # Print the path saved
  cat("The file for collapsing is stored at: ", output.file, "\n")
  return(celltag.obj)
}

#' CellTag Starcode Post Collapsing
#'
#' This function processes the result generated from starcode - https://github.com/gui11aume/starcode.
#' @param celltag.obj A CellTag object with the pre-starcode mapping matrix filled.
#' @param collapsed.rslt.file File path to the collapsed result file
#' @return A CellTag object with collapsed count matrix stored in collapsed.count slot
#' @keywords single-cell RNA-seq data, CellTagging
#' @export
#' @examples
#' CellTagDataPostCollapsing(bam.test.obj, "./collapsing_result.txt")
#' 
CellTagDataPostCollapsing <- function(celltag.obj, collapsed.rslt.file) {
  # Read in the collpased result
  collapsed <- read.table(collapsed.rslt.file, sep = "\t", header = F, stringsAsFactors = F)
  # Read in the file for collapsing
  collapsing <- celltag.obj@pre.starcode[[celltag.obj@curr.version]]
  colnames(collapsing)[c(1:2)] <- c("CellTag", "Cell.Barcode")
  new.collapsing.df <- collapsing
  final.collapsing.df <- data.frame()
  # Process the collapsing data file
  for (i in 1:nrow(collapsed)) {
    curr.row <- collapsed[i,]
    curr.centroid <- curr.row$V1
    curr.count <- curr.row$V2
    curr.collapse.set <- strsplit(curr.row$V3, ",")[[1]]
    
    # ind <- which(collapsing$concat == curr.centroid)
    # 
    # curr.new.row <- new.collapsing.df[ind, ]
    # curr.new.row$value <- curr.count
    
    if (length(curr.collapse.set) > 1) {
      curr.collapse.set <- strsplit(curr.row$V3, ",")[[1]]
      curr.to.collapse <- setdiff(curr.collapse.set, curr.centroid)
      for (j in 1:length(curr.to.collapse)) {
        curr.for.c <- curr.to.collapse[j]
        curr.for.c.ct <- substring(curr.for.c, 1, 8)
        if (curr.for.c.ct != curr.ct) {
          ind <- which(collapsing$concat == curr.to.collapse[j])
          ind.cent <- which(collapsing$concat == curr.centroid)
          new.collapsing.df[ind, "concat"] <- curr.centroid
          new.collapsing.df[ind, "CellTag"] <- collapsing[ind.cent[1], "CellTag"]
          new.collapsing.df[ind, "Cell.Barcode"] <- collapsing[ind.cent[1], "Cell.Barcode"]
        }
      }
      curr.centroid.sub <- new.collapsing.df[which(new.collapsing.df$concat == curr.centroid), ]
      curr.count.new <- sum(curr.centroid.sub$value)
      curr.new.row <- data.frame(concat = curr.centroid, CellTag = unique(curr.centroid.sub$CellTag),
                                 Cell.Barcode = unique(curr.centroid.sub$Cell.Barcode), value = curr.count.new,
                                 stringsAsFactors = F)
    } else {curr.new.row <- new.collapsing.df[which(new.collapsing.df$concat == curr.centroid), ]}

    if (nrow(final.collapsing.df) <= 0){
      final.collapsing.df <- curr.new.row
    } else {
      final.collapsing.df <- rbind(final.collapsing.df, curr.new.row)
    }
  }
  
  final.collapsing.df <- setDT(final.collapsing.df)
  # Regenerate the new matrix
  new.matrix <- dcast(final.collapsing.df, Cell.Barcode~CellTag, fill = 0)
  # Give the matrix rownames
  cell.rnm <- new.matrix$Cell.Barcode
  cnms <- colnames(new.matrix)[2:ncol(new.matrix)]
  new.matrix <- as.matrix(new.matrix[, ..cnms])
  rownames(new.matrix) <- cell.rnm
  # Save the new matrix to the object
  new.obj <- SetCellTagCurrentVersionWorkingMatrix(celltag.obj, "collapsed.count", as(new.matrix, "dgCMatrix"))
  return(new.obj)
}
