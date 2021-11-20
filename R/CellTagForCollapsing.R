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
  
  for.collapse <- as.data.frame(Matrix::summary(umi.matrix))
  for.collapse$i <- rownames(umi.matrix)[for.collapse$i]
  for.collapse$j <- colnames(umi.matrix)[for.collapse$j]

  colnames(for.collapse) <- c("X2", "X1", "value")
  for.collapse$X1 <- as.character(for.collapse$X1)
  for.collapse$X2 <- as.character(for.collapse$X2)
  for.collapse <- for.collapse[which(for.collapse$value > 0), ]
  # Create the contatenation column
  if (length(list.files(celltag.obj@fastq.bam.dir)) > 1) {
    parts.to.paste <- unlist(lapply(strsplit(for.collapse$X2, "_"), function(x) x[2]))
    for.collapse$concat <- paste0(for.collapse$X1, unlist(lapply(strsplit(parts.to.paste, "-"), function(x) x[1])))
    sample.list.prefix <- unique(unlist(lapply(strsplit(for.collapse$X2, "_"), function(x) x[1])))
    r <- apply(as.data.frame(sample.list.prefix), 1, 
               function(x) {
                 for.collapse.sub <- for.collapse[which(startsWith(for.collapse$X2, paste0(x, "_"))), c("concat", "value")]
                 filename.to.save <- paste0(strsplit(output.file, "[.]")[[1]][1], "_", x, ".txt")
                 write.table(for.collapse.sub, filename.to.save, sep = "\t", row.names = F, quote = F, col.names = F)
               })
  } else {
    for.collapse$concat <- paste0(for.collapse$X1, unlist(lapply(strsplit(for.collapse$X2, "-"), function(x) x[1])))
    for.collapse.sub <- for.collapse[, c("concat", "value")]
    write.table(for.collapse.sub, output.file, sep = "\t", row.names = F, quote = F, col.names = F)
  }
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
CellTagDataPostCollapsing <- function(celltag.obj, collapsed.rslt.file, replace.option = FALSE) {
  ultimate.collapsing.df <- data.frame()
  for (i in 1:length(collapsed.rslt.file)) {
    final.collapsing.df <- data.frame()
    # Process this one by one
    curr.file.dir <- collapsed.rslt.file[i]
    print(paste0("Processing ", curr.file.dir))
    # Read in the collpased result
    collapsed <- read.table(curr.file.dir, sep = "\t", header = F, stringsAsFactors = F)
    # Read in the file for collapsing
    if (length(collapsed.rslt.file) > 1) {
      curr.sample.parts <- strsplit(basename(curr.file.dir), "_")[[1]]
      curr.sample <- strsplit(curr.sample.parts[length(curr.sample.parts)], "[.]")[[1]][1]
      collapsing <- celltag.obj@pre.starcode[[celltag.obj@curr.version]]
      collapsing <- collapsing[which(startsWith(collapsing$X2, paste0(curr.sample, "_"))), ]
    } else {
      collapsing <- celltag.obj@pre.starcode[[celltag.obj@curr.version]]
    }
    rownames(collapsing) <- collapsing$concat
    colnames(collapsing)[c(1:2)] <- c("Cell.Barcode", "CellTag")
    new.collapsing.df <- collapsing
    
    cell.bc <- substring(collapsed$V1, 9)
    cell.ct <- substring(collapsed$V1, 1, 8)
    cell.same <- apply(collapsed, 1, 
                       function(x) {
                         cell.bc <- substring(x[1], 9)
                         cell.subset <- strsplit(x[3], ",")[[1]]
                         return(all(endsWith(cell.subset, cell.bc)))
                       })
    
    cell.same.index <- which(cell.same)
    cell.diff.indx <- which(!cell.same)
    
    pb <- txtProgressBar(min = 0, max = length(cell.bc), style = 3)
    pb.count <- 0
    for (csi in cell.same.index) {
      pb.count <- pb.count + 1
      setTxtProgressBar(pb, pb.count)
      
      curr.row <- collapsed[csi,]
      curr.centroid <- curr.row$V1
      curr.count <- curr.row$V2
      curr.ct <- cell.ct[csi]
      
      curr.new.row <- data.frame(row.names = curr.centroid, concat = curr.centroid, CellTag = curr.ct, 
                                 value = curr.count, stringsAsFactors = F)
      
      if (nrow(final.collapsing.df) <= 0){
        final.collapsing.df <- curr.new.row
      } else {
        final.collapsing.df <- rbind(final.collapsing.df, curr.new.row)
      }
    }
    
    for (cdi in cell.diff.indx) {
      pb.count <- pb.count + 1
      setTxtProgressBar(pb, pb.count)
      
      curr.row <- collapsed[cdi,]
      curr.centroid <- curr.row$V1
      curr.count <- curr.row$V2
      curr.collapse.set <- strsplit(curr.row$V3, ",")[[1]]
      curr.ct <- cell.ct[cdi]
      curr.bc <- cell.bc[cdi]
      
      same.concat <- curr.collapse.set[which(endsWith(curr.collapse.set, curr.bc))]
      curr.to.collapse <- setdiff(same.concat, curr.centroid)
      
      if (length(curr.to.collapse) > 0) {
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
        curr.new.row <- data.frame(row.names = curr.centroid, concat = curr.centroid, CellTag = curr.ct,
                                   value = curr.count.new, stringsAsFactors = F)
      }else {
        curr.new.row <- new.collapsing.df[same.concat, c("concat", "CellTag", "value")]
      }
      curr.diff.rows <- new.collapsing.df[setdiff(curr.collapse.set, same.concat), c("concat", "CellTag", "value")]
      
      final.collapsing.df <- rbind(final.collapsing.df, curr.new.row)
      final.collapsing.df <- rbind(final.collapsing.df, curr.diff.rows)
    }

    if (length(which(is.na(final.collapsing.df$concat))) > 0) final.collapsing.df <- final.collapsing.df[-which(is.na(final.collapsing.df$concat)), ]
    rownames(final.collapsing.df) <- final.collapsing.df$concat
    final.collapsing.df <- cbind(final.collapsing.df, collapsing[rownames(final.collapsing.df), c("Cell.Barcode", "concat")])
    
    if (nrow(ultimate.collapsing.df) <= 0) {
      ultimate.collapsing.df <- final.collapsing.df
    } else {
      ultimate.collapsing.df <- rbind(ultimate.collapsing.df, final.collapsing.df)
    }
    rownames(ultimate.collapsing.df) <- NULL
    close(pb)
  }

  df <- transform(ultimate.collapsing.df, Cell.Barcode = factor(Cell.Barcode), CellTag = factor(CellTag))

  celltag.count.sparse <- sparseMatrix(as.integer(df$Cell.Barcode), as.integer(df$CellTag), x = df$value)
  colnames(celltag.count.sparse) <- levels(df$CellTag)
  rownames(celltag.count.sparse) <- levels(df$Cell.Barcode)
  
  # Save the new matrix to the object
  new.obj <- SetCellTagCurrentVersionWorkingMatrix(celltag.obj, "collapsed.count", as(celltag.count.sparse, "dgCMatrix"), replace = replace.option)
  return(new.obj)
}
