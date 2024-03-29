#' Fastq Process Function
#'
#' This function extracts CellTags from the raw fastq sequencing file, provides counts of each CellTag and sorts them in desending order.
#' @param fastq.file The input fastq/bam data directory
#' @param pattern The pattern to seek for
#' @param short.nt.before.tag A short sequence before the 8nt tag to help more specific identification
#' @param short.nt.after.tag A short sequence after the 8nt tag to help more specific identification
#' @return A list contains count table of CellTags. If requested to save fullTag counts, i.e. save.fullTag.counts = TRUE, return a list of both 8nt tags and full sequences count. Otherwise, a list of 8nt tags counts. 
#' @keywords single-cell RNA-seq data, CellTagging
#' @export
#' @examples
#' fastq.process("data.fastq", "CCGGT[ATCG]{8}GAATTC", "CCGGT", "GAATTC")
#' 
fastq.process <- function(fastq.file, pattern, short.nt.before.tag, short.nt.after.tag) {
  con <- file(fastq.file, "r")
  
  # Get the sequences containing the tags (with both full tag region and only 8nt tag)
  seq.list <- c()
  filtered.sequences <- c()
  full.tag.seq <- c()
  only.tag.seq <- c()
  print("Reading File......")
  # Get the size of the bam file
  fq.size <- file.size(fastq.file)
  total <- fq.size/(1000000 * 101)
  # Initialize the progress bar
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  # Initialize the count
  count <- 0
  while(TRUE) {
    curr.lines <- readLines(con, 1000000)
    if (length(curr.lines) == 0) break
    else {
      curr.seqs <- curr.lines[seq(2, 1000000, by = 4)]
      seq.list <- c(seq.list, curr.seqs)
      reg.rslt <- regexpr(pattern, curr.seqs, ignore.case = TRUE, perl = TRUE)
      contain.idx <- which(reg.rslt > 0)
      curr.f.seq <- curr.seqs[contain.idx]
      
      filtered.sequences <- c(filtered.sequences, curr.f.seq)
      start.loc <- reg.rslt[contain.idx]
      end.loc <- start.loc + nchar(short.nt.before.tag) + 8 + nchar(short.nt.after.tag) - 1
      curr.full.tag <- substr(curr.f.seq, start = start.loc, stop = end.loc)
      only.tag <- substr(curr.full.tag, start = (nchar(short.nt.before.tag) + 1), stop = (nchar(short.nt.before.tag) + 8))
      full.tag.seq <- c(full.tag.seq, curr.full.tag)
      only.tag.seq <- c(only.tag.seq, only.tag)
    }
    count <- count + 1
    if (count > total) {
      count <- total
    }
    setTxtProgressBar(pb, count)
  }
  close(con)
  close(pb)
  rslt <- list(full.tag.seq, only.tag.seq)
  return(rslt)
}

#' Bam File Process Function
#'
#' This function extracts CellTags from the bam sequencing file, provides cell barcode, umi and their corresponding celltag information.
#' @param bam.file The input bam data directory
#' @param pattern The pattern to seek for
#' @param short.nt.before.tag A short sequence before the 8nt tag to help more specific identification
#' @param short.nt.after.tag A short sequence after the 8nt tag to help more specific identification
#' @return A data table contains cell barcode, celltag and umi information
#' @keywords single-cell RNA-seq data, CellTagging
#' @export
#' @examples
#' bam.process("data.fastq", "CCGGT[ATCG]{8}GAATTC", "CCGGT", "GAATTC")
#' 
bam.process <- function(bam.file, pattern, short.nt.before.tag, short.nt.after.tag, technique) {
  # Install Rsamtools
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  if (!requireNamespace("Rsamtools", quietly = TRUE)) {
    BiocManager::install("Rsamtools")
  }
  library(Rsamtools)
  # Get the bam file
  bamFile <- BamFile(bam.file)
  # Get the size of the bam file
  bam.size <- file.size(bam.file)
  total <- bam.size/(1000000 * 82.99)
  print(paste0("Reading ", bam.file, " ..."))
  # Initialize the progress bar
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  # Initialize the number of lines to read at once
  yieldSize(bamFile) <- 1000000
  open(bamFile)
  if (tolower(technique) == "10x") {
      parameters <- ScanBamParam(what = scanBamWhat(), tag = c("CB", "GN", "UB", "CR"))
  } else {
      if (tolower(technique) == "dropseq") {
          parameters <- ScanBamParam(what = scanBamWhat(), tag = c("XC", "GN", "XM", "GE"))
      } else {
          if (tolower(technique) == "zumi") {
              parameters <- ScanBamParam(what = scanBamWhat(), tag = c("BC", "GN", "UB", "CR"))
          } else {
              stop("We don't support your current single-cell sequencing technology. Please contact us to add.")
          }
      }
  }
  bam.parsed.df <- data.table()
  count <- 0
  while(TRUE) {
    curr.read <- scanBam(bamFile, param = parameters)[[1]]
#    print(count)
    if (length(curr.read$qname) <= 0) {
      break
    } else {
      # Read in all information
      curr.seqs <- as.character(curr.read$seq)
      # Check if the sequences contain the celltag motif
      reg.rslt <- regexpr(pattern, curr.seqs, ignore.case = TRUE, perl = TRUE)
      contain.idx <- which(reg.rslt > 0)
      if (length(contain.idx) > 0) {
          if (tolower(technique) == "10x") {
            curr.cell.bc <- curr.read$tag$CB
            curr.umi <- curr.read$tag$UB
          } else {
            if (tolower(technique) == "dropseq") {
                curr.cell.bc <- curr.read$tag$XC
                curr.umi <- curr.read$tag$XM
            } else if (tolower(technique) == "zumi") {
                curr.cell.bc <- curr.read$tag$BC
                curr.umi <- curr.read$tag$UB
            }
          }
        curr.cell.tag <- rep(NA, length(curr.read$qname))
        if (!(is.null(curr.cell.bc) | is.null(curr.umi))) {
          # Initialize the current data table
          curr.df <- data.table(Cell.BC = curr.cell.bc, UMI = curr.umi, Cell.Tag = curr.cell.tag)
          curr.f.seq <- curr.seqs[contain.idx]
          start.loc <- reg.rslt[contain.idx]
          end.loc <- start.loc + nchar(short.nt.before.tag) + 8 + nchar(short.nt.after.tag) - 1
          
          curr.full.tag <- substr(curr.f.seq, start = start.loc, stop = end.loc)
          only.tag <- substr(curr.full.tag, start = (nchar(short.nt.before.tag) + 1), stop = (nchar(short.nt.before.tag) + 8))
          
          curr.df$Cell.Tag[contain.idx] <- only.tag
          # Add to the current data frame
          if (nrow(bam.parsed.df) <= 0) {
            bam.parsed.df <- curr.df[contain.idx,]
          } else {
            bam.parsed.df <- rbind(bam.parsed.df, curr.df[contain.idx, ])
          }
        }
      }
    }
    count <- count + 1
    setTxtProgressBar(pb, count)
  }
  close(bamFile)
  close(pb)
  return(bam.parsed.df)
}

#' CellTag Pattern Calling Function
#'
#' This function provides motif patterns corresponding to the input celltag version
#' @param celltag.version Which CellTag version are you investigating?
#' @return A list containing the pattern, nucleotides to look for before/after the motif
#' @keywords single-cell RNA-seq data, CellTagging
#' @export
#' @examples
#' CellTagPatternCalling("v1")
#' 
CellTagPatternCalling <- function(celltag.version) {
  celltag.df <- data.frame(version = c("v1", "v2", "v3"),
                           nt.before.tag = c("GGT", "GTGATG", "TGTACG"),
                           stringsAsFactors = F)
  rownames(celltag.df) <- celltag.df$version
  short.nt.before.tag <- celltag.df[celltag.version, "nt.before.tag"]
  short.nt.after.tag <- "GAATTC"
  
  pattern <- paste0(short.nt.before.tag, "[ATCG]{8}", short.nt.after.tag)
  return(c(pattern, short.nt.before.tag, short.nt.after.tag))
}

#' CellTag Barcode Aggregation function
#'
#' This function allows barcode aggregation of multiple-file processing. 
#' @param file.list files in a list to aggregate in order same as the BAM files
#' @param output.file where to save this aggregated output file. Should be a .tsv file.
#' @return A list containing the pattern, nucleotides to look for before/after the motif
#' @keywords single-cell RNA-seq data, CellTagging
#' @export
#' @examples
#' Barcode.Aggregate(list("barcodes_1.tsv", "barcodes_2.tsv"), output.file = "barcode_aggr.tsv")
#' 
Barcode.Aggregate <- function(file.list, output.file) {
  final.bc <- c()
  for (i in 1:length(file.list)) {
    curr.prefix <- paste0("Sample-", i, "_")
    curr.file <- file.list[[i]]
    curr.bc <- read.table(curr.file, header = F, stringsAsFactors = F)
    bc.to.save <- paste0(curr.prefix, curr.bc[,1])
    final.bc <- c(final.bc, bc.to.save)
  }
  write.table(as.data.frame(final.bc), output.file, sep = "\t", row.names = F, col.names = F, quote = F)
}


GetCellTagCurrentVersionWorkingMatrix <- function(celltag.obj, slot.to.select) {
  curr.mtx <- slot(celltag.obj, slot.to.select)
  if (nrow(curr.mtx) <= 0) {
    return(curr.mtx)
  } else {
    curr.version <- celltag.obj@curr.version
    curr.mtx.sub <- curr.mtx[, which(startsWith(colnames(curr.mtx), curr.version))]
    colnames(curr.mtx.sub) <- gsub(pattern = paste0(curr.version, "."), replacement = "", colnames(curr.mtx.sub))
    full.mtx.sub <- curr.mtx.sub[Matrix::rowSums(is.na(curr.mtx.sub)) != ncol(curr.mtx.sub),]
    
    return(full.mtx.sub)
  }
}

SetCellTagCurrentVersionWorkingMatrix <- function(celltag.obj, slot.to.set, final.to.set, replace = FALSE) {
  cop.final <- final.to.set
  colnames(cop.final) <- paste0(celltag.obj@curr.version, ".", colnames(cop.final))
  curr.version.existing.mtx <- GetCellTagCurrentVersionWorkingMatrix(celltag.obj, slot.to.set)
  
  if (replace) {
    slot(celltag.obj, slot.to.set) <- cop.final
    return(celltag.obj)
  }
  
  if (sum(dim(slot(celltag.obj, slot.to.set))) <= 0) {
    slot(celltag.obj, slot.to.set) <- cop.final
  } else  {
    curr.existing.mtx <- slot(celltag.obj, slot.to.set)
    if (ncol(curr.version.existing.mtx) > 0) {
      curr.ver.exist.colnames <- paste0(celltag.obj@curr.version, ".", colnames(curr.version.existing.mtx))
      indx <- which(colnames(curr.existing.mtx) %in% curr.ver.exist.colnames)
      curr.existing.mtx <- curr.existing.mtx[, -indx]
    }
    new.rownames <- unique(c(rownames(curr.existing.mtx), rownames(cop.final)))

    diff.rnms <- setdiff(new.rownames, rownames(cop.final))
    cop.comp.mtx <- matrix(NA, nrow = length(diff.rnms), ncol = ncol(cop.final))
    rownames(cop.comp.mtx) <- diff.rnms
    colnames(cop.comp.mtx) <- colnames(cop.final)

    diff.rnms.2 <- setdiff(new.rownames, rownames(curr.existing.mtx))
    cem.comp.mtx <- matrix(NA, nrow = length(diff.rnms.2), ncol = ncol(curr.existing.mtx))
    rownames(cem.comp.mtx) <- diff.rnms.2
    colnames(cem.comp.mtx) <- colnames(curr.existing.mtx)

    to.merge.mtx.cop <- rbind(cop.final, cop.comp.mtx)
    to.merge.mtx.cem <- rbind(curr.existing.mtx, cem.comp.mtx)

    if (ncol(to.merge.mtx.cem) <= 0) {
      new.mtx <- to.merge.mtx.cop[,colnames(cop.final)]
    } else {
      new.mtx <- cbind(to.merge.mtx.cop[new.rownames,], to.merge.mtx.cem[new.rownames, ])
    }
    
    slot(celltag.obj, slot.to.set) <- new.mtx
  }
  
  return(celltag.obj)
} 


