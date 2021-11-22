#' CellTag Extraction Function
#'
#' This function extracts CellTags from the raw fastq/bam sequencing file. If it is a fastq file, provides counts of each CellTag and sorts them in desending order. If it is a bam file, returns the barcode, umi, celltag information.
#' @param celltag.obj A CellTag object initialized with path to the fastq/bam file
#' @param celltag.version The CellTag version to extract
#' @param technique The technique used for scRNA-seq, Default to 10x. Currently enabled for 10x and dropseq.
#' @return A CellTag object with attribute (bam.parse.rslt) filled
#' @keywords single-cell RNA-seq data, CellTagging
#' @export
#' @examples
#' CellTagExtraction(bam.test.obj)
#' 
CellTagExtraction <- function(celltag.obj, celltag.version, technique = "10x") {
  celltag.obj@curr.version <- celltag.version
  if (file_test("-f", celltag.obj@fastq.bam.dir)) {
    fastq.bam.input <- celltag.obj@fastq.bam.dir
  } else {
    fastq.bam.input <- list.files(celltag.obj@fastq.bam.dir, full.names = T)
  }
  file.extension.unique <- unique(file_ext(fastq.bam.input))
  
  if (length(celltag.obj@celltag.version) > 0) {
    if (celltag.obj@curr.version %in% celltag.obj@celltag.version) {
      print("This CellTag has already been processed!")
    } else {
      celltag.obj@celltag.version <- c(celltag.obj@celltag.version, celltag.obj@curr.version)
    }
  } else {
    celltag.obj@celltag.version <- celltag.obj@curr.version
  }
  
  p.calling <- CellTagPatternCalling(celltag.version)
  
  if (endsWith(file.extension.unique, "fastq") || endsWith(file.extension.unique, "fq")) {
    if (length(fastq.bam.input) > 1) {
      stop("Please process the whitelist files one at a time!")
    }
    rslt <- fastq.process(fastq.file = fastq.bam.input, pattern = p.calling[1], p.calling[2], p.calling[3])
    celltag.obj@fastq.full.celltag[[celltag.version]] <- rslt[[1]]
    celltag.obj@fastq.only.celltag[[celltag.version]] <- rslt[[2]]
  }
  if (endsWith(file.extension.unique, "bam")) {
    rslt <- NULL
    for (i in 1:length(fastq.bam.input)) {
      curr.rslt <- bam.process(bam.file = fastq.bam.input[i], pattern = p.calling[1], p.calling[2], p.calling[3], technique)
      if (length(fastq.bam.input) > 1) curr.rslt$Cell.BC <- paste0("Sample-", i, "_", curr.rslt$Cell.BC)
      if (is.null(rslt)) {
        rslt <- curr.rslt
      } else {
        rslt <- rbind(rslt, curr.rslt, fill = TRUE)
      }
    }
    celltag.obj@bam.parse.rslt[[celltag.version]] <- rslt
  }
  
  return(celltag.obj)
}

