#' CellTag Extraction Function
#'
#' This function extracts CellTags from the raw fastq/bam sequencing file. If it is a fastq file, provides counts of each CellTag and sorts them in desending order. If it is a bam file, returns the barcode, umi, celltag information.
#' @param celltag.obj A CellTag object initialized with path to the fastq/bam file
#' @return A CellTag object with attribute (bam.parse.rslt) filled
#' @keywords single-cell RNA-seq data, CellTagging
#' @export
#' @examples
#' CellTagExtraction(bam.test.obj)
#' 
CellTagExtraction <- function(celltag.obj, celltag.version) {
  celltag.obj@curr.version <- celltag.version
  fastq.bam.input <- celltag.obj@fastq.bam.dir
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
  
  if (endsWith(fastq.bam.input, "fastq")) {
    rslt <- fastq.process(fastq.file = fastq.bam.input, pattern = p.calling[1], p.calling[2], p.calling[3])
    celltag.obj@fastq.full.celltag[[celltag.version]] <- rslt[[1]]
    celltag.obj@fastq.only.celltag[[celltag.version]] <- rslt[[2]]
  }
  if (endsWith(fastq.bam.input, "bam")) {
    rslt <- bam.process(bam.file = fastq.bam.input, pattern = p.calling[1], p.calling[2], p.calling[3])
    celltag.obj@bam.parse.rslt[[celltag.version]] <- rslt
  }
  
  return(celltag.obj)
}

