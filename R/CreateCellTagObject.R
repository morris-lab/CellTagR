#' Create a New CellTag Object 
#'
#' This function creates a CellTag object that contains the basic information required for the object
#' @param object.name The name of the object
#' @param fastq.bam.input The input fastq/bam data directory
#' @param celltag.version Which version of CellTags are you working with?
#' @return A CellTag Object with open attributes that can be filled as analysis moving along
#' @keywords single-cell RNA-seq data, CellTagging
#' @export
#' @examples
#' CellTagObejct("hf1.d15.test", "hf1.d15.bam", "v1")
#'
CellTagObject <- function(object.name, fastq.bam.directory) {
  ct <- new("CellTag", obj.name = object.name, fastq.bam.dir = fastq.bam.directory)
  return(ct)
}
