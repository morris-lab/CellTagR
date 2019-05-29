#' @export
setClass("CellTag",
         slots = list(obj.name = "character",
                      fastq.bam.dir = "character",
                      curr.version = "character",
                      celltag.version = "character",
                      fastq.full.celltag = "ANY",
                      fastq.only.celltag = "ANY",
                      celltag.freq.stats = "ANY",
                      whitelist = "ANY",
                      bam.parse.rslt = "ANY",
                      celltag.stats = "ANY",
                      pre.starcode = "ANY",
                      raw.count = "dgCMatrix",
                      collapsed.count = "dgCMatrix",
                      whitelisted.count = "dgCMatrix",
                      metric.filtered.count = "dgCMatrix",
                      binary.mtx = "dgCMatrix",
                      jaccard.mtx = "dgCMatrix",
                      clone.composition = "ANY",
                      clone.size.info = "ANY",
                      celltag.aggr.final = "data.frame",
                      network.link.list = "ANY",
                      nodes = "ANY",
                      network = "ANY"))
#' @export
setMethod("show",
          "CellTag",
          function(object) {
            cat("Object name: ", object@obj.name, "\n")
            cat("Raw CellTag Counts = ", (ncol(object@raw.count)), "\n")
            cat("Raw Number of Cells with CellTag = ", nrow(object@raw.count), "\n")
            cat("Collapsed CellTag Counts = ", ncol(object@collapsed.count), "\n")
            cat("Whitelisted CellTag Counts = ", (ncol(object@whitelisted.count)), "\n")
            cat("Whitelisted Number of Cells with CellTag = ", nrow(object@whitelisted.count), "\n")
          })

