library(data.table)
library(parallel)

meta.data <- read.table("/scratch/smlab/CellTag_paper_analysis/permutation_test/final.drop.seq.10x.meta.data.txt", sep = "\t", 
                        stringsAsFactors = F, header = T, row.names = 1)

meta.data.orig <- read.table("/scratch/smlab/CellTag_paper_analysis/permutation_test/qp.comb.clone.meta.data.drop.seq.10X.txt", sep = "\t", 
                             stringsAsFactors = F, header = T, row.names = 1)

clones <- meta.data.orig[,c("hf1.v1", "hf1.v2.1", "hf1.v2.2", "hf2.v1", "hf2.v2.1", "hf2.v2.2")]

hf1 <- clones[, c("hf1.v1", "hf1.v2.1", "hf1.v2.2")]
hf2 <- clones[, c("hf2.v1", "hf2.v2.1", "hf2.v2.2")]

colnames(hf1) <- c("v1.1", "v2.1", "v2.2")
colnames(hf2) <- c("v1.1", "v2.1", "v2.2")

hf2.not.na.v1.1 <- which(!is.na(hf2$v1.1))
hf2.not.na.v2.1 <- which(!is.na(hf2$v2.1))
hf2.not.na.v2.2 <- which(!is.na(hf2$v2.2))
index.v1 <- intersect(which(is.na(hf1$v1.1)), hf2.not.na.v1.1)
hf1$v1.1[index.v1] <- 2000 + hf2$v1.1[index.v1]
hf1$v2.1[hf2.not.na.v2.1] <- 2000 + hf2$v2.1[hf2.not.na.v2.1]
hf1$v2.2[hf2.not.na.v2.2] <- 2000 + hf2$v2.2[hf2.not.na.v2.2]

all.clones <- hf1
# v2.1.unique.clones <- unique(all.clones$v2.1)
# clone.2.1 <- all.clones[,2]
# clone.2.1.count <- as.data.frame(table(clone.2.1))
# over.10.freq <- clone.2.1.count[which(clone.2.1.count$Freq > 10), ]

v1.1.unique.clones <- unique(all.clones$v1.1)
clone.1.1 <- all.clones[,1]
clone.1.1.count <- as.data.frame(table(clone.1.1))
over.10.freq.1.1 <- clone.1.1.count[which(clone.1.1.count$Freq > 10), ]


# high.number.clone.2.1 <- hf1[which(hf1$v2.1 %in% over.10.freq$clone.2.1), ]
high.number.clone.1.1 <- hf1[which(hf1$v1.1 %in% over.10.freq.1.1$clone.1.1), ]

####
# Fast sampling only
sampling <- function(clone.id, clone.info, over.threshold.df, subset.df=NULL) {
  curr.count <- over.threshold.df[as.character(clone.id), "Freq"]
  curr.cell.barcode <- rownames(clone.info)[which(clone.info$v1.1 == clone.id)]
  replicate.num <- ceiling(nrow(clone.info)/curr.count)
  barcode.names <- rownames(clone.info)[which(startsWith(rownames(clone.info), "_10X"))]
  if (!is.null(subset.df)){
    bc.nams <- rownames(subset.df)
    perm.subset <- replicate(replicate.num, sample(bc.nams, curr.count))
  }
  perm <- replicate(replicate.num, sample(barcode.names, curr.count))
  clone.perm <- replicate(replicate.num, sample(curr.cell.barcode, curr.count))
  return(list(perm, clone.perm, perm.subset))
}

# clones.id <- as.numeric(as.character(over.10.freq$clone.2.1))
# over.10.freq$clone.2.1 <- as.integer(as.character(over.10.freq$clone.2.1))
# rownames(over.10.freq) <- over.10.freq$clone.2.1

clones.id <- as.numeric(as.character(over.10.freq.1.1$clone.1.1))
over.10.freq.1.1$clone.1.1 <- as.integer(as.character(over.10.freq.1.1$clone.1.1))
rownames(over.10.freq.1.1) <- over.10.freq.1.1$clone.1.1
num.to.rep <- as.data.frame(seq(1, 50))

hf1.w.tp <- cbind(hf1, timepoints = unlist(lapply(strsplit(rownames(hf1), "-"), function(x) x[length(x)])))
hf1.w.tp.10x.only <- hf1.w.tp[which(startsWith(rownames(hf1.w.tp), "_10X")), ]
hf1.w.tp.10x.only$timepoints <- as.integer(as.character(hf1.w.tp.10x.only$timepoints))
hf1.aft.tp.3 <- hf1.w.tp.10x.only[which(hf1.w.tp.10x.only$timepoints >= 3), ]
rslt <- apply(num.to.rep,1,
              function(x) {
                sampling.ls <- mclapply(over.10.freq.1.1$clone.1.1, sampling, 
                                        clone.info = hf1, over.threshold.df = over.10.freq.1.1, subset.df = hf1.w.tp.10x.only,
                                        mc.cores = 24)
                return(sampling.ls)
              })


save(rslt, file = "/scratch/smlab/CellTag_paper_analysis/permutation_test/sampling_output_v1.RData")

# load("/scratch/smlab/CellTag_paper_analysis/permutation_test/sampling_output.RData")
# Calculate the percentages
percentage.perm.calc <- function(col.num, meta.data.original, x) {
  curr.samp.perm <- x[[1]][,col.num]
  curr.samp.perm.clone <- x[[2]][,col.num]
  curr.samp.perm.subset <- x[[3]][,col.num]

  perm.cluster.0.8 <- meta.data.original[curr.samp.perm, "res.0.8"]
  clone.perm.cluster.0.8 <- meta.data.original[curr.samp.perm.clone, "res.0.8"]
  subset.perm.cluster.0.8 <- meta.data.original[curr.samp.perm.subset, "res.0.8"]

  perm.percent <- length(which(perm.cluster.0.8 == 5)) * 100/length(perm.cluster.0.8)
  clone.perm.percent <- length(which(clone.perm.cluster.0.8 == 5)) * 100/length(clone.perm.cluster.0.8)
  subset.perm.percent <- length(which(subset.perm.cluster.0.8 == 5)) * 100/length(subset.perm.cluster.0.8)

  return(c(perm.percent, clone.perm.percent, subset.perm.percent))
}

percentage.ls <- lapply(rslt,
                                function(x) {
                                  rep <-
                                    lapply(x,
                                      function(x) {
                                         curr.cell.bar <- x[[2]][1,1]
                                         clone.id <- hf1.w.tp.10x.only[curr.cell.bar, ]$v1.1
                                         curr.cell.barcode <- rownames(hf1.w.tp.10x.only)[which(hf1.w.tp.10x.only$v1.1 == clone.id)]
                                         clone.cluster.0.8 <- meta.data.orig[curr.cell.barcode, "res.0.8"]
                                         percent.null <- length(which(clone.cluster.0.8 == 5)) * 100/length(clone.cluster.0.8)
                                         perc.calc.rslt.ls <- mclapply(seq(1, ncol(x[[1]])), percentage.perm.calc,
                                           meta.data.original = meta.data.orig, x = x,
                                           mc.cores = 24)
                                         percent.perm <- unlist(lapply(perc.calc.rslt.ls, function(x) x[1]))
                                         percent.clone.perm <- unlist(lapply(perc.calc.rslt.ls, function(x) x[2]))
                                         percent.subset.perm <- unlist(lapply(perc.calc.rslt.ls, function(x) x[3]))
                                         return(list(clone.id, percent.null, percent.perm, percent.clone.perm, percent.subset.perm))
                                      }
                                    )
                                  return(rep)
                                }
                              )

save(percentage.ls, file = "/scratch/smlab/CellTag_paper_analysis/permutation_test/percentage_over_10_v1.RData")

perm.test.super.ls <- lapply(percentage.ls,
                             function(x) {
                               p.value.ls <- lapply(x, 
                                                    function(y) {
                                                      clone.id <- y[[1]]
                                                      null.percent <- y[[2]]
                                                      real.distribution <- y[[5]]
                                                      curr.p <- sum(real.distribution > null.percent)/length(real.distribution)
                                                      return(data.frame(clone.num = clone.id, p.val = curr.p))
                                                    })
                               p.value.df <- rbindlist(p.value.ls)
                               return(p.value.df)
                             })

perm.df <- data.frame()
for (i in 1:length(perm.test.super.ls)) {
  curr.df <- perm.test.super.ls[[i]]
  if (ncol(perm.df) == 0) {
    perm.df <- curr.df
  } else {
    perm.df <- cbind(perm.df, curr.df[,2])
  }
}

clone.vec <- perm.df$clone.num
perm.df <- as.data.frame(perm.df[,-1])
rownames(perm.df) <- clone.vec
perm.df <- cbind(perm.df, avg = rowMeans(perm.df))

p.val.df <- data.frame(clone.id = rownames(perm.df), avg.p = perm.df[rownames(perm.df),]$avg)

save(p.val.df, file = "/scratch/smlab/CellTag_paper_analysis/permutation_test/p_value_v1.RData")