library(gaston)
library(RAINBOWR)

# create save folder
saveFolder <- "data/genome/"
if (!dir.exists(saveFolder)) {
  dir.create(saveFolder)
}

# read the vcf file
genomeFile <- "raw_data/genome/Gm198_HCDB_190207.fil.snp.remHet.MS0.95_bi_MQ20_DP3-1000.MAF0.025.imputed.v2.chrnum.vcf.gz"
genomeRaw <- read.vcf(genomeFile)
genomeRaw@snps$id <- paste0("Chr", genomeRaw@snps$chr, "_", genomeRaw@snps$pos)

GRM <- function(genomeRaw, threshold, saveFolder) {
  # filtering the SNPs data
  genomeFil <- select.snps(genomeRaw, condition = maf >= 0.05)
  genomeFil <- LD.thin(genomeFil, threshold = threshold)
  
  genomeMat <- as.matrix(genomeFil) - 1
  # sum(apply(genomeMat == 0, 1, sum))
  # genomeMat[1:5, 1:5]
  filename <- paste0("amat", ncol(genomeMat), "SNP.csv")
  amat <- calcGRM(genoMat = genomeMat, methodGRM = "A.mat")
  write.csv(x = amat, file = paste0(saveFolder, filename))
}
# create the normal GRM 
GRM(genomeRaw = genomeRaw, threshold = 0.8, saveFolder = saveFolder)

# # create the marker csv (LD = 0.1 - 0.5)
# for (i in seq(0.1, 0.5, 0.1)) {
#   genomeFil <- select.snps(genomeRaw, condition = maf >= 0.05)
#   genomeFil <- LD.thin(genomeFil, threshold = i)
#   
#   genomeMat <- as.matrix(genomeFil) - 1
#   genomeRaw@snps$pos
#   dim(genomeMat)
#   filename <- paste0("genoMarker_LD", i, "_SNP", ncol(genomeMat), ".csv")
#   write.csv(x = genomeMat, file = paste0(saveFolder, filename))
# }
# 
