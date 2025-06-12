#This is an example of the R code you run in order to prepare the matrices for R in order to finding the celltype and gene expression
#Use R v.4.4.0 or later 

library(Matrix)

refmat <- Matrix::readMM("matrix_og_barcode/cleaned_Lib7_refmat.mtx")
altmat <- Matrix::readMM("matrix_og_barcode/cleaned_Lib7_altmat.mtx")

snps <- read.table("sample7_snps.txt", stringsAsFactors=F)
cellBarcodes <- read.table("og_barcodes/sample7_barcodes.txt", stringsAsFactors = F)$V1
snpnames <- paste0(snps$V1, ':', snps$V2, '-', snps$V2)
rownames(refmat) <- rownames(altmat) <- cellBarcodes
colnames(refmat) <- colnames(altmat) <- snpnames

row.names(refmat) <- trimws(paste0("Lib7_",row.names(refmat)))
row.names(altmat) <- trimws(paste0("Lib7_",row.names(altmat)))

sref <- summary(refmat)
salt <- summary(altmat)

altdf  <- data.frame(
  cell = rownames(altmat)[salt$i],       
  i     = salt$i,                    
  j     = colnames(altmat)[salt$j],      
  alt   = salt$x                      
)


refdf  <- data.frame(
  cell = rownames(refmat)[sref$i],       
  i     = sref$i,                    
  j     = colnames(refmat)[sref$j],      
  ref  = sref$x                      
)

marc <- merge(refdf, altdf, by=c("i", "j"))
marcdf <- marc[order(marc$j),]


write.table(marcdf,"complete_matrice_for_R/Lib7_counts.txt", quote = F)

