# Library
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
library(org.Hs.eg.db)
library(GenomicRanges)
library(dplyr)

# Set working directory
setwd("/Users/karen/mount/chuk/YTHDC1_HyperTRIBE/Project_10325_B/JAX_0381/")

# Import HyperTRIBE data
df <- read.csv("/Users/karen/mount/chuk/YTHDC1_HyperTRIBE/Project_10325_B/JAX_0381/MOLM13_YTHDC1_snp_counts_snp_counts_significance_fpkm.csv")

# Prepare significant edit sites
df.sig <- subset(df, p.adj < 0.05)
df.sig$grange.start <- df.sig$start - 100
df.sig$grange.end <- df.sig$start + 100
hypertribe.cluster.df <- df.sig %>% select(c("seqnames", "grange.start", "grange.end", "strand", "gene.symbol"))

# Prepare non-significant edit sites (background)
df.not.sig <- subset(df, p.adj >= 0.05)
df.not.sig$grange.start <- df.not.sig$start - 100
df.not.sig$grange.end <- df.not.sig$start + 100
hypertribe.cluster.df.not.sig <- df.not.sig %>% select(c("seqnames", "grange.start", "grange.end", "strand", "gene.symbol"))

# Convert dataframe to GenomicRanges object
hypertribe.cluster.gr <- makeGRangesFromDataFrame(hypertribe.cluster.df,
                         keep.extra.columns=TRUE,
                         seqnames.field=c("seqnames"),
                         start.field="grange.start",
                         end.field=c("grange.end"),
                         strand.field="strand" )

hypertribe.cluster.gr.not.sig <- makeGRangesFromDataFrame(hypertribe.cluster.df.not.sig,
                                                          keep.extra.columns=TRUE,
                                                          seqnames.field=c("seqnames"),
                                                          start.field="grange.start",
                                                          end.field=c("grange.end"),
                                                          strand.field="strand" )

# Get genome
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene # annotation database of mouse genome
exons.by.transcript.id <- GenomicFeatures::exonsBy(txdb, by=("tx"), use.names = TRUE) # Group exons by transcript ID
genome <- BSgenome.Hsapiens.UCSC.hg19 # mouse genome

# Get sequence based on grange coordinates
prepare.homer.input <- function(cluster.gr, cluster.df) {
  
  hypertribe.cluster.sequence <- Biostrings::getSeq(genome, cluster.gr)
  hypertribe.cluster.sequence.df <- as.data.frame(hypertribe.cluster.sequence)
  rownames(hypertribe.cluster.sequence.df) <- paste(cluster.df$gene.symbol, 
                                                    cluster.df[,1],
                                                    cluster.df[,2],
                                                    cluster.df[,3],
                                                    rownames(hypertribe.cluster.sequence.df),
                                                    sep=".")
  colnames(hypertribe.cluster.sequence.df) <- "sequence"
  
  return(hypertribe.cluster.sequence.df)
  
}

hypertribe.cluster.sequence <- prepare.homer.input(hypertribe.cluster.gr, hypertribe.cluster.df)
hypertribe.cluster.sequence.not.sig <- prepare.homer.input(hypertribe.cluster.gr.not.sig, hypertribe.cluster.df.not.sig)


# write fasta file using dataframe that has column names "ensembl_gene_id" and "3utr".
writeFasta<-function(data, filename){
  
  data <- as.data.frame(data)
  
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", rownames(data)[rowNum], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"sequence"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

# Write all sequences
writeFasta(hypertribe.cluster.sequence,
           paste0("MOLM13_YTHDC1_HyperTRIBE_significant_hypertribe_clusters_100bp_window.fa"))
writeFasta(hypertribe.cluster.sequence.not.sig,
           paste0("MOLM13_YTHDC1_HyperTRIBE_NOT_significant_hypertribe_clusters_100bp_window.fa"))



