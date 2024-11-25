library(maftools)
library(data.table)
library(tidyverse)
library(DESeq2)
library(RColorBrewer)
library(gplots)
library(EnhancedVolcano)


tcga_avail = tcgaAvailable()

args <- commandArgs(trailingOnly = TRUE)

tumor <- args[1]
gene <- toupper(args[2])
#gene <- args[2]
outdir <- args[3]
input_file <- args[4]


# tumor="PAAD"
# gene="TP53"
# outdir="/Users/gerardo/Desktop/"

cx = tcgaLoad(study = tumor, source = "Firehose")

#define what clinical feature you want to split on
tipo <- gene

first.cond <- subsetMaf(cx, genes = tipo[1])


P_mut <- first.cond@data$Tumor_Sample_Barcode
P_mut <- Reduce(intersect,list(cx@data$Tumor_Sample_Barcode, P_mut))


P <- setdiff(cx@data$Tumor_Sample_Barcode, P_mut)

#apro il dataframe che contiene l'espressione dei geni
#"/Users/chiaranotturnogranieri/Downloads/notturno/TCGA/split_count_tumor_ensg/",tumor,".tsv"
df <- fread(paste(input_file, sep = ""), sep = "\t", index = "Genes")
cts <- as.data.frame(df)
cts <- cts %>% column_to_rownames(., var = 'Genes')

#selezioni dal dataframe solo i campioni per cui ho l'annotazione della mutazione della proteina
df_col <- Reduce(intersect, list(colnames(df), c(P_mut, P)))
cts <- cts[, df_col]

#seleziono dai pazienti wt tutti quelli per cui ho i dati di trascrittomica
P <- Reduce(intersect,list(colnames(df), P))

#seleziono dai pazienti mut tutti quelli per cui ho i dati di trascrittomica
P_mut <- Reduce(intersect,list(colnames(df), P_mut))



#exit with error if sample sizes are too small
if(length(P_mut)<5){ #non ci sono abbastanza campioni mutati
    quit(status=1) 
} else if (length(P)<5){
    quit(status=2) #non ci sono campioni abbastanza campioni non mutati
}


sample_type <- c(rep(c("MUT"),times=c(length(P_mut))),rep(c("WT"),times=c(length(P))))
sample_id <- c(colnames(cts))

coldata <- data.frame(Sample_ID = sample_id, Type = sample_type)
coldata <- coldata %>% column_to_rownames(., var = 'Sample_ID')

dds <- DESeqDataSetFromMatrix(countData = cts, colData= coldata, design = ~Type)

#inserire almeno il numero della metà dei campioni
keep <- rowSums(counts(dds) >= 10) >= length(sample_id)/2
dds <- dds[keep,]
dds2 <- DESeq(dds)


res <- results(dds2, contrast = c("Type", 'MUT', 'WT'),alpha=0.05)

write.table(res, file = paste(outdir,'/res_',tumor,"_",gene,'.txt',sep=""), sep="\t", row.names = TRUE, col.names = TRUE)

vsd <- varianceStabilizingTransformation(dds2)

sampleDist <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDist)
rownames (sampleDistMatrix) <- paste(vsd$Type, vsd$patient, sep="-")
colnames(sampleDistMatrix) <- paste(vsd$Type, vsd$patient, sep="-")

colours <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

output_file <- paste0(outdir, "/", "heatmap_", tumor, "_", gene, ".png")
jpeg(output_file, width = 1200, height = 1000, quality = 100)
heatmap.2(sampleDistMatrix, trace="none", col=colours, cexRow = 0.3, cexCol = 0.3)
dev.off()

output_file <- paste0(outdir, "/", "PCA_", tumor, "_", gene, ".png")
jpeg(output_file, width = 1200, height = 1000, quality = 100)
plotPCA(vsd, intgroup = c("Type"))
dev.off()


topVarGeni <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20) 

write.table(assay(vsd), file = paste(outdir,'/top_50_vsd_',tumor,'_',gene,'.tsv',sep=""), sep="\t", row.names = TRUE, col.names = TRUE)


output_file <- paste0(outdir, "/", "Top50genes_", tumor, "_", gene, ".png")
jpeg(output_file, width = 1200, height = 1000, quality = 100)
heatmap.2(assay(vsd)[topVarGeni, ], cexRow = 0.35, cexCol = 0.36, scale = "row", trace = "none", dendrogram = "column", col = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255))
dev.off()


output_file <- paste0(outdir, "/", "EnhancedVolcano_", tumor, "_", gene, ".png")
jpeg(output_file, width = 1200, height = 1000, quality = 100)
EnhancedVolcano(res,
                  lab = rownames(res),
                  x = 'log2FoldChange',
                  y = 'padj')
dev.off()
