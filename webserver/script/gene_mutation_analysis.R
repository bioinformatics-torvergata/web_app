library(maftools)



# Recupera gli argomenti dal terminale
args <- commandArgs(trailingOnly = TRUE)
tumore <- args[1]
print(tumore)
gene_input <- args[2]
output_dir <- args[3]


tcga_avail <- tcgaAvailable()
#cx <- tcgaLoad(study = tumore)

cx <- tcgaLoad(study = tumore, source = "Firehose")

#somaticInteractions
output_file <- paste0( output_dir, "/", tumore, "_", gene_input, "_somaticInteractions.png")
png(output_file, width = 1200, height = 800, res = 150)
results <- somaticInteractions(maf = cx, top = 10, pvalue = c(0.05, 0.01))
dev.off()
output_csv <- paste0(output_dir, "/", tumore, "_", gene_input, "_results.csv")
write.csv(results, output_csv)


#PlotVAF
output_file <- paste0(output_dir, "/", tumore, "_", gene_input, "_TumorVAF.png")
png(output_file, width = 1200, height = 800, res = 150)
plotVaf(maf = cx, vafCol = 'i_TumorVAF_WU', top = 10)
dev.off()


#LollipopPlot
output_file <- paste0(output_dir, "/", tumore, "_", gene_input, "_lollipopPlot.png")
print(output_file)
png(output_file, width = 1200, height = 800, res = 150)
lollipopPlot(
  maf = cx,
  gene = gene_input,
  AACol = 'Protein_Change',
  showMutationRate = TRUE
)
dev.off()