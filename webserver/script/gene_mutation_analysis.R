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

#LollipopPlot
output_file <- paste0(output_dir, "/", tumore, "_", gene_input, "_lollipopPlot.png")
print(output_file)

png(output_file, width = 1200, height = 800, res = 150)
lollipopPlot(
  maf = cx,
  gene = gene_input,
  AACol = 'Protein_Change',
  showMutationRate = TRUE,
  showDomainLabel=FALSE,
)
dev.off()