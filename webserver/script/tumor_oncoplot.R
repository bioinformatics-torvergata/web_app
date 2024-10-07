library(maftools)

tcga_avail <- tcgaAvailable()

# Recupera gli argomenti dal terminale
args <- commandArgs(trailingOnly = TRUE)

tumore <- args[1]
output_dir <- args[2]
number <- args[3]

cx <- tcgaLoad(study = tumore, source = "Firehose")

#oncoplot for top ten mutated genes.
output_file <- paste0(output_dir, "/", tumore, "_oncoplot.png")
png(output_file, width = 1200, height = 1000, res = 150)
oncoplot(maf = cx, top = number)
dev.off()