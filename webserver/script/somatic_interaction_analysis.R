library(maftools)

tcga_avail <- tcgaAvailable()

# Recupera gli argomenti dal terminale
args <- commandArgs(trailingOnly = TRUE)

tumore <- args[1]
output_dir <- args[2]

cx <- tcgaLoad(study = tumore, source = "Firehose")


#somaticInteractions
output_file <- paste0( output_dir, "/", tumore, "_somaticInteractions.png")
png(output_file, width = 1200, height = 800, res = 150)
results <- somaticInteractions(maf = cx, top = 10, pvalue = c(0.05, 0.01))
dev.off()
output_csv <- paste0(output_dir, "/", tumore, "_results.csv")
write.csv(results, output_csv, row.names = FALSE)

