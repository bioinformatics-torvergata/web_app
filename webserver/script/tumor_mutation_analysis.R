library(maftools)

tcga_avail <- tcgaAvailable()

#head(tcga_avail, 35)


# Recupera gli argomenti dal terminale
args <- commandArgs(trailingOnly = TRUE)


# Prendi il primo argomento (il nome del tumore)
tumore <- args[1]
output_dir <- args[2]

cx <- tcgaLoad(study = tumore, source = "Firehose")

#save plot in output_dir
output_file <- paste0(output_dir,"/", tumore, "_maf_summary.png")
print(output_file)
png(output_file, width = 1200, height = 800, res = 150)

# Genera il plot
plotmafSummary(maf = cx, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

