library(maftools)

tcga_avail <- tcgaAvailable()

# Recupera gli argomenti dal terminale
args <- commandArgs(trailingOnly = TRUE)

tumore <- args[1]
output_dir <- args[2]

cx <- tcgaLoad(study = tumore, source = "Firehose")

#Plotsummarymaf
output_file <- paste0(output_dir, "/", tumore, "_maf_summary.png")
png(output_file, width = 1200, height = 800, res = 150)
plotmafSummary(maf = cx, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

#plot titv summary
output_file <- paste0(output_dir, "/", tumore, "_Titv.png")
png(output_file, width = 1200, height = 800, res = 150)
cx.titv = titv(maf = cx, plot = FALSE, useSyn = TRUE)
plotTiTv(res = cx.titv)
dev.off()


# #oncoplot for top ten mutated genes.
# output_file <- paste0(output_dir, "/", tumore, "_oncoplot.png")
# png(output_file, width = 1200, height = 800, res = 150)
# oncoplot(maf = cx, top = 15)
# dev.off()


# #somaticInteractions
# output_file <- paste0( output_dir, "/", tumore, "_somaticInteractions.png")
# png(output_file, width = 1200, height = 800, res = 150)
# results <- somaticInteractions(maf = cx, top = 10, pvalue = c(0.05, 0.01))
# dev.off()
# output_csv <- paste0(output_dir, "/", tumore, "_results.csv")
# write.csv(results, output_csv)


# #PlotVAF
# output_file <- paste0(output_dir, "/", tumore, "_TumorVAF.png")
# png(output_file, width = 1200, height = 800, res = 150)
# plotVaf(maf = cx, vafCol = 'i_TumorVAF_WU', top = 10)
# dev.off()
