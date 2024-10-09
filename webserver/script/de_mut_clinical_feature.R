library(maftools)

tcga_avail <- tcgaAvailable()

# Recupera gli argomenti dal terminale
args <- commandArgs(trailingOnly = TRUE)

tumore <- args[1]
clinical_feature <- args[2]
output_dir <- args[3]

cx <- tcgaLoad(study = tumore, source = "Firehose")

#define what clinical feature you want to split on
features <- cx@clinical.data[[clinical_feature]]

features <- features[!is.na(features)]
features <- unique(features)

tipo <- c(features[1],features[2])

tum1 <- tipo[1]
tum2 <- tipo[2]


first.cond <- subsetMaf(cx,clinQuery = paste0(clinical_feature, " == tipo[1]"))
second.cond <- subsetMaf(cx,clinQuery = paste0(clinical_feature, " == tipo[2]"))

#Considering only genes which are mutated in at-least in 5 samples in one of the cohort to avoid bias due to genes mutated in single sample.
pt.vs.rt <- mafCompare(m1 = first.cond, m2 = second.cond, m1Name = tum1, m2Name = tum2, minMut = 5)




# Controlla se ci sono più di 20 geni nei risultati
if (nrow(pt.vs.rt$results) > 5) {
  # Se ci sono più di 20 geni, seleziona i top 20 basati sul p-value
  topGenes <- pt.vs.rt$results[order(pt.vs.rt$results$pval), ][1:10, ]
  # Mantieni l'oggetto mafCompare ma sostituisci i risultati con i top 10 geni
  pt.vs.rt.top10 <- pt.vs.rt
  pt.vs.rt.top10$results <- topGenes
  # Crea il forest plot solo per i top 10 geni
  output_file <- paste0(output_dir, "/", tumore, "_forestPlot.png")
  png(output_file, width = 1200, height = 800, res = 150)
  forestPlot(mafCompareRes = pt.vs.rt.top10, pVal = 0.1)
  dev.off()
  output_csv <- paste0(output_dir, "/", tumore, "_", clinical_feature, "_results.csv")
  write.csv(pt.vs.rt$results, output_csv, row.names = FALSE)
} else {
  # Se ci sono 10 o meno geni, crea il forest plot con tutti i geni
    output_file <- paste0(output_dir, "/", tumore, "_forestPlot.png")
    png(output_file, width = 1200, height = 800, res = 150)
    forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.1)
    dev.off()

    output_csv <- paste0(output_dir, "/", tumore, "_", clinical_feature, "_results.csv")
    write.csv(pt.vs.rt$results, output_csv, row.names = FALSE)   
}



output_file <- paste0(output_dir, "/", tumore, "_coBarplot.png")
png(output_file, width = 1200, height = 800, res = 150)
coBarplot(m1 = first.cond, m2 = second.cond, m1Name = tum1, m2Name = tum2)
dev.off()