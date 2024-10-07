library(maftools)

tcga_avail <- tcgaAvailable()

# Recupera gli argomenti dal terminale
args <- commandArgs(trailingOnly = TRUE)

tumore <- args[1]
clinical_feature <- args[2]
output_dir <- args[3]

cx <- tcgaLoad(study = tumore, source = "Firehose")

#define what clinical feature you want to split on
#clinical_feature <- "history_of_diabetes"
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
#print(pt.vs.rt)

output_file <- paste0(output_dir, "/", tumore, "_forestPlot.png")
png(output_file, width = 1200, height = 800, res = 150)
forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.1)
dev.off()


output_file <- paste0(output_dir, "/", tumore, "_coBarplot.png")
png(output_file, width = 1200, height = 800, res = 150)
coBarplot(m1 = first.cond, m2 = second.cond, m1Name = tum1, m2Name = tum2)
dev.off()