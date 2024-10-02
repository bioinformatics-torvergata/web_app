library(maftools)



# Recupera gli argomenti dal terminale
args <- commandArgs(trailingOnly = TRUE)
tumore <- args[1]
gene_input <- args[2]
output_dir <- args[3]




tcga_avail <- tcgaAvailable()


cx <- tcgaLoad(study = tumore, source = "Firehose")

#Survival analysis based on grouping of "geneX" mutation status
output_file <- paste0(output_dir, "/", tumore, "_", gene_input, "_survival.png")
png(output_file, width = 1200, height = 800, res = 150)
mafSurvival(maf = cx, genes =gene_input, time = 'days_to_last_followup', Status = 'vital_status', isTCGA = TRUE)
dev.off()


#for more gene
# prog_geneset = survGroup(maf = cx, top = 20, geneSetSize = 3, time = "days_to_last_followup", Status = "vital_status", verbose = FALSE)
# print(prog_geneset)

# g.set<-strsplit(prog_geneset$Gene_combination[1], "_")[[1]]
# mafSurvGroup(maf = cx, geneSet = c(g.set), time = "days_to_last_followup", Status = "vital_status")
