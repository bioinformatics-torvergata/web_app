library(maftools)



# Recupera gli argomenti dal terminale
args <- commandArgs(trailingOnly = TRUE)
tumore <- args[1]
#print(tumore)
gene_input <- toupper(args[2])
#gene_input <- args[2]
output_dir <- args[3]



tcga_avail <- tcgaAvailable()
#cx <- tcgaLoad(study = tumore)

cx <- tcgaLoad(study = tumore, source = "Firehose")

#lista geni per quel tumore:
gene_list <- unique(cx@data$Hugo_Symbol)

if (gene_input %in% gene_list) {
  #LollipopPlot
  output_file <- paste0(output_dir, "/", tumore, "_", gene_input, "_lollipopPlot.png")
  #print(output_file)
 
  png(output_file, width = 1200, height = 800, res = 150)
  lollipopPlot(
    maf = cx,
    gene = gene_input,
    AACol = 'Protein_Change',
    showMutationRate = TRUE,
    showDomainLabel=FALSE,
  )
  dev.off()
  
  gene_data <- cx@data[cx@data$Hugo_Symbol == gene_input, ]
  output_txt<- paste0(output_dir,"/result.txt" )
  write.table(gene_data, file = output_txt, sep = "\t", row.names = FALSE, quote = FALSE)
  
} else{
  return(0)
}