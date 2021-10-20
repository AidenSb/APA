library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
library(Sierra)
library(presto)
library(ggplot2)
library(msigdbr)
library(foreach)
library(doParallel)
library(fgsea)
library(Seurat)
library(tidyverse)
setwd('/data/APAproject/post_qual/data/Dimitry_Velmeshev_ASD/')

neurons_peaks_so <- readRDS('rdata/neurons_peaks_so.rds')
View(PlotUTRLengthShift)
tst <- trace(PlotUTRLengthShift, edit=TRUE)          ## edit the title of the plots


reference.file <- '/home/aiden/data/refgenome/refdata-gex-GRCh38-2020-A/genes//genes.gtf'
gtf_gr <- rtracklayer::import(reference.file)
gtf_TxDb <- GenomicFeatures::makeTxDbFromGFF(reference.file, format="gtf")

## test unit
IE_neruons_UTR_shift <- DetectUTRLengthShift(peaks.object = neurons_peaks_so, 
                                             gtf_gr = gtf_gr,
                                             gtf_TxDb = gtf_TxDb,
                                             population.1 = "Excitatory_Neurons", 
                                             population.2 = "Inhibitory_Neurons")
results_root <- '/data/APAproject/post_qual/data/Dimitry_Velmeshev_ASD/results/'
plotname <- paste0(results_root, 'tst.pdf')
pdf(file=plotname, width = 8.27, height = 11.69, paper='A4r', onefile=T)
PlotUTRLengthShift(IE_neruons_UTR_shift)
dev.off()
write.table(IE_neruons_UTR_shift, file=paste0(results_root,'Inhibitory_Excitatory_3UTR_shift.tsv'), sep='\t')

IE_neruons_APA_usage <- DetectAEU(peaks.object = neurons_peaks_so, 
                                 gtf_gr = gtf_gr,
                                 gtf_TxDb = gtf_TxDb,
                                 do.MAPlot = T,
                                 population.1 = "Neu-NRGN-II", 
                                 population.2 = "L5/6-CC")
## now lets do this for all the combinations of the neuron types


reference.file <- '/home/aiden/data/refgenome/refdata-gex-GRCh38-2020-A/genes//genes.gtf'
gtf_gr <- rtracklayer::import(reference.file)
gtf_TxDb <- GenomicFeatures::makeTxDbFromGFF(reference.file, format="gtf")



for(i in 1:ncol(combinations)) {       # for-loop over columns
  cell1 = combinations[1, i]
  cell2 = combinations[2, i]
  print(cell1)
  print(cell2)
  p <- NULL
  res.table <- NULL
  res.table <- DetectUTRLengthShift(peaks.object = neurons_peaks_so, 
                                    gtf_gr = gtf_gr,
                                    gtf_TxDb = gtf_TxDb,
                                    population.1 = cell1, 
                                    population.2 = cell2)
  
  print("Global 3UTR shift analysis done and now doing APAusage analysis")
  apa.res <- NULL
  tryCatch({apa.res <- DetectAEU(peaks.object = neurons_peaks_so, 
                                 gtf_gr = gtf_gr,
                                 gtf_TxDb = gtf_TxDb,
                                 do.MAPlot = F,
                                 population.1 = cell1,
                                 population.2 = cell2)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  cell1 = gsub('/','-',cell1)
  cell2 = gsub('/','-',cell2)
  plotoutname <- paste0(results_root,'neuron_subtypes/',cell1,"_vs_",cell2,'_3UTR_shift_plot.pdf')
  tablename <- paste0(results_root,'neuron_subtypes/',cell1,"_vs_",cell2,'_3UTR_shift_table.tsv')
  APAoutename <- paste0(results_root,'neuron_subtypes/',cell1,"_vs_",cell2,'_APAusage_table.tsv')
  print(plotoutname)
  tryCatch({
  write.table(res.table, file=tablename, sep='\t')
  p <- PlotUTRLengthShift(res.table)
  print('second check point')
  write.table(apa.res, file=APAoutename, sep='\t')}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  pdf(file=plotoutname, width = 8.27, height = 11.69, paper='A4r', onefile=T)
  print(p)
  dev.off()
  
}

### repeat the operation above for PFC only neurons
PFC_neurons_peaks_so <- subset(neurons_peaks_so, subset=region=='PFC')
saveRDS(PFC_neurons_peaks_so,'rdata/PFC_neurons_peaks_so.rds')



for(i in 1:ncol(combinations)) {       # for-loop over columns
  cell1 = combinations[1, i]
  cell2 = combinations[2, i]
  print(cell1)
  print(cell2)
  p <- NULL
  res.table <- NULL
  res.table <- DetectUTRLengthShift(peaks.object = PFC_neurons_peaks_so, 
                                    gtf_gr = gtf_gr,
                                    gtf_TxDb = gtf_TxDb,
                                    population.1 = cell1, 
                                    population.2 = cell2)
  
  print("Global 3UTR shift analysis done and now doing APAusage analysis")
  apa.res <- NULL
  tryCatch({apa.res <- DetectAEU(peaks.object = PFC_neurons_peaks_so, 
                                 gtf_gr = gtf_gr,
                                 gtf_TxDb = gtf_TxDb,
                                 do.MAPlot = F,
                                 population.1 = cell1,
                                 population.2 = cell2)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  cell1 = gsub('/','-',cell1)
  cell2 = gsub('/','-',cell2)
  plotoutname <- paste0(results_root,'neuron_subtypes/PFC_only/',cell1,"_vs_",cell2,'_3UTR_shift_plot.pdf')
  tablename <- paste0(results_root,'neuron_subtypes/PFC_only/',cell1,"_vs_",cell2,'_3UTR_shift_table.tsv')
  APAoutename <- paste0(results_root,'neuron_subtypes/PFC_only/',cell1,"_vs_",cell2,'_APAusage_table.tsv')
  print(plotoutname)
  tryCatch({
    write.table(res.table, file=tablename, sep='\t')
    p <- PlotUTRLengthShift(res.table)
    print('second check point')
    write.table(apa.res, file=APAoutename, sep='\t')}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  pdf(file=plotoutname, width = 8.27, height = 11.69, paper='A4r', onefile=T)
  print(p)
  dev.off()
  
}




### repeat the operation above for All the celltypes in PFC dataset
PFC_all_peaks_so <- readRDS('rdata/Control_celltypes_peaks.rds')
PFC_all_peaks_so <-  subset(PFC_all_peaks_so, subset=region=='PFC')
PFC_all_peaks_so <- SetIdent(PFC_all_peaks_so, value = 'cluster')
combinations <- combn(unique(PFC_all_peaks_so$cluster),2)

for(i in 1:ncol(combinations)) {       # for-loop over columns
  cell1 = combinations[1, i]
  cell2 = combinations[2, i]
  print(cell1)
  print(cell2)
  p <- NULL
  res.table <- NULL
  res.table <- DetectUTRLengthShift(peaks.object = PFC_all_peaks_so, 
                                    gtf_gr = gtf_gr,
                                    gtf_TxDb = gtf_TxDb,
                                    population.1 = cell1, 
                                    population.2 = cell2)
  
  print("Global 3UTR shift analysis done and now doing APAusage analysis")
  apa.res <- NULL
  tryCatch({apa.res <- DetectAEU(peaks.object = PFC_all_peaks_so, 
                                 gtf_gr = gtf_gr,
                                 gtf_TxDb = gtf_TxDb,
                                 do.MAPlot = F,
                                 population.1 = cell1,
                                 population.2 = cell2)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  cell1 = gsub('/','-',cell1)
  cell2 = gsub('/','-',cell2)
  plotoutname <- paste0(results_root,'all_vs_all_celltypes/',cell1,"_vs_",cell2,'_3UTR_shift_plot.pdf')
  tablename <- paste0(results_root,'all_vs_all_celltypes/',cell1,"_vs_",cell2,'_3UTR_shift_table.tsv')
  APAoutename <- paste0(results_root,'all_vs_all_celltypes/',cell1,"_vs_",cell2,'_APAusage_table.tsv')
  print(plotoutname)
  tryCatch({
    write.table(res.table, file=tablename, sep='\t')
    p <- PlotUTRLengthShift(res.table)
    print('second check point')
    write.table(apa.res, file=APAoutename, sep='\t')}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  pdf(file=plotoutname, width = 8.27, height = 11.69, paper='A4r', onefile=T)
  print(p)
  dev.off()
  
}


