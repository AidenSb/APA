library(Sierra)
library(presto)
library(ggplot2)
library(msigdbr)
library(foreach)
library(doParallel)
library(fgsea)
library(Seurat)
library(tidyverse)


setwd('/data/APAproject/post_qual/data/control_vs_pathology')


## read in the Exc_CTandAD peaks so
Exc_peaks_AD_CT <- readRDS('rdata/Exc_peaks_CTandAD_merged.rds')

#lets get the combinations of celltype vs assigned subclusters
table(Exc_peaks_AD_CT$subclusterAssignment)

ADsc <- c("SFG:Exc.s0","SFG:Exc.s2" ,"SFG:Exc.s1","SFG:Exc.s10","SFG:Exc.s8","SFG:Exc.s9" ,"SFG:Exc.s4" ,
        "SFG:Exc.s7","SFG:Exc.s6","SFG:Exc.s3","SFG:Exc.s5" )
CTsc <- c("L5/6","L2/3","L5/6-CC","Neu-NRGN-I","Neu-mat","L4","Neu-NRGN-II")


results_root <- 'results/Exc_AD_vs_DimCT/'

for (ct_sc in CTsc){
  for (ad_sc in ADsc){
    cell1 = ct_sc
    cell2 = ad_sc
    sanity <- paste0(cell1, '--vs--' , cell2)
    print(sanity)
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
    plotoutname <- paste0(results_root,cell1,"_vs_",cell2,'_3UTR_shift_plot.pdf')
    tablename <- paste0(results_root,cell1,"_vs_",cell2,'_3UTR_shift_table.tsv')
    APAoutename <- paste0(results_root,cell1,"_vs_",cell2,'_APAusage_table.tsv')
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
}
######
##############################################################################################################
#################################### Dim neurons vs Kapmann Neurons ###############
# read in the Dimitery neurosns