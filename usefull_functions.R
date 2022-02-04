### Useful functions!



## get the relative locoation and you can decide on the distal vs proximal 
get_relative_pos <- function(pos,sitNum){
    relative_loc <- (pos - 1) / (sitNum -1)
    return(relative_loc)
}

map_to_ref_UTR <- function(peak_id, utr3.ref, all_expressed_peaks_df){
    ## 2- make the query granges to intersect with the reference 3UTR
    strand = sub(".*:.*:.*-.*:(.*)", "\\1", peak_id)
    strand = plyr::mapvalues(x = strand, from = c("1", "-1"), to = c("+", "-"))
    peak.remainder = sub(".*:(.*:.*-.*):.*", "\\1", peak_id)

    peaks.expressed.granges = paste0(peak.remainder, ":", strand)
    expressed.peaks.gr <- GenomicRanges::GRanges(peaks.expressed.granges)
    granges_peaks_mapping_table <- data.frame(PeakID = peak_id,
                                                row.names = peaks.expressed.granges, 
                                                stringsAsFactors = FALSE)
    all_UTR_3_hits <- GenomicRanges::findOverlaps(expressed.peaks.gr , utr3.ref, type = "any")
    utr3.mappings <- as.data.frame(all_UTR_3_hits)  ## this returns the location of the matched ref 3UTR range
    ##
    query.hit.df <- as.data.frame(expressed.peaks.gr[utr3.mappings$queryHits, ])
    subject.hit.df <- as.data.frame(utr3.ref[utr3.mappings$subjectHits, ],
                                    row.names = as.character(1:nrow(utr3.mappings)))  ## this is basically the UTRs region that your peaks map :) and this is what we want
    query.hit.df %>% dplyr::mutate(granges_peak = paste0(seqnames,":",start,"-",end,":",strand)) -> query.hit.df
    peak.ids <- granges_peaks_mapping_table[as.character(query.hit.df$granges_peak), 'PeakID']


    peaks_mapped_to_ref_UTR <- data.frame(subject.hit.df,
                                          data.frame(peak_ID=peak.ids, peak_width=query.hit.df$width))
    
    ## get the peaks that expressed in AD
    idx = which(peaks_mapped_to_ref_UTR$peak_ID %in% all_expressed_peaks_df$peaks) 
    peaks_mapped_to_ref_UTR = peaks_mapped_to_ref_UTR[idx,]
    return(peaks_mapped_to_ref_UTR)
    
}

################
## this function takes in the distal peak and all the expressed UTR peaks for the gene and returns most 
## proximal non DE peak coordinates
get_proximal_peak_id <- function(distal_peak, mapped_utrs, df){
    distance = 1e15
    proximal_id <- 'Not-found'
    for (peak in unlist(unique(mapped_utrs['peak_ID']))){
        if (!peak %in% df$peak_ids){
            strand = sub(".*:.*:.*-.*:(.*)", "\\1", distal_peak)
            distal_p_start = sub(".*:.*:(.*)-.*:.*", "\\1", distal_peak)
            distal_p_end = sub(".*:.*:.*-(.*):.*", "\\1", distal_peak)
            proximal_p_start = sub(".*:.*:(.*)-.*:.*", "\\1", peak)
            proximal_p_end = sub(".*:.*:.*-(.*):.*", "\\1", peak)
            if (strand=='1'){
                current_dist <- as.numeric(distal_p_end) - as.numeric(proximal_p_end)
		if (current_dist <= distance & current_dist > 0){
                    distance = current_dist
                    proximal_id <- peak
                }

            } else if (strand == '-1') {
                    current_dist <- as.numeric(proximal_p_start) - as.numeric(distal_p_start) 
                    if (current_dist <= distance & current_dist > 0){
                        distance = current_dist
                        proximal_id <- peak
                    }
            }
        }

    }
    return(proximal_id)
}

################
########################################
# this function returns the peak id of most distal and most proximal( to the distal) of the input gene.
# arguments:  gene name,   df: the modified version of DetectUTRLengthShift output.
## to modify the DetectUTRLengthShift out put use the add_gene_info function
get_utr <- function(gene, lengthened_utrs_df, all_expressed_peaks_df, utr3.ref) {
    tmp_df <- lengthened_utrs_df %>% filter(gene_id == gene) #ATP1B1
    tmp_df = tmp_df[tmp_df['SiteLocation']== max(unique(tmp_df['SiteLocation'])),] ## this finds the most distal peak for a gene
    most_distal_peak_id <- as.character(tmp_df$Row.names)
    peaks_to_map_to_ref_UTR <- as.data.frame(all_expressed_peaks_df[which(sub("(.*):.*:.*-.*:.*", "\\1", all_expressed_peaks_df$peaks) == gene),])
    colnames(peaks_to_map_to_ref_UTR) <- 'peak_ids'
    # map the expressed peaks of the gene to reference UTR locations
    if (nrow(peaks_to_map_to_ref_UTR) >= 1) {
    mapped_utrs <- map_to_ref_UTR(peaks_to_map_to_ref_UTR$peak_ids, utr3.ref, all_expressed_peaks_df)
    most_proximal_peak_id <- get_proximal_peak_id(most_distal_peak_id,
                                                  mapped_utrs, lengthened_utrs_df)

    res <- list(result= c(gene,most_distal_peak_id, most_proximal_peak_id),
                all_peaks=peaks_to_map_to_ref_UTR, mapped_utrs=mapped_utrs)
    dis_prox <- c(most_distal_peak_id, most_proximal_peak_id)
    } else {
	    dis_prox <- c(most_distal_peak_id, 'Not-found')
    }
    return(dis_prox)
    
}

################################################### you need the genes to operate and get their utr region (region between proximal and distal peak sites) 


#tst_res <- invisible(lapply(X = genes_to_operate$gene_id,FUN = get_utr, df=sig_upregulated_long_utrs))
#res_df <- t(as.data.frame(tst_res, col.names=F))
#colnames(res_df) <- c('distal_peak', 'proximal_peak')
#Genes_dis_prox_df <- data.frame(Genes_dis_prox_df, res_df)

## get the bed file from df that is outout lines above

get_utr_bed_file <- function(df){

    df['bed_out'] <- 'NAN'
    for (row in 1:nrow(df)){
        gene_id <- df[row, 'gene_id']
        dist_p <- df[row, 'distal_peak']
        prox_p <- df[row, 'proximal_peak']
        strand = sub(".*:.*:.*-.*:(.*)", "\\1", dist_p)
        chr = sub(".*:(.*):.*-.*:.*", "\\1", dist_p)
        if (strand == '1'){
            prox_p_start <- sub(".*:.*:(.*)-.*:.*", "\\1", prox_p)
            dist_p_end <- sub(".*:.*:.*-(.*):.*", "\\1", dist_p)
            utr_start_coordinate <- as.numeric(prox_p_start) - 1000  ## get 1kb inward of proximal peak
            utr_end_coordinate <- as.numeric(dist_p_end) + 200   ## add 200bp to the end of peak just to be safe
        } else if (strand == '-1'){
            dist_p_start <- sub(".*:.*:(.*)-.*:.*", "\\1", dist_p)
            prox_p_end <- sub(".*:.*:.*-(.*):.*", "\\1", prox_p)
            utr_start_coordinate <- as.numeric(dist_p_start) - 200 ## add 200bp to the end of peak just to be safe
            utr_end_coordinate <- as.numeric(prox_p_end) + 1000    ## get 1kb inward of proximal peak
        }
        fasta_header_name <- paste0(chr,':',gene_id,':',as.character(utr_start_coordinate),":",
                                  as.character(utr_end_coordinate))
        bed_out <- paste0(chr, '\t' ,as.character(utr_start_coordinate), '\t',
                          as.character(utr_end_coordinate), '\t', fasta_header_name)
        df[row, 'bed_out'] = bed_out

    }
    return(df)
}


####################################################################################
##### a function to add gene name and extra peak info the output of DetectUTRLengthShift
## arguments:  the df output of DetectUTRLengthShift,  peak.annotations: the peak annotation df from Sierra output
#peak.annotations <- read.table("sierra_AD_CT_mix/AD_CT_peaks_annotated.txt",
#                               header = TRUE,
#                               sep = "\t",
#                               row.names = 1,
#                               stringsAsFactors = FALSE)

add_gene_peak_info <- function(test.res, peak.annotations){
    tmp_df <- data.frame(chr = peak.annotations[rownames(peak.annotations) %in% rownames(test.res),]['seqnames'],
                     geneID = peak.annotations[rownames(peak.annotations) %in% rownames(test.res),]['gene_id'],
                     width = peak.annotations[rownames(peak.annotations) %in% rownames(test.res),]['width'],
                     start_pos = peak.annotations[rownames(peak.annotations) %in% rownames(test.res),]['start'],
                     end_pos = peak.annotations[rownames(peak.annotations) %in% rownames(test.res),]['end'])
    main_table <- test.res[1:9]
    final_table <- merge(main_table,tmp_df, by = 0)

    return(final_table)

}

############################################################################################
### a function to operate on the output of DetectUTRshift function,  it adds the information to the 
### df and then extract the significant longer UTRs and then return those genes names
### arguments:  dataframe from DetectUTRshift,   output: list of genes with significant UTR lengtheneing
get_sig_longer_UTRs <- function(df, peak.annotations){
    df2 <- add_gene_peak_info(df, peak.annotations)
    df2 <- df2 %>% filter(Log2_fold_change >= .5)
    df2['peak_relative_postion'] <- mapply(get_relative_pos, df2$SiteLocation, df2$NumSites)
    df2 <- df2 %>% filter(peak_relative_postion >= 0.5)
    df2$peak_ids <- df2$Row.names
    genes_to_operate <- data.frame(unique(df2['gene_id']))
    res <- list(genes=genes_to_operate, df=df2)
    return(res)
}


#######################################################################################
#######################################################################################
######### wrapper########################33

get_UTR_seqs <- function(results, so, peaks_ann, gtf_TxDb){
    ### 2-  get the reference UTR regions for all the genes
    utr3.ref <- GenomicFeatures::threeUTRsByTranscript(gtf_TxDb)
    utr3.ref <- unlist(utr3.ref)
    for (case in results){
        tryCatch(
            {
            name <- case[[1]]
            df <- get(case[[2]])
            print(paste('processing the ', name, 'results'))



            ### 4- get all the expressed peaks in the cell population you are studying
            all.peaks.expressed <- as.data.frame(GetExpressedPeaks(so,
                                                                   population.1 = name))
            colnames(all.peaks.expressed) <- 'peaks'

            ### 5- get the significantly lengthened and expressed UTRs with the designed function
            ###    this will add more info to the input data frame
            res <- get_sig_longer_UTRs(df,peaks_ann)
            tmp_genes_to_operate <- res[['genes']]
            tmp_lenghtened_UTRs <- res[['df']]

            ### 6- get the utr region using the distal and proximal peaks locations.
            ###    there are a lot going on in this step to identify the right expressed
            ###    most proximal peak location and then retruning proximal and distal peak locations for genes
	    tmp_lengthened_genes_UTR_reg <- invisible(lapply(X = tmp_genes_to_operate$gene_id,FUN = get_utr,
                                                             lengthened_utrs_df=tmp_lenghtened_UTRs,
                                                             all_expressed_peaks_df=all.peaks.expressed,
                                                             utr3.ref=utr3.ref))
	    for (i in 1:length(tmp_lengthened_genes_UTR_reg)){
    		if (length(tmp_lengthened_genes_UTR_reg[[i]]) > 2){
        	tmp_lengthened_genes_UTR_reg[[i]] <- tmp_lengthened_genes_UTR_reg[[i]][1:2]
    		}
		}
	    tmp_lengthened_genes_UTR_reg <- data.frame(tmp_lengthened_genes_UTR_reg)
	    tmp_lengthened_genes_UTR_reg <- t(as.data.frame(tmp_lengthened_genes_UTR_reg, col.names=F))
            colnames(tmp_lengthened_genes_UTR_reg) <- c('distal_peak', 'proximal_peak')
            tmp_Genes_dis_prox_df <- data.frame(tmp_genes_to_operate, tmp_lengthened_genes_UTR_reg)
            tmp_Genes_dis_prox_df = tmp_Genes_dis_prox_df[!grepl("Not-found", tmp_Genes_dis_prox_df$proximal_peak),]


            ### 7- get the bed file for the UTR region between most proximal and distal peaks and save the file
	    tmp_lengthened_UTR_region_bed <- get_utr_bed_file(tmp_Genes_dis_prox_df)
            outname <- paste0('results/',name,'_upregulated_longer_utr_region.bed')
            write.table(tmp_lengthened_UTR_region_bed['bed_out'], file=outname,
                       row.names=F, col.names=F, quote=F)
                },
            error=function(cond) {
                message(paste("couldnt get any UTR seqs for the population:", name))
                message("Here's the original error message:")
                message(cond)
                # Choose a return value in case of error
                return(NA)
            })
    }
}
