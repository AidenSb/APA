### Useful functions!



## get the relative locoation and you can decide on the distal vs proximal 
get_relative_pos <- function(pos,sitNum){
    relative_loc <- (pos - 1) / (sitNum -1)
    return(relative_loc)
}

map_to_ref_UTR <- function(peak_id){
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
    peaks_mapped_to_ref_UTR ## the final table with all the peaks mapped to their ref 3UTR
    
    ## get the peaks that expressed in AD
    idx = which(peaks_mapped_to_ref_UTR$peak_ID %in% tst.peaks.expressed) 
    peaks_mapped_to_ref_UTR = peaks_mapped_to_ref_UTR[idx,]
    return(peaks_mapped_to_ref_UTR)
    
}

################
## this function takes in the distal peak and all the expressed UTR peaks for the gene and returns most 
## proximal non DE peak coordinates
get_proximal_peak_id <- function(distal_peak, mapped_utrs, df){
    distance = 1e15
    for (peak in unlist(unique(mapped_utrs['peak_ID']))){
        if (!peak %in% df$peak_ids){
            strand = sub(".*:.*:.*-.*:(.*)", "\\1", distal_peak)
            distal_p_start = sub(".*:.*:(.*)-.*:.*", "\\1", distal_peak)
            distal_p_end = sub(".*:.*:.*-(.*):.*", "\\1", distal_peak)
            proximal_p_start = sub(".*:.*:(.*)-.*:.*", "\\1", peak)
            proximal_p_end = sub(".*:.*:.*-(.*):.*", "\\1", peak)
            if (strand=='1'){
                current_dist <- as.numeric(distal_p_end) - as.numeric(proximal_p_start)
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
get_utr <- function(gene, df) {
    tmp_df <- df %>% filter(gene_id == gene) #ATP1B1
    tmp_df = tmp_df[tmp_df['SiteLocation']== max(unique(tmp_df['SiteLocation'])),] ## this finds the most distal peak for a gene
    most_distal_peak_id <- as.character(tmp_df$Row.names)
    peaks_to_map_to_ref_UTR <- as.data.frame(AD.peaks.expressed[which(startsWith(AD.peaks.expressed$AD.peaks.expressed,gene)),])
    colnames(peaks_to_map_to_ref_UTR) <- 'peak_ids'
    # map the expressed peaks of the gene to reference UTR locations
    mapped_utrs <- map_to_ref_UTR(peaks_to_map_to_ref_UTR$peak_ids)
    most_proximal_peak_id <- get_proximal_peak_id(most_distal_peak_id,
                                                  mapped_utrs, df)

    res <- list(result= c(gene,most_distal_peak_id, most_proximal_peak_id),
                all_peaks=peaks_to_map_to_ref_UTR, mapped_utrs=mapped_utrs)
    dis_prox <- c(most_distal_peak_id, most_proximal_peak_id)
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





