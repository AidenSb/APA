c9als_exc = read.csv('/home/aiden/codes/APA_stuff/post_qual/APA/For_ALS_atlas_paper/C9ALSFTLD_Excitatory_utr_shift_table.tsv',
sep = '\t', header=T)
head(c9als_exc)


library(dplyr)

c9als_exc <- c9als_exc %>% filter(FC_direction=='Up') %>% filter(SiteLocation > NumSites/2)

library(stringr)

ids = rownames(c9als_exc)

gene_names <- str_split(ids, ':', simplify = T)[0:length(ids),1]

c9als_exc$geneIDs <- gene_names
head(c9als_exc)


c9als_exc <- c9als_exc %>% arrange(desc(Log2_fold_change))
head(c9als_exc)

out_df <- data.frame(c9als_exc$geneIDs, c9als_exc$Log2_fold_change)
colnames(out_df) <- c('geneID', "LFC")

length(out_df$geneID)
length(unique(out_df$geneID))


out_df <- out_df[!duplicated(out_df$geneID),]
head(out_df)

tail(out_df)

write.table(out_df, file='/home/aiden/codes/APA_stuff/post_qual/APA/For_ALS_atlas_paper/C9ALS_Exc_lengthened_transcripts.tsv',
            sep='\t',quote=F, row.names = F)


