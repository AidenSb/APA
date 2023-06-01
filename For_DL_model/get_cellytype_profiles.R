install.packages('Seurat')
library(Seurat)
als_atlas <- readRDS('/data1/APA/Paul_ALS_Data/ALS_snRNA_final_with_scVI.RDS')


## lets save the 
sub_cells <- c()

for (ct in unique(als_atlas$diagnosis_Velm_cellsubtype)){
  if (grepl( 'control', ct, fixed = TRUE)){
    sub_cells <- c(sub_cells, ct)
  }
}

main_cells <- c()

for (ct in unique(als_atlas$diagnosis_celltype)){
  if (grepl( 'CTRL', ct, fixed = TRUE)){
    main_cells <- c(main_cells, ct)
  }
}



all_avg_embeddings = list()
for (sub_ct in sub_cells){
  idx = which(als_atlas$diagnosis_Velm_cellsubtype == sub_ct)
  all_avg_embeddings[[sub_ct]] <- colMeans(Embeddings(object = als_atlas[["scvi"]])[idx,])
}

for (sub_ct in main_cells){
  idx = which(als_atlas$diagnosis_celltype == sub_ct)
  all_avg_embeddings[[sub_ct]] <- colMeans(Embeddings(object = als_atlas[["scvi"]])[idx,])
}

ctrl_profile <- data.frame(all_avg_embeddings)
######################################################################
rbps <- read.table('/data1/APA/Paul_ALS_Data/RBP_Information.txt', sep='\t', header=T)
######
RBPs_to_add <- rownames(als_atlas)[(rownames(als_atlas) %in% rbps$RBP_Name)]
##
table(als_atlas$diagnosis_celltype)
als_atlas <- SetIdent(als_atlas, value = als_atlas$diagnosis_celltype)
# make C9ALS profiles:
main_rbp_profiles <- list()
for (ct in unique(als_atlas$celltype)){
  diff <- FindMarkers(
    als_atlas,
    ident.1 = paste0('C9ALS_', ct),
    ident.2 = paste0('CTRL_', ct),
    features = RBPs_to_add,
    slot = "data",logfc.threshold = 0.05,
    test.use = "wilcox",
    min.pct = 0.05
    )
  main_rbp_profiles[[ct]] <- diff
}
##
table(als_atlas$diagnosis_Velm_cellsubtype)
als_atlas <- SetIdent(als_atlas, value = als_atlas$diagnosis_Velm_cellsubtype)
# make C9ALS profiles:
sub_rbp_profiles <- list()
for (ct in unique(als_atlas$Velm_cellsubtype)){
  diff <- FindMarkers(
    als_atlas,
    ident.1 = paste0('C9ALSFTLD_', ct),
    ident.2 = paste0('control_', ct),
    features = RBPs_to_add,
    slot = "data",logfc.threshold = 0.05,
    test.use = "wilcox",
    min.pct = 0.05
  )
  sub_rbp_profiles[[ct]] <- diff
}

for (ct in unique(als_atlas$Velm_cellsubtype)){
  print(dim(sub_rbp_profiles[[ct]] ))
}

## combind 
cmb = c(sub_rbp_profiles, main_rbp_profiles)

### 

tst = merge(main_rbp_profiles[['Inhibitory']]['avg_log2FC'], sub_rbp_profiles[[ "L4" ]]['avg_log2FC'],
            all = T, by = 'row.names',suffixes =c('_Inh','_L4') )

##
ct_names = unique(c(names(sub_rbp_profiles), names(main_rbp_profiles)))
merged_dfs = main_rbp_profiles[[ct_names[1]]]['avg_log2FC']

for (i in 2:length(ct_names)-1){
  print(i)
  merged_dfs = merge(main_rbp_profiles[[ct_names[i]]]['avg_log2FC'], sub_rbp_profiles[[ct_names[i+1]]]['avg_log2FC'],
              all = T, by = 'row.names',suffixes =c(paste0('_', ct_names[i]),paste0('_', ct_names[i+1])) )
}

r = c()
for (name in names(sub_rbp_profiles)){
  r = c(r, rownames(sub_rbp_profiles[[name]]))
}
r = unique(r)
vals = rep(0, length(r))
temp_df = data.frame(vals)
rownames(temp_df) = r
cnames <- c()
for (i in 2:length(sub_rbp_profiles)-1){
  print(i)
  temp_df = merge(temp_df, sub_rbp_profiles[[names(sub_rbp_profiles)[i]]]['avg_log2FC'],
                     all = T, by = 'row.names')
  cnames = c(cnames, names(sub_rbp_profiles)[i] )
}

###
library(purrr)
library(dplyr)

# Extract the tst column from each dataframe in main_dfs and convert row.names to a regular column
tst_columns <- map2(sub_rbp_profiles, names(sub_rbp_profiles), function(x, name) {
  data.frame(rn = row.names(x), tst = x[, "avg_log2FC"]) %>%
    setNames(c("rn", name))
})

# Bind the tst columns together into a single dataframe
big_df <- reduce(tst_columns, full_join, by = "rn")


####
tst_columns <- map2(main_rbp_profiles, names(main_rbp_profiles), function(x, name) {
  data.frame(rn = row.names(x), tst = x[, "avg_log2FC"]) %>%
    setNames(c("rn", name))
})

# Bind the tst columns together into a single dataframe
big_df_main <- reduce(tst_columns, full_join, by = "rn")
##

final_df <- merge(big_df,big_df_main, by='rn')
### 
clns <- c("rn","Excitatory","Astrocytes","Inhibitory","L5/6-CC", "AST-PP","IN-SST","Oligodendrocytes.x",
          "Microglia.y", "OPC.y", "L2/3","IN-VIP","IN-PV","L4","AST-FB", "IN-SV2C", "L5/6")
#
final_df <- final_df[clns]
colnames(final_df) <- c("rn","Excitatory","Astrocytes","Inhibitory","L5/6-CC", "AST-PP","IN-SST","Oligodendrocytes",
                        "Microglia", "OPC", "L2/3","IN-VIP","IN-PV","L4","AST-FB", "IN-SV2C", "L5/6")
#
final_df[is.na(final_df)] <- 0
#########
rownames(final_df) <- final_df$rn
##
clns <- c("Excitatory","Astrocytes","Inhibitory","L5/6-CC", "AST-PP","IN-SST","Oligodendrocytes",
          "Microglia", "OPC", "L2/3","IN-VIP","IN-PV","L4","AST-FB", "IN-SV2C", "L5/6")

final_df <- select(final_df, -1)
##

profile <- select(ctrl_profile, 1:18, -16)
colnames(profile) <- c('L5-6-CC','Oligodendrocytes','AST-PP','IN-SST','Microglia','OPC','L2-3',
                            'IN-VIP','IN-PV','L4','AST-FB','IN-SV2C','L5-6','Endothelial', 'Excitatory',
                            'Astrocytes','Inhibitory')

###
colnames(final_df) <- c( "Excitatory","Astrocytes","Inhibitory", "L5-6-CC", "AST-PP","IN-SST","Oligodendrocytes",
                         "Microglia", "OPC","L2-3","IN-VIP","IN-PV","L4", "AST-FB","IN-SV2C", "L5-6")

order <-   c('Excitatory','Astrocytes','Inhibitory','Oligodendrocytes','Microglia','OPC','L5-6-CC','L2-3',
             'AST-PP','IN-SST','IN-VIP','IN-PV','L4','AST-FB','IN-SV2C','L5-6')
final_df <- select(final_df, order)

profile <- select(profile, order)

final_profile <- rbind(profile, final_df)

####
write.table(final_profile, file='/data1/APA/Paul_ALS_Data/bams_in/subscelltype_bamfiles/Mapper_outs/V2/celltype_profiles.tsv', 
            sep='\t',quote = F, row.names = T)
