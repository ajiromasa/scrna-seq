# pseudobulk analysis, Uba1 mutant hematocytes.
# 2026/08/25

# library ----
lapply(c("tidyverse", "edgeR", "fgsea", "org.Hs.eg.db",
         "ComplexHeatmap", "circlize", "grid"),
       library, character.only = TRUE)
# DEG ----
for(i in c("HSC","LMPP","MPP","CMP","GMP","MEP")){
  df <- read.table(paste0(i,"_counts.tsv"),
                   sep = "\t",
                   row.names = 1,
                   header = T)
  df <- arrange(df, row.names(df))
  name <- paste0(i,"_count")
  assign(name,df)
}

HSC_count <- HSC_count[,c(1,5,2,4,6,3)]
MPP_count <- MPP_count[,c(1,6,2,3,4,5)]
LMPP_count <- LMPP_count[,c(2,1,3,4,5,6)]
CMP_count <- CMP_count[,c(2,1,6,3,5,4)]
MEP_count <- MEP_count[,c(3,2,6,1,5,4)]
GMP_count <- GMP_count[,c(1,6,2,4,5,3)]

LK_count <- HSC_count + MPP_count + LMPP_count
nLK_count <- CMP_count + MEP_count + GMP_count

group<- c("mt","mt", "mt","wt","wt","wt")

for (i in c("HSC","LMPP","MPP",
            "CMP","GMP","MEP",
            "LK","nLK")){
  df <- get(paste0(i,"_count"))
  C <- DGEList(df, group = group)
  C <- calcNormFactors(C) 
  C <- estimateCommonDisp(C)
  C <- estimateTagwiseDisp(C)
  nc <- as.data.frame(C$pseudo.counts) %>% arrange(row.names(.))
  comp <- exactTest(C, pair = c("wt", "mt"))
  tmp <- topTags(comp, n=nrow(comp$table))
  deg <- arrange(tmp$table, rownames(tmp$table))
  deg.cpm <- subset(deg, deg$logCPM >= log2(3))
  deg.cpm.p <- subset(deg.cpm.fc, deg.cpm.fc$logFC >= log2(1.25))
  deg.cpm.n <- subset(deg.cpm.fc, deg.cpm.fc$logFC <= -log2(1.25))
  df2 <- deg.cpm %>% mutate(metrix = -log10(PValue)*sign(logFC))
  ranks <- df2$metrix
  names(ranks) <- row.names(deg.cpm)
  ranks <- na.omit(replace(ranks, is.infinite(ranks), NA))
  
  # fgsea
  pathways <- gmtPathways("../GSEA/msigdb_v2023.2.Mm_GMTs/mh.all.v2023.2.Mm.symbols.gmt")
  fgseaRes <- fgseaMultilevel(pathways = pathways, stats = ranks)
  fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES))
  
  df.gsea <- data.frame(pathway = fgseaResTidy$pathway,
                        pval = fgseaResTidy$pval,
                        padj = fgseaResTidy$padj,
                        ES = fgseaResTidy$ES,
                        NES = fgseaResTidy$NES,
                        size = fgseaResTidy$size,
                        leadingEdge = sapply(fgseaResTidy$leadingEdge,
                                             function(x) paste(x, collapse = ", "))
  )
  assign(paste0(i,"_nc"), nc)
  assign(paste0(i,"_fc1"), deg)
  assign(paste0(i,"_fc"), deg.cpm)
  assign(paste0(i,"_fc_p"), deg.cpm.p)
  assign(paste0(i,"_fc_n"), deg.cpm.n)
  assign(paste0(i,"_gsea"), df.gsea)
}

path <- c(
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING"
)

type <- c("HSC", "LMPP", "MPP","CMP", "GMP", "MEP","LK", "nLK")

## matrix
padj_mat <- matrix(
  NA_real_,
  nrow = length(path),
  ncol = length(type),
  dimnames = list(path, type)
)

nes_mat <- matrix(
  NA_real_,
  nrow = length(path),
  ncol = length(type),
  dimnames = list(path, type)
)

for (i in type) {
  gsea <- get(paste0(i, "_gsea"))
  idx <- match(path, gsea$pathway)
  padj_mat[, i] <- gsea$padj[idx]
  nes_mat[, i]  <- gsea$NES[idx]
}

logp_mat <- -log10(padj_mat)
threshold <- -log10(0.25)
finite_max <- max(logp_mat[is.finite(logp_mat)], na.rm = TRUE)

if (any(is.infinite(logp_mat))) {
  logp_mat[is.infinite(logp_mat)] <- finite_max
}

max_logp <- max(logp_mat, na.rm = TRUE)
red_fun <- colorRamp2(
  c(threshold, max_logp),
  c("#F2F2F2", "red")
)

blue_fun <- colorRamp2(
  c(threshold, max_logp),
  c("#F2F2F2", "blue")
)

ht <- Heatmap(
  logp_mat,
  name = "-log10(padj)",
  col = c("white" = "white"),
  show_heatmap_legend = F,
  rect_gp = gpar(type = "none"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  column_names_side = "top",
  cell_fun = function(j, i, x, y, width, height, fill) {
    val <- logp_mat[i, j]
    nes <- nes_mat[i, j]
    if (is.na(val) || is.na(nes)) {
      cell_col <- "black"
    } else if (val < threshold) {
      cell_col <- "gray70"
    } else if (nes > 0) {
      cell_col <- red_fun(val)
    } else if (nes < 0) {
      cell_col <- blue_fun(val)
    } else {
      cell_col <- "gray70"
    }
    grid.rect(
      x = x,
      y = y,
      width = width,
      height = height,
      gp = gpar(
        fill = cell_col,
        col = "black",
        lwd = 1
      )
    )
  }
)
lgd <- Legend(
  title = expression(-log[10](padj)),
  at = c(threshold, max_logp),
  labels = c(
    sprintf("%.2f", threshold),
    sprintf("%.2f", max_logp)
  ),
  col_fun = colorRamp2(
    c(-max_logp, -threshold, threshold, max_logp),
    c("blue", "#F2F2F2", "#F2F2F2", "red")
  )
)
lgd_ns <- Legend(
  labels = "padj > 0.25",
  legend_gp = gpar(fill = "gray60"),
  title = NULL
)
draw(
  ht,
  annotation_legend_list = list(lgd,lgd_ns),
  padding = unit(c(5,40,5,5),"mm")
)

# Heatmap projection by pathways ----

target_pathways <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_INFLAMMATORY_RESPONSE")
target_pathways_2 <- c("HALLMARK_INTERFERON_ALPHA_RESPONSE","HALLMARK_INTERFERON_GAMMA_RESPONSE")

get_leading_edge <- function(gsea_df){
  df <- subset(
    gsea_df,
    pathway %in% target_pathways & pval < 0.05 # no cutoff set for target_pathways_2/LMPP
  )
  unique(trimws(unlist(strsplit(df$leadingEdge, ","))))
}

genes <- get_leading_edge(__gsea_output__)
colnames(__cell_types_nc__)[1:6] <- c(paste0("__colname_1__", 1:3),paste0("__colname_2__", 1:3))
mat_all <- HSC_nc[, 1:6]
genes_use <- intersect(genes, rownames(mat_all))
mat_heat <- as.matrix(mat_all[genes_use, ])
mat_z <- t(scale(t(mat_heat)))
mat_z[is.na(mat_z)] <- 0
mat_z[is.infinite(mat_z)] <- 0
group <- factor(
  c(rep("__group_name_1__", 3), rep("__group_name_2__", 3)),
  levels = c("__group_name_1__", "__group_name_2__")
)

ha <- HeatmapAnnotation(
  group = group,
  col = list(
    group = c(
      __group_name_1__ = "#ff00ff",
      __group_name_2__ = "#00cccc"
    )
  )
)

col_fun <- colorRamp2(
  c(-2, 0, 2),
  c("white", "#ffb266", "red")
)

Heatmap(
  mat_z,
  name = "row Z-score",
  col = col_fun,
  top_annotation = ha,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = T,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 7),
  column_names_gp = gpar(fontsize = 8),
  column_names_rot = 45,
  border = TRUE,
  rect_gp = gpar(col = "grey80", lwd = 0.3)
)

fc2 <- __cell_type_fc__ %>% mutate(group = "ns")
fc_p2 <- __cell_type_fc_p__ %>% mutate(group = "positive")
fc_n2 <- __cell_type_fc_n__ %>% mutate(group = "negative")
df_long <- rbind(__cell_type_fc2__,
                 __cell_type_fc_p2__,
                 __cell_type_fc_n2__)

ggplot(df_long, aes(x = logFC, y = -log10(FDR), fill = group)) +
  geom_jitter(
    shape = 21,
    size = 2.5,
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      "ns" = "grey80",
      "positive" = "#ff0000",
      "negative" = "#6666ff"
    )
  ) +
  geom_text_repel(
    data = subset(
      df_long,
      group != "ns" &
        grepl(genes_label, gene)
    ),
    aes(label = row.names(subset(df_long, group != "ns" & grepl(genes_label, row.names(df_long))))),
    nudge_y = 0.5,
    min.segment.length = 0,
    segment.color = "black",
    segment.size = 0.3,
    segment.linetype = 3,
    size = 4
  ) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    panel.grid = element_blank()
  ) 
dev.off()

sessionInfo()
R version 4.4.2 (2024-10-31)
Platform: aarch64-apple-darwin20
Running under: macOS Sequoia 15.7.7

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

locale:
  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: Asia/Tokyo
tzcode source: internal

attached base packages:
  [1] stats4    grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
  [1] fgsea_1.32.4           ggrepel_0.9.6          clusterProfiler_4.14.6 org.Hs.eg.db_3.20.0    AnnotationDbi_1.68.0   Biobase_2.66.0        
[7] gplots_3.2.0           edgeR_4.4.2            limma_3.62.2           rstatix_0.7.3          ggpubr_0.6.3           rtracklayer_1.66.0    
[13] GenomicRanges_1.58.0   GenomeInfoDb_1.42.3    IRanges_2.40.1         S4Vectors_0.44.0       BiocGenerics_0.52.0    matrixStats_1.5.0     
[19] circlize_0.4.16        ComplexHeatmap_2.22.0  lubridate_1.9.4        forcats_1.0.1          stringr_1.6.0          dplyr_1.1.4           
[25] purrr_1.2.0            readr_2.1.6            tidyr_1.3.1            tibble_3.3.0           ggplot2_4.0.1          tidyverse_2.0.0       

loaded via a namespace (and not attached):
  [1] RColorBrewer_1.1-3          rstudioapi_0.17.1           jsonlite_2.0.0              shape_1.4.6.1               magrittr_2.0.4             
[6] ggtangle_0.0.8              farver_2.1.2                GlobalOptions_0.1.2         fs_1.6.6                    BiocIO_1.16.0              
[11] zlibbioc_1.52.0             vctrs_0.6.5                 memoise_2.0.1               Rsamtools_2.22.0            RCurl_1.98-1.17            
[16] ggtree_3.14.0               S4Arrays_1.6.0              curl_7.0.0                  broom_1.0.10                gridGraphics_0.5-1         
[21] SparseArray_1.6.2           Formula_1.2-5               KernSmooth_2.23-26          plyr_1.8.9                  cachem_1.1.0               
[26] GenomicAlignments_1.42.0    igraph_2.2.1                lifecycle_1.0.4             iterators_1.0.14            pkgconfig_2.0.3            
[31] gson_0.1.0                  Matrix_1.7-4                R6_2.6.1                    fastmap_1.2.0               GenomeInfoDbData_1.2.13    
[36] MatrixGenerics_1.18.1       clue_0.3-66                 aplot_0.2.9                 enrichplot_1.26.6           digest_0.6.39              
[41] colorspace_2.1-2            patchwork_1.3.2             RSQLite_2.4.4               labeling_0.4.3              timechange_0.3.0           
[46] httr_1.4.7                  abind_1.4-8                 mgcv_1.9-4                  compiler_4.4.2              bit64_4.6.0-1              
[51] withr_3.0.2                 doParallel_1.0.17           S7_0.2.1                    backports_1.5.0             BiocParallel_1.40.2        
[56] carData_3.0-5               DBI_1.2.3                   R.utils_2.13.0              ggsignif_0.6.4              rappdirs_0.3.3             
[61] DelayedArray_0.32.0         rjson_0.2.23                gtools_3.9.5                caTools_1.18.3              tools_4.4.2                
[66] ape_5.8-1                   R.oo_1.27.1                 glue_1.8.0                  restfulr_0.0.16             nlme_3.1-168               
[71] GOSemSim_2.32.0             reshape2_1.4.5              cluster_2.1.8.1             generics_0.1.4              gtable_0.3.6               
[76] tzdb_0.5.0                  R.methodsS3_1.8.2           data.table_1.17.8           hms_1.1.4                   car_3.1-3                  
[81] utf8_1.2.6                  XVector_0.46.0              foreach_1.5.2               pillar_1.11.1               yulab.utils_0.2.1          
[86] splines_4.4.2               treeio_1.30.0               lattice_0.22-7              bit_4.6.0                   tidyselect_1.2.1           
[91] GO.db_3.20.0                locfit_1.5-9.12             Biostrings_2.74.1           SummarizedExperiment_1.36.0 statmod_1.5.1              
[96] stringi_1.8.7               UCSC.utils_1.2.0            lazyeval_0.2.2              ggfun_0.2.0                 yaml_2.3.10                
[101] codetools_0.2-20            qvalue_2.38.0               ggplotify_0.1.3             cli_3.6.5                   dichromat_2.0-0.1          
[106] Rcpp_1.1.0                  png_0.1-8                   XML_3.99-0.20               parallel_4.4.2              blob_1.2.4                 
[111] DOSE_4.0.1                  bitops_1.0-9                tidytree_0.4.6              scales_1.4.0                crayon_1.5.3               
[116] GetoptLong_1.0.5            rlang_1.1.6                 fastmatch_1.1-6             cowplot_1.2.0               KEGGREST_1.46.0           
