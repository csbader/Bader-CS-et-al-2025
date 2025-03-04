#######################################
####### Bader CS, et al. 2025.  #######
#######   Fig. 4, S2, S5, S6    #######
#######   scRNA-seq Analysis    #######
#######################################


## Load packages ##

library(Seurat)
library(BPCells)
library(tidyverse)
library(data.table)
library(scRepertoire)
library(clipr)
library(RColorBrewer)

#### mRNA ####

# Load rna.rds file
rna <- load("data/rna.rds")

# Load vdj.rds file
vdj <- load("data/vdj.rds")


#### RNA Quality Control ####

# Add Log10GenesPerUMI
rna$log10GenesPerUMI <- log10(rna$nFeature_RNA) / log10(rna$nCount_RNA)

# Visualize Data

# Cells per sample
metadata <- rna@meta.data
metadata %>%
  ggplot(aes(x=sample.name, fill=sample.name)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

# Transcripts per cell
metadata %>%
  ggplot(aes(color=sample.name, x=nCount_RNA, fill=sample.name)) +
  geom_density(alpha=0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 60)

# Genes detected per cell
metadata %>%
  ggplot(aes(color=sample.name, x=nFeature_RNA, fill=sample.name)) +
  geom_density(alpha=0.2) +
  theme_classic() +
  scale_x_log10() +
  geom_vline(xintercept=30)

# Visualize the distribution of genes detected per cell via boxplot
metadata %>% 
  ggplot(aes(x=sample.name, y=log10(nFeature_RNA), fill=sample.name)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")

# Correlation between genes detected and number of UMIs
metadata %>%
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA)) +
  geom_point() +
  stat_smooth(method=lm) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  geom_vline(xintercept = 4000, color = "black") +
  geom_hline(yintercept = 100, color = "black") +
  geom_hline(yintercept = 10, color = "black") +
  facet_wrap(~sample.name)

# Visualize the overall complexity of gene expression
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color=sample.name, fill=sample.name)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept=0.6)

# Filter outliers
filtered_rna <- rna %>%
  subset(subset = (nFeature_RNA >= 30) &
           (nCount_RNA >= 60) &
           (log10GenesPerUMI > 0.6))
# removed 48187 cells

# Remove genes that are expressed in less than 10 cells
counts <- GetAssayData(object = filtered_rna, slot = "counts")
nonzero <- counts > 0
keep_genes <- rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]
filtered_rna <- CreateSeuratObject(filtered_counts, meta.data = filtered_rna@meta.data)
# removed 7 genes

#### Write to bpcells to run on SSD instead of memory ####

# Write the counts layer to a directory
filtered_rna[["RNA"]]$counts <- as(object = filtered_rna[["RNA"]]$counts, Class = "dgCMatrix")
write_matrix_dir(mat = filtered_rna[["RNA"]]$counts, dir = "data/bpcells/matrix")

# Save metadata
metadata <- filtered_rna@meta.data
saveRDS(metadata, file = "data/metadata.rds")

# Write settings, clear environment and reload data
options(future.globals.maxSize = 3e+09)
options(Seurat.object.assay.version = "v5")

counts.mat <- open_matrix_dir(dir = "data/bpcells/matrix")   
metadata <- readRDS("data/metadata.rds")

# Create Seurat object
filtered_rna <- CreateSeuratObject(counts = counts.mat, meta.data = metadata)

#### Examine Unintegrated Data ####

# Run standard Seurat pipeline
filtered_rna <- SCTransform(filtered_rna)
filtered_rna <- RunPCA(filtered_rna, verbose = TRUE)
filtered_rna <- FindNeighbors(filtered_rna, dims = 1:30, reduction = "pca")
filtered_rna <- FindClusters(filtered_rna, resolution = 2, cluster.name = "unintegrated_clusters")
filtered_rna <- RunUMAP(filtered_rna, dims = 1:30, reduction = "pca", reduction.name = "umap unintegrated")

# View unintegrated UMAP by sample
DimPlot(filtered_rna, reduction = "umap.unintegrated", group.by = "sample.name", raster = FALSE)

#### Integrate data ####

# Split Seurat object by sample name as integration variable
split_rna <- SplitObject(filtered_rna, split.by = "sample.name")

# Normalize data by sample name
split_rna <- lapply(X = split_rna, FUN = function (x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})

# Select integration features
features <- SelectIntegrationFeatures(object.list = split_rna)

# Scale and run PCA by sample name
split_rna <- lapply(X = split_rna, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

# Find integration anchors
anchors <- FindIntegrationAnchors(object.list = split_rna, reduction = "rpca",
                                  dims = 1:30)

# Integrate data
integrated_rna <- IntegrateData(anchorset = anchors, dims = 1:30)

# Scale integrated data and run PCA
integrated_rna <- ScaleData(object = integrated_rna)
integrated_rna <- RunPCA(object = integrated_rna)

# Choose highly variable PCs
ElbowPlot(integrated_rna, ndims = 50) # 18 dims

# View integrated UMAP
integrated_rna <- RunUMAP(integrated_rna,
                          dims = 1:18,
                          reduction = "pca")

DimPlot(integrated_rna, group.by = "sample.name", raster = FALSE)

#### Subset Tcon sample.tags ####

tcon_seurat <- subset(integrated_rna, cell.type == "Tcon") #120269 cells

#### Tcon Clustering ####

# Re-cluster subsetted T cells
tcon_seurat <- RunPCA(object = tcon_seurat)
ElbowPlot(tcon_seurat, ndims = 50) # 18 dims

# UMAP
tcon_seurat <- RunUMAP(tcon_seurat,
                       dims = 1:18,
                       reduction = "pca")

# Find clusters at resolution 0.6
tcon_seurat <- FindNeighbors(object = tcon_seurat,
                             dims = 1:18)

tcon_seurat <- FindClusters(object = tcon_seurat,
                            resolution = 0.6)

# Plot UMAP
Idents(object = tcon_seurat) <- "integrated_snn_res.0.6" # 20 clusters
DimPlot(tcon_seurat,
        reduction = "umap",
        label = TRUE,
        label.size = 6,
        raster = FALSE)

# Find all markers
DefaultAssay(tcon_seurat) <- "SCT"
tcon_seurat <- JoinLayers(tcon_seurat)
tcon_seurat <- PrepSCTFindMarkers(tcon_seurat, assay = "SCT", verbose = TRUE)
rna_markers <- FindAllMarkers(object = tcon_seurat,
                              only.pos = TRUE,
                              logfc.threshold = 0.25,
                              assay = "SCT")

view(rna_markers)

# Remove myeloid and B cells (first time)

tcon_seurat <- subset(tcon_seurat, idents = c("4", "13", "15"), invert = TRUE) # removed 14223 cells

# Re-cluster remaining cells

# UMAP
DefaultAssay(tcon_seurat) <- "integrated"

tcon_seurat <- RunUMAP(tcon_seurat,
                       dims = 1:18,
                       reduction = "pca")

# Find clusters at resolution 0.4
tcon_seurat <- FindNeighbors(object = tcon_seurat,
                             dims = 1:18)

tcon_seurat <- FindClusters(object = tcon_seurat,
                            resolution = c(0.4, 1.2))

# Resolution = 0.4
Idents(object = tcon_seurat) <- "integrated_snn_res.0.4"

DimPlot(tcon_seurat,
        label = TRUE,
        raster = FALSE,
        split.by = "sample.name") + NoLegend()

DefaultAssay(tcon_seurat) <- "SCT"
tcon_seurat <- PrepSCTFindMarkers(tcon_seurat, assay = "SCT", verbose = TRUE)
rna_markers <- FindAllMarkers(object = tcon_seurat,
                              only.pos = TRUE,
                              logfc.threshold = 0.25,
                              assay = "SCT")

view(rna_markers)

# Remove myeloid and B cells again (second time)

tcon_seurat <- subset(tcon_seurat, idents = c("13"), invert = TRUE) # removed 736 cells

# Remove B cells (third round after increasing resolution to 1.2)

# Resolution = 1.2
Idents(object = tcon_seurat) <- "integrated_snn_res.1.2"

DimPlot(tcon_seurat,
        label = TRUE,
        raster = FALSE,
        split.by = "sample.name") + NoLegend()

DefaultAssay(tcon_seurat) <- "SCT"
tcon_seurat <- PrepSCTFindMarkers(tcon_seurat, assay = "SCT", verbose = TRUE)
rna_markers <- FindAllMarkers(object = tcon_seurat,
                              only.pos = TRUE,
                              logfc.threshold = 0.25,
                              assay = "SCT")

View(rna_markers)

tcon_seurat <- subset(tcon_seurat, idents = c("4"), invert = TRUE) # removed 7743 cells

# Re-cluster final remaining cells

# UMAP
DefaultAssay(tcon_seurat) <- "integrated"

tcon_seurat <- RunUMAP(tcon_seurat,
                       dims = 1:18,
                       reduction = "pca")

# Find clusters
tcon_seurat <- FindNeighbors(object = tcon_seurat,
                             dims = 1:18)

tcon_seurat <- FindClusters(object = tcon_seurat,
                            resolution = 1.2)

# Resolution = 1.2
Idents(object = tcon_seurat) <- "integrated_snn_res.1.2"

DimPlot(tcon_seurat,
        label = TRUE,
        raster = FALSE,
        split.by = "sample.name") + NoLegend()

DefaultAssay(tcon_seurat) <- "SCT"
tcon_seurat <- PrepSCTFindMarkers(tcon_seurat, assay = "SCT", verbose = TRUE)
rna_markers <- FindAllMarkers(object = tcon_seurat,
                              only.pos = TRUE,
                              logfc.threshold = 0.25,
                              assay = "SCT")

view(rna_markers)

#### Tcon Annotation ####

# Fig. 4B
DimPlot(tcon_seurat,
        label = TRUE, raster = FALSE, label.size = 8) +
  scale_color_hue(labels = clusters) + NoLegend()

# Fig. 4C
DotPlot(tcon_seurat, features = c("CCR7", "FOXP3", "IKZF2", "PCNA", "CHI3L2",
                                  "KLRG1", "IFNG", "PDCD1", "GZMK", "KCNE3"),
        dot.scale = 10,
        dot.min = 0.1,
        col.min = 0.1) +
  theme(text = element_text(size = 30),
        axis.text.x = element_text(size = 18, face = "italic", angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 18)) +
  labs(x = "Gene", y = "Cluster", color = "Normalized Expression") +
  scale_color_gradientn(colours = rev(brewer.pal(name = "RdYlBu", n = 10)),
                        limits = c(0, 3), oob = scales::squish, name = "log2 (count + 1)")

# Fig. S1
FeaturePlot(tcon_seurat, features = c("IL2RA", "FOXP3", "IKZF2"),
            label = TRUE, raster = FALSE, pt.size = 1) +
  scale_color_viridis_c()

# Fig. 4D

# export cell number in each cluster by sample
md2 <- tcon_seurat@meta.data %>% 
  as.data.table

md2 <- md2[, .N, by = c("sample.name", "integrated_snn_res.1.2")] %>%
  dcast(., sample.name ~ integrated_snn_res.1.2, value.var = "N")

# covert NA to 0
md2[is.na(md2)] <- 0

# move sample name to rownames
md2_raw <- md2 %>%
  column_to_rownames("sample.name")

md2_pro <- md2 %>% 
  column_to_rownames("sample.name")

# calculate frequency
for (row in 1:nrow(md2_pro)){
  md2_pro[row,] <- round((md2_pro[row,]/sum(md2_pro[row,])), 4)
}

# Write to clipboard for easy copy to excel

# raw cell number
colnames(md2_raw) <- levels(Idents(tcon_seurat))
rownames(md2_raw) <- c("A", "B", "C", "D", "E", "F", "G", "I", "J", "K", "L")
write_clip(md2_raw)

# cell proportion
colnames(md2_pro) <- clusters
rownames(md2_pro) <- c("A", "B", "C", "D", "E", "F", "G", "I", "J", "K", "L")
write_clip(md2_pro)

#### Tcon VDJ ####

# Take only cells that were kept in Tcon scRNA-seq analysis
tcon_cell_id <- which(rownames(vdj) %in% rownames(tcon_seurat@meta.data))
filtered_tcr <- dplyr::slice(vdj, tcon_cell_id)

# Remove cells that have both TCR and BCR
filtered_tcr <- filtered_tcr %>% #336 cells removed
  rownames_to_column("cell") %>%
  filter(TCR_Paired_Chains == "FALSE" | 
           BCR_Paired_Chains == "FALSE") %>%
  column_to_rownames("cell")

#### Wrangle Tcon VDJ data for scRepertoire ####

# Remove B cells and modify format to match expected columns
filtered_rep_t <- filtered_tcr %>% # cells
  filter(TCR_Paired_Chains == "TRUE" |
           (filtered_tcr$TCR_Alpha_Gamma_Read_Count > 0 &
              filtered_tcr$TCR_Beta_Delta_Read_Count == 0 &
              filtered_tcr$BCR_Heavy_Read_Count == 0 &
              filtered_tcr$BCR_Light_Read_Count == 0) |
           (filtered_tcr$TCR_Alpha_Gamma_Read_Count == 0 &
              filtered_tcr$TCR_Beta_Delta_Read_Count > 0 &
              filtered_tcr$BCR_Heavy_Read_Count == 0 &
              filtered_tcr$BCR_Light_Read_Count == 0)) %>%
  select(starts_with("TCR"), "Sample_Name") %>%
  add_column("barcode" = "none", 
             "is_cell" = TRUE, 
             "contig_id" = "none", 
             "high_confidence" = "TRUE", 
             "length" = "None", 
             "chain" = "none", 
             "v_gene" = "none", 
             "d_gene" = "None", 
             "j_gene" = "none", 
             "c_gene" = "none", 
             "full_length" = "TRUE", 
             "productive" = "TRUE", 
             "cdr3" = "none", 
             "cdr3_nt" = "none", 
             "reads" = "none", 
             "umis" = "None", 
             "raw_clonotype_id" = "None", 
             "raw_consensus_id" = "None")

filtered_rep_t$barcode <- rownames(filtered_rep_t)

tcr_type <- substr(filtered_rep_t$TCR_Alpha_Gamma_V_gene_Dominant, start = 1, stop = 3)
filtered_rep_t$chain <- tcr_type

tcr_agv <- as.data.frame(do.call(rbind, str_split(filtered_rep_t[, "TCR_Alpha_Gamma_V_gene_Dominant"], "\\*")))
filtered_rep_t$v_gene <- tcr_agv$V1

tcr_agj <- as.data.frame(do.call(rbind, str_split(filtered_rep_t[, "TCR_Alpha_Gamma_J_gene_Dominant"], "\\*")))
filtered_rep_t$j_gene <- tcr_agj$V1 

tcr_agc <- as.data.frame(do.call(rbind, str_split(filtered_rep_t[, "TCR_Alpha_Gamma_C_gene_Dominant"], "\\*")))
filtered_rep_t$c_gene <- tcr_agc$V1

filtered_rep_t$cdr3 <- filtered_rep_t[, "TCR_Alpha_Gamma_CDR3_Translation_Dominant"]

filtered_rep_t$cdr3_nt <- filtered_rep_t[, "TCR_Alpha_Gamma_CDR3_Nucleotide_Dominant"]

# Add id number for sorting later, duplicate rows to have one chain per row, overwrite second chain data
filtered_rep_t_dup <- filtered_rep_t
filtered_rep_t_dup <- add_column(filtered_rep_t_dup, id = 1:length(rownames(filtered_rep_t_dup)))
filtered_rep_t_dup <- bind_rows(filtered_rep_t_dup, filtered_rep_t_dup)

tcr_type2 <- substr(filtered_rep_t$TCR_Beta_Delta_V_gene_Dominant, start = 1, stop = 3)
filtered_rep_t_dup$chain[(nrow(filtered_rep_t)+1):nrow(filtered_rep_t_dup)] <- tcr_type2

tcr_bdv <- as.data.frame(do.call(rbind, str_split(filtered_rep_t_dup[(nrow(filtered_rep_t)+1):nrow(filtered_rep_t_dup), 
                                                                     "TCR_Beta_Delta_V_gene_Dominant"], "\\*")))
filtered_rep_t_dup$v_gene[(nrow(filtered_rep_t)+1):nrow(filtered_rep_t_dup)] <- tcr_bdv$V1

tcr_bdd <- as.data.frame(do.call(rbind, str_split(filtered_rep_t_dup[(nrow(filtered_rep_t)+1):nrow(filtered_rep_t_dup), 
                                                                     "TCR_Beta_Delta_D_gene_Dominant"], "\\*")))
filtered_rep_t_dup$d_gene[(nrow(filtered_rep_t)+1):nrow(filtered_rep_t_dup)] <- tcr_bdd$V1

tcr_bdj <- as.data.frame(do.call(rbind, str_split(filtered_rep_t_dup[(nrow(filtered_rep_t)+1):nrow(filtered_rep_t_dup), 
                                                                     "TCR_Beta_Delta_J_gene_Dominant"], "\\*")))
filtered_rep_t_dup$j_gene[(nrow(filtered_rep_t)+1):nrow(filtered_rep_t_dup)] <- tcr_bdj$V1

tcr_bdc <- as.data.frame(do.call(rbind, str_split(filtered_rep_t_dup[(nrow(filtered_rep_t)+1):nrow(filtered_rep_t_dup), 
                                                                     "TCR_Beta_Delta_C_gene_Dominant"], "\\*")))
filtered_rep_t_dup$c_gene[(nrow(filtered_rep_t)+1):nrow(filtered_rep_t_dup)] <- tcr_bdc$V1

filtered_rep_t_dup$cdr3[(nrow(filtered_rep_t)+1):nrow(filtered_rep_t_dup)] <- 
  filtered_rep_t_dup[(nrow(filtered_rep_t)+1):nrow(filtered_rep_t_dup),
                     "TCR_Beta_Delta_CDR3_Translation_Dominant"]

filtered_rep_t_dup$cdr3_nt[(nrow(filtered_rep_t)+1):nrow(filtered_rep_t_dup)] <- 
  filtered_rep_t_dup[(nrow(filtered_rep_t)+1):nrow(filtered_rep_t_dup),
                     "TCR_Beta_Delta_CDR3_Nucleotide_Dominant"]

# Re-order the cells so that contigs are next to each other, assign contig ids, and collect wanted columns
filtered_rep_t_ord <- filtered_rep_t_dup[order(filtered_rep_t_dup$id),]
filtered_rep_t_sub <- filtered_rep_t_ord[, 17:36]

#### Import Tcon TCR data into scRepertoire ####

# Alpha Beta T cells
t1 <- filter(filtered_rep_t_sub, Sample_Name == "vdj_183a" & (chain == "TRA" | chain == "TRB"))
t2 <- filter(filtered_rep_t_sub, Sample_Name == "vdj_183b" & (chain == "TRA" | chain == "TRB"))
t3 <- filter(filtered_rep_t_sub, Sample_Name == "vdj_183c" & (chain == "TRA" | chain == "TRB"))
t4 <- filter(filtered_rep_t_sub, Sample_Name == "vdj_183d" & (chain == "TRA" | chain == "TRB"))
t5 <- filter(filtered_rep_t_sub, Sample_Name == "vdj_183e" & (chain == "TRA" | chain == "TRB"))
t6 <- filter(filtered_rep_t_sub, Sample_Name == "vdj_183f" & (chain == "TRA" | chain == "TRB"))
t7 <- filter(filtered_rep_t_sub, Sample_Name == "vdj_183g" & (chain == "TRA" | chain == "TRB"))
t8 <- filter(filtered_rep_t_sub, Sample_Name == "vdj_183i" & (chain == "TRA" | chain == "TRB"))
t9 <- filter(filtered_rep_t_sub, Sample_Name == "vdj_183j" & (chain == "TRA" | chain == "TRB"))
t10 <- filter(filtered_rep_t_sub, Sample_Name == "vdj_183k" & (chain == "TRA" | chain == "TRB"))
t11 <- filter(filtered_rep_t_sub, Sample_Name == "vdj_183l" & (chain == "TRA" | chain == "TRB"))

t_list <- list(t1[2:19], t2[2:19], t3[2:19], t4[2:19],
               t5[2:19], t6[2:19], t7[2:19], t8[2:19],
               t9[2:19], t10[2:19], t11[2:19])

# Combine the contigs for T cells
ab_tcr <- combineTCR(t_list,
                     samples = c("a", "b", "c", "d",
                                 "e", "f", "g", "i",
                                 "j", "k", "l"),
                     removeNA = TRUE)

# Add treatment group
ab_tcr <- addVariable(ab_tcr, variable.name = "group",
                      variables= c("PBSC", "OrcaT", "PBSC", "PBSC",
                                   "OrcaT", "PBSC", "OrcaT", "OrcaT",
                                   "PTC", "PTC", "PTC"))

# Fig. 4E
clonalDiversity(ab_tcr, cloneCall = "aa", group.by = "sample", x.axis = "group", exportTable = TRUE)

#### Combine Tcon VDJ with Tcon seurat ####

# Make row names match
rnames <- rnames[tcon_cell_id]
tcon_seurat <- RenameCells(tcon_seurat, new.names = rnames)

tcon_seurat <- combineExpression(ab_tcr, tcon_seurat,
                                 cloneCall = "aa",
                                 group.by = "sample",
                                 cloneSize = c(Rare = 0.001, Small = 0.005, Medium = 0.01,
                                               Large = 0.05, Hyperexpanded = 1))


# Set color palette
colorblind_vector <- colorRampPalette(c("firebrick2", "darkorange", 
                                        "dodgerblue", "purple", "forestgreen"))

# Change cloneType to factors
tcon_seurat@meta.data$cloneSize <- factor(tcon_seurat@meta.data$cloneSize, 
                                          levels = c("Hyperexpanded (0.05 < X <= 1)",
                                                     "Large (0.01 < X <= 0.05)", 
                                                     "Medium (0.005 < X <= 0.01)", 
                                                     "Small (0.001 < X <= 0.005)", 
                                                     "Rare (0 < X <= 0.001)"))
# Fig. 4F
clonalOverlay(subset(tcon_seurat, (group == "PBSC")),
              reduction = "umap",
              freq.cutpoint = 2,
              bins = 10,
              facet = "group") +
  guides(color = "none") +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        axis.title = element_text(size = 20))

clonalOverlay(subset(tcon_seurat, (group == "OrcaT")),
              reduction = "umap",
              freq.cutpoint = 2,
              bins = 10,
              facet = "group") +
  guides(color = "none") +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        axis.title = element_text(size = 20))

tcon_seurat@meta.data[(which(tcon_seurat@meta.data$Frequency > 0.0005)),] %>% 
  group_by(group) %>% 
  count(integrated_snn_res.0.4) %>%
  print(n=42)

# Frequency of clones by sample in each cluster (only counting clones with clonalFreq >2)
clone_freq <- tcon_seurat@meta.data[(which(tcon_seurat@meta.data$clonalFrequency >= 2)),] %>% 
  group_by(sample.name, group, integrated_snn_res.1.2) %>% 
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(group)

# Assign sample names
clone_freq$sample.name <- factor(clone_freq$sample.name,
                                 levels = c("B", "E", "G", "I",
                                            "J", "K", "L",
                                            "A", "C", "D", "F"))

# Export for easy pasting into excel
write_clip(pivot_wider(clone_freq, id_cols = "sample.name", names_from = "integrated_snn_res.1.2", values_from = c("freq")) %>%
             replace(is.na(.), 0))



#### Subset Treg sample.tags ####

treg_seurat <- subset(integrated_rna, cell.type == "Treg") #73563 cells

# Re-cluster subsetted T cells
treg_seurat <- RunPCA(object = treg_seurat)
ElbowPlot(treg_seurat, ndims = 50) # 30 dims

#### Treg Clustering ####

# UMAP
treg_seurat <- RunUMAP(treg_seurat,
                       dims = 1:30,
                       reduction = "pca")

# Find clusters
treg_seurat <- FindNeighbors(object = treg_seurat,
                             dims = 1:30)

treg_seurat <- FindClusters(object = treg_seurat,
                            resolution = 0.4)

# Plot UMAP
Idents(object = treg_seurat) <- "integrated_snn_res.0.4"
DimPlot(treg_seurat,
        label = TRUE,
        group.by = "ident",
        raster = FALSE) + NoLegend()

DefaultAssay(treg_seurat) <- "SCT"
treg_seurat <- JoinLayers(treg_seurat)
treg_seurat <- PrepSCTFindMarkers(treg_seurat, assay = "SCT", verbose = TRUE)
rna_markers <- FindAllMarkers(object = treg_seurat,
                              only.pos = TRUE,
                              logfc.threshold = 0.25,
                              assay = "SCT")

view(rna_markers)

# Remove myeloid and B cells (first time)

treg_seurat <- subset(treg_seurat, idents = c("2", "16"), invert = TRUE) # removed 8287 cells

# Re-cluster remaining cells

# UMAP
DefaultAssay(treg_seurat) <- "integrated"

treg_seurat <- RunUMAP(treg_seurat,
                       dims = 1:30,
                       reduction = "pca")

# Find clusters
treg_seurat <- FindNeighbors(object = treg_seurat,
                             dims = 1:30)

treg_seurat <- FindClusters(object = treg_seurat,
                            resolution = 0.4)

Idents(object = treg_seurat) <- "integrated_snn_res.0.4"

DefaultAssay(treg_seurat) <- "SCT"
treg_seurat <- PrepSCTFindMarkers(treg_seurat, assay = "SCT", verbose = TRUE)
rna_markers <- FindAllMarkers(object = treg_seurat,
                              only.pos = TRUE,
                              logfc.threshold = 0.25,
                              assay = "SCT")

view(rna_markers)

# Remove myeloid and B cells again (second time)

treg_seurat <- subset(treg_seurat, idents = c("1"), invert = TRUE) # removed 10168 cells

# Re-cluster remaining cells

# UMAP
DefaultAssay(treg_seurat) <- "integrated"

treg_seurat <- RunUMAP(treg_seurat,
                       dims = 1:30,
                       reduction = "pca")

# Find clusters
treg_seurat <- FindNeighbors(object = treg_seurat,
                             dims = 1:30)

treg_seurat <- FindClusters(object = treg_seurat,
                            resolution = 0.4)

# Resolution = 0.4
Idents(object = treg_seurat) <- "integrated_snn_res.0.4"

DimPlot(treg_seurat,
        label = TRUE, raster= FALSE) + NoLegend()

DefaultAssay(treg_seurat) <- "SCT"
treg_seurat <- PrepSCTFindMarkers(treg_seurat, assay = "SCT", verbose = TRUE)
rna_markers <- FindAllMarkers(object = treg_seurat,
                              only.pos = TRUE,
                              logfc.threshold = 0.25,
                              assay = "SCT")

view(rna_markers)

#### Treg Annotation ####

# Fig. S2A
DimPlot(treg_seurat,
        label = TRUE, raster = FALSE, label.size = 8) +
  scale_color_hue(labels = clusters) + NoLegend()

# Fig. S2B
DotPlot(treg_seurat, features = c("HLA-DRA", "IKZF2", "PCNA", "F5", "CCR4", "NCR3", #1000x1000
                                  "LRRC32", "TIGIT", "NINJ2", "ENTPD1", "PMCH", "LAIR2", "CXCR6", "LAG3",
                                  "IL12RB2", "TNFRSF4", "RORC", "CCR3", "CCR8", "ARG1", "CCR9", "CCR10"),
        dot.scale = 10,
        dot.min = 0.1,
        col.min = 0.1) +
  theme(text = element_text(size = 30),
        axis.text.x = element_text(size = 18, face = "italic", angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 18)) +
  labs(x = "Gene", y = "Cluster", color = "Normalized Expression") +
  scale_color_gradientn(colours = rev(brewer.pal(name = "RdYlBu", n = 10)),
                        limits = c(0, 3), oob = scales::squish, name = "log2 (count + 1)")

# Fig. S2C

# export cell number in each cluster by sample
md2 <- treg_seurat@meta.data %>% 
  as.data.table
md2 <- md2[, .N, by = c("sample.name", "integrated_snn_res.0.4")] %>% 
  dcast(., sample.name ~ integrated_snn_res.0.4, value.var = "N")

# covert NA to 0
md2[is.na(md2)] <- 0

# move sample name to rownames
md2_raw <- md2 %>%
  column_to_rownames("sample.name")

md2_pro <- md2 %>% 
  column_to_rownames("sample.name")

# calculate frequency
for (row in 1:nrow(md2_pro)){
  md2_pro[row,] <- round((md2_pro[row,]/sum(md2_pro[row,])), 4)
}

# Write to clipboard for easy copy to excel

# raw cell number
colnames(md2_raw) <- levels(Idents(tcon_seurat))
rownames(md2_raw) <- c("A", "B", "C", "D", "E", "F", "G", "I", "J", "K", "L")
write_clip(md2_raw)

# cell proportion
colnames(md2_pro) <- clusters
rownames(md2_pro) <- c("A", "B", "C", "D", "E", "F", "G", "I", "J", "K", "L")
write_clip(md2_pro)

#### Treg VDJ ####

# Take only cells that were kept in scRNA-seq analysis
treg_cell_id <- which(rownames(t_tcr) %in% rownames(treg_seurat@meta.data))
filtered_tcr <- dplyr::slice(t_tcr, treg_cell_id)

# Remove cells that have both TCR and BCR
filtered_tcr <- filtered_tcr %>% #336 cells removed
  rownames_to_column("cell") %>%
  filter(TCR_Paired_Chains == "FALSE" | 
           BCR_Paired_Chains == "FALSE") %>%
  column_to_rownames("cell")

#### Wrangle Treg VDJ data for scRepertoire ####

# T cells
filtered_rep_t <- filtered_tcr %>% # cells
  filter(TCR_Paired_Chains == "TRUE" |
           (filtered_tcr$TCR_Alpha_Gamma_Read_Count > 0 &
              filtered_tcr$TCR_Beta_Delta_Read_Count == 0 &
              filtered_tcr$BCR_Heavy_Read_Count == 0 &
              filtered_tcr$BCR_Light_Read_Count == 0) |
           (filtered_tcr$TCR_Alpha_Gamma_Read_Count == 0 &
              filtered_tcr$TCR_Beta_Delta_Read_Count > 0 &
              filtered_tcr$BCR_Heavy_Read_Count == 0 &
              filtered_tcr$BCR_Light_Read_Count == 0)) %>%
  select(starts_with("TCR"), "Sample_Name") %>%
  add_column("barcode" = "none", 
             "is_cell" = TRUE, 
             "contig_id" = "none", 
             "high_confidence" = "TRUE", 
             "length" = "None", 
             "chain" = "none", 
             "v_gene" = "none", 
             "d_gene" = "None", 
             "j_gene" = "none", 
             "c_gene" = "none", 
             "full_length" = "TRUE", 
             "productive" = "TRUE", 
             "cdr3" = "none", 
             "cdr3_nt" = "none", 
             "reads" = "none", 
             "umis" = "None", 
             "raw_clonotype_id" = "None", 
             "raw_consensus_id" = "None")

filtered_rep_t$barcode <- rownames(filtered_rep_t)

tcr_type <- substr(filtered_rep_t$TCR_Alpha_Gamma_V_gene_Dominant, start = 1, stop = 3)
filtered_rep_t$chain <- tcr_type

tcr_agv <- as.data.frame(do.call(rbind, str_split(filtered_rep_t[, "TCR_Alpha_Gamma_V_gene_Dominant"], "\\*")))
filtered_rep_t$v_gene <- tcr_agv$V1

tcr_agj <- as.data.frame(do.call(rbind, str_split(filtered_rep_t[, "TCR_Alpha_Gamma_J_gene_Dominant"], "\\*")))
filtered_rep_t$j_gene <- tcr_agj$V1 

tcr_agc <- as.data.frame(do.call(rbind, str_split(filtered_rep_t[, "TCR_Alpha_Gamma_C_gene_Dominant"], "\\*")))
filtered_rep_t$c_gene <- tcr_agc$V1

filtered_rep_t$cdr3 <- filtered_rep_t[, "TCR_Alpha_Gamma_CDR3_Translation_Dominant"]

filtered_rep_t$cdr3_nt <- filtered_rep_t[, "TCR_Alpha_Gamma_CDR3_Nucleotide_Dominant"]

# Add id number for sorting later, duplicate rows to have one chain per row, overwrite second chain data
filtered_rep_t_dup <- filtered_rep_t
filtered_rep_t_dup <- add_column(filtered_rep_t_dup, id = 1:length(rownames(filtered_rep_t_dup)))
filtered_rep_t_dup <- bind_rows(filtered_rep_t_dup, filtered_rep_t_dup)

tcr_type2 <- substr(filtered_rep_t$TCR_Beta_Delta_V_gene_Dominant, start = 1, stop = 3)
filtered_rep_t_dup$chain[(nrow(filtered_rep_t)+1):nrow(filtered_rep_t_dup)] <- tcr_type2

tcr_bdv <- as.data.frame(do.call(rbind, str_split(filtered_rep_t_dup[(nrow(filtered_rep_t)+1):nrow(filtered_rep_t_dup), 
                                                                     "TCR_Beta_Delta_V_gene_Dominant"], "\\*")))
filtered_rep_t_dup$v_gene[(nrow(filtered_rep_t)+1):nrow(filtered_rep_t_dup)] <- tcr_bdv$V1

tcr_bdd <- as.data.frame(do.call(rbind, str_split(filtered_rep_t_dup[(nrow(filtered_rep_t)+1):nrow(filtered_rep_t_dup), 
                                                                     "TCR_Beta_Delta_D_gene_Dominant"], "\\*")))
filtered_rep_t_dup$d_gene[(nrow(filtered_rep_t)+1):nrow(filtered_rep_t_dup)] <- tcr_bdd$V1

tcr_bdj <- as.data.frame(do.call(rbind, str_split(filtered_rep_t_dup[(nrow(filtered_rep_t)+1):nrow(filtered_rep_t_dup), 
                                                                     "TCR_Beta_Delta_J_gene_Dominant"], "\\*")))
filtered_rep_t_dup$j_gene[(nrow(filtered_rep_t)+1):nrow(filtered_rep_t_dup)] <- tcr_bdj$V1

tcr_bdc <- as.data.frame(do.call(rbind, str_split(filtered_rep_t_dup[(nrow(filtered_rep_t)+1):nrow(filtered_rep_t_dup), 
                                                                     "TCR_Beta_Delta_C_gene_Dominant"], "\\*")))
filtered_rep_t_dup$c_gene[(nrow(filtered_rep_t)+1):nrow(filtered_rep_t_dup)] <- tcr_bdc$V1

filtered_rep_t_dup$cdr3[(nrow(filtered_rep_t)+1):nrow(filtered_rep_t_dup)] <- 
  filtered_rep_t_dup[(nrow(filtered_rep_t)+1):nrow(filtered_rep_t_dup),
                     "TCR_Beta_Delta_CDR3_Translation_Dominant"]

filtered_rep_t_dup$cdr3_nt[(nrow(filtered_rep_t)+1):nrow(filtered_rep_t_dup)] <- 
  filtered_rep_t_dup[(nrow(filtered_rep_t)+1):nrow(filtered_rep_t_dup),
                     "TCR_Beta_Delta_CDR3_Nucleotide_Dominant"]

# Re-order the cells so that contigs are next to each other, assign contig ids, and collect wanted columns
filtered_rep_t_ord <- filtered_rep_t_dup[order(filtered_rep_t_dup$id),]
filtered_rep_t_sub <- filtered_rep_t_ord[, 17:36]

#### Import Treg TCR data into scRepertoire ####

# Alpha Beta T cells
t1 <- filter(filtered_rep_t_sub, Sample_Name == "vdj_183a" & (chain == "TRA" | chain == "TRB"))
t2 <- filter(filtered_rep_t_sub, Sample_Name == "vdj_183b" & (chain == "TRA" | chain == "TRB"))
t3 <- filter(filtered_rep_t_sub, Sample_Name == "vdj_183c" & (chain == "TRA" | chain == "TRB"))
t4 <- filter(filtered_rep_t_sub, Sample_Name == "vdj_183d" & (chain == "TRA" | chain == "TRB"))
t5 <- filter(filtered_rep_t_sub, Sample_Name == "vdj_183e" & (chain == "TRA" | chain == "TRB"))
t6 <- filter(filtered_rep_t_sub, Sample_Name == "vdj_183f" & (chain == "TRA" | chain == "TRB"))
t7 <- filter(filtered_rep_t_sub, Sample_Name == "vdj_183g" & (chain == "TRA" | chain == "TRB"))
t8 <- filter(filtered_rep_t_sub, Sample_Name == "vdj_183i" & (chain == "TRA" | chain == "TRB"))
t9 <- filter(filtered_rep_t_sub, Sample_Name == "vdj_183j" & (chain == "TRA" | chain == "TRB"))
t10 <- filter(filtered_rep_t_sub, Sample_Name == "vdj_183k" & (chain == "TRA" | chain == "TRB"))
t11 <- filter(filtered_rep_t_sub, Sample_Name == "vdj_183l" & (chain == "TRA" | chain == "TRB"))

t_list <- list(t1[2:19], t2[2:19], t3[2:19], t4[2:19],
               t5[2:19], t6[2:19], t7[2:19], t8[2:19],
               t9[2:19], t10[2:19], t11[2:19])

# Combine the contigs for Tregs
ab_tcr <- combineTCR(t_list,
                     samples = c("a", "b", "c", "d",
                                 "e", "f", "g", "i",
                                 "j", "k", "l"),
                     removeNA = TRUE)

# Add treatment group variable
ab_tcr <- addVariable(ab_tcr, variable.name = "group",
                      variables= c("PBSC", "OrcaT", "PBSC", "PBSC",
                                   "OrcaT", "PBSC", "OrcaT", "OrcaT",
                                   "PTC", "PTC", "PTC"))

# Fig. S2D
clonalDiversity(ab_tcr, cloneCall = "aa", group.by = "sample", x.axis = "group", exportTable = TRUE)

#### Combine Treg VDJ with Treg seurat ####

# Make rownames match
rnames <- rnames[treg_cell_id]
treg_seurat <- RenameCells(treg_seurat, new.names = rnames)

treg_seurat <- combineExpression(ab_tcr, treg_seurat,
                                 cloneCall = "aa",
                                 group.by = "sample",
                                 cloneSize = c(Rare = 0.001, Small = 0.005, Medium = 0.01,
                                               Large = 0.05, Hyperexpanded = 1))
# Fig. S2E
clonalOverlay(subset(treg_seurat, (group == "SOC")),
              reduction = "umap",
              freq.cutpoint = 2,
              bins = 10,
              facet = "group") +
  guides(color = "none") +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        axis.title = element_text(size = 20))

clonalOverlay(subset(treg_seurat, (group == "OrcaT")),
              reduction = "umap",
              freq.cutpoint = 2,
              bins = 10,
              facet = "group") +
  guides(color = "none") +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        axis.title = element_text(size = 20))

# Plot frequency of clones by sample in each cluster (only counting clones with clonefreq > 2)
clone_freq <- treg_seurat@meta.data[(which(treg_seurat@meta.data$clonalFrequency >= 2)),] %>% 
  group_by(sample.name, group, integrated_snn_res.0.4) %>% 
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(group)

# Add sample names
clone_freq$sample.name <- factor(clone_freq$sample.name,
                                 levels = c("B", "E", "G", "I",
                                            "J", "K", "L",
                                            "A", "C", "D", "F"))
# Export for easy pasting into excel
write_clip(pivot_wider(clone_freq, id_cols = "sample.name", names_from = "integrated_snn_res.0.4", values_from = c("n", "freq")))

sessionInfo()

# R version 4.4.1 (2024-06-14 ucrt)
# Platform: x86_64-w64-mingw32/x64
# Running under: Windows 11 x64 (build 26100)
# 
# Matrix products: default
# 
# 
# locale:
#   [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8
# [4] LC_NUMERIC=C                           LC_TIME=English_United States.utf8    
# 
# time zone: America/Los_Angeles
# tzcode source: internal
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] RColorBrewer_1.1-3 clipr_0.8.0        scRepertoire_2.0.0 data.table_1.16.0  BPCells_0.2.0      Seurat_5.1.0       SeuratObject_5.0.2
# [8] sp_2.1-4           lubridate_1.9.3    forcats_1.0.0      stringr_1.5.1      dplyr_1.1.4        purrr_1.0.2        readr_2.1.5       
# [15] tidyr_1.3.1        tibble_3.2.1       ggplot2_3.5.1      tidyverse_2.0.0   
# 
# loaded via a namespace (and not attached):
#   [1] cubature_2.1.1              RcppAnnoy_0.0.22            splines_4.4.1               later_1.3.2                 polyclip_1.10-7            
# [6] fastDummies_1.7.4           lifecycle_1.0.4             globals_0.16.3              lattice_0.22-6              MASS_7.3-60.2              
# [11] magrittr_2.0.3              plotly_4.10.4               httpuv_1.6.15               sctransform_0.4.1           spam_2.10-0                
# [16] spatstat.sparse_3.1-0       reticulate_1.38.0           cowplot_1.1.3               pbapply_1.7-2               abind_1.4-5                
# [21] zlibbioc_1.50.0             Rtsne_0.17                  GenomicRanges_1.56.1        ggraph_2.2.1                BiocGenerics_0.50.0        
# [26] tweenr_2.0.3                evmix_2.12                  GenomeInfoDbData_1.2.12     IRanges_2.38.1              S4Vectors_0.42.1           
# [31] ggrepel_0.9.5               irlba_2.3.5.1               listenv_0.9.1               spatstat.utils_3.1-0        iNEXT_3.0.1                
# [36] MatrixModels_0.5-3          goftest_1.2-3               RSpectra_0.16-2             spatstat.random_3.3-1       fitdistrplus_1.2-1         
# [41] parallelly_1.38.0           leiden_0.4.3.1              codetools_0.2-20            DelayedArray_0.30.1         ggforce_0.4.2              
# [46] tidyselect_1.2.1            UCSC.utils_1.0.0            farver_2.1.2                viridis_0.6.5               matrixStats_1.3.0          
# [51] stats4_4.4.1                spatstat.explore_3.3-2      jsonlite_1.8.8              tidygraph_1.3.1             progressr_0.14.0           
# [56] ggridges_0.5.6              ggalluvial_0.12.5           survival_3.8-3              tools_4.4.1                 stringdist_0.9.12          
# [61] ica_1.0-3                   Rcpp_1.0.13                 glue_1.8.0                  gridExtra_2.3               SparseArray_1.4.8          
# [66] MatrixGenerics_1.16.0       GenomeInfoDb_1.40.1         withr_3.0.1                 fastmap_1.2.0               fansi_1.0.6                
# [71] SparseM_1.84-2              digest_0.6.37               timechange_0.3.0            R6_2.5.1                    mime_0.12                  
# [76] colorspace_2.1-1            scattermore_1.2             tensor_1.5                  spatstat.data_3.1-2         utf8_1.2.4                 
# [81] generics_0.1.3              graphlayouts_1.1.1          httr_1.4.7                  htmlwidgets_1.6.4           S4Arrays_1.4.1             
# [86] uwot_0.2.2                  pkgconfig_2.0.3             gtable_0.3.5                lmtest_0.9-40               SingleCellExperiment_1.26.0
# [91] XVector_0.44.0              htmltools_0.5.8.1           dotCall64_1.1-1             scales_1.3.0                Biobase_2.64.0             
# [96] png_0.1-8                   spatstat.univar_3.0-0       ggdendro_0.2.0              rstudioapi_0.16.0           tzdb_0.4.0                 
# [101] reshape2_1.4.4              rjson_0.2.22                nlme_3.1-164                zoo_1.8-12                  cachem_1.1.0               
# [106] KernSmooth_2.23-24          parallel_4.4.1              miniUI_0.1.1.1              pillar_1.9.0                grid_4.4.1                 
# [111] vctrs_0.6.5                 RANN_2.6.2                  VGAM_1.1-11                 promises_1.3.0              xtable_1.8-4               
# [116] cluster_2.1.6               truncdist_1.0-2             cli_3.6.3                   compiler_4.4.1              rlang_1.1.4                
# [121] crayon_1.5.3                future.apply_1.11.2         plyr_1.8.9                  stringi_1.8.4               viridisLite_0.4.2          
# [126] deldir_2.0-4                munsell_0.5.1               gsl_2.1-8                   lazyeval_0.2.2              spatstat.geom_3.3-2        
# [131] quantreg_5.98               Matrix_1.7-0                RcppHNSW_0.6.0              hms_1.1.3                   patchwork_1.2.0            
# [136] future_1.34.0               shiny_1.9.1                 SummarizedExperiment_1.34.0 evd_2.3-7                   ROCR_1.0-11                
# [141] igraph_2.0.3                memoise_2.0.1   