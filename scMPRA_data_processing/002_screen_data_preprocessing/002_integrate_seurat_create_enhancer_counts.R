#In this script we merge multiple cellranger outputs/10x reactions and plot qc measurements of the different runs

library(Seurat)
library(tidyverse)
library(viridis)
library(Matrix)
library(data.table)
library(actuar)
library(ComplexHeatmap)

set.seed(456)

setwd("path_to_script")

assoc <- read_csv("./example_feature_ref.csv")
assoc <- assoc %>% as.data.frame %>% 
  select(CRS, guide_id, gfp_id, reads)

#read in TAP-Seq genes: 
tap_genes <- read.csv("~/cluster/project/SCG4SYN/Experiments/220923_5TAP_panel_scMPRA/HSC_panel/HSC_panel_selected_genes.csv") %>% unlist() %>% as.vector()

#these are the TAP genes that are not captured in the data (so i exclude them)
tap_not_there <- c("BC051077", "Ebf1", "Klf4", "Ogn", "Spdef", "Tfap2a", "Vpreb3", "Prom1")
tap_measured <- subset(tap_genes, !(tap_genes %in% tap_not_there))

#read in data from all samples 
sam7.data <- Read10X(data.dir = "./out/scMPRA_sam07_filtered_bcs_test/outs/filtered_feature_bc_matrix/")
sam7 <- CreateSeuratObject(counts =  sam7.data$`Gene Expression`, project = "scMPRA007")

sam8.data <- Read10X(data.dir = "./out/scMPRA_sam08_filtered_bcs_test/outs/filtered_feature_bc_matrix/")
sam8 <- CreateSeuratObject(counts =  sam8.data$`Gene Expression`, project = "scMPRA008")

# create a new seurat object for each of the samples
sam7.custom <- CreateSeuratObject(counts = sam7.data$`Custom`)
sam8.custom <- CreateSeuratObject(counts = sam8.data$`Custom`)

#create 3 different seurat objects for the 3 different library types
gfp.assay <- subset(sam7.custom, features = grep("reporter", rownames(sam7.custom)))
guide.assay <- subset(sam7.custom, features = grep("guide", rownames(sam7.custom)))
guide.assay <- guide.assay[rowSums(guide.assay, slot = "counts") > 0,]

#i think this takes ages, so if we dont need it lets not do it 
sam7[["guide"]] <- guide.assay@assays$RNA
sam7[["gfp"]] <- gfp.assay@assays$RNA

#create 3 different seurat objects for the 3 different library types
gfp.assay <- subset(sam8.custom, features = grep("reporter", rownames(sam8.custom)))
guide.assay <- subset(sam8.custom, features = grep("guide", rownames(sam8.custom)))
guide.assay <- guide.assay[rowSums(guide.assay, slot = "counts") > 0,]

#i think this takes ages, so if we dont need it lets not do it 
sam8[["guide"]] <- guide.assay@assays$RNA
sam8[["gfp"]] <- gfp.assay@assays$RNA
#______________________merged all the data together in separate seurat objects 

#save the merged seurat object 
#Combine the seurat objexts 
#
#Combine the seurat objexts 
gex.combined <- merge(sam7, y = c(sam8), add.cell.ids = c("007", "008"), project = "scMPRA")

gfp_hist <- ggplot(qc_combined@meta.data, aes(x = nFeature_gfp)) +
  geom_histogram(fill = "#608c83", binwidth  = 0.5) + 
  ggtitle("Number of GFP features per Cell") + 
  xlab("GFP Features per Cell") + 
  ylab("Counts") + 
  xlim(c(-1,10)) + 
  theme_bw() + 
  facet_wrap(~ orig.ident)


guide_hist <- ggplot(qc_combined@meta.data, aes(x = nFeature_guide)) + 
  geom_histogram(fill = "#513c59", binwidth = 0.5) + 
  ggtitle("Number of Guide features per Cell") +
  xlab("Guide Features per Cell") + 
  ylab("Counts") + 
  xlim(c(-1,10)) + 
  theme_bw() + 
  facet_wrap(~ orig.ident)

gfp_hist + guide_hist


#______________________Lets do some QC________________________________________

#How many cells do we have before qc? 
table(gex.combined$orig.ident)

VlnPlot(gex.combined, features = c("nFeature_RNA", "nCount_RNA"), pt.size = 0.1, ncol = 3)  

#QC in numbers 
qc_measures_sum <- gex.combined@meta.data %>% 
  dplyr::group_by(orig.ident) %>% 
  dplyr::summarise(min_nCountRNA = min(nCount_RNA), 
                   max_nCountRNA = max(nCount_RNA), 
                   mean_nCountRNA = mean(nCount_RNA), 
                   min_nFeatures = min(nFeature_RNA),
                   max_nFeatures = max(nFeature_RNA), 
                   mean_nFeatures = mean(nFeature_RNA),
                   mean_nFeatures_gfp = mean(nFeature_gfp), 
                   mean_nFeatures_guide = mean(nFeature_guide)
  )


#perform basic qc & filter cells 
cells_to_keep <- WhichCells(gex.combined, expression = nFeature_RNA > 40 & nFeature_RNA < 200 & nCount_RNA > 150 & nCount_RNA < 6000) 
qc_combined <- subset(gex.combined, cells = cells_to_keep)

table(qc_combined$orig.ident)

qc_measures_sum_after <- qc_combined@meta.data %>% 
  dplyr::group_by(orig.ident) %>% 
  dplyr::summarise(min_nCountRNA = min(nCount_RNA), 
                   max_nCountRNA = max(nCount_RNA), 
                   mean_nCountRNA = mean(nCount_RNA), 
                   min_nFeatures = min(nFeature_RNA),
                   max_nFeatures = max(nFeature_RNA), 
                   mean_nFeatures = mean(nFeature_RNA),
                   mean_nFeatures_gfp = mean(nFeature_gfp), 
                   mean_nFeatures_guide = mean(nFeature_guide))

VlnPlot(qc_combined, features = c("nFeature_RNA", "nCount_RNA"), pt.size = 0.1, ncol = 3)  

#______this filtering can be adjusted_____________
#__________________________________________________

#normalize the data after joining the count layers 
#Variable features before running SCT transform bc otherwise it doesnt work 

joined_layers <- JoinLayers(qc_combined)
joined_layers <- NormalizeData(joined_layers)

joined_layers <- FindVariableFeatures(joined_layers, selection.method = "vst", nfeatures = 500)
variable_features <- VariableFeatures(joined_layers)

variable_tap <- intersect(variable_features, tap_genes) #these are 100 genes 

#save the list of variable tap genes here 
#write_csv(as.data.frame(variable_tap), "./list_of_highly_variable_tap_genes.csv")

qc_combined <- SCTransform(qc_combined, verbose = FALSE)
DefaultAssay(qc_combined) <- "SCT"

#specificy the TAP genes as highly variable ones 
VariableFeatures(qc_combined) <- tap_genes

#Scale the data 
#scaling on only variable features, bc all the genes scaled values take up a looot of memory 
qc_combined <- ScaleData(qc_combined, features = tap_genes) #or tap measured 

#Compute PCA
qc_combined <- RunPCA(qc_combined, features = variable_tap) #or tap measured 

DimPlot(qc_combined, reduction = "pca", group.by = 'orig.ident') + NoLegend()
ElbowPlot(qc_combined)

#Find the nearest neighbors 
qc_combined <- FindNeighbors(qc_combined, dims = 1:12) #30 works quite well already, everything > 30 is too much.I guesssomething between 20 & 30 is ideal 
qc_combined <- FindClusters(qc_combined, resolution = 0.35)
qc_combined <- RunUMAP(qc_combined, dims = 1:12)



cols <- c("#B0799AFF", "#17154FFF","#D2D2C3FF","#5B859EFF", "#2F357CFF", "#6C5D9EFF", "#9D9CD5FF", "#F6B3B0FF", "#E48171FF", "#BF3729FF", "#F5BB50FF", "#ADA43BFF", "#355828FF", "#E69B00FF")
another_pal<- c("#17154FFF","#3E5496FF", "#8290BBFF", "#008ECEFF", "#59C7EBFF", "#077187FF", "#6AAAB7FF", "#0A9086FF", "#54BFB7FF", "#8E2043FF", "#BC7A8FFF", "#E0607EFF", "#ECA0B2FF", "#FEA090FF", "#FECFC7FF", "#B8BCC1FF", "#E1E2E5FF", "#ADA43BFF", "#355828FF") 


DimPlot(qc_combined, reduction = "umap", shuffle = TRUE)
DimPlot(qc_combined, reduction = "umap", group.by = "orig.ident", shuffle = TRUE)
DimPlot(qc_combined, reduction = "umap", group.by = 'SCT_snn_res.0.35', shuffle = TRUE, cols = another_pal, label = TRUE)

FeaturePlot(qc_combined, features = c("nFeature_RNA", "nCount_RNA"), reduction = "umap") & 
  scale_color_viridis(option = "magma")

FeaturePlot(qc_combined, features = c("nFeature_gfp", "nCount_gfp"), reduction = "umap") & 
  scale_color_viridis(option = "magma")

FeaturePlot(qc_combined, features = c("nFeature_guide", "nCount_guide"), reduction = "umap") & 
  scale_color_viridis(option = "magma")

#plot some marker genes important for hematopoiesis 
FeaturePlot(qc_combined, features = c("Ctss", "Hbb-bt", "Gata2", "Mpo", "Gzmb", "Kit", "Itga2b", "S100a8", "Vwf", "Hbb-bs", "Ccl4", "Dlk1", "Mmp12", "Cd7"), max.cutoff = 3, cols = c("grey", "red"))


#projection of cell types 
#or just directly project int onto reference of Roberts and Alejo and everything 
load("~/cluster/lvelten/Analysis/SCG4SYN/Citeseq_Lenti/February_all_integrated.rda")
DimPlot(combined, reduction = "umap", label = TRUE)

#save seurat object 
#saveRDS(qc_combined, "./scMPRA4_11qced_all_layers.RDS")

#___________look at the features now___________________________________________________
qc_combined <- readRDS("./scMPRA4_11qced_all_layers.RDS")

#looking at the feature and count distribution across runs 
gfp_features <- ggplot(qc_combined@meta.data, aes(y=nFeature_gfp, x=orig.ident)) + 
  geom_violin() + 
  geom_boxplot(width=0.1) + 
  theme_bw()

gfp_counts <- ggplot(qc_combined@meta.data, aes(y=nCount_gfp, x=orig.ident)) + 
  geom_violin() + 
  geom_boxplot(width=0.1) + 
  theme_bw()

guide_features <- ggplot(qc_combined@meta.data, aes(y=nFeature_guide, x=orig.ident)) + 
  geom_violin() + 
  geom_boxplot(width=0.1) + 
  theme_bw()

guide_counts <- ggplot(qc_combined@meta.data, aes(y=nCount_guide, x=orig.ident)) + 
  geom_violin() + 
  geom_boxplot(width=0.1) + 
  theme_bw()

(gfp_features + gfp_counts) /
  (guide_features + guide_counts)


#now we perform: guide calling and then later associate the corresponding GFP count 

#create binary matrix for guide_1,...,guide_n with TRUE or FLASE values if guide count > threshold. 
qc_combined <- JoinLayers(qc_combined, assay = "guide")
qc_combined <- JoinLayers(qc_combined, assay = "gfp")

# convert to triplet form, 
#Where i represents the row indices, j represents the column indices, and x represents the values of the non-zero elements of the sparse matrix.
sparse_guide_long <- summary(qc_combined@assays$guide$counts)
sparse_gfp_long <- summary(qc_combined@assays$gfp$counts)

#assign the appropriate guide/gfp and cell name 
sparse_guide_long$guide <- rownames(qc_combined@assays$guide)[sparse_guide_long$i]
sparse_gfp_long$gfp <- rownames(qc_combined@assays$gfp)[sparse_gfp_long$i]
  
sparse_guide_long$cell <- colnames(qc_combined)[sparse_guide_long$j]
sparse_gfp_long$cell <- colnames(qc_combined)[sparse_gfp_long$j]

sparse_guide_long$guide_id <- str_split_i(sparse_guide_long$guide, pattern = "-", 1)
sparse_gfp_long$gfp_id <- str_split_i(sparse_gfp_long$gfp, pattern = "-", 1) 

nrow(sparse_guide_long)
nrow(sparse_gfp_long)

#we have 133095 guide observations and 94515 GFP observations


#______unique guide_gfp_combo filtering, based on association counts 
# in association we have 1,869,829 mio reads 
suma_per_gfp <- assoc %>% 
  group_by(guide_id) %>% 
  dplyr::mutate(prop_reads_GFP_per_guide = reads/sum(reads))

#pick the combinations with > 80% reads
selected_combos <- suma_per_gfp %>% filter(prop_reads_GFP_per_guide > 0.8)
dim(selected_combos) # left with 1,737087

#now we have to get rid of the rest in the oservations table
sparse_guide_long_filtered_unique_gfps <- merge(sparse_guide_long, selected_combos[,c("guide_id", "gfp_id")])
nrow(sparse_guide_long_filtered_unique_gfps)

#you end up with 126,276 observations, where you have a unique GFP and guide assignment (this is an average of 3 observations per cell)
# a bit better than before (before with unique/non-unique filtering we had 124,046 obs). 

#___________________________________________________________

#lets now map the counts 
#first i aggregate the guide longtable with_gfp_reporter and the GRE
sparse_guide_long_filtered_unique_gfps <- merge(sparse_guide_long_filtered_unique_gfps, assoc[,c("guide_id", "CRS", 'gfp_id')]) %>% 
  select(guide_id, x, cell, CRS, gfp_id) %>% 
  dplyr::rename(guide_count = x) 

sparse_gfp_long <- sparse_gfp_long %>% 
  select(gfp_id, cell, x) %>% 
  dplyr::rename(gfp_count = x)

final_GRE_counts <- left_join(sparse_guide_long_filtered_unique_gfps, sparse_gfp_long)

#adding meta data 
meta_file <- read.table("~/cluster/project/SCG4SYN/LibraryDesign/scMPRA/Library_Beta/final_meta_libBeta.csv", sep = ";", header = TRUE)

final_GRE_counts <- final_GRE_counts %>% dplyr::rename(Seqname = CRS)
final_GRE_counts <- left_join(final_GRE_counts, meta_file) 

#here im deciding the cut-off, wich cell is presence barcode positive 

guide_threshold <- 2

final_GRE_counts <- final_GRE_counts %>% 
  filter(guide_count > guide_threshold)

dim(final_GRE_counts)  

#saveRDS(final_GRE_counts, "./libX_matched_gfp_guide_enhancer_counts.rds") #these count tables can be found on figshare




