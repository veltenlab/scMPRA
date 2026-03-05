# Lentiviral single-cell MPRA of synthetic enhancers

Code and instructions to reproduce the analyses and figures for: 
Lentiviral single-cell MPRA of synthetic enhancers reveals motif affinity-based encoding of cell type specificity (Julia Rühle @ Velten Lab / CRG, 2026).
[preprint](https://www.biorxiv.org/content/10.64898/2026.02.27.708495v1)

## Repository layout

```text
.
├── scMPRA_data_processing/
│   ├── 001_CRS_barcode_association/      # scripts to process barcode association library
│   └── 002_screen_data_preprocessing/    # scripts to process GFP, guide & GEX libraries
└── README.md
```

## Data availability 

- **Preprocessed objects used for plotting**: 
   - Seurat objects containing raw and normalized gene expression 
   - Preprocessed enhancer count data with synthetic enhancer annotation (for Library PoC, Library Alpha, Library Beta)
   
   are available at [https://doi.org/10.6084/m9.figshare.31420724](https://doi.org/10.6084/m9.figshare.31420724)
   
## Reproducibility 

For the experimental protocol see:

### Association Library processing: 
Generate look-up table that connects presence barcode, CRE and quantitative barcode (GFP). 

  1. Create a reference fasta file with all sequences in CRE library (such as example_fasta_CRE_reference.fa)
  2. `scMPRA_data_processing/001_CRS_barcode_association/001_pipeline_singlecell.sh` 
  Aligns reads to CREs & readsout and counts barcodes. 
  3. `scMPRA_data_processing/001_CRS_barcode_association/003_barcode_assoc_report.Rmd` 
  Plots basic Library QC features. 
  4. `scMPRA_data_processing/001_CRS_barcode_association/004_create_feature_ref_filtered.R` 
  Produces look-up table (example_feature_ref_filtered.csv) that is input file for the step 1 in the Screen Libraries processing. 
  
### Screen Libraries processing: 
Generate count tables for guideRNA library (presence barcodes), GFP library (presence barcodes) & targeted transcriptomes:

  1. `002_screen_data_preprocessing/cellranger_scMPRA_count.sh`
  Create filtered_feature_bc_matrix using Cellranger (using example_feature_ref_filtered.csv & library_scMPRA.csv)
  2. `002_integrate_seurat_create_enhancer_counts.sh`
  Create seurat object & enhancer quantification table.
  

## Citation

If you use this code, please cite:

Lentiviral single-cell MPRA of synthetic enhancers reveals motif affinity-based encoding of cell type specificity. 
Rühle, J., Frömel, R., Bernal Martínez, A., Szu-Tu, C., Bowness, J., Velten, L. 
bioRxiv 2026.02.27.708495; doi: https://doi.org/10.64898/2026.02.27.708495







  

  

   



