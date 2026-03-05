#This script extracts the barcodes of GFP and the guide RNAs from the CRS library 
#we need to pass these barcodes to cellranger, to be able to map the BC counts to the respective GRE 
#at the same time, we need a file that later on lets us back track the GREs to the barcodes in the seurat object, after we counted them using cellranger from the seurat

#For cellranger we need a .csv file in the following format (ids or sequence cant be duplicated!)
#id, name, read, pattern, sequence, feature_type 

# i used to do it like this: 
#for GFP: reporter1, reporter1_inert, R2, CTTGCTCACCATGGTGGCGACCGGT(BC) this should still work, ACGGTTAGTCGTCGG, Custom 
#for guide: guide1, guide1_inert,R2,TTCTAGCTCTAAAAC(BC),CCCCCCTAACAAGTG,Custom

#the barcode of the mapped file: 0:15 bp-> GFP bc 0:15, guide BC 16:30
#from 5' to 3' of transcript 

#we have to reverse complement all the sequences,
#because in the scMPRA screen we are sequencing them with R2 

setwd("~/cluster/project/SCG4SYN/Experiments/240524_scMPRA_LibBeta_HSC_iseq/")
bc_crs_ass_Beta <- read.csv(gzfile("./CRS_BC_ASS/outs/libBeta_novaseq_no_bc_corr.csv.gz"), sep = "\t") 

#according to the log-log plot i would choose 10 as a cutoff for correct barcodes 
bc_crs_ass_Beta_filt <- subset(bc_crs_ass_Beta, READS > 10)

#first we define a function to compute the reverse complemetn 
reverse_complement <- function(dna_sequence) {
  # Define a dictionary to map nucleotides to their complements
  complement_dict <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")
  
  # Convert the input sequence to uppercase to handle both upper and lower case
  dna_sequence <- toupper(dna_sequence)
  
  # Reverse the sequence and replace each nucleotide with its complement
  complemented_sequence <- sapply(strsplit(dna_sequence, ""), function(nucleotide) {
    paste0(rev(sapply(nucleotide, function(n) complement_dict[n])), collapse = "")
  })
  
  return(complemented_sequence)
}

#in the CRS barcode mapping file, we merged the guide and GFP barcode to one single barcode, here im splitting them 
barcodes <- bc_crs_ass_Beta_filt$BARCODE

gfp_bcs <- unlist(lapply(as.vector(barcodes), function(barcode) {
  sub_gfp <- substr(barcode, 0, 15)
  reverse_complement(sub_gfp)
}))

guide_bcs <- unlist(lapply(as.vector(barcodes), function(barcode) {
  sub_guide <- substr(barcode, 16, 30)
  reverse_complement(sub_guide)
}))

#this is the file we create not for cellranger, but later on for the seurat object to map the bc to the GRE 
mapping_file <- tibble(gfp_bc_fwd = substr(barcodes, 0, 15), 
                       guide_bc_fwd = substr(barcodes, 16, 30), 
                       gfp_bc_rc = gfp_bcs, 
                       guide_bc_rc = guide_bcs,
                       CRS = bc_crs_ass_Beta_filt$CRS,
                       reads = bc_crs_ass_Beta_filt$READS)

#make the list of barcodes unique, because they cant be duplicated in the cellranger reference file 
gfp_bcs <- cbind(gfp_bcs, bc_crs_ass_Beta_filt$CRS) %>%  as.data.frame()
guide_bcs <- cbind(guide_bcs, bc_crs_ass_Beta_filt$CRS) %>% as.data.frame()

gfps_unique <- gfp_bcs %>% distinct(gfp_bcs, .keep_all = TRUE)
guides_unique <- guide_bcs %>% distinct(guide_bcs, .keep_all = TRUE)

#set up the .csv table for cellranger 
#creating the first column 
num_elements <- length(barcodes)

num_gfp <- nrow(gfps_unique)
num_gfp_range <- seq(1,num_gfp)
num_gfp_range_with_leading_zeros <- str_pad(num_gfp_range, 8, pad = "0")

num_guide <- nrow(guides_unique)
num_guide_range <- seq(1,num_guide)
num_guide_range_with_leading_zeros <- str_pad(num_guide_range, 8, pad = "0")

#generate list to name the elements, the id column  
reporter_list <- sapply(num_gfp_range_with_leading_zeros, function(i) {
  paste0("reporter", i)
})

guide_list <- sapply(num_guide_range_with_leading_zeros, function(i) {
  paste0("guide", i)
})

id <- c(reporter_list, guide_list) 
 
#create the second column (i think we should put GRE name here (if we can))
name <- c(gfps_unique$V2, guides_unique$V2)
name <- gsub(",", "_", name)

#create the third column 
read <- rep("R2", num_gfp + num_guide)

#create the 4th column 
pattern <- c(rep("CTTGCTCACCATGGTGGCGACCGGT(BC)", num_gfp), rep("TTCTAGCTCTAAAAC(BC)", num_guide)) 

#create the 5th column 
sequence <- c(gfps_unique$gfp_bcs, guides_unique$guide_bcs)

#create the 6th column 
feature_type <- rep("Custom", num_gfp + num_guide) 

# create this kind of name that includes id + name GRE
full_id <- paste0(id, "_", name)
name <- full_id

#merge all the columns to one table 
final_feature_ref_file <- cbind(id, name, read, pattern, sequence, feature_type) 

#write.table(final_feature_ref_file, file="./screen_HSC/001_cellranger/001_create_seurat_CRS_level/feature_ref_filtered.csv", sep=",", row.names = FALSE, quote=FALSE)

#now im finishing the association table for the seurat object 
#split final_feature_ref_file in guide and GFP

#final_guides <- final_feature_ref_file[grep('guide', final_feature_ref_file[,1]),]
#final_gfps <- final_feature_ref_file[grep('reporter', final_feature_ref_file[,1]),]

#colnames(final_guides) <- c("id", "name", "read", "pattern", "guide_bc_rc", "feature_type")
#colnames(final_gfps) <- c("id", "name", "read", "pattern", "gfp_bc_rc", "feature_type")

#always needs to be double-checked here 
#mapping_file <- merge(x = mapping_file, y = final_guides[,c("guide_bc_rc", "id")])
#mapping_file <- mapping_file %>% dplyr::rename(guide_id = id)

#mapping_file <- merge(x = mapping_file, y = final_gfps[,c("gfp_bc_rc", "id")])
#mapping_file <- mapping_file %>% dplyr::rename(gfp_id = id)

#fix the GRE names (always make problems because of their commas in the names)
#mapping_file$CRS <- gsub(",", "_", mapping_file$CRS)

#write.table(mapping_file, file="/scMPRA_data_processing/002_screen_data_preprocessing/feature_ref_filtered_for_seurat.csv", sep=",", row.names = FALSE, quote=FALSE)

