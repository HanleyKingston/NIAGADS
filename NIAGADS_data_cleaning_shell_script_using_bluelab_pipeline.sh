data.dir <- "/nfs/beluga0_home/DATA/NIAGADS/NG00020/data/"
initial_qc.dir <- "/nfs/beluga0_home/DATA/NIAGADS/NG00020/initial_qc/out/"



# 1. Convert plink files from build 36 to build 37
## produces:
### plink files (ignore): /nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/data
### updated-map gds (to be used for variant QC, PC and GRM generation, and association testing): /nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/data/D_LOAD_families_010413_aut.gds

```{r}
map_37 <- read.table("/nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/data/Human610-Quadv1_H.bed", header = FALSE, skip = 1)

#Load Plink map file to compair against
bim <- read.table("/nfs/beluga0_home/DATA/NIAGADS/NG00020/NG00020_genotype_files/LOAD_families_010413.bim")

head(bim)
#  V1            V2 V3    V4 V5 V6
#1 10  rs35418599DB  0 59083  0  0
#2 10 cnvi0015449DB  0 78874  0  0
#Since variant IDs end in DB, and don't want to alter this file, will add DB to end of variant IDs in map_37

colnames(bim) <- c("chr", "variant.name", "cM", "pos", "allele_1", "allele_2")
colnames(map_37) <- c("chr", "start_pos", "end_pos", "variant.name")
#Note: variant.ids that we want (rsIDs) are further down in file
#Note: SNP position corresponds with end pos. Plink only wants a single position

head(map_37)
#  chr start_pos   end_pos variant.name
#1   9 139906358 139906359       200003
#2   9 139926401 139926402       200006
#3   2 220084901 220084902       200047

map_37$chr <- gsub("chr", "",  map_37$chr)
map_37$variant.name <- paste0(map_37$variant.name, "DB")

sum(map_37$variant.name %in% bim$variant.name)/nrow(map_37)
#[1] 1 #All variants from map file are in bim file

sum(bim$variant.name %in% map_37$variant.name)/nrow(bim)
#[1] 0.9897101 #99% of variants from bim file are in map file

#What about autosomes only?
sum(bim[bim$chr != "X" & bim$chr != "XY" & bim$chr != "Y",]$variant.name %in% map_37[map_37$chr != "X" & map_37$chr != "Y",]$variant.name)/sum(bim$chr != "X" & bim$chr != "XY" & bim$chr != "Y")
#[1] 0.9902216

#Note: variants in gds that aren't in map file are unchanged when converting in plink, so
#Make a vector of variants not in build 37 to remove from plink files
exclude_SNPs <- bim[!(bim$variant.name %in% map_37$variant.name), "variant.name", drop = TRUE]
length(exclude_SNPs)
#[1] 6389

write(exclude_SNPs, "/nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/data/exclude_SNPs_plink.txt", ncolumns = 1)

#Generate plink friendly files
map_37_pos <- map_37[, c("variant.name", "end_pos")]
map_37_chr <- map_37[, c("variant.name", "chr")]


write.table(map_37_pos, "/nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/data/map_37_pos.txt", col.names = FALSE, quote = FALSE, row.names = FALSE)

write.table(map_37_chr, "/nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/data/map_37_chr.txt", col.names = FALSE, quote = FALSE, row.names = FALSE)
```

```{bash}
cd /nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/data
```


```{bash}

## Remove SNPs not in map 37 file
~/programs/plinkV1.9/plink --noweb --bfile /nfs/beluga0_home/DATA/NIAGADS/NG00020/NG00020_genotype_files/LOAD_families_010413 --exclude exclude_SNPs_plink.txt --make-bed --out LOAD_families_010413_only_build37_SNPs
#620901 variants loaded from .bim file.
#614512 variants remaining.

## update pos
~/programs/plinkV1.9/plink --noweb --bfile A_LOAD_families_010413_only_build37_SNPs --update-map map_37_pos.txt --make-bed --out B_LOAD_families_010413_new_pos
#614512 variants and 9396 people pass filters and QC.
#Note: variants will no be ordered by new positions until make-bed is run again!

## update chr
~/programs/plinkV1.9/plink --noweb --bfile B_LOAD_families_010413_new_pos --update-map map_37_chr.txt --update-chr --make-bed --out C_LOAD_families_010413_new_pos_chr
#614512 variants and 9396 people pass filters and QC.

## exclude non-autosomes
~/programs/plinkV1.9/plink --noweb --bfile C_LOAD_families_010413_new_pos_chr --autosome --make-bed --out D_LOAD_families_010413_aut
#594742 variants and 9396 people pass filters and QC.

```

```{bash}
#Convert to GDS
/acct/hkings/git/blue_pipeline/assoc_pipeline/plink_to_gds.R D_LOAD_families_010413_aut
```


# 2. Combine phenotype files and perform QC 
## Produces:
### annotated phenotype dataframe (PCs will be added and will be used for assoc testing): /nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/pheno/out/pheno_adf.rds
### sample filter (to be used for generatign PCs and GRM): /nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/pheno/out/pass_pheno_qc.rds
### summary figures and tables

```{r}
#Note: Walk through script in R, don't execute as a singel script
/nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/pheno/pheno_qc.Rmd
#Output files are here: /nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/pheno/out
```


# 3. Variant and sample QC - check for missingness
## produces:
### summary of sample missingness: /nfs/beluga _home/ANALYSIS/NIAGADS/NG00020/hkings/initial_qc/out/missing_by_sample.rds
### sample missingness figures
### Note there are no samples that fail based on our missingngess filters. If there were, these would need to be intersected with those that fail phenotype QC
### summary of variant qc measures: /nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/initial_qc/out/variant_metrics.rds
### variant filter (to be fed to LD-pruning; doe not include HWE filtering, this will be added before association testing): /nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/initial_qc/out/vars_qc_keep.rds
### variant missingness figures


```{bash}
## sample QC
/acct/hkings/git/blue_pipeline/assoc_pipeline/sample_qc.R /nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/data/D_LOAD_families_010413_aut.gds --out_prefix /nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/initial_qc/out/
#The samples with 100% missingness are the "dummy" founder samples included for plink files


## variant QC
### Note: must include sampel filter because there are a lot of dummy samples and otherwise this will lead to high missingness
/acct/hkings/git/blue_pipeline/assoc_pipeline/variant_qc.R /nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/data/D_LOAD_families_010413_aut.gds --sample_id /nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/pheno/out/pass_pheno_qc.rds --out_prefix /nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/initial_qc/out/ 
```


#Extract variant and sample filters from dataframes
```{r samp_qc}
library(dplyr)

samp_qc <- readRDS("/nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/initial_qc/out/missing_by_sample.rds")
summary(samp_qc$missing.rate)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.03986 0.04020 0.04201 0.46705 1.00000 1.00000
```
There are variants with high missingness...
```{r varmet}
varmet <- readRDS("/nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/initial_qc/out/variant_metrics.rds")
summary(varmet$missing.rate)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#0.0000000 0.0002267 0.0004533 0.0409082 0.0009066 1.0000000
```

...as well as rare variants...

```{r maf}
summary(varmet$maf)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.0000  0.1044  0.2176  0.2286  0.3526  0.5000
nrow(varmet)
#[1] 594742

...so we'll create a variant keep list:


#Get counts:
nrow(varmet[!(varmet$missing.rate < 0.05),])
#[1] 24745 #excluded for missingness
nrow(varmet[varmet$maf < 0.05 & varmet$maf > 0,])
#[1] 49050 #excluded for MAF > 0.05

keep_vars <- filter(varmet, missing.rate < 0.05, !(maf < 0.05 & maf > 0))$variant.id

length(keep_vars)
#[1] 521088

saveRDS(keep_vars, "/nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/initial_qc/out/vars_qc_keep.rds")

```



# 4. Develop PCs and GRM
## produces: (all files in /nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/pcs_grm_pheno_only/out)
### pruned SNP list:pruned_snps.rds
### intermediate PC and GRM files (king, and A PC and GRM files)
### final PC files:b_pcair.rds (complete PC-Air object), b_pcair_pcs.rds (just the PCs to be merged with phenotype file for assoc. testing), b_pcair_unrels.rds, b_pcair_rels.rds
### final GRM files:b_pcr_mat.rds (relatedness matrix to be passed to association testing), b_pcrelate.rds (complete PC-Relate object)
### Plots:
#from PC-Relate: b_pcrelate_kinship.png (kinship plot)
#from PC-AiR:  b_paracord.png, b_pairs.png (pairwise PC scatterplots), b_pc12.png (PC1xPC2), b_pc34.png (PC3xPC4), b_pc_scree.png (% variance explained by each PC)

```{bash}
cd /nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/pcs_grm_pheno_only/out

## run LD-pruning:
/nfs/beluga2_home/hkings/git/blue_pipeline/assoc_pipeline/ld_pruning.R /nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/data/D_LOAD_families_010413_aut.gds --variant_id /nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/initial_qc/out/vars_qc_keep.rds --sample_id /nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/pheno/out/pass_pheno_qc.rds > LD_pruning.log &
#Used defaults for:
#maf: 0.05
#missingness: 0.05
#r-threshold: sqrt(0.1)
#Did not excldue "problematic" PC regions like HLA or LCT
#95,639 markers are selected in total.

## Run King
/nfs/beluga2_home/hkings/git/blue_pipeline/assoc_pipeline/king_grm.R /nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/data/D_LOAD_families_010413_aut.gds --variant_id pruned_snps.rds --sample_id /nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/pheno/out/pass_pheno_qc.rds > King.log &

## Plot King
/nfs/beluga2_home/hkings/git/blue_pipeline/assoc_pipeline/kinship_plot.R king_out.rds --is_king --out_file king_kinship.png > King_plot.log &

## Run 1st iteration PC air
/nfs/beluga2_home/hkings/git/blue_pipeline/assoc_pipeline/pcair.R /nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/data/D_LOAD_families_010413_aut.gds king_grm.rds king_grm.rds --out_prefix a_ --variant_id pruned_snps.rds --sample_id /nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/pheno/out/pass_pheno_qc.rds --num_core 4 > a_PC_AiR.log &
#The first file is the kin_file and the secodn file is the div_file
#used defaults for:
#kin_thresh: 2 ^ (-9 / 2)
#div_thresh: 2 ^ (-9 / 2)


## Run 1st iteration PC-Relate
/nfs/beluga2_home/hkings/git/blue_pipeline/assoc_pipeline/pcrelate.R /nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/data/D_LOAD_families_010413_aut.gds a_pcair.rds --sample_block_size 6000 --variant_block 10000 --out_prefix a_ --n_pcs 4 --variant_id pruned_snps.rds --sample_id /nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/pheno/out/pass_pheno_qc.rds > a_PC_relate.log &
#used defaults for:
#scale_kin: 1
#sparse_thresh: 0, "Threshold for sparsifying GRM (will be multiplied by scale_kin)"
  add_argument("--small_samp_correct",
               "Flag to implement small-sample correction",
               flag = TRUE) %>%



## Run 2nd iteraction PC-AiR
/nfs/beluga2_home/hkings/git/blue_pipeline/assoc_pipeline/pcair.R /nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/data/D_LOAD_families_010413_aut.gds a_pcr_mat.rds king_grm.rds --out_prefix b_ --variant_id pruned_snps.rds --sample_id /nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/pheno/out/pass_pheno_qc.rds --num_core 4 > b_PC_AiR.log &
#The first file is the kin_file and the second file is the div_file (note: that now the PC-Relate output is used as the kin_file)
#used defaults for:
#kin_thresh: 2 ^ (-9 / 2)
#div_thresh: 2 ^ (-9 / 2)

## Plot PCs
/nfs/beluga2_home/hkings/git/blue_pipeline/assoc_pipeline/pca_plots.R b_pcair.rds --phenotype_file /nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/pheno/out/pheno_adf.rds --group race_ethnicity --n_pairs 6 --out_prefix b_ > PC_plots.log &


## Run 2nd iteration PC-Relate
/nfs/beluga2_home/hkings/git/blue_pipeline/assoc_pipeline/pcrelate.R /nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/data/D_LOAD_families_010413_aut.gds b_pcair.rds --out_prefix b_ --sample_block_size 6000 --variant_block 10000 --n_pcs 4 --scale_kin 2 --variant_id pruned_snps.rds --sample_id /nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/pheno/out/pass_pheno_qc.rds > b_PC_relate.log &
#Note: scale_kin 2 argument added

## Plot PC-Relate
/nfs/beluga2_home/hkings/git/blue_pipeline/assoc_pipeline/kinship_plot.R b_pcrelate.rds --out_file b_pcrelate_kinship.png > pcrelate_plot.log &

```

# 5. Set up files for association tests: Exclude identical twins from sample keep vector and add PCs to annotated phenotype dataframe


```{r}
## 5a. Identify identical twins & generate keep_samples vector that excludes them

"%notin%" <- Negate("%in%")

#Read in PCrelate object and existing keep_samples vector
library(GENESIS)
mypcrel <- readRDS("/nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/pcs_grm_pheno_only/out/b_pcrelate.rds")
kinship <- mypcrel$kinBtwn

#Read in keep_samples file
keep_samples <- readRDS("/nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/pheno/out/pass_pheno_qc.rds")

#Find likely identical twins
identical_twins <- kinship[kinship$kin > 0.45, c("kin", "k2", "ID1", "ID2")]
identical_twins
saveRDS(identical_twins, "/nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/pheno/out/identical_twins.rds")

#Remove identical twins from keep_samples vector
exclude_twins <- kinship[kinship$kin > 0.45, "ID1"] #ID1 is always the lowest of the 2 number IDs
sum(keep_samples %in% exclude_twins)
#[1] 13

keep_samples_noTwins <- keep_samples[keep_samples %notin% exclude_twins]

saveRDS(keep_samples_noTwins, file = "/nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/pheno/out/pass_pheno_qc_noTwins.rds")

length(keep_samples_noTwins)
#[1] 4399

exclude_twins2 <- append(exclude_twins, kinship[kinship$kin > 0.45, 



## 5b. Add PCs to annotated phenotype file
library(Biobase)
library(SeqArray)

#Read in phenotype file
annot <- readRDS("/nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/pheno/out/pheno_adf.rds")
#extract just phenotype info:
phenotype <- pData(annot)

# Get PCs
pca <- readRDS("/nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/pcs_grm_pheno_only/out/b_pcair.rds")

#get gds and extract sample IDs
gds <- seqOpen("/nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/data/D_LOAD_families_010413_aut.gds")
gds.id <- seqGetData(gds, "sample.id")


##The null model requires a sample ID column:
phenotype$sample.id <- phenotype$SUBJ_NO

##Read in PCA covariates -will read in 8 but expect to use 4 (based on: file:///C:/Users/hking/OneDrive/my%20documents/UW/thesis%20and%20research/NIAGADS/NG00020_setup.html)
pcs.df <- as.data.frame(pca$vectors[,1:8])
colnames(pcs.df) <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8")

pcs.df$sample.id <- row.names(pcs.df)

#Add PCA covariates to phenotype data of interest by subject number
#Only select columns of interest
merged_phen <- merge(phenotype[, c("sample.id", "SEX", "Case_Control", "ConType", "BirthYr_pheno", "AgeDeath", "Race", "Hispanic", "race_ethnicity", "E2", "E4", "E3", "E4_carrier", "E4_het_hom", "E2_carrier", "E2_het_hom", "age", "age_source", "surv", "devres", "martingale", "pass_qc")], pcs.df, by = "sample.id", all.x=TRUE)

#order phenotype df to match order of gds ids
merged_phen <- merged_phen[match(gds.id, merged_phen$sample.id),]
identical(as.character(merged_phen$sample.id), as.character(gds.id))

colnames(merged_phen)

##annotate data:
metadata <- data.frame(labelDescription = c(
  "sample.id",
  "sex (1= male, 0 = female)",
  "case-control status",
  "ConType",
  "birth year",
  "age at death",
  "reported race",
  "reported hispanic status",
  "combined descriptor of race and Hispanic status",
  "count of E2 allele",
  "count of E4 allele",
  "count of E3 allele",
  "whether they are a carrier of the E4 allele (1 = heterozygous or homozygous for E4)",
  "whether they are a heterozygote (0) or homozygote (1) for E4 allele; homozygote for non-E4 = NA",
  "whether they are a carrier of the E2 allele (1 = heterozygous or homozygous for E2)",
  "whether they are a heterozygote (0) or homozygote (1) for E2 allele; homozygote for non-E2 = NA",
  "age censored; estimated age of symptom onset (AgeDem) for cases; AgeAtLastEval or (if not available) Sampled.Age for controls",
  "first principal component",
  "second principle component",
  "third principle component",
  "fourth principle component",
  "fifth principle component",
  "sixth principle component",
  "seventh principle component",
  "eighth principle component",
  "survival time / censoring age",
  "the source of the censoring age variable",
  "deviance residuals from survival data",
  "martingale residuals from survival data",
  "Whether they passed phenotype QC filters, does not account for identical twins"
))


#Re-generate annotated dataframe
#(This is a bit round about... it would be better to just add PCs to already made annotated dataframe, but I couldn't figure out how to do that)
annot <- AnnotatedDataFrame(merged_phen, metadata)

#Save annot as an R object
saveRDS(annot, file = "/nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/pheno/out/pheno_for_assoc.rds")
```



# 6. Calculate and plot PC correlations
## produces: (all files are in /nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/pcs_grm_pheno_only/out)
### pc-correlation file: /nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/pcs_grm_pheno_only/out/snp_corr.rds
### pc-correlation plots: /nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/pcs_grm_pheno_only/out/pc_corr.png
### zoomed-in plots from PC2 correlation: /nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/pcs_grm_pheno_only/out/PC2_chr<>_corr_peak.png

```{BASH}
cd /nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/pcs_grm_pheno_only/out

R -q --vanilla --args --pca_file "b_pcair.rds" --gds_file "/nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/data/D_LOAD_families_010413_aut.gds" --block-size 60000 --outfile "snp_corr.rds" < /nfs/beluga2_home/hkings/NIAGADS/8_30_scripts/calculate_SNP_PC_corr.R > calculate_snp_pc_corr.log &


#Note: Walk through script in R, don't execute as a single script
#MUST UPDATE WHEN HWE COMPLETED
/nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/pcs_grm_pheno_only/plot_snp_pc_corr.R
```



# 7. Test HWE in set of unrelated controls, controling for ancestry through PCs





# 8. Association testing:
See README's in AAO, AD, E4_carrier, and E2_carrier
#Shared files for assoc testing:
#phenotype file: /nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/pheno/out/pheno_for_assoc.rds
#sample filter: /nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/pheno/out/pass_pheno_qc_noTwins.rds
#gds: /nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/data/D_LOAD_families_010413_aut.gds
#GRM: /nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/pcs_grm_pheno_only/out/b_pcr_mat.rds
#TEMPORARY variant filter (no HWE applied, so need to be updated):/nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/initial_qc/out/vars_qc_keep.rds







# 9. Additional data summary for paper:
```{r}
library(SeqArray)
pheno <- readRDS("/nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/pheno/out/pheno_for_assoc.rds")
pheno2 <- pData(pheno)


## 9a. Sample summaries
### Filter phenotype data based on sample filter
pheno2 <- pheno2[pheno2$sample.id %in% samps,]
samps <- readRDS("/nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/pheno/out/pass_pheno_qc_noTwins.rds")
length(samps)
#[1] 4399


table(pheno2$Case_Control)
#   0    1
#2221 2178



## 9b. Variant summaries - NOT DONE, need HWE
vars <- readRDS("vars_for_assoc_test")
```


## 9c. calculate type I error rate

cd /nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/calculate_typeI_error
```{r}
# Determine Bonferoni-corrected p-value threshold using genetic type I error calculator (GEC)

#What are variant IDs look like in bim file (column 2)?
head /home/DATA/NIAGADS/NG00020/NG00020_genotype_files/LOAD_families_010413.bim
#10      rs35418599DB    0       59083   0       0
#10      cnvi0015449DB   0       78874   0       0

#What do sample ids look like in fam file (column 2)?
head /home/DATA/NIAGADS/NG00020/NG00020_genotype_files/LOAD_families_010413.fam
#1600 160001 0 0 1 -9
#1600 160003 0 0 2 -9


#Convert variant filter into a line-seperated list by rsid and cnvi to match bin file
#--extract accepts list of variant IDs as 1 per line or space-seperated
library(SeqArray)
library(SeqVarTools)
gds <- seqOpen("/nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/data/D_LOAD_families_010413_aut.gds")

head(seqGetData(gds, "annotation/id"))
#[1] "rs35418599DB"  "cnvi0015449DB"

keep_vars <- readRDS("vars_for_assoc_test")

seqSetFilter(gds, variant.id = keep_vars)

#Save a vector of filtered variant IDs
keep_vars_id <- seqGetData(gds, "annotation/id")

write(keep_vars_id, file = "vars_qc_keep_GEC.txt", ncolumns = 1)


#Save sample filter in format for plink
#--keep accepts a space/tab-delimited text file with family IDs in the first column and within-family IDs in the second column

keep_samps <- readRDS("/nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/hkings/pheno/out/pass_pheno_qc_noTwins.rds")

#For some reason, this seems to also require family IDs, so will add those from the fam file
fam <- read.table("/home/DATA/NIAGADS/NG00020/NG00020_genotype_files/LOAD_families_010413.fam")

IDs_keep <- fam[fam$V2 %in% keep_samps, c("V1", "V2")]
nrow(IDs_keep)
#[1] 4385

write.table(IDs_keep, file = "pass_pheno_qc_noTwins_GEC.txt", sep = " ",  col.names = F, row.names = F)
```

# Remove non-QC passing variants and excluded samples from Plink binary files
~/programs/plinkV1.9/plink --bfile /home/DATA/NIAGADS/NG00020/NG00020_genotype_files/LOAD_families_010413 --extract vars_qc_keep_GEC.txt --keep pass_pheno_qc_noTwins_GEC.txt --make-bed --out plink_filt
```

# Calculate Type I error rate

java -jar gec.jar --effect-number --plink-binary plink_filt --genome --error 0.05 --out type_1_error_thresh



```

