# NIAGADS
Note: This is a summary of where to find files and analyses. README's with actual code are listed and described after each number

# Initial data sources:
## Phenotype
>--Phenotype file: /nfs/beluga0_home/DATA/NIAGADS/NG00020/NG00020_LOAD_phenotype_files/NG00020_LOAD_phenotype/phenotype/2009_08_24_phenotype_data.csv
>--APOE genotypes: /nfs/beluga0_home/DATA/NIAGADS/NG00020/NG00020_LOAD_phenotype_files/NG00020_LOAD_phenotype/phenotype/2009_12_17 APOE results.csv
>--Age and sample source: /nfs/beluga0_home/DATA/NIAGADS/NG00020/NG00020_LOAD_phenotype_files/NG00020_LOAD_phenotype/phenotype/2009_11_16\ Age\ and\ Source\ File.csv	(note: source of Sampled.Age and AgeAtLastEval variable... the phenotype file also has an AgeAtLastEval variable, but it is less complete)

## PCs/GRM
--PC Relate GRM (by Jai): /nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/pcs_grm_pheno_only/out/b_pcrelate.rds
--PC AiR PCs (by Jai): /nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/pcs_grm_pheno_only/out/b_pcair.rds

## Variant
--gds: /nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/LOAD_families_010413_genotyped.gds


# 1. PHENOTYPES
## /acct/hkings/NIAGADS/DATA/README_NIAGADS_phenotype_harmonization.txt
-Combines phenotype file, APOE genotypes file, and Age and sample sources file

-creates variables:
	-E2 (count)
	-E3 (count)
	-E4 (count)
	-E2_carrier
	-E4_carrier
	-E4_het_hom
	-censored_age (AAO for cases, AgeAtLastEval for controls)
	-surv (survival time)
	-Dev_Res (deviance residuals based on hazzard regression of "surv"... for AAO analysis)

-Extracts PCs 1-8 from b_pcair.rds and adds as variables (the 1st 4 PCs are use in analyses)

-generates annotated dataframe (/acct/hkings/NIAGADS/DATA/annot.rds)

-generates sample filters:
	-/acct/hkings/NIAGADS/DATA/keep_samples_noTwins.rds... samples excluded for:
		-no APOE genotype
		-1 individual per twin pair excluded based on PC Relate object	(Note: file without twins excluded is called /acct/hkings/NIAGADS/DATA/keep_samples_e4_count)
		-No Case_Control status
	-/acct/hkings/NIAGADS/DATA/keep_samples_censored_analyses_noTwins.rds 	(This is the sample filter used in ALL analyses)... samples excluded for all above and also:
		-no AAO for cases
		-AAO after AgeDxDem for cases
		-AAO > 10 years before AgeDxDem for controls
		-no AgeAtLastEval for controls
		-AgeAtLastEval differs between phenotype and age_and_sample_source file


# 2. VARIANTS
## /acct/hkings/NIAGADS/DATA/README_NIAGADS_variant_filter.txt
-Generates a variant sumamry file /acct/hkings/NIAGADS/DATA/variant_metrics.rds:
	-variant.id
	-missing.rate
	-ref.freq
	-maf

-Checks missingness rate for samples (none are < 0.05%)

-plots variant MAF, missingness, and sample missingness

-Generates a variant filter ~/NIAGADS/DATA/vars_qc_pre_HWE.rds:
~/NIAGADS/DATA/vars_qc_pre_HWE.rds (note: this is the exact same file as /acct/hkings/NIAGADS/DATA/vars_qc_keep.rds, which is used in assoc tests until HWE done)
	-maf < 0.05
	-maf > 0
	-missing.rate < 0.05
	-HWE not calculated yet!!!

-calculated type I error... outputs saved to: /acct/hkings/NIAGADS/calculate_typeI_error/



# 3. ANALYSES
## files used in all analyses can be found here: /acct/hkings/NIAGADS/DATA/ and here: /nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/pcs_grm_pheno_only/out/
	-annotated pehnotype file (with PCs): /acct/hkings/NIAGADS/DATA/annot.rds
	-gds: /nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/LOAD_families_010413_genotyped.gds
	-GRM: /nfs/beluga0_home/ANALYSIS/NIAGADS/NG00020/pcs_grm_pheno_only/out/b_pcr_mat.rds
	-sample filter: ~/NIAGADS/DATA/keep_samples_censored_analyses_noTwins.rds
	-variant filter (WILL BE UPDATED AFTER HWE DONE): ~/NIAGADS/DATA/vars_qc_keep.rds

## config files, analyses, and follow-up annotation can be found in subfolders here: /acct/hkings/NIAGADS/analyses
-Analyses included in paper are (note: covars for all includes GRM):
	-AD (binary case/control); covars: PCs 1-4: /acct/hkings/NIAGADS/analyses/AD/AD_4PCs/assoc_readme.txt
	-AD (binary case/control); covars: PCs 1-4, E2 (count), E4 (count): /acct/hkings/NIAGADS/analyses/AD/AD_xE4E2_4PCs/assoc_readme.txt
	-AAO (DevRes continuous); covars: PCs 1-4: /acct/hkings/NIAGADS/analyses/AAO/AAO_4PCs/assoc_readme.txt
	-AAO (DevRes continuous); covars: PCs 1-4, E2 (count), E4 (count): /acct/hkings/NIAGADS/analyses/AAO/AAO_xE4E2_4PCs/assoc_readme.txt
	-E2_carrier; covars: PCs 1-4: /acct/hkings/NIAGADS/analyses/APOE/E2_carrier/E2_carrier_4PCs/assoc_readme.txt
	-E2_carrier; covars: none: /acct/hkings/NIAGADS/analyses/APOE/E2_carrier/E2_carrier_noPCs/assoc_readme.txt
	-E4_carrier; covars: PCs 1-4: /acct/hkings/NIAGADS/analyses/APOE/E4_carrier/E4_carrier_4PCs/assoc_readme.txt
	-E4_carrier; covars: none: /acct/hkings/NIAGADS/analyses/APOE/E4_carrier/E4_carrier_noPCs/assoc_readme.txt
	-E4_het_hom (NOT DONE); covars: PCs 1-4: /acct/hkings/NIAGADS/analyses/APOE/E4_het_hom/E4_het_hom_4PCs/assoc_readme.txt
	-E4_het_hom (NOT DONE); covars: none: /acct/hkings/NIAGADS/analyses/APOE/E4_het_hom/E4_het_hom_noPCs/assoc_readme.txt
	-I also referenced this in the paper:
		-sex; autosome-only; covars: PCs 1-4: /acct/hkings/NIAGADS/analyses/sex/sex_4PCs/assoc_readme.txt
		-sex; autosome-only; covars: none: /acct/hkings/NIAGADS/analyses/sex/sex_noPCs/assoc_readme.txt
(note: analyses not listed above have not been kept up to date and may not utilize the same input files)



# 4. Analysis of GWAS hits and comparisons across analyses and to other papers: see assoc_readme.txt in each analysis subfolder
-Identifies lead SNPs below a P-value threshold (5x10^-5) for each GWAS
	-ex. output: /acct/hkings/NIAGADS/analyses/APOE/E2_carrier/E2_carrier_noPCs/SNP_peaks.txt

-Identifies and plots ranges surrounding each lead SNP (note: ranges are an estimate that should be checked visually)
	-ex. output: /acct/hkings/NIAGADS/analyses/APOE/E2_carrier/E2_carrier_noPCs/E2_carrier_noPCs_peak_ranges.txt

-Looks for overlap between E2/E4 (no PCs) GWAS hits (from peak_ranges files) and hits from other papers:
	-Andrews meta analysis: /acct/hkings/NIAGADS/DATA/andrews36_ranges.txt
		-see: /acct/hkings/NIAGADS/DATA/notes_on_Andrews.txt
		-ex. output: /acct/hkings/NIAGADS/analyses/APOE/E2_carrier/E2_carrier_noPCs/E2_carrier_noPCs_Andrews_intersect_ranges.txt
	-Grinde long-range LD regions: /acct/hkings/NIAGADS/DATA/Grinde36.txt
		-see: /acct/hkings/NIAGADS/DATA/notes_on_Grinde.txt
		ex. output: /acct/hkings/NIAGADS/analyses/APOE/E2_carrier/E2_carrier_noPCs/E2_carrier_noPCs_Grinde_intersect_ranges.txt

-Cross-analysis checks: /acct/hkings/NIAGADS/analyses/cross_analyses_tests_readme.txt
	-correlation between AD GWAS with and without APOE genotype control & correaltion between AAO GWAS with and without APOE genotype control
	-Look up lead SNPs from E4_carrier GWAS in AD and AAO (with and without APOE genotype control) GWAS
	-Look up lead SNPs from E2_carrier GWAS in AD and AAO (with and without APOE genotype control) GWAS 




# 5. Data summary (generate table 1 / other data checks): /acct/hkings/NIAGADS/DATA/README_additional_data_summary_for_paper.txt

## Past_AD_GWAS_and_summary
Summarizes results from prior GWAS of AD (that also control for APOE E2 and/or E4 in some way)  
Tests correlation between these studies and plots Miami plots 
