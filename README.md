# Protein-prediction-models
Plamsa protein prediction models using cis and trans SNPs

In subcohort 1 (N=2,481), protein levels were log transformed and adjusted for age, sex, duration between blood draw and processing, and the first three principal components of ancestry. For the rank-inverse normalized residuals of each protein of interest, we followed the TWAS/FUSION framework [13] to develop genetic prediction models, using nearby SNPs (within 100kb) of potentially associated SNPs as potential predictors. A false discovery rate (FDR) < 0.05 and  *P* -value ≤ 5×10^-8^ were used to determine potentiallly associated SNPs in *cis-* and *trans-* region. We defined  *cis-* region as a region within 1 Mb of the transcriptional start site (TSS) of the gene encoding the target protein of interest. Then, we extracted all SNPs located within 100 kb of the aforementioned potentially associated SNPs to serve as predictors for building protein prediction models, excluding the ambiguous SNPs. Four methods, namely, best linear unbiased predictor, elastic net, LASSO, and top1 were used for establishing the models.

For developed protein prediction models with prediction performance (R ^2^ ) of at least 0.01, we further conducted external validation using subcohort 2 (N=820) data.



example codes:
Rscript FUSION.assoc_test.R --sumstats Alzheimer_GWAS_summary_example.txt --weights_dir models/ --weights pos_file/ATP1A1.11993.227.3.pos --ref_ld_chr LD_reference/ATP1A1.11993.227.3.ld.1000g. --chr Z --out ATP1A1.11993.227.3_on_AD_association.txt

Rscript FUSION.assoc_test.R --sumstats Alzheimer_GWAS_summary_example.txt --weights_dir models/ --weights pos_file/DCK.9836.20.3.pos --ref_ld_chr LD_reference/DCK.9836.20.3.ld.1000g. --chr Z --out DCK.9836.20.3_on_AD_association.txt
