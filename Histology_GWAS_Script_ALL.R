
library(tidyr)
library(writexl)


1.#make bed
system("plink --vcf Histology_samples_SNPs_alone_with_Artificial_ID.vcf.gz --chr-set 29 --nonfounders --allow-no-sex --allow-extra-chr --pheno Histology_Phenotypes.txt --all-pheno --autosome --make-bed --out Histology_Genotype")
system("plink --bfile Histology_Genotype --chr-set 29 --nonfounders --allow-no-sex --allow-extra-chr --make-bed")

#Do QC 
system("plink --bfile Histology_Genotype --chr-set 29 --allow-extra-chr --nonfounders --keep-allele-order --allow-no-sex --pheno Histology_Phenotypes.txt --all-pheno --maf 0.1 --hwe 0.1 --geno 0.0005 --make-bed --out After_QC_Histology_Genotype_SNPs_Alone_Genotype_MAF_GENO_AND_HWE")


#running PCA in PLINK
system ("plink --bfile After_QC_Histology_Genotype_SNPs_Alone_Genotype_MAF_GENO_AND_HWE --chr-set 29 --nonfounders --allow-no-sex --pca --out Histology_Genotype_PCA")

#Reading the result filesfrom PCA
eigenvalues<- read.delim("Histology_Genotype_PCA.eigenval", header = F)
eigenvectors<- read.table("Histology_Genotype_PCA.eigenvec", header = F)
data.frame(eigenvectors)

#Proportion of variation captured by each vector
eigen_percent<- round((eigenvalues/ (sum(eigenvalues))*100),2)

#
library(tidyverse)
#plotting the PCA 
ggplot(data = eigenvectors) +
  geom_point(mapping = aes(x= V3, y= V4, color= V2), size = 3, show.legend = TRUE)+
  geom_hline(yintercept = 0, linetype= "dotted")+
  geom_vline(xintercept = 0, linetype= "dotted")+
  labs(title = "PCA",
       x = paste0("Principal component 1 (",eigen_percent[1,1]," %)"),
       y = paste0("Principal component 2 (",eigen_percent[2,1]," %)")) 
#convert to data frame ()
df_PCA<- data.frame(eigenvectors)


#write the dataframe into text file and save it.
write.table(df_PCA, file = "df_PCA.txt", append = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)


#from here, opened the df_PCA.txt and select 20PCs + Family 
## Load required package
library(dplyr)

## Read PCA eigenvectors from PLINK 
pcs <- read.table(
  "Histology_Genotype_PCA.eigenvec",
  header = FALSE
)

## Name columns: FID, IID, PC1, PC2, ...
colnames(pcs) <- c(
  "FID",
  "IID",
  paste0("PC", 1:(ncol(pcs) - 2))
)


## keep PC1–PC20
pcs <- pcs %>% select(FID, IID, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, PC11, PC12, PC13, PC14, PC15, PC16, PC17, PC18, PC19, PC20)

## 3 Read  phenotype file and extract Family
pheno <- read.table(
  "Histology_Phenotypes.txt",
  header = TRUE
)

fam <- pheno %>%
  select(FID, IID, Family)

## Merge PCs + Family by FID/IID
covar <- pcs %>%
  left_join(fam, by = c("FID", "IID"))

## Recode Family as numeric 
covar$Family <- as.numeric(as.factor(covar$Family))

## 6️⃣ Write final PLINK covariate file
write.table(
  covar,
  file = "PC1-PC5_Family.covar",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)



#perform association analysis with the linear function (that included Principal components PC as cofactors) and --adjust for multiple testing
system("plink --bfile After_QC_Histology_Genotype_SNPs_Alone_Genotype_MAF_GENO_AND_HWE --chr-set 29 --nonfounders --allow-no-sex --linear hide-covar --pheno Histology_Phenotypes.txt --all-pheno --adjust --ci 0.95 --covar PC1-PC5_Family.covar --out Association_Analysis_Histology_with_PCA_and_Family")