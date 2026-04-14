library("qqman")
library("svglite")
library("ggplot2")

setwd("C://Users/orsoo/OneDrive/Documents/Aquaculture/Muscle Histology Work/GWAS")


input <- read.table("117_sig_snps_snpeff_out/Association_Analysis_Histology_with_PCA_and_Family.Filled_Area.assoc.linear.adjusted.txt" , header =  T)


pvals <- input$UNADJ[!is.na(input$UNADJ)]

pval_num <- as.numeric(pvals)


chisq <- qchisq(1 - pval_num , df = 1 )

lambda <- median(chisq) / qchisq(0.5 ,1 )


qqnorm(input$UNADJ)
qqline(input$UNADJ)




