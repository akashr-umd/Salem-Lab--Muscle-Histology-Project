library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(patchwork)


setwd("C:/Users/orsoo/OneDrive/Documents/Aquaculture/Muscle Histology Work/GWAS")
data <- read.table("117_sig_snps_snpeff_out/Association_Analysis_Histology_with_PCA_and_Family.Filled_Area.assoc.linear.adjusted.txt" , header = T)




# Load your data
# It must contain: CHR, GENOME_POS, SNP, BONF
# For example:
# data <- read.delim("your_file.tsv", header = TRUE)

# Define expected chromosome order
all_chr_levels <- as.character(1:28)

# Prepare data
data <- data %>%
  mutate(
    CHR = as.character(CHR),
    CHR = factor(CHR, levels = all_chr_levels),
    GENOME_POS = as.numeric(GENOME_POS)
  )

# Compute chromosome lengths and offsets
chromosome_offsets <- data %>%
  group_by(CHR) %>%
  summarize(chr_len = max(GENOME_POS, na.rm = TRUE), .groups = "drop") %>%
  complete(CHR = all_chr_levels, fill = list(chr_len = 0)) %>%
  arrange(as.numeric(as.character(CHR))) %>%
  mutate(chr_offset = lag(cumsum(chr_len), default = 0))

# Join offsets and compute cumulative position
data <- data %>%
  left_join(chromosome_offsets, by = "CHR") %>%
  mutate(
    cum_pos = GENOME_POS + chr_offset,
    logBONF = -log10(BONF) ,
    logFDR_BY = -log10(FDR_BY) , 
    logUNADJ = -log10(UNADJ)

  )



# Prepare x-axis labels at chromosome centers
axis_df <- chromosome_offsets %>%
  mutate(center = chr_offset + chr_len / 2)

# Build plot
manhattan <- ggplot(data, aes(x = cum_pos, y = logUNADJ, color = CHR)) +
  geom_point(size = 1.5, alpha = 0.8) +
  geom_hline(yintercept = 6.1, linetype = "dashed", color = "red") +
  annotate("text" , x = 900000000 , y = 5.9 , label = "BONFERRONI FDR ADJUSTED P < 0.05 CUTOFF") + 
  scale_x_continuous(
    label = axis_df$CHR,
    breaks = axis_df$center,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_color_manual(values = rep(brewer.pal(8, "Dark2"), length.out = 28)) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14), 
    axis.title = element_text(size = 14)
  ) +
  labs(x = "Chromosome", y = "-log10(unadjusted p-value)")


manhattan



data_chrom_2 <- data[data$CHR == "2" , ]



format_millions <- function(x) {
  paste0(x / 1e6, "Mb")
}






data_chrom_2_subset <- data_chrom_2[data_chrom_2$FDR_BY <= 0.05 ,  ]
data_chrom_2_subset <- data_chrom_2_subset[data_chrom_2_subset$logUNADJ  > 6.3 ,  ]


data_chrom_2_subset_adenylate_cyclase_type_3  <- data_chrom_2_subset[data_chrom_2_subset$GENOME_POS  < 43131104 ,  ]




manhattan_adenylate_cyclase_type_3 <- ggplot(data_chrom_2_subset_adenylate_cyclase_type_3, aes(GENOME_POS, logUNADJ, color = CHR)) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_manual(values =   c("#E7298AFF")) + 
  geom_hline(yintercept = 6.1, linetype = "dashed", color = "red") +
  theme_bw() +
  theme(axis.title.x = element_blank() , legend.position = "none") + 
  labs(y = "-log10(unadjusted p-value)") +
  scale_x_continuous(position = "top", labels = format_millions) + 
  scale_y_reverse()
  
  


manhattan_adenylate_cyclase_type_3




data_chrom_2_subset_TGF2Beta <- data_chrom_2_subset[data_chrom_2_subset$GENOME_POS  > 43131104 & data_chrom_2_subset$GENOME_POS  < 47638526,  ]




manhattan_TGF2Beta <- ggplot(data_chrom_2_subset_TGF2Beta, aes(GENOME_POS, logUNADJ, color = CHR)) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_manual(values =   c("#E7298AFF")) + 
  geom_hline(yintercept = 6.1, linetype = "dashed", color = "red") +
  theme_bw() +
  theme(axis.title.x = element_blank() , legend.position = "none") + 
  labs(y = "-log10(unadjusted p-value)") +
  scale_x_continuous(position = "top", labels = format_millions) + 
  scale_y_reverse() 



manhattan_TGF2Beta





data_chrom_2_subset_47Mb <- data_chrom_2_subset[data_chrom_2_subset$GENOME_POS  > 47001104 ,  ]




manhattan_47Mb <- ggplot(data_chrom_2_subset_47Mb, aes(GENOME_POS, logUNADJ, color = CHR)) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_manual(values =   c("#E7298AFF")) + 
  geom_hline(yintercept = 6.1, linetype = "dashed", color = "red") +
  theme_bw() +
  theme(axis.title.x = element_blank() , legend.position = "none") + 
  labs(y = "-log10(unadjusted p-value)") +
  scale_x_continuous(position = "top", labels = format_millions) + 
  scale_y_reverse() 



manhattan_47Mb






Patch <- manhattan / ( manhattan_adenylate_cyclase_type_3 | manhattan_TGF2Beta | manhattan_47Mb ) + 
  plot_annotation(title = 'Signficaint QTL associated with fiber density size and growth performance traits' , tag_levels = 'A')


ggsave("Manhattan_filled_area_patch_3.svg" , Patch , height = 10 , width = 14 , dpi = 300)



