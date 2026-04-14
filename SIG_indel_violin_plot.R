library(readxl)
library(tidyr)
library(ggplot2)
library(svglite)


setwd("C:/Users/orsoo/OneDrive/Documents/Aquaculture/Muscle Histology Work/Haplotype Analysis")


metadata <- read_excel("Metadata_Table.xlsx", sheet = "Sheet1", col_names = FALSE)

cat_cell_count <- metadata[3:nrow(metadata), 1:4] %>%
  setNames(c("Sample_ID", "Phenotype_Cell", "Class_Cell", "Family_Cell")) %>%
  filter(!is.na(Sample_ID))

cat_filled_area <- metadata[3:nrow(metadata), 7:10] %>%
  setNames(c("Sample_ID", "Phenotype_Filled", "Class_Filled", "Family_Filled")) %>%
  filter(!is.na(Sample_ID))

quant <- metadata[3:nrow(metadata), 13:16] %>%
  setNames(c("Sample_ID", "Family", "Filled_Area", "Cell_Count")) %>%
  filter(!is.na(Sample_ID)) %>%
  mutate(Filled_Area = as.numeric(Filled_Area),
         Cell_Count = as.numeric(Cell_Count))

metadata_merged <- quant %>%
  left_join(cat_cell_count, by = "Sample_ID") %>%
  left_join(cat_filled_area, by = "Sample_ID") %>%
  mutate(
    Class_Cell = ifelse(Class_Cell == "High", 1, ifelse(Class_Cell == "Low", 0, NA)),
    Class_Filled = ifelse(Class_Filled == "High", 1, ifelse(Class_Filled == "Low", 0, NA))
  )

# Read genotypes (check.names = FALSE is crucial)
genotypes <- read.table("chr2_indels_full_genotypes.txt", 
                        header = TRUE, sep = "\t", 
                        stringsAsFactors = FALSE, check.names = FALSE)

genotypes_long <- genotypes %>%
  pivot_longer(cols = -c(ID, CHROM, POS, REF, ALT), 
               names_to = "Sample_ID", values_to = "GT") %>%
  mutate(
    Sample_ID = as.character(Sample_ID),
    Dose = case_when(
      GT %in% c("0/0", "0|0") ~ 0,
      GT %in% c("0/1", "1/0", "0|1", "1|0") ~ 1,
      GT %in% c("1/1", "1|1") ~ 2,
      TRUE ~ NA_real_
    )
  ) %>%
  select(ID, Sample_ID, Dose) %>%
  pivot_wider(names_from = ID, values_from = Dose)

df <- metadata_merged %>%
  left_join(genotypes_long, by = "Sample_ID")



#INTERESTED INDELS
indels <- c(
  "2_43212512_C_CAACACAGCGAGACACAACAAG",
  "2_43212632_T_TGAGAGAGAGAGA"
)




for (indel in indels) {
  cat("\n=== Tests for INDEL:", indel, "===\n")
  
  # Quantitative: Cell Count
  quant_cell <- df %>% filter(!is.na(Cell_Count) & !is.na(.data[[indel]]))
  if (nrow(quant_cell) > 0) {
    # Safe: use the filtered data's column directly
    #lm_cell <- lm(Cell_Count ~ quant_cell[[indel]], data = quant_cell)
    lm_cell <- lm(Cell_Count ~ quant_cell[[indel]] + factor(quant_cell$Family), data = quant_cell)
    print(summary(lm_cell))
  } else {
    cat("No data for Cell Count quantitative.\n")
  }
  
  # Quantitative: Filled Area
  quant_filled <- df %>% filter(!is.na(Filled_Area) & !is.na(.data[[indel]]))
  if (nrow(quant_filled) > 0) {
    #lm_filled <- lm(Filled_Area ~ quant_filled[[indel]], data = quant_filled)
    lm_filled <- lm(Filled_Area ~ quant_filled[[indel]] + factor(quant_filled$Family), data = quant_filled)
    print(summary(lm_filled))
  } else {
    cat("No data for Filled Area quantitative.\n")
  }
}



#Find the one sample outlier
# For INDEL 2
indel2 <- "2_43212632_T_TGAGAGAGAGAGA"
low_outlier <- quant_filled %>%
  filter(.data[[indel2]] == 2) %>%
  select(Sample_ID, Filled_Area, Family, everything())
print(low_outlier)

#The outlier is sample 269 and it has one of the lowest filled area phenotype


#Re-run without that one outlier sample (269)
quant_filled_no_outlier <- quant_filled %>% filter(.data[[indel2]] != 2 | is.na(.data[[indel2]]))
lm_no_out <- lm(Filled_Area ~ quant_filled_no_outlier[[indel2]] + factor(quant_filled_no_outlier$Family),
                data = quant_filled_no_outlier)
summary(lm_no_out)






#PLOT
quant_filled_no_outlier$adjusted_Filled <- lm_no_out$fitted.values  

ggplot(quant_filled_no_outlier, aes(x = factor(quant_filled_no_outlier[[indel2]]), y = adjusted_Filled)) +
  geom_boxplot() + geom_jitter(width = 0.2) +
  labs(title = "Adjusted Filled Area by Dose (Family-Adjusted Fitted Values)")



colnames(quant_filled_no_outlier)





out_plot <- ggplot(quant_filled_no_outlier,
  aes(x = factor(quant_filled_no_outlier[[indel2]]),
    y = adjusted_Filled , fill = factor(quant_filled_no_outlier[[indel2]]))) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_jitter(width = 0.1) + 
  labs(
    title = " Plot of Fitted myofiber CSA by INDEL genotype",
    y=expression(paste( "Fitted myofiber Cross Sectional Area ", "(px"^"2",")"))  ,
    x = "Genotypes Chr2 - 43212632") +
  theme_bw() +
  scale_x_discrete(labels = c("0" = "T / T", "1" = "T / TGAGAGAGAGAGA")) +
  theme(legend.position="none")
  
  


ggsave("Sig_indel_violin_plot.svg" , out_plot , width = 6 , height = 6)

ggsave("Sig_indel_violin_plot.png" , out_plot , width = 6 , height = 6 , dpi = 300)



