setwd("~/Desktop/data")
library(OlinkAnalyze)
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(tidyverse)
library(limma)
library(lme4)
library(readr)
library(janitor)
library(ggrepel)
library(clusterProfiler)
require(DOSE)
load("WGCNA_analysis.RData")

data <- read_NPX("~/Desktop/data/Q-17124_Data Delivery/Q-17124_NPX_2025-03-24.parquet")
data <-  data[!str_detect(data$AssayType, "ctrl"), ]
metadata <- read.csv("~/Desktop/data/sampling_18022024/olink_sampled_patient_information_all.csv")

merged_data <- data %>%
  inner_join(metadata, by = "SampleID") 
merged_data$m3_stroke_followup <- as.factor(merged_data$m3_stroke_followup)

data_pivot <- merged_data %>%
  dplyr::select(SampleID, Assay, NPX) %>%
  pivot_wider(
    names_from = Assay,
    values_from = NPX
  ) %>% # remove NA rows
  filter(rowSums(is.na(.)) == 0)

# Get expression matrix
expr_matrix <- data_pivot %>%
  column_to_rownames("SampleID") %>%
  as.matrix() %>%
  t()

col_data <- metadata %>%
  filter(SampleID %in% colnames(expr_matrix)) %>%
  # bin NIHSS score into categories
  mutate(nih_admission = cut(nih_admission, breaks = c(-Inf, 0, 4, 15, 20, Inf), labels = c("0", "1-4", "5-15", "16-20", ">20"))) %>%
  mutate(
    m3_stroke_followup = factor(m3_stroke_followup),
    sex = factor(sex),        # Convert sex to factor
    age = as.numeric(age),     # Ensure age is numeric
    samplesource = factor(Sample_Source),
    nih_admission = as.factor(nih_admission), # FIXME: keep as numeric or categorical? 
    etiology_toast = as.factor(etiology_toast), 
  ) %>%
  column_to_rownames("SampleID")
# add in information about plate number
plate_info <- data[,c("SampleID", "PlateID")] %>% distinct()
col_data$PlateID <- plate_info$PlateID[match(rownames(col_data), plate_info$SampleID)]
col_data <- col_data %>%
  mutate(
    PlateID = as.factor(PlateID)
  )

expr_matrix <- expr_matrix[, rownames(col_data)]

col_data$etiology_toast[is.na(col_data$etiology_toast)] <- 4

#
expr_matrix_expr <- expr_matrix[, rownames(wgcna_protein_data)] # can use them if you want to eliminate the "outliers" based on hierarchical clustering from WGCNA
col_data_expr <- col_data[rownames(wgcna_protein_data), ]




col_data_expr_5 <- col_data_expr[col_data_expr$etiology_toast==5,]
col_data_expr_1 <- col_data_expr[col_data_expr$etiology_toast==1,]
col_data_expr_2 <- col_data_expr[col_data_expr$etiology_toast==2,]
col_data_expr_3 <- col_data_expr[col_data_expr$etiology_toast==3,]
col_data_expr_4 <- col_data_expr[col_data_expr$etiology_toast==4,]

design <- model.matrix(~ m3_stroke_followup + age + sex + nih_admission + etiology_toast, data = col_data_expr) # left out sample, data = col_data) # left out sample source
design <- model.matrix(~ m3_stroke_followup + age, data = col_data_expr_1)

fit <- lmFit(expr_matrix_expr[, rownames(col_data_expr_1)], design, method="robust") # ,rownames(col_data_expr_5)
fit <- eBayes(fit)

# Extract results for m3_stroke_followup
res <- topTable(fit, coef = "m3_stroke_followup1", number = Inf, sort.by = "P",
                adjust.method="BH")

res_df_limma <- res %>%
  rownames_to_column(var = "Assay") %>%
  mutate(
    neg_log10_padj = -log10(adj.P.Val),
    neg_log10_p = -log10(P.Value),
    sig = case_when(
      P.Value < 0.05 & abs(logFC) > 1.0 ~ "Significant",
      TRUE ~ "Not significant"
    )
  )

# Compute average NPX by group
mean_npx_by_group <- merged_data %>%
  group_by(Assay, m3_stroke_followup) %>%
  summarise(mean_NPX = mean(NPX, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = m3_stroke_followup, values_from = mean_NPX,
              names_prefix = "m3_") %>%
  mutate(
    delta_NPX = m3_1 - m3_0
  )
res_df_limma <- res_df_limma %>%
  left_join(mean_npx_by_group %>% dplyr::select(Assay, delta_NPX), by = "Assay")
top_hits <- res_df_limma %>%
  filter(sig == "Significant") %>%
  arrange(desc(logFC)) %>%
  slice_head(n = 5) %>%
  bind_rows(
    res_df_limma %>%
      filter(sig == "Significant") %>%
      arrange(logFC) %>%
      slice_head(n = 5)
  )
# drop duplicates in top_hits
top_hits <- top_hits %>%
  distinct(Assay, .keep_all = TRUE)

# Volcano plot
ggplot(res_df_limma, aes(x = logFC, y = neg_log10_p)) +
  geom_point(aes(color = sig), alpha = 0.6) +
  #geom_text(data = top_hits, aes(label = Assay), vjust = 1, hjust = 1, size = 3) +
  geom_text_repel(data = top_hits, aes(label = Assay), size = 5, max.overlaps = 20) +
  scale_color_manual(values = c("Significant" = "red", "Not significant" = "grey")) +
  # geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  labs(
    title = "Month 3 Recurrent Patients vs. Non-Recurrent",
    x = "Log2 Fold Change (NPX)",
    # x = "Delta NPX",
    y = "-Log10 Unadjusted P-value",
    color = "Significance"
  ) +
  theme_minimal() + xlim(-5, 5) + theme(axis.text=element_text(size=12))

write.csv(res_df_limma, file="limma_DE_results.csv")