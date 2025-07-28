# utils/r_runner.py
import subprocess
import tempfile
import os
import threading
import queue
import time
from pathlib import Path
import logging

class RScriptRunner:
    def __init__(self, data_dir="~/Desktop/data", timeout=600):
        self.data_dir = Path(data_dir).expanduser()
        self.timeout = timeout
        self.temp_dir = Path(tempfile.mkdtemp())
        # Use absolute path for results directory relative to current working directory
        self.results_dir = Path.cwd() / "results"
        self.results_dir.mkdir(exist_ok=True)
        
        # Setup logging
        logging.basicConfig(level=logging.INFO)
        self.logger = logging.getLogger(__name__)
        
        # Validate R installation
        self._validate_r_installation()
        
    def _validate_r_installation(self):
        """Check if R is installed and accessible"""
        try:
            result = subprocess.run(['Rscript', '--version'], 
                                  capture_output=True, text=True, timeout=10)
            if result.returncode != 0:
                raise RuntimeError("R is not properly installed or not in PATH")
            self.logger.info("R installation validated successfully")
        except (FileNotFoundError, subprocess.TimeoutExpired) as e:
            raise RuntimeError(f"R is not available: {e}")
        except Exception as e:
            self.logger.warning(f"Could not validate R installation: {e}")
        
    def create_limma_script(self, etiology_group, analysis_name, p_value_threshold=0.05, logfc_threshold=1.0):
        """Create customized limma script"""
        # Get absolute paths
        results_limma_dir = self.results_dir / "limma_results"
        
        template = f'''
# Set working directory to absolute path
setwd("{str(self.data_dir)}")

# Load required libraries
suppressPackageStartupMessages({{
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
}})

# Load data
tryCatch({{
    load("WGCNA_analysis.RData")
    data <- read_NPX("Q-17124_Data Delivery/Q-17124_NPX_2025-03-24.parquet")
    data <- data[!str_detect(data$AssayType, "ctrl"), ]
    metadata <- read.csv("sampling_18022024/olink_sampled_patient_information_all.csv")
}}, error = function(e) {{
    cat("Error loading data:", e$message, "\\n")
    quit(status = 1)
}})

# Data processing
merged_data <- data %>%
    inner_join(metadata, by = "SampleID") 
merged_data$m3_stroke_followup <- as.factor(merged_data$m3_stroke_followup)

data_pivot <- merged_data %>%
    dplyr::select(SampleID, Assay, NPX) %>%
    pivot_wider(names_from = Assay, values_from = NPX) %>%
    filter(rowSums(is.na(.)) == 0)

expr_matrix <- data_pivot %>%
    column_to_rownames("SampleID") %>%
    as.matrix() %>%
    t()

col_data <- metadata %>%
    filter(SampleID %in% colnames(expr_matrix)) %>%
    mutate(
        nih_admission = cut(nih_admission, breaks = c(-Inf, 0, 4, 15, 20, Inf), 
                           labels = c("0", "1-4", "5-15", "16-20", ">20")),
        m3_stroke_followup = factor(m3_stroke_followup),
        sex = factor(sex),
        age = as.numeric(age),
        samplesource = factor(Sample_Source),
        nih_admission = as.factor(nih_admission),
        etiology_toast = as.factor(etiology_toast)
    ) %>%
    column_to_rownames("SampleID")

# Add plate information
plate_info <- data[,c("SampleID", "PlateID")] %>% distinct()
col_data$PlateID <- plate_info$PlateID[match(rownames(col_data), plate_info$SampleID)]
col_data <- col_data %>% mutate(PlateID = as.factor(PlateID))

expr_matrix <- expr_matrix[, rownames(col_data)]
col_data$etiology_toast[is.na(col_data$etiology_toast)] <- 4

# Use WGCNA filtered data
expr_matrix_expr <- expr_matrix[, rownames(wgcna_protein_data)]
col_data_expr <- col_data[rownames(wgcna_protein_data), ]

# Group-specific filtering and design matrix
'''
        
        if etiology_group == "all":
            template += '''
col_data_selected <- col_data_expr
design <- model.matrix(~ m3_stroke_followup + age + sex + nih_admission + etiology_toast, 
                      data = col_data_selected)
'''
        else:
            template += f'''
col_data_selected <- col_data_expr[col_data_expr$etiology_toast == {etiology_group}, ]
if (nrow(col_data_selected) < 10) {{
    cat("Warning: Very few samples in etiology group {etiology_group}\\n")
}}
design <- model.matrix(~ m3_stroke_followup + age, data = col_data_selected)
'''

        template += f'''
# Limma analysis
tryCatch({{
    fit <- lmFit(expr_matrix_expr[, rownames(col_data_selected)], design, method="robust")
    fit <- eBayes(fit)
    
    res <- topTable(fit, coef = "m3_stroke_followup1", number = Inf, sort.by = "P",
                    adjust.method = "BH")
    
    # Process results
    res_df_limma <- res %>%
        rownames_to_column(var = "Assay") %>%
        mutate(
            neg_log10_padj = -log10(adj.P.Val),
            neg_log10_p = -log10(P.Value),
            sig = case_when(
                P.Value < {p_value_threshold} & abs(logFC) > {logfc_threshold} ~ "Significant",
                TRUE ~ "Not significant"
            )
        )
    
    # Add UniProt IDs for enrichment analysis
    assay_uniprot_map <- data %>%
        dplyr::select(Assay, UniProt) %>%
        distinct()
    
    res_df_limma <- res_df_limma %>%
        left_join(assay_uniprot_map, by = "Assay")
    
    # Calculate mean NPX by group
    mean_npx_by_group <- merged_data %>%
        filter(SampleID %in% rownames(col_data_selected)) %>%
        group_by(Assay, m3_stroke_followup) %>%
        summarise(mean_NPX = mean(NPX, na.rm = TRUE), .groups = "drop") %>%
        pivot_wider(names_from = m3_stroke_followup, values_from = mean_NPX,
                    names_prefix = "m3_") %>%
        mutate(delta_NPX = m3_1 - m3_0)
    
    res_df_limma <- res_df_limma %>%
        left_join(mean_npx_by_group %>% dplyr::select(Assay, delta_NPX), by = "Assay")
    
    # Add sample size information
    res_df_limma$n_samples <- nrow(col_data_selected)
    res_df_limma$analysis_group <- "{analysis_name}"
    
    # Ensure results directories exist
    dir.create("{results_limma_dir}", recursive = TRUE, showWarnings = FALSE)
    
    # Save results
    write.csv(res_df_limma, file = "{results_limma_dir}/limma_results_{analysis_name}.csv", row.names = FALSE)
    
    # Save summary statistics
    summary_stats <- list(
        total_proteins = nrow(res_df_limma),
        significant_proteins = sum(res_df_limma$sig == "Significant"),
        upregulated = sum(res_df_limma$sig == "Significant" & res_df_limma$logFC > 0),
        downregulated = sum(res_df_limma$sig == "Significant" & res_df_limma$logFC < 0),
        sample_size = nrow(col_data_selected)
    )
    
    cat("Analysis completed successfully\\n")
    cat("Total proteins analyzed:", summary_stats$total_proteins, "\\n")
    cat("Significant proteins:", summary_stats$significant_proteins, "\\n")
    cat("Sample size:", summary_stats$sample_size, "\\n")
    
}}, error = function(e) {{
    cat("Error in limma analysis:", e$message, "\\n")
    quit(status = 1)
}})
'''
        
        return template
    
    def create_enrichment_script(self, analysis_name, go_pvalue_cutoff=0.05):
        """Create enrichment analysis script with GO, KEGG, and REACTOME"""
        # Get absolute paths
        results_limma_dir = self.results_dir / "limma_results"
        results_enrichment_dir = self.results_dir / "enrichment_results"
        
        return f'''
# Set working directory  
setwd("{self.data_dir}")

# Load required libraries
suppressPackageStartupMessages({{
    library("org.Hs.eg.db", character.only = TRUE)
    library(AnnotationDbi)
    library(clusterProfiler)
    library(enrichplot)
    library(ggplot2)
    library(dplyr)
    library(OlinkAnalyze)
    library(ReactomePA)
}})

# Load limma results
tryCatch({{
    res_df_limma <- read.csv("{results_limma_dir}/limma_results_{analysis_name}.csv", 
                             stringsAsFactors = FALSE)
    
    # Check if UniProt column exists
    if (!"UniProt" %in% colnames(res_df_limma)) {{
        cat("Error: UniProt column not found in limma results\\n")
        quit(status = 1)
    }}
    
    cat("Loaded limma results with", nrow(res_df_limma), "proteins\\n")
    
}}, error = function(e) {{
    cat("Error loading limma results:", e$message, "\\n")
    quit(status = 1)
}})

# Remove rows without UniProt IDs
res_df_limma <- res_df_limma[!is.na(res_df_limma$UniProt) & res_df_limma$UniProt != "", ]
cat("After filtering, have", nrow(res_df_limma), "proteins with UniProt IDs\\n")

# Create ranked gene list for UniProt IDs
protein_list <- res_df_limma$logFC
names(protein_list) <- res_df_limma$UniProt
protein_list <- na.omit(protein_list)
protein_list <- sort(protein_list, decreasing = TRUE)

cat("Prepared UniProt gene list with", length(protein_list), "proteins\\n")

# Ensure results directories exist
dir.create("{results_enrichment_dir}", recursive = TRUE, showWarnings = FALSE)

# Function to safely perform GO GSEA
safe_go_gsea <- function(ont_type, pval_cutoff = 0.05) {{
    cat("Starting GO", ont_type, "enrichment analysis with p-value cutoff:", pval_cutoff, "\\n")
    tryCatch({{
        result <- gseGO(
            geneList = protein_list,
            ont = ont_type,
            keyType = "UNIPROT", 
            pAdjustMethod = "hochberg",
            verbose = FALSE,
            pvalueCutoff = pval_cutoff,
            OrgDb = "org.Hs.eg.db"
        )
        
        cat("GO GSEA completed for", ont_type, "with", nrow(as.data.frame(result)), "results\\n")
        
        if (nrow(as.data.frame(result)) > 0) {{
            df_result <- as.data.frame(result)
            df_result$analysis_group <- "{analysis_name}"
            write.csv(df_result, 
                     file = paste0("{results_enrichment_dir}/enrichment_", 
                                  tolower(ont_type), "_{analysis_name}.csv"), 
                     row.names = FALSE)
            cat("Saved GO", ont_type, "results with", nrow(df_result), "terms\\n")
            return(TRUE)
        }} else {{
            cat("No significant GO", ont_type, "terms found\\n")
            return(FALSE)
        }}
    }}, error = function(e) {{
        cat("Error in GO", ont_type, "enrichment:", e$message, "\\n")
        return(FALSE)
    }})
}}

# Function to safely perform KEGG GSEA
safe_kegg_gsea <- function(pval_cutoff = 0.05) {{
    cat("Starting KEGG enrichment analysis with p-value cutoff:", pval_cutoff, "\\n")
    
    tryCatch({{
        result <- gseKEGG(
            geneList = protein_list,
            organism = "hsa",
            pAdjustMethod = "hochberg",
            verbose = FALSE,
            pvalueCutoff = pval_cutoff,
            keyType = "uniprot"
        )
        
        cat("KEGG GSEA completed with", nrow(as.data.frame(result)), "results\\n")
        
        if (nrow(as.data.frame(result)) > 0) {{
            df_result <- as.data.frame(result)
            df_result$analysis_group <- "{analysis_name}"
            write.csv(df_result, 
                     file = "{results_enrichment_dir}/enrichment_kegg_{analysis_name}.csv",
                     row.names = FALSE)
            cat("Saved KEGG results with", nrow(df_result), "pathways\\n")
            return(TRUE)
        }} else {{
            cat("No significant KEGG pathways found\\n")
            return(FALSE)
        }}
    }}, error = function(e) {{
        cat("Error in KEGG enrichment:", e$message, "\\n")
        return(FALSE)
    }})
}}

# Function to safely perform REACTOME GSEA
safe_reactome_gsea <- function(pval_cutoff = 0.05) {{
    cat("Starting REACTOME enrichment analysis with p-value cutoff:", pval_cutoff, "\\n")
    
    tryCatch({{
        # Convert UniProt to Entrez IDs using bitr
        protein_list_reactome <- protein_list
        entrez_conversion <- bitr(names(protein_list_reactome), 
                                 fromType = "UNIPROT", 
                                 toType = "ENTREZID", 
                                 OrgDb = "org.Hs.eg.db")
        
        if (nrow(entrez_conversion) == 0) {{
            cat("No Entrez IDs available for REACTOME analysis\\n")
            return(FALSE)
        }}
        
        # Map the names and remove NAs/duplicates
        names(protein_list_reactome) <- entrez_conversion$ENTREZID[match(names(protein_list_reactome), entrez_conversion$UNIPROT)]
        protein_list_reactome <- protein_list_reactome[!is.na(names(protein_list_reactome))]
        protein_list_reactome <- protein_list_reactome[!duplicated(names(protein_list_reactome))]
        
        cat("Prepared", length(protein_list_reactome), "genes for REACTOME analysis\\n")
        
        result <- gsePathway(
            geneList = protein_list_reactome,
            pAdjustMethod = "hochberg",
            verbose = FALSE,
            pvalueCutoff = pval_cutoff
        )
        
        cat("REACTOME GSEA completed with", nrow(as.data.frame(result)), "results\\n")
        
        if (nrow(as.data.frame(result)) > 0) {{
            df_result <- as.data.frame(result)
            df_result$analysis_group <- "{analysis_name}"
            write.csv(df_result, 
                     file = "{results_enrichment_dir}/enrichment_reactome_{analysis_name}.csv",
                     row.names = FALSE)
            cat("Saved REACTOME results with", nrow(df_result), "pathways\\n")
            return(TRUE)
        }} else {{
            cat("No significant REACTOME pathways found\\n")
            return(FALSE)
        }}
    }}, error = function(e) {{
        cat("Error in REACTOME enrichment:", e$message, "\\n")
        return(FALSE)
    }})
}}

# Perform enrichment analyses
bp_success <- safe_go_gsea("BP", {go_pvalue_cutoff})
mf_success <- safe_go_gsea("MF", {go_pvalue_cutoff})  
cc_success <- safe_go_gsea("CC", {go_pvalue_cutoff})
kegg_success <- safe_kegg_gsea({go_pvalue_cutoff})
reactome_success <- safe_reactome_gsea({go_pvalue_cutoff})

# Create summary
summary_enrichment <- data.frame(
    analysis_group = "{analysis_name}",
    bp_terms = ifelse(bp_success, "Available", "No significant terms"),
    mf_terms = ifelse(mf_success, "Available", "No significant terms"), 
    cc_terms = ifelse(cc_success, "Available", "No significant terms"),
    kegg_pathways = ifelse(kegg_success, "Available", "No significant pathways"),
    reactome_pathways = ifelse(reactome_success, "Available", "No significant pathways"),
    total_proteins = length(protein_list)
)

write.csv(summary_enrichment, 
          file = "{results_enrichment_dir}/enrichment_summary_{analysis_name}.csv",
          row.names = FALSE)

cat("Enrichment analysis completed for {analysis_name}\\n")
cat("Results saved: GO-BP =", bp_success, ", GO-MF =", mf_success, ", GO-CC =", cc_success, "\\n")
cat("KEGG =", kegg_success, ", REACTOME =", reactome_success, "\\n")
'''

    def run_script(self, script_content, script_name, timeout=None):
        """Execute R script with proper error handling"""
        if timeout is None:
            timeout = self.timeout
            
        script_path = self.temp_dir / script_name
        
        # Write script to temporary file
        with open(script_path, 'w') as f:
            f.write(script_content)
        
        # Create results directories
        (self.results_dir / "limma_results").mkdir(exist_ok=True)
        (self.results_dir / "enrichment_results").mkdir(exist_ok=True)
        
        try:
            # Run R script
            result = subprocess.run(
                ['Rscript', str(script_path)],
                cwd=str(self.data_dir),
                capture_output=True,
                text=True,
                timeout=timeout
            )
            
            self.logger.info(f"Script {script_name} completed with return code: {result.returncode}")
            
            if result.returncode == 0:
                return True, result.stdout
            else:
                self.logger.error(f"Script failed: {result.stderr}")
                return False, result.stderr
                
        except subprocess.TimeoutExpired:
            error_msg = f"Script {script_name} timed out after {timeout} seconds"
            self.logger.error(error_msg)
            return False, error_msg
        except Exception as e:
            error_msg = f"Unexpected error running {script_name}: {str(e)}"
            self.logger.error(error_msg)
            return False, error_msg
    
    def run_limma_analysis(self, etiology_group, analysis_name, p_value_threshold=0.05, logfc_threshold=1.0):
        """Run limma analysis for specific group"""
        script_content = self.create_limma_script(etiology_group, analysis_name, p_value_threshold, logfc_threshold)
        return self.run_script(script_content, f"limma_{analysis_name}.R")
    
    def run_enrichment_analysis(self, analysis_name, go_pvalue_cutoff=0.05):
        """Run enrichment analysis"""
        script_content = self.create_enrichment_script(analysis_name, go_pvalue_cutoff)
        return self.run_script(script_content, f"enrichment_{analysis_name}.R")
    
    def cleanup(self):
        """Clean up temporary files"""
        import shutil
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)


