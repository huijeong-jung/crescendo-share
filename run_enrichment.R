



library("org.Hs.eg.db", character.only = TRUE)

# add the uniprot ids to res_df_limma
for (assay in res_df_limma$Assay) {
  # get the uniprot id from the assay name
  uniprot_id <- data[data$Assay == assay, "UniProt"] %>%
    unique() 
  res_df_limma$UniProt[res_df_limma$Assay == assay] <- uniprot_id
}
res_df_limma <- as.data.frame(res_df_limma, stringsAsFactors = FALSE)
res_df_limma_ <- apply(res_df_limma, 2, as.character)

protein_list <- res_df_limma$logFC
names(protein_list) <- res_df_limma$UniProt
protein_list <- na.omit(protein_list)
protein_list <- sort(protein_list, decreasing = TRUE) 
gse_all <- gseGO(geneList = protein_list,
                 ont = "ALL", 
                 keyType = "UNIPROT",
                 pAdjustMethod = "hochberg",
                 verbose=TRUE,
                 pvalueCutoff = 0.05,
                 OrgDb = "org.Hs.eg.db")


gse_cc <- gseGO(geneList = protein_list,
                ont = "ALL", 
                keyType = "UNIPROT",
                pAdjustMethod = "hochberg",
                verbose=TRUE,
                pvalueCutoff = 0.05,
                OrgDb = "org.Hs.eg.db")
gse_bp <- gseGO(geneList = protein_list,
                ont = "BP",
                keyType = "UNIPROT",
                pAdjustMethod = "hochberg",
                verbose=TRUE,
                pvalueCutoff = 0.3,
                OrgDb = "org.Hs.eg.db")
gse_mf <- gseGO(geneList = protein_list,
                ont = "MF",
                keyType = "UNIPROT",
                pAdjustMethod = "hochberg",
                verbose=TRUE,
                pvalueCutoff = 0.05,
                OrgDb = "org.Hs.eg.db")


dotplot(gse_cc, showCategory=25) +
  theme_minimal() +theme(axis.text.y = element_text(size=10, color="black"), axis.text.x = element_text(size=13, color="black")) +
  ggtitle("Top Enriched GO Terms")  
dotplot(gse_bp, showCategory=25) +
  theme_minimal() +theme(axis.text.y = element_text(size=10, color="black"), axis.text.x = element_text(size=13, color="black")) +
  ggtitle("Top Enriched GO Terms") 
dotplot(gse_mf, showCategory=25) +
  theme_minimal() +theme(axis.text.y = element_text(size=10, color="black"), axis.text.x = element_text(size=13, color="black")) +
  ggtitle("Top Enriched GO Terms") 


gse_bp <- as.data.frame(gse_bp)
term <- gse_bp[gse_bp$Description == "protein-RNA complex organization",]
coregenes <- strsplit(term$core_enrichment, "/")[[1]]
# select out these genes from res_df_limma_, convert them to "Assay" column
res_df_limma_ <- as.data.frame(res_df_limma_)
coregenes <- res_df_limma_[res_df_limma_$UniProt %in% coregenes, "Assay"]


#ggsave("dotplot_GO_terms.png", width = 8, height = 9)
# make plot bigger
ridgeplot(gse, showCategory = 10) + labs(x = "enrichment distribution") + 
  theme(axis.text.y = element_text(size = 9)) +
  ggtitle("Distribution of Enrichment Scores") + theme_minimal()
library(enrichplot)
# gse_simple <- simplify(gse, by = "p.adjust", select_fun = min) # not using the gse_simple for now
emapplot(pairwise_termsim(gse), layout="kk",
         showCategory = 30, 
         node_label = "category") + theme(axis.text.r = element_text(size=))