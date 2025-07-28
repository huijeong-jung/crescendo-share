# Proteomics Dashboard - Stroke Etiology Analysis

An interactive Streamlit dashboard for exploring differential expression and pathway enrichment results from proteomics analysis of stroke etiologies.

## 🎯 Features

- **Interactive Volcano Plots**: Dynamic P-value and logFC thresholds
- **Pathway Enrichment**: GO (BP/MF/CC), KEGG, and REACTOME analysis
- **Multiple Stroke Types**: 6 different etiology classifications
- **Real-time Filtering**: Adjust significance thresholds instantly
- **Downloadable Results**: Export filtered results as CSV

## 📊 Pre-computed Results Included

This repository contains pre-computed analysis results, so users can explore findings without needing to install R or run computationally intensive analyses locally.

### Available Analyses:
- All Groups (Holistic)
- Cardioembolism
- Large Artery Atherosclerosis  
- Small Vessel Occlusion
- Stroke of Other Etiology
- Stroke of Undetermined Etiology

## 🚀 Quick Access

**Deployed App**: [Link will be available after deployment]

## 🛠️ Local Development

```bash
# Clone repository
git clone https://github.com/yourusername/proteomics-dashboard.git
cd proteomics-dashboard

# Install dependencies
pip install -r requirements.txt

# Run app
streamlit run deploy_app.py
```

## 📈 Analysis Methods

- **Differential Expression**: Limma with robust fitting
- **Significance Thresholds**: P-value < 0.05, |logFC| > 1.0
- **Pathway Analysis**: clusterProfiler (GO), KEGG.db, ReactomePA
- **Multiple Testing**: Benjamini-Hochberg correction

## 📁 Repository Structure

```
├── deploy_app.py              # Main Streamlit application
├── results/
│   ├── limma_results/         # Differential expression results
│   └── enrichment_results/    # Pathway enrichment results
├── utils/                     # Helper functions
├── requirements.txt           # Python dependencies
└── README.md                 # This file
```

## 🔬 Citation

If you use these results in your research, please cite:
[Your publication details here]

## 📞 Contact

For questions about the analysis or methodology, please contact:
[Your contact information]
