# 🚀 Deployment Guide for Proteomics Dashboard

## Overview
This Streamlit app provides interactive visualization of proteomics differential expression and pathway enrichment results for stroke etiology analysis.

## 📊 Pre-computed Results Included
- **Limma Results**: Differential expression analysis for 6 stroke etiology groups
- **GO Enrichment**: Biological Process, Molecular Function, Cellular Component
- **KEGG Pathways**: Metabolic and signaling pathway enrichment
- **REACTOME**: Pathway analysis results
- **Interactive Thresholds**: P-value < 0.05, |logFC| > 1.0

## 🌐 Deployment Options

### Option 1: Streamlit Community Cloud (Recommended)

1. **Create GitHub Repository**:
   ```bash
   git commit -m "Add proteomics dashboard with pre-computed results"
   git remote add origin https://github.com/yourusername/proteomics-dashboard.git
   git push -u origin main
   ```

2. **Deploy**:
   - Visit [share.streamlit.io](https://share.streamlit.io)
   - Connect your GitHub account
   - Select your repository
   - Set main file to `deploy_app.py`
   - Deploy!

3. **Features**:
   - ✅ Free hosting
   - ✅ Auto-deploys on git push
   - ✅ Custom domain available
   - ✅ HTTPS included

### Option 2: Docker Deployment

1. **Build and Run**:
   ```bash
   docker-compose up --build
   ```

2. **Access**:
   - Local: http://localhost:8501
   - Production: Configure reverse proxy

### Option 3: Local Sharing

```bash
streamlit run deploy_app.py --server.port 8501
# Then use ngrok or similar for public access
```

## 📁 Repository Structure
```
├── deploy_app.py          # Deployment-ready app (no R dependencies)
├── app.py                 # Full app (requires R)
├── requirements.txt       # Python dependencies
├── results/              # Pre-computed analysis results
│   ├── limma_results/    # Differential expression
│   └── enrichment_results/ # Pathway enrichment
├── utils/                # Helper modules
└── .streamlit/           # Streamlit configuration
```

## 🔧 Configuration

### Streamlit Settings (`.streamlit/config.toml`):
```toml
[server]
maxUploadSize = 200

[theme]
primaryColor = "#FF6B35"
backgroundColor = "#FFFFFF"
secondaryBackgroundColor = "#F0F2F6"
```

## 🚀 Quick Start for Users

1. **No Installation Required**: Access the deployed app directly
2. **Select Stroke Etiology**: Choose from dropdown menu
3. **Adjust Thresholds**: Interactive P-value and logFC controls
4. **Explore Results**: 
   - Volcano plots with dynamic thresholds
   - GO/KEGG/REACTOME enrichment tables
   - Downloadable results

## 💡 Key Features

- **Read-Only Mode**: No R computation required for deployment
- **Interactive Visualization**: Real-time threshold adjustment
- **Pre-computed Results**: All analyses already completed
- **Mobile Friendly**: Responsive design
- **Fast Loading**: CSV-based data loading

## 🔒 Data Privacy

- Raw data (`data/` folder) excluded from repository
- Only aggregated analysis results included
- No patient-identifiable information

## 📞 Support

For questions about the analysis methods or results, please contact the research team.
