#!/bin/bash

# Local development script
export STREAMLIT_SERVER_HEADLESS=false
export DATA_DIR="../data"

# Create directories if they don't exist
mkdir -p results/limma_results
mkdir -p results/enrichment_results

# Install Python dependencies
pip install -r requirements.txt

# Run Streamlit app
streamlit run app.py --server.port 8501
