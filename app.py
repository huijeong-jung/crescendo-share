import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import os
import json
from pathlib import Path
import time
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
import logging

# Import custom utilities (assuming they're in utils/ directory)
try:
    from utils.r_runner import RScriptRunner
    from utils.plotting import PlotGenerator
except ImportError:
    st.error("Could not import utilities. Please ensure utils/ directory contains r_runner.py and plotting.py")
    st.stop()

# Configure page
st.set_page_config(
    page_title="🧬 Proteomics Analysis Dashboard",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for better styling
st.markdown("""
<style>
    .main-header {
        font-size: 2.5rem;
        font-weight: 700;
        color: #1f77b4;
        text-align: center;
        margin-bottom: 2rem;
    }
    .analysis-card {
        background-color: #f8f9fa;
        padding: 1rem;
        border-radius: 10px;
        border-left: 4px solid #1f77b4;
        margin: 1rem 0;
    }
    .metric-card {
        background-color: #e3f2fd;
        padding: 1rem;
        border-radius: 8px;
        text-align: center;
    }
    .success-box {
        background-color: #d4edda;
        border: 1px solid #c3e6cb;
        color: #155724;
        padding: 1rem;
        border-radius: 5px;
        margin: 0.5rem 0;
    }
    .error-box {
        background-color: #f8d7da;
        border: 1px solid #f5c6cb;
        color: #721c24;
        padding: 1rem;
        border-radius: 5px;
        margin: 0.5rem 0;
    }
</style>
""", unsafe_allow_html=True)

class AnalysisDashboard:
    def __init__(self):
        self.analysis_groups = {
            "All Groups (Holistic)": "all",
            "Large Artery Atherosclerosis": "1", 
            "Cardioembolism": "2",
            "Small Vessel Occlusion": "3", 
            "Stroke of Other Etiology": "4",
            "Stroke of Undetermined Etiology": "5"
        }
        self.plotter = PlotGenerator()
        
        # Initialize session state
        if 'analysis_results' not in st.session_state:
            st.session_state.analysis_results = {}
        if 'analysis_running' not in st.session_state:
            st.session_state.analysis_running = False
            
    def create_sidebar(self):
        """Create sidebar for configuration"""
        st.sidebar.header("🔧 Analysis Configuration")
        
        # Data directory
        data_dir = st.sidebar.text_input(
            "📁 Data Directory Path:",
            value="~/Desktop/data/proteomics/share_results/data",
            help="Path to directory containing your data files"
        )
        
        # Analysis selection
        selected_analyses = st.sidebar.multiselect(
            "📊 Select Analysis Groups:",
            options=list(self.analysis_groups.keys()),
            default=["All Groups (Holistic)"],
            help="Choose which etiology groups to analyze"
        )
        
        # GO/Enrichment analysis settings
        st.sidebar.markdown("### 🧬 GO Analysis Settings")
        enable_go_analysis = st.sidebar.checkbox(
            "Run GO Term Enrichment", 
            value=True,
            help="Perform Gene Ontology enrichment analysis"
        )
        
        go_pvalue_cutoff = st.sidebar.number_input(
            "GO P-value cutoff:",
            min_value=0.01, max_value=0.3, value=0.05, step=0.01,
            help="P-value cutoff for GO enrichment"
        )
        
        # Advanced options
        with st.sidebar.expander("⚙️ Advanced Options"):
            timeout = st.number_input("Script Timeout (seconds)", 
                                    min_value=60, max_value=1800, value=300)
            parallel_execution = st.checkbox("Run analyses in parallel", value=True)
            
            st.markdown("#### 🗂️ Cache Management")
            clear_session_cache = st.button("🗑️ Clear Session Cache")
            clear_file_cache = st.button("🗂️ Clear All Result Files", 
                                       help="Delete all saved analysis result files")
            force_rerun = st.checkbox("🔄 Force re-run analyses (ignore cached results)", 
                                    help="Re-run analyses even if results already exist")
            
            if clear_session_cache:
                st.session_state.analysis_results = {}
                st.success("Session cache cleared!")
                
            if clear_file_cache:
                import shutil
                try:
                    if Path("results/limma_results").exists():
                        for file in Path("results/limma_results").glob("*.csv"):
                            file.unlink()
                    if Path("results/enrichment_results").exists():
                        for file in Path("results/enrichment_results").glob("*.csv"):
                            file.unlink()
                    st.session_state.analysis_results = {}
                    st.success("All result files cleared!")
                except Exception as e:
                    st.error(f"Error clearing files: {e}")
                
            # Show existing result files
            if st.button("📁 Show Cached Files"):
                limma_files = list(Path("results/limma_results").glob("*.csv")) if Path("results/limma_results").exists() else []
                enrichment_files = list(Path("results/enrichment_results").glob("*.csv")) if Path("results/enrichment_results").exists() else []
                
                if limma_files or enrichment_files:
                    st.info(f"Found {len(limma_files)} limma files and {len(enrichment_files)} enrichment files")
                else:
                    st.info("No cached result files found")
        
        return data_dir, selected_analyses, timeout, parallel_execution, enable_go_analysis, go_pvalue_cutoff, force_rerun
    
    def validate_data_files(self, data_dir):
        """Validate that required data files exist"""
        data_path = Path(data_dir).expanduser()
        
        required_files = [
            "WGCNA_analysis.RData",
            "Q-17124_Data Delivery/Q-17124_NPX_2025-03-24.parquet",
            "sampling_18022024/olink_sampled_patient_information_all.csv"
        ]
        
        missing_files = []
        for file in required_files:
            if not (data_path / file).exists():
                missing_files.append(file)
        
        return missing_files
    
    def check_existing_results(self, analysis_key):
        """Check if analysis results already exist"""
        limma_path = Path("results") / "limma_results" / f"limma_results_{analysis_key}.csv"
        enrichment_paths = {
            'bp': Path("results") / "enrichment_results" / f"enrichment_bp_{analysis_key}.csv",
            'mf': Path("results") / "enrichment_results" / f"enrichment_mf_{analysis_key}.csv", 
            'cc': Path("results") / "enrichment_results" / f"enrichment_cc_{analysis_key}.csv",
            'kegg': Path("results") / "enrichment_results" / f"enrichment_kegg_{analysis_key}.csv",
            'reactome': Path("results") / "enrichment_results" / f"enrichment_reactome_{analysis_key}.csv",
            'summary': Path("results") / "enrichment_results" / f"enrichment_summary_{analysis_key}.csv"
        }
        
        limma_exists = limma_path.exists()
        
        # Check if specific enrichment types exist
        go_exists = all(enrichment_paths[ont].exists() for ont in ['bp', 'mf', 'cc'])
        kegg_exists = enrichment_paths['kegg'].exists()
        reactome_exists = enrichment_paths['reactome'].exists()
        
        # Consider enrichment complete only if ALL expected files exist
        enrichment_exists = go_exists and kegg_exists and reactome_exists
        
        return {
            'limma_exists': limma_exists,
            'enrichment_exists': enrichment_exists,
            'go_exists': go_exists,
            'kegg_exists': kegg_exists, 
            'reactome_exists': reactome_exists,
            'limma_path': limma_path,
            'enrichment_paths': enrichment_paths
        }
    
    def get_results_cache_key(self, analysis_key, p_value_threshold, logfc_threshold, go_pvalue_cutoff):
        """Generate cache key based on analysis parameters"""
        return f"{analysis_key}_p{p_value_threshold}_lfc{logfc_threshold}_go{go_pvalue_cutoff}"
    
    def display_r_status(self):
        """Display R environment status"""
        st.sidebar.markdown("---")
        st.sidebar.markdown("### 🔧 Environment Status")
        
        try:
            # Test R runner creation
            test_runner = RScriptRunner(data_dir="./data", timeout=10)
            st.sidebar.success("✅ R Environment Ready")
            st.sidebar.info(f"📁 Results Directory: {test_runner.results_dir}")
            test_runner.cleanup()
        except Exception as e:
            st.sidebar.error("❌ R Environment Issue")
            st.sidebar.error(f"Error: {str(e)}")
            with st.sidebar.expander("Troubleshooting"):
                st.markdown("""
                **Common Issues:**
                1. R is not installed or not in PATH
                2. Required R packages are missing
                3. Data directory permissions
                
                **Quick Fix:**
                ```bash
                # Install R packages
                R -e "install.packages(c('OlinkAnalyze', 'dplyr', 'ggplot2'))"
                ```
                """)
    
    def run_single_analysis(self, analysis_name, etiology_group, data_dir, timeout, 
                           p_value_threshold=0.05, logfc_threshold=1.0, 
                           enable_go_analysis=True, go_pvalue_cutoff=0.05, force_rerun=False):
        """Run analysis for a single group"""
        analysis_key = analysis_name.lower().replace(' ', '_').replace('(', '').replace(')', '')
        
        try:
            # Check if results already exist (but respect force_rerun flag)
            existing_results = self.check_existing_results(analysis_key)
            
            # Check if we can use cached results
            use_cached_limma = existing_results['limma_exists'] and not force_rerun
            use_cached_enrichment = existing_results['enrichment_exists'] and enable_go_analysis and not force_rerun
            
            if use_cached_limma:
                st.info(f"📂 Found existing limma results for {analysis_name} - loading from cache")
            elif existing_results['limma_exists'] and force_rerun:
                st.info(f"🔄 Found existing limma results for {analysis_name} - re-running due to force flag")
            
            if use_cached_enrichment:
                st.info(f"📂 Found existing GO enrichment results for {analysis_name} - loading from cache")
            elif existing_results['enrichment_exists'] and force_rerun and enable_go_analysis:
                st.info(f"🔄 Found existing GO enrichment results for {analysis_name} - re-running due to force flag")
            
            # Initialize runner only if we need to run analyses
            runner = None
            if not use_cached_limma or (enable_go_analysis and not use_cached_enrichment):
                # Expand the data directory path
                data_dir_expanded = str(Path(data_dir).expanduser().resolve())
                runner = RScriptRunner(data_dir=data_dir_expanded, timeout=timeout)
                
                # Debug output
                st.info(f"🔍 Running analysis: {analysis_name} (key: {analysis_key})")
                st.info(f"📁 Data directory: {data_dir_expanded}")
                st.info(f"💾 Results will be saved to: {runner.results_dir}")
                st.info(f"🎯 Thresholds: P-value < {p_value_threshold}, |logFC| > {logfc_threshold}")
            
            # Run or load Limma analysis
            limma_success = False
            limma_output = ""
            
            if use_cached_limma:
                limma_success = True
                limma_output = "Loaded from existing results"
                st.success(f"✅ Limma results loaded from cache for {analysis_name}")
            else:
                # Run Limma analysis with custom thresholds
                limma_success, limma_output = runner.run_limma_analysis(
                    etiology_group, analysis_key, p_value_threshold, logfc_threshold
                )
                
                if not limma_success:
                    st.error(f"❌ Limma analysis failed for {analysis_name}")
                    with st.expander("View Error Details"):
                        st.code(limma_output, language="text")
                    return {
                        'name': analysis_name,
                        'key': analysis_key,
                        'status': 'failed',
                        'error': f"Limma analysis failed: {limma_output}",
                        'limma_success': False,
                        'enrichment_success': False
                    }
            
            # Run or load Enrichment analysis (GO terms) if enabled
            enrichment_success = False
            enrichment_output = "GO analysis disabled"
            
            if enable_go_analysis:
                if use_cached_enrichment:
                    enrichment_success = True
                    enrichment_output = "Loaded from existing results"
                    st.success(f"✅ GO enrichment results loaded from cache for {analysis_name}")
                else:
                    # Run GO enrichment analysis
                    enrichment_success, enrichment_output = runner.run_enrichment_analysis(
                        analysis_key, go_pvalue_cutoff
                    )
                    
                    # Show enrichment status
                    if not enrichment_success:
                        st.warning(f"⚠️ GO enrichment analysis failed for {analysis_name}")
                        with st.expander("View GO Enrichment Error Details"):
                            st.code(enrichment_output, language="text")
            else:
                st.info(f"🔄 GO analysis skipped for {analysis_name} (disabled in settings)")
            
            result = {
                'name': analysis_name,
                'key': analysis_key,
                'status': 'completed' if limma_success and (enrichment_success or not enable_go_analysis) else 'partial',
                'limma_success': limma_success,
                'enrichment_success': enrichment_success,
                'limma_output': limma_output,
                'enrichment_output': enrichment_output,
                'go_analysis_enabled': enable_go_analysis
            }
            
            # Load results if successful
            if limma_success:
                result['limma_data'] = self.load_limma_results(analysis_key)
            if enrichment_success:
                result['enrichment_data'] = self.load_enrichment_results(analysis_key)
            
            # Cleanup runner if it was created
            if runner is not None:
                runner.cleanup()
            
            return result
            
        except Exception as e:
            return {
                'name': analysis_name,
                'key': analysis_key,
                'status': 'error',
                'error': str(e),
                'limma_success': False,
                'enrichment_success': False
            }
    
    def load_limma_results(self, analysis_key):
        """Load limma results from CSV"""
        try:
            results_path = Path("results") / "limma_results" / f"limma_results_{analysis_key}.csv"
            if results_path.exists():
                return pd.read_csv(results_path)
        except Exception as e:
            st.error(f"Error loading limma results for {analysis_key}: {e}")
        return None
    
    def load_enrichment_results(self, analysis_key):
        """Load enrichment results from CSV files"""
        results = {}
        try:
            enrichment_dir = Path("results") / "enrichment_results"
            
            # Load GO terms
            for ont in ['bp', 'mf', 'cc']:
                file_path = enrichment_dir / f"enrichment_{ont}_{analysis_key}.csv"
                if file_path.exists():
                    df = pd.read_csv(file_path)
                    results[ont] = df
            
            # Load KEGG pathways
            kegg_path = enrichment_dir / f"enrichment_kegg_{analysis_key}.csv"
            if kegg_path.exists():
                df = pd.read_csv(kegg_path)
                results['kegg'] = df
            
            # Load REACTOME pathways
            reactome_path = enrichment_dir / f"enrichment_reactome_{analysis_key}.csv"
            if reactome_path.exists():
                df = pd.read_csv(reactome_path)
                results['reactome'] = df
            
            # Load summary
            summary_path = enrichment_dir / f"enrichment_summary_{analysis_key}.csv"
            if summary_path.exists():
                results['summary'] = pd.read_csv(summary_path)
                
        except Exception as e:
            st.error(f"Error loading enrichment results for {analysis_key}: {e}")
        
        return results if results else None
    
    def display_limma_results(self, analysis_name, limma_data, p_value_threshold=0.05, logfc_threshold=1.0):
        """Display limma analysis results with dynamic thresholds"""
        if limma_data is None or limma_data.empty:
            st.warning("No limma results available")
            return
        
        st.markdown(f"### 📈 Differential Expression Results - {analysis_name}")
        
        # Interactive threshold controls
        st.markdown("#### 🎯 Interactive Significance Thresholds")
        col_p, col_lfc = st.columns(2)
        
        with col_p:
            dynamic_p_threshold = st.number_input(
                "P-value threshold (raw):",
                min_value=0.001, max_value=0.1, value=p_value_threshold, step=0.005,
                format="%.3f",
                help="Raw p-value cutoff (not adjusted) - adjust and see results update instantly",
                key=f"p_thresh_{analysis_name}"
            )
        
        with col_lfc:
            dynamic_logfc_threshold = st.number_input(
                "LogFC threshold:",
                min_value=0.0, max_value=3.0, value=logfc_threshold, step=0.1,
                help="Adjust log fold-change cutoff and see results update instantly",
                key=f"lfc_thresh_{analysis_name}"
            )
        
        # Apply dynamic significance classification
        limma_data_dynamic = limma_data.copy()
        
        # Ensure numeric columns are properly typed
        limma_data_dynamic['adj.P.Val'] = pd.to_numeric(limma_data_dynamic['adj.P.Val'], errors='coerce')
        limma_data_dynamic['logFC'] = pd.to_numeric(limma_data_dynamic['logFC'], errors='coerce')
        limma_data_dynamic['P.Value'] = pd.to_numeric(limma_data_dynamic['P.Value'], errors='coerce')
        
        # Recalculate significance based on current thresholds
        limma_data_dynamic['sig'] = 'Not significant'
        significant_mask = (
            (limma_data_dynamic['P.Value'] < dynamic_p_threshold) & 
            (abs(limma_data_dynamic['logFC']) > dynamic_logfc_threshold)
        )
        limma_data_dynamic.loc[significant_mask, 'sig'] = 'Significant'
        
        # Recalculate -log10(p-value) for plotting
        limma_data_dynamic['neg_log10_p'] = -np.log10(limma_data_dynamic['P.Value'])
        
        # Summary metrics with dynamic thresholds
        col1, col2, col3, col4 = st.columns(4)
        
        total_proteins = len(limma_data_dynamic)
        significant = len(limma_data_dynamic[limma_data_dynamic['sig'] == 'Significant'])
        upregulated = len(limma_data_dynamic[(limma_data_dynamic['sig'] == 'Significant') & (limma_data_dynamic['logFC'] > 0)])
        downregulated = len(limma_data_dynamic[(limma_data_dynamic['sig'] == 'Significant') & (limma_data_dynamic['logFC'] < 0)])
        
        with col1:
            st.metric("Total Proteins", total_proteins)
        with col2:
            st.metric("Significant", significant, f"{significant/total_proteins*100:.1f}%")
        with col3:
            st.metric("Upregulated", upregulated)
        with col4:
            st.metric("Downregulated", downregulated)
        
        # Show threshold info
        st.info(f"Using thresholds: P-value < {dynamic_p_threshold}, |logFC| > {dynamic_logfc_threshold}")
        
        # Volcano plot with dynamic thresholds
        volcano_plot = self.plotter.create_volcano_plot(
            limma_data_dynamic, 
            f"({analysis_name})",
            p_threshold=dynamic_p_threshold,
            logfc_threshold=dynamic_logfc_threshold
        )
        st.plotly_chart(volcano_plot, use_container_width=True)
        
        # All proteins table (ordered by logFC, with significant ones highlighted)
        st.markdown("#### 📋 All Proteins (Ordered by LogFC)")
        
        # Sort by logFC descending, using dynamic data
        limma_data_sorted = limma_data_dynamic.sort_values(['logFC'], ascending=False)
        
        # Create a color-coded display using dynamic significance
        def highlight_significant(row):
            if row['sig'] == 'Significant':
                return ['background-color: #ffeb3b; font-weight: bold'] * len(row)
            else:
                return [''] * len(row)
        
        display_cols = ['Assay', 'logFC', 'P.Value', 'adj.P.Val', 'sig', 'neg_log10_p']
        available_display_cols = [col for col in display_cols if col in limma_data_sorted.columns]
        
        styled_df = limma_data_sorted[available_display_cols].style.apply(highlight_significant, axis=1)
        st.dataframe(styled_df, use_container_width=True, height=500)
        
        # Download all results (with dynamic significance)
        csv_all = limma_data_sorted.to_csv(index=False)
        st.download_button(
            "📥 Download All Results",
            data=csv_all,
            file_name=f"all_proteins_{analysis_name.lower().replace(' ', '_')}_p{dynamic_p_threshold}_lfc{dynamic_logfc_threshold}.csv",
            mime="text/csv"
        )
        
        # Top significant proteins table (using dynamic thresholds)
        if significant > 0:
            st.markdown("#### 🎯 Top Significant Proteins")
            sig_proteins = limma_data_dynamic[limma_data_dynamic['sig'] == 'Significant'].sort_values('P.Value')
            display_cols = ['Assay', 'logFC', 'P.Value', 'adj.P.Val', 'neg_log10_p']
            st.dataframe(
                sig_proteins[display_cols].head(20),
                use_container_width=True,
                height=400
            )
            
            # Download button for significant proteins
            csv = sig_proteins.to_csv(index=False)
            st.download_button(
                "📥 Download Significant Proteins",
                data=csv,
                file_name=f"significant_proteins_{analysis_name.lower().replace(' ', '_')}_p{dynamic_p_threshold}_lfc{dynamic_logfc_threshold}.csv",
                mime="text/csv"
            )
        else:
            st.info(f"No proteins meet the current significance criteria (P < {dynamic_p_threshold}, |logFC| > {dynamic_logfc_threshold})")
    
    def display_enrichment_results(self, analysis_name, enrichment_data):
        """Display enrichment analysis results"""
        if enrichment_data is None:
            st.warning("No enrichment results available")
            return
        
        st.markdown(f"### 🎯 Pathway Enrichment Results - {analysis_name}")
        
        # Create tabs for different ontologies
        available_ontologies = []
        if 'bp' in enrichment_data and not enrichment_data['bp'].empty:
            available_ontologies.append("Biological Process")
        if 'mf' in enrichment_data and not enrichment_data['mf'].empty:
            available_ontologies.append("Molecular Function")
        if 'cc' in enrichment_data and not enrichment_data['cc'].empty:
            available_ontologies.append("Cellular Component")
        if 'kegg' in enrichment_data and not enrichment_data['kegg'].empty:
            available_ontologies.append("KEGG Pathways")
        if 'reactome' in enrichment_data and not enrichment_data['reactome'].empty:
            available_ontologies.append("REACTOME Pathways")
        
        if not available_ontologies:
            st.info("No significant enrichment terms found")
            return
        
        tabs = st.tabs(available_ontologies)
        
        for i, ont_name in enumerate(available_ontologies):
            with tabs[i]:
                # Map ontology names to keys
                ont_key_map = {
                    "Biological Process": "bp",
                    "Molecular Function": "mf", 
                    "Cellular Component": "cc",
                    "KEGG Pathways": "kegg",
                    "REACTOME Pathways": "reactome"
                }
                ont_key = ont_key_map[ont_name]
                
                if ont_key in enrichment_data:
                    ont_data = enrichment_data[ont_key]
                    
                    # Summary metrics
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        st.metric("Enriched Terms", len(ont_data))
                    with col2:
                        if 'NES' in ont_data.columns:
                            avg_nes = ont_data['NES'].mean()
                            st.metric("Avg NES", f"{avg_nes:.2f}")
                    with col3:
                        if 'p.adjust' in ont_data.columns:
                            min_padj = ont_data['p.adjust'].min()
                            st.metric("Min Adj P-val", f"{min_padj:.2e}")
                    
                    # Enrichment plot
                    enrich_plot = self.plotter.create_enrichment_dotplot(
                        ont_data, f"- {ont_name} ({analysis_name})"
                    )
                    if enrich_plot:
                        st.plotly_chart(enrich_plot, use_container_width=True)
                    
                    # Results table
                    st.markdown(f"#### 📋 {ont_name} Terms")
                    
                    # Define display columns based on database type
                    if ont_key in ['bp', 'mf', 'cc']:
                        # GO terms columns
                        display_cols = ['Description', 'NES', 'pvalue', 'p.adjust', 'setSize']
                    elif ont_key == 'kegg':
                        # KEGG pathways columns
                        display_cols = ['Description', 'NES', 'pvalue', 'p.adjust', 'setSize']
                    elif ont_key == 'reactome':
                        # REACTOME pathways columns  
                        display_cols = ['Description', 'NES', 'pvalue', 'p.adjust', 'setSize']
                    else:
                        # Default columns
                        display_cols = ['Description', 'NES', 'pvalue', 'p.adjust', 'setSize']
                    
                    available_cols = [col for col in display_cols if col in ont_data.columns]
                    
                    # Show all available columns if none of the expected ones are found
                    if not available_cols:
                        available_cols = list(ont_data.columns)[:6]  # Show first 6 columns
                        st.info(f"Showing available columns for {ont_name}: {available_cols}")
                    
                    st.dataframe(
                        ont_data[available_cols].head(20),
                        use_container_width=True,
                        height=400
                    )
                    
                    # Download button
                    csv = ont_data.to_csv(index=False)
                    st.download_button(
                        f"📥 Download {ont_name} Results",
                        data=csv,
                        file_name=f"enrichment_{ont_key}_{analysis_name.lower().replace(' ', '_')}.csv",
                        mime="text/csv",
                        key=f"download_{ont_key}_{analysis_name}"
                    )
    
    def display_comparative_analysis(self):
        """Display comparative analysis across all groups"""
        if not st.session_state.analysis_results:
            return
        
        st.markdown("## 🔍 Comparative Analysis")
        
        # Collect summary data
        summary_data = []
        for result in st.session_state.analysis_results.values():
            if result.get('limma_data') is not None:
                limma_data = result['limma_data']
                total = len(limma_data)
                significant = len(limma_data[limma_data['sig'] == 'Significant'])
                upregulated = len(limma_data[(limma_data['sig'] == 'Significant') & (limma_data['logFC'] > 0)])
                downregulated = len(limma_data[(limma_data['sig'] == 'Significant') & (limma_data['logFC'] < 0)])
                
                summary_data.append({
                    'name': result['name'],
                    'total_proteins': total,
                    'significant_proteins': significant,
                    'upregulated': upregulated,
                    'downregulated': downregulated
                })
        
        if summary_data:
            # Summary bar plot
            summary_plot = self.plotter.create_summary_barplot(summary_data)
            st.plotly_chart(summary_plot, use_container_width=True)
            
            # Summary table
            st.markdown("### 📊 Summary Statistics")
            summary_df = pd.DataFrame(summary_data)
            summary_df['% Significant'] = (summary_df['significant_proteins'] / summary_df['total_proteins'] * 100).round(1)
            st.dataframe(summary_df, use_container_width=True)
            
            # Top proteins across analyses
            self.display_common_proteins()
    
    def display_common_proteins(self):
        """Display proteins that are significant across multiple analyses"""
        st.markdown("### 🎯 Common Significant Proteins")
        
        # Collect significant proteins from all analyses
        all_significant = {}
        
        for result in st.session_state.analysis_results.values():
            if result.get('limma_data') is not None:
                limma_data = result['limma_data']
                significant = limma_data[limma_data['sig'] == 'Significant']
                
                for _, protein in significant.iterrows():
                    protein_name = protein['Assay']
                    if protein_name not in all_significant:
                        all_significant[protein_name] = []
                    
                    all_significant[protein_name].append({
                        'analysis': result['name'],
                        'logFC': protein['logFC'],
                        'p_value': protein['P.Value'],
                        'adj_p_value': protein['adj.P.Val']
                    })
        
        # Find proteins significant in multiple analyses
        multi_significant = {k: v for k, v in all_significant.items() if len(v) > 1}
        
        if multi_significant:
            st.markdown(f"Found {len(multi_significant)} proteins significant in multiple analyses:")
            
            # Create comparison table
            comparison_data = []
            for protein, analyses in multi_significant.items():
                for analysis_data in analyses:
                    comparison_data.append({
                        'Protein': protein,
                        'Analysis': analysis_data['analysis'],
                        'logFC': analysis_data['logFC'],
                        'P-value': analysis_data['p_value'],
                        'Adj P-value': analysis_data['adj_p_value']
                    })
            
            comparison_df = pd.DataFrame(comparison_data)
            st.dataframe(comparison_df, use_container_width=True, height=400)
            
            # Download button
            csv = comparison_df.to_csv(index=False)
            st.download_button(
                "📥 Download Common Proteins",
                data=csv,
                file_name="common_significant_proteins.csv",
                mime="text/csv"
            )
        else:
            st.info("No proteins were found to be significant across multiple analyses.")
    
    def run_analyses(self, selected_analyses, data_dir, timeout, parallel_execution,
                    p_value_threshold, logfc_threshold, enable_go_analysis, go_pvalue_cutoff, force_rerun):
        """Run selected analyses"""
        if not selected_analyses:
            st.error("Please select at least one analysis group.")
            return
        
        # Validate data files
        missing_files = self.validate_data_files(data_dir)
        if missing_files:
            st.error("❌ Missing required data files:")
            for file in missing_files:
                st.write(f"- {file}")
            return
        
        st.session_state.analysis_running = True
        
        progress_container = st.container()
        with progress_container:
            st.info("🚀 Starting analyses...")
            progress_bar = st.progress(0)
            status_text = st.empty()
            
            results = {}
            
            if parallel_execution and len(selected_analyses) > 1:
                # Run analyses in parallel
                with ThreadPoolExecutor(max_workers=min(3, len(selected_analyses))) as executor:
                    future_to_analysis = {
                        executor.submit(
                            self.run_single_analysis,
                            analysis_name,
                            self.analysis_groups[analysis_name],
                            data_dir,
                            timeout,
                            p_value_threshold,
                            logfc_threshold,
                            enable_go_analysis,
                            go_pvalue_cutoff,
                            force_rerun
                        ): analysis_name for analysis_name in selected_analyses
                    }
                    
                    completed = 0
                    for future in as_completed(future_to_analysis):
                        analysis_name = future_to_analysis[future]
                        try:
                            result = future.result()
                            results[result['key']] = result
                            
                            completed += 1
                            progress = completed / len(selected_analyses)
                            progress_bar.progress(progress)
                            status_text.text(f"Completed: {analysis_name} ({completed}/{len(selected_analyses)})")
                            
                            # Show immediate status
                            if result['status'] == 'completed':
                                st.success(f"✅ {analysis_name}: Analysis completed successfully")
                            elif result['status'] == 'partial':
                                st.warning(f"⚠️ {analysis_name}: Partial completion (check individual steps)")
                            else:
                                st.error(f"❌ {analysis_name}: {result.get('error', 'Unknown error')}")
                                
                        except Exception as e:
                            st.error(f"❌ {analysis_name}: Unexpected error - {str(e)}")
            else:
                # Run analyses sequentially
                for i, analysis_name in enumerate(selected_analyses):
                    status_text.text(f"Running: {analysis_name}...")
                    
                    result = self.run_single_analysis(
                        analysis_name,
                        self.analysis_groups[analysis_name],
                        data_dir,
                        timeout,
                        p_value_threshold,
                        logfc_threshold,
                        enable_go_analysis,
                        go_pvalue_cutoff,
                        force_rerun
                    )
                    
                    results[result['key']] = result
                    
                    progress = (i + 1) / len(selected_analyses)
                    progress_bar.progress(progress)
                    
                    # Show status
                    if result['status'] == 'completed':
                        st.success(f"✅ {analysis_name}: Analysis completed successfully")
                    elif result['status'] == 'partial':
                        st.warning(f"⚠️ {analysis_name}: Partial completion")
                    else:
                        st.error(f"❌ {analysis_name}: {result.get('error', 'Unknown error')}")
            
            # Update session state
            st.session_state.analysis_results.update(results)
            st.session_state.analysis_running = False
            
            # Final status
            successful = sum(1 for r in results.values() if r['status'] == 'completed')
            st.success(f"🎉 Analysis complete! {successful}/{len(selected_analyses)} analyses successful.")
    
    def display_results_tabs(self):
        """Display results in organized tabs"""
        if not st.session_state.analysis_results:
            st.info("No analysis results available. Run analyses to see results.")
            return
        
        st.markdown("## 📊 Analysis Results")
        
        # Create tabs for each completed analysis
        completed_analyses = [
            (key, result) for key, result in st.session_state.analysis_results.items()
            if result.get('status') in ['completed', 'partial']
        ]
        
        if not completed_analyses:
            st.warning("No completed analyses found.")
            return
        
        # Add comparative analysis tab
        tab_names = [result['name'] for _, result in completed_analyses] + ["Comparative Analysis"]
        tabs = st.tabs(tab_names)
        
        # Individual analysis tabs
        for i, (key, result) in enumerate(completed_analyses):
            with tabs[i]:
                st.markdown(f"### Results for {result['name']}")
                
                # Create sub-tabs for limma and enrichment
                if result.get('limma_success') and result.get('enrichment_success') and result.get('go_analysis_enabled', True):
                    limma_tab, enrichment_tab = st.tabs(["📈 Differential Expression", "🎯 GO Term Enrichment"])
                    
                    with limma_tab:
                        self.display_limma_results(result['name'], result.get('limma_data'), 0.05, 1.0)
                    
                    with enrichment_tab:
                        self.display_enrichment_results(result['name'], result.get('enrichment_data'))
                        
                elif result.get('limma_success'):
                    self.display_limma_results(result['name'], result.get('limma_data'), 0.05, 1.0)
                    if not result.get('go_analysis_enabled', True):
                        st.info("🔄 GO Term Enrichment was disabled for this analysis")
                    
                elif result.get('enrichment_success'):
                    self.display_enrichment_results(result['name'], result.get('enrichment_data'))
                    
                else:
                    st.error("No successful analyses to display")
                    if 'error' in result:
                        st.error(f"Error: {result['error']}")
        
        # Comparative analysis tab
        with tabs[-1]:
            self.display_comparative_analysis()
    
    def run(self):
        """Main dashboard function"""
        # Header
        st.markdown('<h1 class="main-header">🧬 CRESCENDO Proteomics Analysis Dashboard</h1>', 
                   unsafe_allow_html=True)
        st.markdown("**Differential Expression & Pathway Enrichment Analysis**")
        st.markdown("---")
        
        # Sidebar configuration
        data_dir, selected_analyses, timeout, parallel_execution, enable_go_analysis, go_pvalue_cutoff, force_rerun = self.create_sidebar()
        
        # R Environment status check
        self.display_r_status()
        
        # Main content
        col1, col2 = st.columns([3, 1])
        
        with col1:
            st.markdown("### 🎯 Analysis Overview")
            st.markdown("""
            This dashboard performs comprehensive proteomics analysis including:
            - **Differential Expression Analysis** using limma with customizable thresholds
            - **Gene Ontology (GO) Term Enrichment Analysis** using clusterProfiler
            - **Interactive Visualizations** with volcano plots and enrichment plots
            - **Comparative Analysis** across etiology groups
            - **All Proteins View** ordered by logFC with significant proteins highlighted
            - **Customizable significance thresholds** for p-values and log fold-changes
            """)
        
        with col2:
            # Run analysis button
            run_button = st.button(
                "🚀 Run Analysis",
                type="primary",
                disabled=st.session_state.analysis_running,
                use_container_width=True
            )
            
            if st.session_state.analysis_running:
                st.info("⏳ Analysis in progress...")
            
            # Clear results button
            if st.button("🗑️ Clear All Results", use_container_width=True):
                st.session_state.analysis_results = {}
                st.success("Results cleared!")
                st.rerun()
        
        # Run analyses if button clicked
        if run_button:
            # Set default thresholds for R script analysis (display uses dynamic thresholds)
            p_value_threshold = 0.05  # Default for R script
            logfc_threshold = 1.0     # Default for R script
            
            self.run_analyses(selected_analyses, data_dir, timeout, parallel_execution,
                            p_value_threshold, logfc_threshold, enable_go_analysis, go_pvalue_cutoff, force_rerun)
            st.rerun()  # Refresh to show results
        
        # Display results
        self.display_results_tabs()
        
        # Instructions
        with st.expander("📋 Instructions & Requirements"):
            st.markdown("""
            ### 📁 Required Data Files
            Ensure these files exist in your data directory:
            - `WGCNA_analysis.RData` - Pre-processed WGCNA results
            - `Q-17124_Data Delivery/Q-17124_NPX_2025-03-24.parquet` - NPX expression data
            - `sampling_18022024/olink_sampled_patient_information_all.csv` - Sample metadata
            
            ### 🔧 R Dependencies
            Required R packages:
            ```r
            install.packages(c("OlinkAnalyze", "dplyr", "ggplot2", "stringr", 
                              "tidyr", "tidyverse", "limma", "lme4", "readr", 
                              "janitor", "ggrepel"))
            BiocManager::install(c("clusterProfiler", "DOSE", "org.Hs.eg.db", "enrichplot"))
            ```
            
            ### 📊 Analysis Groups
            - **All Groups**: Combined analysis across all etiology groups
            - **Large Artery Atherosclerosis (Group 1)**: Large vessel stroke
            - **Cardioembolism (Group 2)**: Cardiac source embolism
            - **Small Vessel Occlusion (Group 3)**: Lacunar stroke
            - **Stroke of Other Etiology (Group 4)**: Other determined etiology
            - **Stroke of Undetermined Etiology (Group 5)**: Cryptogenic stroke
            
            ### ⚙️ Advanced Features
            - **Customizable Thresholds**: Adjust p-value and logFC cutoffs for significance
            - **GO Term Analysis**: Toggle Gene Ontology enrichment analysis on/off
            - **All Proteins View**: See all proteins ordered by logFC with significant ones highlighted
            - **Parallel Execution**: Run multiple analyses simultaneously
            - **Interactive Plots**: Hover for detailed information
            - **Data Export**: Download results as CSV files
            - **Comparative Analysis**: Cross-group protein comparison
            
            ### 🎯 New Features Added
            1. **Significance Thresholds**: Customize p-value (default: 0.05) and logFC (default: 1.0) thresholds
            2. **GO Analysis Control**: Enable/disable GO term enrichment with custom p-value cutoffs
            3. **Enhanced Protein Table**: View all proteins sorted by logFC with highlighting for significant ones
            4. **Better Error Reporting**: Detailed error messages with expandable sections
            5. **Result Caching**: Automatically loads existing results to avoid re-running analyses
            6. **Force Re-run Option**: Override caching to re-run analyses with new parameters
            """)


def main():
    """Main application entry point"""
    dashboard = AnalysisDashboard()
    dashboard.run()


if __name__ == "__main__":
    main()