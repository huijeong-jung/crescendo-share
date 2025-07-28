# deploy_app.py - Streamlit Cloud deployment version
import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path
import sys
import os
import numpy as np

# Handle imports more robustly for deployment
try:
    # Add current directory to path
    current_dir = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, current_dir)
    
    from utils.plotting import PlotGenerator
    plot_generator = PlotGenerator()
    
    def create_volcano_plot(data, p_threshold=0.05, logfc_threshold=1.0, analysis_name="Analysis"):
        return plot_generator.create_volcano_plot(data, f"({analysis_name})", p_threshold, logfc_threshold)
    
    # Make plot_generator globally accessible    
    globals()['plot_generator'] = plot_generator
        
except ImportError:
    # Fallback: define the function inline if import fails
    import numpy as np
    plot_generator = None  # No plot generator available
    
    def create_volcano_plot(data, p_threshold=0.05, logfc_threshold=1.0, analysis_name="Analysis"):
        """Fallback volcano plot function"""
        import plotly.graph_objects as go
        
        # Calculate significance
        data = data.copy()
        data['significant'] = (
            (data['P.Value'] < p_threshold) & 
            (abs(data['logFC']) > logfc_threshold)
        )
        
        # Create colors
        data['color'] = data['significant'].map({True: 'red', False: 'lightgray'})
        
        # Create plot
        fig = go.Figure()
        
        # Add points
        fig.add_trace(go.Scatter(
            x=data['logFC'],
            y=-np.log10(data['P.Value']),
            mode='markers',
            marker=dict(color=data['color'], size=6, opacity=0.7),
            text=data['Assay'],
            hovertemplate='<b>%{text}</b><br>logFC: %{x:.3f}<br>-log10(P): %{y:.3f}<extra></extra>',
            showlegend=False
        ))
        
        # Add threshold lines
        fig.add_hline(y=-np.log10(p_threshold), 
                     line_dash="dash", line_color="blue", opacity=0.7)
        fig.add_vline(x=logfc_threshold, 
                     line_dash="dash", line_color="blue", opacity=0.7)
        fig.add_vline(x=-logfc_threshold, 
                     line_dash="dash", line_color="blue", opacity=0.7)
        
        # Update layout
        fig.update_layout(
            title=f"Volcano Plot - {analysis_name}",
            xaxis_title="Log2 Fold Change",
            yaxis_title="-Log10(P-value)",
            template="plotly_white",
            width=800,
            height=600
        )
        
        return fig

st.set_page_config(
    page_title="Proteomics Analysis Results",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

class ProteomicsResultsViewer:
    """Read-only viewer for pre-computed proteomics results"""
    
    def __init__(self):
        self.results_dir = Path("results")
        self.limma_dir = self.results_dir / "limma_results"
        self.enrichment_dir = self.results_dir / "enrichment_results"
        
        # Available analyses
        self.available_analyses = {
            "All Groups (Holistic)": "all_groups_holistic",
            "Cardioembolism": "cardioembolism", 
            "Large Artery Atherosclerosis": "large_artery_atherosclerosis",
            "Small Vessel Occlusion": "small_vessel_occlusion",
            "Stroke of Other Etiology": "stroke_of_other_etiology",
            "Stroke of Undetermined Etiology": "stroke_of_undetermined_etiology"
        }
    
    def load_limma_results(self, analysis_key):
        """Load limma results for a specific analysis"""
        limma_file = self.limma_dir / f"limma_results_{analysis_key}.csv"
        if limma_file.exists():
            return pd.read_csv(limma_file)
        return None
    
    def load_enrichment_results(self, analysis_key, ontology):
        """Load enrichment results for a specific analysis and ontology"""
        enrichment_file = self.enrichment_dir / f"enrichment_{ontology}_{analysis_key}.csv"
        if enrichment_file.exists():
            return pd.read_csv(enrichment_file)
        return None
    
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
        
        # Apply dynamic filtering
        limma_data_filtered = limma_data.copy()
        limma_data_filtered['dynamic_sig'] = (
            (limma_data_filtered['P.Value'] < dynamic_p_threshold) & 
            (limma_data_filtered['logFC'].abs() > dynamic_logfc_threshold)
        )
        
        # Summary statistics
        total_proteins = len(limma_data_filtered)
        significant_proteins = len(limma_data_filtered[limma_data_filtered['dynamic_sig']])
        upregulated = len(limma_data_filtered[
            (limma_data_filtered['dynamic_sig']) & (limma_data_filtered['logFC'] > 0)
        ])
        downregulated = len(limma_data_filtered[
            (limma_data_filtered['dynamic_sig']) & (limma_data_filtered['logFC'] < 0)
        ])
        
        # Display summary
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("Total Proteins", total_proteins)
        with col2:
            st.metric("Significant", significant_proteins)
        with col3:
            st.metric("Upregulated", upregulated)
        with col4:
            st.metric("Downregulated", downregulated)
        
        # Volcano plot
        st.markdown("#### 🌋 Volcano Plot")
        if 'create_volcano_plot' in globals():
            volcano_fig = create_volcano_plot(
                limma_data_filtered, 
                p_threshold=dynamic_p_threshold,
                logfc_threshold=dynamic_logfc_threshold
            )
            st.plotly_chart(volcano_fig, use_container_width=True)
        else:
            # Fallback volcano plot
            fig = px.scatter(
                limma_data_filtered,
                x='logFC',
                y='neg_log10_p',
                color='dynamic_sig',
                hover_data=['Assay', 'P.Value'],
                title="Volcano Plot",
                labels={'logFC': 'Log2 Fold Change', 'neg_log10_p': '-Log10(P-value)'}
            )
            fig.add_hline(y=-np.log10(dynamic_p_threshold), line_dash="dash", line_color="red")
            fig.add_vline(x=dynamic_logfc_threshold, line_dash="dash", line_color="red")
            fig.add_vline(x=-dynamic_logfc_threshold, line_dash="dash", line_color="red")
            st.plotly_chart(fig, use_container_width=True)
        
        # Results table
        st.markdown("#### 📊 Significant Proteins")
        significant_results = limma_data_filtered[limma_data_filtered['dynamic_sig']].sort_values('P.Value')
        
        if not significant_results.empty:
            display_columns = ['Assay', 'logFC', 'P.Value', 'adj.P.Val', 'UniProt']
            if 'delta_NPX' in significant_results.columns:
                display_columns.append('delta_NPX')
                
            st.dataframe(
                significant_results[display_columns].head(50),
                use_container_width=True
            )
            
            # Download button for significant results
            csv = significant_results.to_csv(index=False)
            st.download_button(
                label="📥 Download Significant Results (CSV)",
                data=csv,
                file_name=f"significant_proteins_{analysis_name.lower().replace(' ', '_')}.csv",
                mime='text/csv'
            )
        else:
            st.info("No proteins meet the current significance criteria.")
        
        # Full protein list with highlighting
        st.markdown("#### 📋 All Proteins (Significant Highlighted)")
        
        # Prepare display data with highlighting
        full_results = limma_data_filtered.copy().sort_values('P.Value')
        display_columns = ['Assay', 'logFC', 'P.Value', 'adj.P.Val', 'UniProt']
        if 'delta_NPX' in full_results.columns:
            display_columns.append('delta_NPX')
        
        # Function to highlight significant rows
        def highlight_significant_rows(row):
            if row['dynamic_sig']:
                return ['background-color: #ffeb3b; font-weight: bold'] * len(row)
            else:
                return [''] * len(row)
        
        # Apply styling and display
        styled_df = full_results[display_columns].style.apply(highlight_significant_rows, axis=1)
        st.dataframe(
            styled_df,
            use_container_width=True,
            height=400
        )
        
        # Download button for all results
        csv_all = full_results.to_csv(index=False)
        st.download_button(
            label="📥 Download All Results (CSV)",
            data=csv_all,
            file_name=f"all_proteins_{analysis_name.lower().replace(' ', '_')}.csv",
            mime='text/csv',
            key=f"download_all_{analysis_name}"
        )
    
    def display_enrichment_results(self, analysis_name, analysis_key):
        """Display enrichment analysis results"""
        st.markdown(f"### 🧬 Pathway Enrichment Results - {analysis_name}")
        
        # Create tabs for different ontologies
        tab1, tab2, tab3, tab4, tab5 = st.tabs([
            "GO: Biological Process", 
            "GO: Molecular Function", 
            "GO: Cellular Component",
            "KEGG Pathways",
            "REACTOME Pathways"
        ])
        
        ontologies = {
            "GO: Biological Process": ("bp", tab1),
            "GO: Molecular Function": ("mf", tab2), 
            "GO: Cellular Component": ("cc", tab3),
            "KEGG Pathways": ("kegg", tab4),
            "REACTOME Pathways": ("reactome", tab5)
        }
        
        for ont_name, (ont_key, tab) in ontologies.items():
            with tab:
                enrichment_data = self.load_enrichment_results(analysis_key, ont_key)
                
                if enrichment_data is not None and not enrichment_data.empty:
                    st.markdown(f"#### {ont_name} Results")
                    
                    # Display summary
                    st.metric("Enriched Terms", len(enrichment_data))
                    
                    # Dot plot of top terms using the utility function
                    if len(enrichment_data) > 0:
                        try:
                            if plot_generator:
                                # Use the plotting utility function
                                enrich_plot = plot_generator.create_enrichment_dotplot(
                                    enrichment_data, f"- {ont_name}"
                                )
                                if enrich_plot:
                                    st.plotly_chart(enrich_plot, use_container_width=True)
                                else:
                                    raise Exception("Plot generation failed")
                            else:
                                raise Exception("No plot generator available")
                        except Exception as e:
                            # Fallback dot plot implementation
                            st.warning(f"Using fallback plot due to: {str(e)[:100]}...")
                            top_terms = enrichment_data.head(15)
                            
                            # Prepare data for dot plot
                            y_values = top_terms['Description']
                            x_values = top_terms['NES'] if 'NES' in top_terms.columns else top_terms['enrichmentScore']
                            
                            # Use p-value for color if available, otherwise use enrichment score
                            if 'pvalue' in top_terms.columns:
                                color_values = -np.log10(top_terms['pvalue'] + 1e-10)  # Add small value to avoid log(0)
                                color_label = '-Log10(P-value)'
                                color_scale = 'Viridis'
                            else:
                                color_values = x_values
                                color_label = 'Enrichment Score'
                                color_scale = 'Viridis'
                            
                            # Use gene ratio or set size for dot size if available
                            if 'setSize' in top_terms.columns:
                                size_values = top_terms['setSize']
                                size_label = 'Gene Set Size'
                            else:
                                size_values = [20] * len(top_terms)  # Default size
                                size_label = None
                            
                            fig = px.scatter(
                                x=x_values,
                                y=y_values,
                                color=color_values,
                                size=size_values,
                                title=f"Top {ont_name} Terms",
                                labels={
                                    'x': 'Normalized Enrichment Score' if 'NES' in top_terms.columns else 'Enrichment Score',
                                    'y': 'Pathway/Term',
                                    'color': color_label,
                                    'size': size_label
                                },
                                color_continuous_scale=color_scale,
                                size_max=15
                            )
                            
                            # Customize layout
                            fig.update_layout(
                                height=500,
                                yaxis={'categoryorder': 'total ascending'},  # Order by enrichment score
                                showlegend=True
                            )
                            
                            # Add vertical line at x=0 if using NES
                            if 'NES' in top_terms.columns:
                                fig.add_vline(x=0, line_dash="dash", line_color="gray", opacity=0.5)
                            
                            st.plotly_chart(fig, use_container_width=True)
                    
                    # Results table
                    display_columns = ['Description', 'setSize', 'enrichmentScore']
                    if 'NES' in enrichment_data.columns:
                        display_columns.append('NES')
                    if 'pvalue' in enrichment_data.columns:
                        display_columns.extend(['pvalue', 'p.adjust'])
                    
                    st.dataframe(
                        enrichment_data[display_columns].head(20),
                        use_container_width=True
                    )
                    
                    # Download button
                    csv = enrichment_data.to_csv(index=False)
                    st.download_button(
                        label=f"📥 Download {ont_name} Results",
                        data=csv,
                        file_name=f"{ont_key}_enrichment_{analysis_key}.csv",
                        mime='text/csv',
                        key=f"download_{ont_key}_{analysis_key}"
                    )
                else:
                    st.info(f"No {ont_name.lower()} enrichment results available.")
    
    def run(self):
        """Main app interface"""
        st.title("🧬 Proteomics Analysis Results Dashboard")
        st.markdown("""
        **Explore pre-computed differential expression and pathway enrichment results.**
        
        - 📊 **Interactive Thresholds**: Adjust significance criteria in real-time
        - 🌋 **Volcano Plots**: Visualize differential expression
        - 🧬 **Pathway Analysis**: GO, KEGG, and REACTOME enrichment results
        - 📥 **Download Options**: Export filtered results
        """)
        
        # Sidebar for analysis selection
        with st.sidebar:
            st.header("🎯 Analysis Selection")
            selected_analysis = st.selectbox(
                "Choose Analysis:",
                list(self.available_analyses.keys()),
                help="Select the stroke etiology group or holistic analysis"
            )
            
            analysis_key = self.available_analyses[selected_analysis]
            
            st.markdown("---")
            st.markdown("### 📁 Available Results")
            
            # Check what results are available
            limma_available = (self.limma_dir / f"limma_results_{analysis_key}.csv").exists()
            enrichment_available = any([
                (self.enrichment_dir / f"enrichment_{ont}_{analysis_key}.csv").exists()
                for ont in ['bp', 'mf', 'cc', 'kegg', 'reactome']
            ])
            
            st.write(f"✅ Differential Expression" if limma_available else "❌ Differential Expression")
            st.write(f"✅ Pathway Enrichment" if enrichment_available else "❌ Pathway Enrichment")
        
        # Main content
        analysis_key = self.available_analyses[selected_analysis]
        
        # Load and display limma results
        limma_data = self.load_limma_results(analysis_key)
        if limma_data is not None:
            self.display_limma_results(selected_analysis, limma_data)
            st.markdown("---")
            
        # Display enrichment results
        self.display_enrichment_results(selected_analysis, analysis_key)

# Run the app
if __name__ == "__main__":
    import numpy as np  # Import here to avoid issues
    viewer = ProteomicsResultsViewer()
    viewer.run()
