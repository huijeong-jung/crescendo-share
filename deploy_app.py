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
    
    def create_volcano_plot(data, analysis_name="Analysis", p_threshold=0.05, logfc_threshold=1.0):
        """Enhanced wrapper with error handling"""
        if data.empty:
            return None
        # Check if required columns exist
        required_cols = ['logFC', 'P.Value', 'sig', 'neg_log10_p', 'Assay']
        missing_cols = [col for col in required_cols if col not in data.columns]
        if missing_cols:
            print(f"Warning: Missing columns {missing_cols} for volcano plot")
            return None
        try:
            return plot_generator.create_volcano_plot(data, f"({analysis_name})", p_threshold, logfc_threshold)
        except Exception as e:
            print(f"Error creating volcano plot: {e}")
            return None
    
    def create_enrichment_dotplot(data, title_suffix="", max_terms=20):
        """Enhanced wrapper with error handling"""
        # Check if required columns exist
        required_cols = ['NES', 'p.adjust', 'setSize', 'Description']
        missing_cols = [col for col in required_cols if col not in data.columns]
        if missing_cols:
            print(f"Warning: Missing columns {missing_cols} for enrichment plot")
            return None
        if data.empty:
            return None
        try:
            # Try the original plot generator first
            return plot_generator.create_enrichment_dotplot(data, title_suffix, max_terms)
        except Exception as e:
            print(f"Error creating enrichment plot with original function: {e}")
            # Fall back to inline implementation
            try:
                import plotly.graph_objects as go
                import numpy as np
                
                # Sort by NES and take top terms
                df_plot = data.sort_values('NES', ascending=True).tail(max_terms)
                df_plot = df_plot.dropna(subset=['NES', 'p.adjust', 'setSize'])
                
                if df_plot.empty:
                    return None
                
                fig = go.Figure()
                
                # Simple dot plot without complex colorbar
                p_adjust_vals = df_plot['p.adjust'].clip(lower=1e-10)
                color_vals = [-np.log10(p) for p in p_adjust_vals]
                
                fig.add_trace(go.Scatter(
                    x=df_plot['NES'],
                    y=list(range(len(df_plot))),
                    mode='markers',
                    marker=dict(
                        size=[max(8, min(20, np.sqrt(s) * 2)) for s in df_plot['setSize']],
                        color=color_vals,
                        colorscale='Viridis',
                        line=dict(width=1, color='white')
                    ),
                    text=df_plot['Description'],
                    customdata=df_plot[['p.adjust', 'setSize', 'NES']].values,
                    hovertemplate='<b>%{text}</b><br>' +
                                 'NES: %{customdata[2]:.3f}<br>' +
                                 'Set Size: %{customdata[1]}<br>' +
                                 'Adj p-value: %{customdata[0]:.2e}<extra></extra>',
                    showlegend=False
                ))
                
                # Update layout
                fig.update_layout(
                    title=f"Gene Set Enrichment Analysis {title_suffix}",
                    xaxis_title="Normalized Enrichment Score (NES)",
                    yaxis=dict(
                        tickmode='array',
                        tickvals=list(range(len(df_plot))),
                        ticktext=[desc[:60] + '...' if len(desc) > 60 else desc 
                                 for desc in df_plot['Description']]
                    ),
                    template="plotly_white",
                    height=max(400, len(df_plot) * 30),
                    width=1000,
                    margin=dict(l=300)
                )
                
                return fig
                
            except Exception as e2:
                print(f"Fallback enrichment plot also failed: {e2}")
                return None
        
except ImportError:
    # Fallback: define the functions inline if import fails
    import numpy as np
    
    def create_volcano_plot(data, p_threshold=0.05, logfc_threshold=1.0, analysis_name="Analysis"):
        """Enhanced volcano plot function matching utils/plotting.py"""
        import plotly.graph_objects as go
        import numpy as np
        
        fig = go.Figure()
        
        # Color scheme
        colors = {
            'significant': '#FF6B6B',
            'not_significant': '#A8A8A8'
        }
        
        # Add non-significant points
        non_sig = data[data['sig'] == 'Not significant']
        if not non_sig.empty:
            fig.add_trace(go.Scatter(
                x=non_sig['logFC'],
                y=non_sig['neg_log10_p'],
                mode='markers',
                marker=dict(
                    color=colors['not_significant'], 
                    size=4, 
                    opacity=0.6,
                    line=dict(width=0)
                ),
                name='Not significant',
                text=non_sig['Assay'],
                customdata=non_sig[['P.Value', 'adj.P.Val']],
                hovertemplate='<b>%{text}</b><br>' +
                             'logFC: %{x:.3f}<br>' +
                             '-log10(p): %{y:.3f}<br>' +
                             'p-value: %{customdata[0]:.2e}<br>' +
                             'adj p-value: %{customdata[1]:.2e}<extra></extra>'
            ))
        
        # Add significant points
        sig = data[data['sig'] == 'Significant']
        if not sig.empty:
            fig.add_trace(go.Scatter(
                x=sig['logFC'],
                y=sig['neg_log10_p'],
                mode='markers',
                marker=dict(
                    color=colors['significant'], 
                    size=6, 
                    opacity=0.8,
                    line=dict(width=1, color='white')
                ),
                name='Significant',
                text=sig['Assay'],
                customdata=sig[['P.Value', 'adj.P.Val']],
                hovertemplate='<b>%{text}</b><br>' +
                             'logFC: %{x:.3f}<br>' +
                             '-log10(p): %{y:.3f}<br>' +
                             'p-value: %{customdata[0]:.2e}<br>' +
                             'adj p-value: %{customdata[1]:.2e}<extra></extra>'
            ))
        
        # Add dynamic threshold lines
        fig.add_hline(y=-np.log10(p_threshold), line_dash="dash", line_color="blue", 
                      annotation_text=f"p = {p_threshold}")
        fig.add_vline(x=logfc_threshold, line_dash="dash", line_color="red", opacity=0.7,
                      annotation_text=f"logFC = {logfc_threshold}")
        fig.add_vline(x=-logfc_threshold, line_dash="dash", line_color="red", opacity=0.7,
                      annotation_text=f"logFC = -{logfc_threshold}")
        fig.add_vline(x=0, line_dash="dot", line_color="gray", opacity=0.5)
        
        # Update layout
        fig.update_layout(
            title=f"Volcano Plot: Recurrent vs Non-Recurrent {analysis_name}",
            xaxis_title="Log2 Fold Change",
            yaxis_title="-Log10(p-value)",
            template="plotly_white",
            width=800,
            height=600,
            showlegend=True,
            legend=dict(
                yanchor="top",
                y=0.99,
                xanchor="left", 
                x=0.01
            )
        )
        
        return fig

    def create_enrichment_dotplot(data, title_suffix="", max_terms=20):
        """Create enrichment dotplot matching utils/plotting.py"""
        import plotly.graph_objects as go
        import numpy as np
        
        if data.empty:
            return None
            
        # Check if required columns exist
        required_cols = ['NES', 'p.adjust', 'setSize', 'Description']
        missing_cols = [col for col in required_cols if col not in data.columns]
        if missing_cols:
            print(f"Warning: Missing columns {missing_cols} for enrichment plot")
            return None
            
        # Sort by NES and take top terms - handle potential sorting issues
        try:
            df_plot = data.sort_values('NES', ascending=True).tail(max_terms)
        except Exception as e:
            print(f"Error sorting data: {e}")
            df_plot = data.head(max_terms)  # Fallback to first N terms
        
        # Ensure no NaN values
        df_plot = df_plot.dropna(subset=['NES', 'p.adjust', 'setSize'])
        
        if df_plot.empty:
            return None
        
        fig = go.Figure()
        
        # Create scatter plot with safer colorbar configuration
        try:
            # Ensure p.adjust values are valid and > 0
            p_adjust_vals = df_plot['p.adjust'].clip(lower=1e-10)  # Prevent log(0)
            color_vals = [-np.log10(p) for p in p_adjust_vals]
            
            # First try with colorbar
            fig.add_trace(go.Scatter(
                x=df_plot['NES'],
                y=list(range(len(df_plot))),
                mode='markers',
                marker=dict(
                    size=[max(8, min(20, np.sqrt(s) * 2)) for s in df_plot['setSize']],  # Size by gene set size with bounds
                    color=color_vals,
                    colorscale='Viridis',
                    showscale=True,
                    colorbar=dict(
                        title="-log10(adj p-value)",
                        x=1.02,
                        len=0.8
                    ),
                    line=dict(width=1, color='white')
                ),
                text=df_plot['Description'],
                customdata=df_plot[['p.adjust', 'setSize', 'NES']].values,
                hovertemplate='<b>%{text}</b><br>' +
                             'NES: %{customdata[2]:.3f}<br>' +
                             'Set Size: %{customdata[1]}<br>' +
                             'Adj p-value: %{customdata[0]:.2e}<extra></extra>',
                showlegend=False
            ))
        except Exception as e:
            print(f"Colorbar plot failed: {e}")
            try:
                # Try with simpler colorbar
                fig.add_trace(go.Scatter(
                    x=df_plot['NES'],
                    y=list(range(len(df_plot))),
                    mode='markers',
                    marker=dict(
                        size=[max(8, min(20, np.sqrt(s) * 2)) for s in df_plot['setSize']],
                        color=color_vals,
                        colorscale='Viridis',
                        line=dict(width=1, color='white')
                    ),
                    text=df_plot['Description'],
                    customdata=df_plot[['p.adjust', 'setSize', 'NES']].values,
                    hovertemplate='<b>%{text}</b><br>' +
                                 'NES: %{customdata[2]:.3f}<br>' +
                                 'Set Size: %{customdata[1]}<br>' +
                                 'Adj p-value: %{customdata[0]:.2e}<extra></extra>',
                    showlegend=False
                ))
            except Exception as e2:
                print(f"Simpler colorbar also failed: {e2}")
                # Final fallback to solid color
                fig.add_trace(go.Scatter(
                    x=df_plot['NES'],
                    y=list(range(len(df_plot))),
                    mode='markers',
                    marker=dict(
                        size=[max(8, min(20, np.sqrt(s) * 2)) for s in df_plot['setSize']],
                        color='blue',
                        line=dict(width=1, color='white')
                    ),
                    text=df_plot['Description'],
                    customdata=df_plot[['p.adjust', 'setSize', 'NES']].values,
                    hovertemplate='<b>%{text}</b><br>' +
                                 'NES: %{customdata[2]:.3f}<br>' +
                                 'Set Size: %{customdata[1]}<br>' +
                                 'Adj p-value: %{customdata[0]:.2e}<extra></extra>',
                    showlegend=False
                ))
        
        # Update layout
        fig.update_layout(
            title=f"Gene Set Enrichment Analysis {title_suffix}",
            xaxis_title="Normalized Enrichment Score (NES)",
            yaxis=dict(
                tickmode='array',
                tickvals=list(range(len(df_plot))),
                ticktext=[desc[:60] + '...' if len(desc) > 60 else desc 
                         for desc in df_plot['Description']]
            ),
            template="plotly_white",
            height=max(400, len(df_plot) * 30),
            width=1000,
            margin=dict(l=300)  # More space for labels
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
        try:
            volcano_plot = create_volcano_plot(
                limma_data_dynamic, 
                analysis_name=f"({analysis_name})",
                p_threshold=dynamic_p_threshold,
                logfc_threshold=dynamic_logfc_threshold
            )
            if volcano_plot:
                st.plotly_chart(volcano_plot, use_container_width=True)
            else:
                st.warning("Could not create volcano plot - please check data format")
        except Exception as e:
            st.error(f"Error creating volcano plot: {str(e)}")
        
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
    
    def display_enrichment_results(self, analysis_name, analysis_key):
        """Display enrichment analysis results matching app.py"""
        st.markdown(f"### 🎯 Pathway Enrichment Results - {analysis_name}")
        
        # Load enrichment data for all ontologies
        enrichment_data = {}
        ontologies = ['bp', 'mf', 'cc', 'kegg', 'reactome']
        
        for ont in ontologies:
            data = self.load_enrichment_results(analysis_key, ont)
            if data is not None and not data.empty:
                enrichment_data[ont] = data
        
        if not enrichment_data:
            st.warning("No enrichment results available")
            return
        
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
                    
                    # Enrichment dot plot
                    try:
                        enrich_plot = create_enrichment_dotplot(
                            ont_data, f"- {ont_name} ({analysis_name})"
                        )
                        if enrich_plot:
                            st.plotly_chart(enrich_plot, use_container_width=True)
                    except Exception as e:
                        st.warning(f"Could not create enrichment plot for {ont_name}: {str(e)}")
                        # Fallback to simple bar chart
                        if len(ont_data) > 0:
                            top_terms = ont_data.head(15)
                            fig = px.bar(
                                top_terms,
                                x='NES' if 'NES' in top_terms.columns else 'enrichmentScore',
                                y='Description',
                                orientation='h',
                                title=f"Top {ont_name} Terms",
                                labels={'x': 'Enrichment Score', 'y': 'Pathway/Term'}
                            )
                            fig.update_layout(height=500)
                            st.plotly_chart(fig, use_container_width=True)
                    
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
