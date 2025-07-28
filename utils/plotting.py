# utils/plotting.py
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np

class PlotGenerator:
    def __init__(self):
        self.colors = {
            'significant': '#FF6B6B',
            'not_significant': '#A8A8A8',
            'upregulated': '#FF4444', 
            'downregulated': '#4444FF',
            'neutral': '#888888'
        }
    
    def create_volcano_plot(self, df, title_suffix="", p_threshold=0.05, logfc_threshold=1.0):
        """Create interactive volcano plot with dynamic threshold lines"""
        fig = go.Figure()
        
        # Add non-significant points
        non_sig = df[df['sig'] == 'Not significant']
        if not non_sig.empty:
            fig.add_trace(go.Scatter(
                x=non_sig['logFC'],
                y=non_sig['neg_log10_p'],
                mode='markers',
                marker=dict(
                    color=self.colors['not_significant'], 
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
        sig = df[df['sig'] == 'Significant']
        if not sig.empty:
            fig.add_trace(go.Scatter(
                x=sig['logFC'],
                y=sig['neg_log10_p'],
                mode='markers',
                marker=dict(
                    color=self.colors['significant'], 
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
            title=f"Volcano Plot: Recurrent vs Non-Recurrent {title_suffix}",
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
    
    def create_enrichment_dotplot(self, df, title_suffix="", max_terms=20):
        """Create enrichment dotplot"""
        if df.empty:
            return None
            
        # Sort by NES and take top terms
        df_plot = df.sort_values('NES', ascending=True).tail(max_terms)
        
        fig = go.Figure()
        
        # Create scatter plot
        fig.add_trace(go.Scatter(
            x=df_plot['NES'],
            y=list(range(len(df_plot))),
            mode='markers',
            marker=dict(
                size=np.sqrt(df_plot['setSize']) * 2,  # Size by gene set size
                color=-np.log10(df_plot['p.adjust']),
                colorscale='Viridis',
                showscale=True,
                colorbar=dict(
                    title="-log10(adj p-value)",
                    titleside="right"
                ),
                line=dict(width=1, color='white')
            ),
            text=df_plot['Description'],
            customdata=df_plot[['p.adjust', 'setSize', 'NES']],
            hovertemplate='<b>%{text}</b><br>' +
                         'NES: %{customdata[2]:.3f}<br>' +
                         'Set Size: %{customdata[1]}<br>' +
                         'Adj p-value: %{customdata[0]:.2e}<extra></extra>'
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
    
    def create_summary_barplot(self, summary_data):
        """Create summary bar plot across analyses"""
        fig = go.Figure()
        
        categories = ['Total Proteins', 'Significant', 'Upregulated', 'Downregulated']
        
        for i, analysis in enumerate(summary_data):
            fig.add_trace(go.Bar(
                name=analysis['name'],
                x=categories,
                y=[
                    analysis.get('total_proteins', 0),
                    analysis.get('significant_proteins', 0), 
                    analysis.get('upregulated', 0),
                    analysis.get('downregulated', 0)
                ],
                text=[
                    analysis.get('total_proteins', 0),
                    analysis.get('significant_proteins', 0),
                    analysis.get('upregulated', 0), 
                    analysis.get('downregulated', 0)
                ],
                textposition='auto'
            ))
        
        fig.update_layout(
            title="Analysis Summary Across Groups",
            xaxis_title="Protein Categories",
            yaxis_title="Count",
            barmode='group',
            template="plotly_white",
            height=500
        )
        
        return fig

