#!/usr/bin/env python3
"""
CDR-Focused GROMACS Analysis Dashboard
Specialized visualization tool for antibody CDR dynamics analysis
Based on 4G6K gevokizumab structure with comprehensive CDR insights
"""

import os
import re
import glob
import json
import math
from datetime import datetime
from pathlib import Path
from statistics import mean, stdev, median


class CDRAnalysisParser:
    """Enhanced parser specifically designed for CDR analysis"""
    
    def __init__(self, filename):
        self.filename = filename
        self.title = ""
        self.xlabel = "X"
        self.ylabel = "Y"
        self.subtitle = ""
        self.legends = []
        self.data = []
        self.metadata = {}
        self.metrics = {}
        self.analysis_type = self._identify_analysis_type()
        self._parse()
        self._calculate_cdr_metrics()
    
    def _identify_analysis_type(self):
        """Identify the type of analysis based on filename"""
        filename_lower = self.filename.lower()
        if 'cdr' in filename_lower:
            if 'rmsd' in filename_lower:
                return 'cdr_rmsd'
            elif 'rmsf' in filename_lower:
                return 'cdr_rmsf'
            elif 'sasa' in filename_lower:
                return 'cdr_sasa'
            elif 'gyrate' in filename_lower:
                return 'cdr_gyration'
        elif 'rmsd' in filename_lower:
            if 'comparison' in filename_lower:
                return 'rmsd_comparison'
            elif 'protein' in filename_lower:
                return 'protein_rmsd'
        elif 'rmsf' in filename_lower and 'protein' in filename_lower:
            return 'protein_rmsf'
        elif 'sasa' in filename_lower and 'protein' in filename_lower:
            return 'protein_sasa'
        elif 'gyrate' in filename_lower and 'protein' in filename_lower:
            return 'protein_gyration'
        elif 'hbond' in filename_lower:
            return 'hydrogen_bonds'
        elif 'temperature' in filename_lower:
            return 'temperature'
        elif 'pressure' in filename_lower:
            return 'pressure'
        elif 'potential' in filename_lower:
            return 'potential_energy'
        
        return 'unknown'
    
    def _parse(self):
        """Parse XVG file with CDR-specific enhancements"""
        with open(self.filename, 'r') as f:
            lines = f.readlines()
        
        # Extract metadata from header
        for i, line in enumerate(lines):
            if line.startswith('#'):
                if 'Created by:' in line:
                    self.metadata['created_by'] = line.split('Created by:')[1].strip()
                elif 'Command line:' in line and i+1 < len(lines):
                    self.metadata['command'] = lines[i+1].strip('#').strip()
            elif line.startswith('@'):
                if 'title' in line:
                    self.title = self._extract_quoted(line)
                elif 'xaxis' in line and 'label' in line:
                    self.xlabel = self._extract_quoted(line)
                elif 'yaxis' in line and 'label' in line:
                    self.ylabel = self._extract_quoted(line)
                elif 'subtitle' in line:
                    self.subtitle = self._extract_quoted(line)
                elif re.match(r'@\s*s\d+\s+legend', line):
                    self.legends.append(self._extract_quoted(line))
            elif not line.startswith(('@', '#', '&')):
                try:
                    values = [float(x) for x in line.strip().split()]
                    if values:
                        self.data.append(values)
                except ValueError:
                    continue
        
        # Enhanced title for CDR analyses
        if not self.title:
            self.title = self._generate_cdr_title()
    
    def _generate_cdr_title(self):
        """Generate CDR-specific titles"""
        titles = {
            'cdr_rmsd': 'CDR Regions RMSD Analysis',
            'cdr_rmsf': 'CDR Residue Flexibility (RMSF)',
            'cdr_sasa': 'CDR Solvent Accessible Surface Area',
            'cdr_gyration': 'CDR Compactness (Radius of Gyration)',
            'rmsd_comparison': 'Protein vs CDR RMSD Comparison',
            'protein_rmsd': 'Protein Backbone RMSD',
            'protein_rmsf': 'Protein Backbone Flexibility',
            'hydrogen_bonds': 'Intermolecular Hydrogen Bonds',
            'temperature': 'System Temperature',
            'pressure': 'System Pressure',
            'potential_energy': 'Potential Energy'
        }
        return titles.get(self.analysis_type, Path(self.filename).stem.replace('_', ' ').title())
    
    def _extract_quoted(self, line):
        """Extract text between quotes"""
        match = re.search(r'"([^"]*)"', line)
        return match.group(1) if match else ""
    
    def _calculate_cdr_metrics(self):
        """Calculate CDR-specific metrics"""
        if not self.data:
            return
        
        columns = list(zip(*self.data))
        
        self.metrics = {
            'data_points': len(self.data),
            'time_range': [columns[0][0], columns[0][-1]] if columns else [0, 0],
            'analysis_type': self.analysis_type,
            'columns': []
        }
        
        for i, col in enumerate(columns):
            if i == 0:  # Time column
                continue
                
            col_data = list(col)
            col_metrics = {
                'index': i,
                'name': self.legends[i-1] if i-1 < len(self.legends) else f'Column {i}',
                'mean': mean(col_data),
                'std': stdev(col_data) if len(col_data) > 1 else 0,
                'min': min(col_data),
                'max': max(col_data),
                'median': median(col_data),
                'final_value': col_data[-1],
                'initial_value': col_data[0],
                'change': col_data[-1] - col_data[0],
                'stability_score': self._calculate_cdr_stability(col_data),
                'flexibility_assessment': self._assess_cdr_flexibility(col_data),
                'trend': self._calculate_trend(col_data)
            }
            self.metrics['columns'].append(col_metrics)
    
    def _calculate_cdr_stability(self, data):
        """Calculate CDR-specific stability score"""
        if not data or len(data) < 2:
            return 0
        
        # Use last 30% of data for CDR stability assessment
        stable_window = max(20, len(data) // 3)
        recent_data = data[-stable_window:]
        
        if len(recent_data) < 2:
            return 0
        
        cv = stdev(recent_data) / abs(mean(recent_data)) if mean(recent_data) != 0 else float('inf')
        
        # CDR-specific scoring
        if self.analysis_type == 'cdr_rmsd':
            # For RMSD: lower values and lower CV = better
            mean_val = mean(recent_data)
            if mean_val < 0.1 and cv < 0.1:
                return 95
            elif mean_val < 0.2 and cv < 0.15:
                return 85
            elif mean_val < 0.3 and cv < 0.2:
                return 70
            else:
                return max(0, 70 - min(mean_val * 100, 50))
        
        elif self.analysis_type == 'cdr_rmsf':
            # For RMSF: moderate values (0.1-0.3 nm) are good for CDRs
            mean_val = mean(recent_data)
            if 0.1 <= mean_val <= 0.3 and cv < 0.3:
                return 90
            elif 0.05 <= mean_val <= 0.5 and cv < 0.4:
                return 75
            else:
                return max(30, 75 - min(abs(mean_val - 0.2) * 200, 45))
        
        else:
            # General stability score
            if cv < 0.05:
                return 95
            elif cv < 0.1:
                return 85
            elif cv < 0.2:
                return 70
            else:
                return max(0, 70 - min(cv * 100, 70))
    
    def _assess_cdr_flexibility(self, data):
        """Assess CDR flexibility characteristics"""
        if not data or self.analysis_type != 'cdr_rmsf':
            return "N/A"
        
        mean_rmsf = mean(data)
        max_rmsf = max(data)
        
        if mean_rmsf < 0.1:
            return "Rigid"
        elif mean_rmsf < 0.2:
            return "Moderately Flexible"
        elif mean_rmsf < 0.4:
            return "Highly Flexible"
        else:
            return "Extremely Flexible"
    
    def _calculate_trend(self, data):
        """Calculate trend with CDR-specific interpretation"""
        if len(data) < 10:
            return "insufficient_data"
        
        q1_size = len(data) // 4
        first_q = mean(data[:q1_size])
        last_q = mean(data[-q1_size:])
        
        change_percent = ((last_q - first_q) / abs(first_q)) * 100 if first_q != 0 else 0
        
        if abs(change_percent) < 2:
            return "stable"
        elif change_percent > 5:
            return "increasing"
        elif change_percent < -5:
            return "decreasing"
        elif change_percent > 0:
            return "slightly_increasing"
        else:
            return "slightly_decreasing"
    
    def to_plotly_format(self):
        """Convert to Plotly format with CDR-specific styling"""
        if not self.data:
            return None
        
        columns = list(zip(*self.data))
        traces = []
        x_data = list(columns[0])
        
        # CDR-specific color schemes
        cdr_colors = {
            'cdr_rmsd': '#e74c3c',
            'cdr_rmsf': '#9b59b6',
            'cdr_sasa': '#3498db',
            'cdr_gyration': '#f39c12'
        }
        
        base_color = cdr_colors.get(self.analysis_type, '#2c3e50')
        
        if len(columns) == 2:
            trace = {
                'x': x_data,
                'y': list(columns[1]),
                'type': 'scatter',
                'mode': 'lines',
                'name': self.legends[0] if self.legends else 'CDR Data',
                'line': {'width': 3, 'color': base_color}
            }
            traces.append(trace)
        else:
            colors = ['#e74c3c', '#3498db', '#2ecc71', '#f39c12', '#9b59b6', '#1abc9c']
            for i in range(1, len(columns)):
                name = self.legends[i-1] if i-1 < len(self.legends) else f'Data {i}'
                trace = {
                    'x': x_data,
                    'y': list(columns[i]),
                    'type': 'scatter',
                    'mode': 'lines',
                    'name': name,
                    'line': {'width': 2, 'color': colors[(i-1) % len(colors)]}
                }
                traces.append(trace)
        
        layout = {
            'title': {'text': self.title, 'font': {'size': 20, 'color': '#2c3e50'}},
            'xaxis': {'title': self.xlabel, 'gridcolor': '#ecf0f1'},
            'yaxis': {'title': self.ylabel, 'gridcolor': '#ecf0f1'},
            'hovermode': 'x unified',
            'showlegend': len(traces) > 1,
            'height': 500,
            'template': 'plotly_white',
            'plot_bgcolor': 'rgba(0,0,0,0)',
            'paper_bgcolor': 'rgba(0,0,0,0)'
        }
        
        return {'data': traces, 'layout': layout}


class CDRDashboard:
    """Specialized dashboard for CDR analysis"""
    
    def __init__(self, parsers):
        self.parsers = parsers
        self.cdr_residues = [31,32,33,34,35,52,54,55,56,100,101,102,103,104,105,106,151,152,169,170,173,211,212,213,214,216]
        self.metrics = {}
        self.scores = {}
        self._calculate_cdr_metrics()
    
    def _calculate_cdr_metrics(self):
        """Calculate comprehensive CDR-specific metrics"""
        
        self.metrics = {
            'simulation_time': 0,
            'cdr_residue_count': len(self.cdr_residues),
            'analyses_count': len(self.parsers),
            'cdr_stability_scores': {},
            'cdr_values': {},
            'protein_comparison': {},
            'binding_indicators': {}
        }
        
        # Process each analysis
        for parser in self.parsers:
            filename = Path(parser.filename).stem
            analysis_type = parser.analysis_type
            
            if parser.metrics:
                # Update simulation time
                time_range = parser.metrics['time_range']
                sim_time = time_range[1] - time_range[0]
                if sim_time > self.metrics['simulation_time']:
                    self.metrics['simulation_time'] = sim_time
                
                # Store CDR-specific metrics
                for col_metric in parser.metrics['columns']:
                    key = f"{analysis_type}_{col_metric['name']}"
                    self.metrics['cdr_stability_scores'][key] = col_metric['stability_score']
                    self.metrics['cdr_values'][key] = {
                        'mean': col_metric['mean'],
                        'final': col_metric['final_value'],
                        'initial': col_metric['initial_value'],
                        'change': col_metric['change'],
                        'trend': col_metric['trend'],
                        'flexibility': col_metric.get('flexibility_assessment', 'N/A')
                    }
        
        self._calculate_cdr_scores()
    
    def _calculate_cdr_scores(self):
        """Calculate CDR-specific quality scores"""
        
        scores = {}
        
        # 1. CDR Structural Stability (based on CDR RMSD)
        cdr_rmsd_stability = self._get_stability_score('cdr_rmsd')
        cdr_rmsd_value = self._get_final_value('cdr_rmsd')
        
        if cdr_rmsd_stability is not None:
            scores['cdr_structural_stability'] = cdr_rmsd_stability
        else:
            scores['cdr_structural_stability'] = 0
        
        # 2. CDR Flexibility Assessment (based on CDR RMSF)
        cdr_rmsf_stability = self._get_stability_score('cdr_rmsf')
        cdr_rmsf_mean = self._get_mean_value('cdr_rmsf')
        
        if cdr_rmsf_stability is not None and cdr_rmsf_mean is not None:
            # Optimal RMSF for CDRs is 0.1-0.3 nm
            flexibility_score = cdr_rmsf_stability
            if 0.1 <= cdr_rmsf_mean <= 0.3:
                flexibility_score += 20
            scores['cdr_flexibility'] = min(100, flexibility_score)
        else:
            scores['cdr_flexibility'] = 0
        
        # 3. CDR Accessibility (based on CDR SASA)
        cdr_sasa_stability = self._get_stability_score('cdr_sasa')
        scores['cdr_accessibility'] = cdr_sasa_stability if cdr_sasa_stability is not None else 0
        
        # 4. CDR Compactness (based on CDR gyration)
        cdr_gyration_stability = self._get_stability_score('cdr_gyration')
        scores['cdr_compactness'] = cdr_gyration_stability if cdr_gyration_stability is not None else 0
        
        # 5. Relative CDR Performance (comparison with protein)
        protein_rmsd = self._get_final_value('protein_rmsd')
        cdr_rmsd = self._get_final_value('cdr_rmsd')
        
        if protein_rmsd is not None and cdr_rmsd is not None:
            # CDR should be somewhat more flexible than overall protein
            relative_flexibility = cdr_rmsd / protein_rmsd if protein_rmsd > 0 else 1
            if 1.2 <= relative_flexibility <= 2.0:
                scores['relative_performance'] = 90
            elif 1.0 <= relative_flexibility <= 2.5:
                scores['relative_performance'] = 75
            else:
                scores['relative_performance'] = max(50, 75 - abs(relative_flexibility - 1.5) * 30)
        else:
            scores['relative_performance'] = 0
        
        # 6. Binding Site Dynamics (based on hydrogen bonds and overall stability)
        hbond_stability = self._get_stability_score('hydrogen_bonds')
        hbond_mean = self._get_mean_value('hydrogen_bonds')
        
        if hbond_stability is not None and hbond_mean is not None:
            binding_score = hbond_stability * 0.7
            if hbond_mean > 800:
                binding_score += 30
            elif hbond_mean > 600:
                binding_score += 20
            scores['binding_dynamics'] = min(100, binding_score)
        else:
            scores['binding_dynamics'] = 0
        
        # 7. Overall CDR Quality Score
        cdr_scores = [scores['cdr_structural_stability'], scores['cdr_flexibility'], 
                     scores['cdr_accessibility'], scores['cdr_compactness']]
        cdr_scores = [s for s in cdr_scores if s > 0]
        
        if cdr_scores:
            scores['overall_cdr_quality'] = mean(cdr_scores)
        else:
            scores['overall_cdr_quality'] = 0
        
        self.scores = scores
    
    def _get_stability_score(self, analysis_type):
        """Get stability score for a specific analysis type"""
        for key, score in self.metrics['cdr_stability_scores'].items():
            if analysis_type.lower() in key.lower():
                return score
        return None
    
    def _get_final_value(self, analysis_type):
        """Get final value for a specific analysis type"""
        for key, values in self.metrics['cdr_values'].items():
            if analysis_type.lower() in key.lower():
                return values['final']
        return None
    
    def _get_mean_value(self, analysis_type):
        """Get mean value for a specific analysis type"""
        for key, values in self.metrics['cdr_values'].items():
            if analysis_type.lower() in key.lower():
                return values['mean']
        return None


def generate_cdr_dashboard_html(dashboard):
    """Generate CDR-specific dashboard HTML"""
    
    scores = dashboard.scores
    metrics = dashboard.metrics
    
    # CDR-specific score cards
    score_cards = []
    
    cdr_score_definitions = {
        'overall_cdr_quality': {'title': 'Overall CDR Quality', 'icon': '🎯', 'desc': 'Comprehensive CDR dynamics assessment'},
        'cdr_structural_stability': {'title': 'CDR Structural Stability', 'icon': '🏗️', 'desc': 'CDR backbone preservation during simulation'},
        'cdr_flexibility': {'title': 'CDR Flexibility', 'icon': '🌊', 'desc': 'Optimal flexibility for antigen binding'},
        'cdr_accessibility': {'title': 'CDR Accessibility', 'icon': '🔓', 'desc': 'Solvent accessibility of binding regions'},
        'cdr_compactness': {'title': 'CDR Compactness', 'icon': '📦', 'desc': 'Structural compactness of CDR regions'},
        'relative_performance': {'title': 'Relative Performance', 'icon': '⚖️', 'desc': 'CDR dynamics relative to protein backbone'},
        'binding_dynamics': {'title': 'Binding Dynamics', 'icon': '🔗', 'desc': 'Hydrogen bonding and interaction stability'}
    }
    
    for key, score in scores.items():
        if key in cdr_score_definitions:
            info = cdr_score_definitions[key]
            color_class = get_score_color_class(score)
            grade = get_score_grade(score)
            
            card_html = f"""
            <div class="score-card {color_class}">
                <div class="score-header">
                    <span class="score-icon">{info['icon']}</span>
                    <span class="score-title">{info['title']}</span>
                </div>
                <div class="score-value">{score:.1f}</div>
                <div class="score-grade">{grade}</div>
                <div class="score-desc">{info['desc']}</div>
                <div class="progress-bar">
                    <div class="progress-fill" style="width: {score}%"></div>
                </div>
            </div>
            """
            score_cards.append(card_html)
    
    # CDR residue information
    cdr_info = f"""
    <div class="cdr-residue-info">
        <h3>🧬 CDR Residue Information (4G6K Gevokizumab)</h3>
        <div class="residue-details">
            <div class="residue-stat">
                <strong>Total CDR Residues:</strong> {metrics['cdr_residue_count']}
            </div>
            <div class="residue-list">
                <strong>CDR Residues:</strong> 
                <span class="residue-numbers">{', '.join(map(str, dashboard.cdr_residues))}</span>
            </div>
        </div>
    </div>
    """
    
    # Summary statistics
    summary_stats = f"""
    <div class="summary-stats">
        <div class="stat-item">
            <div class="stat-number">{metrics['analyses_count']}</div>
            <div class="stat-label">CDR Analyses</div>
        </div>
        <div class="stat-item">
            <div class="stat-number">{metrics['simulation_time']:.1f}</div>
            <div class="stat-label">Simulation Time (ns)</div>
        </div>
        <div class="stat-item">
            <div class="stat-number">{metrics['cdr_residue_count']}</div>
            <div class="stat-label">CDR Residues</div>
        </div>
    </div>
    """
    
    # Radar chart data
    radar_data = {
        'labels': ['Structural', 'Flexibility', 'Accessibility', 'Compactness', 'Relative', 'Binding'],
        'values': [
            scores.get('cdr_structural_stability', 0),
            scores.get('cdr_flexibility', 0),
            scores.get('cdr_accessibility', 0),
            scores.get('cdr_compactness', 0),
            scores.get('relative_performance', 0),
            scores.get('binding_dynamics', 0)
        ]
    }
    
    dashboard_html = f"""
    <div class="dashboard">
        <div class="dashboard-header">
            <h1>🧬 CDR Dynamics Analysis Dashboard</h1>
            <p>Specialized analysis for antibody CDR regions based on 4G6K gevokizumab structure</p>
        </div>
        
        {cdr_info}
        {summary_stats}
        
        <div class="score-grid">
            {"".join(score_cards)}
        </div>
        
        <div class="radar-container">
            <h3>CDR Quality Assessment Radar</h3>
            <div id="radarChart"></div>
        </div>
        
        <div class="insights">
            <h3>💡 CDR Analysis Insights</h3>
            {generate_cdr_insights(dashboard)}
        </div>
    </div>
    
    <script>
        var radarData = {json.dumps(radar_data)};
        var radarTrace = {{
            type: 'scatterpolar',
            r: radarData.values,
            theta: radarData.labels,
            fill: 'toself',
            name: 'CDR Quality Scores',
            line: {{ color: '#e74c3c' }},
            fillcolor: 'rgba(231, 76, 60, 0.3)'
        }};
        
        var radarLayout = {{
            polar: {{
                radialaxis: {{
                    visible: true,
                    range: [0, 100],
                    tickfont: {{ size: 12 }},
                    gridcolor: 'rgba(0,0,0,0.1)'
                }},
                angularaxis: {{
                    tickfont: {{ size: 14, color: '#2c3e50' }}
                }}
            }},
            showlegend: false,
            height: 400,
            font: {{ family: 'Arial, sans-serif' }}
        }};
        
        Plotly.newPlot('radarChart', [radarTrace], radarLayout, {{responsive: true}});
    </script>
    """
    
    return dashboard_html


def generate_cdr_insights(dashboard):
    """Generate CDR-specific insights"""
    insights = []
    scores = dashboard.scores
    metrics = dashboard.metrics
    
    # Overall CDR assessment
    overall = scores.get('overall_cdr_quality', 0)
    if overall >= 85:
        insights.append("✅ <strong>Excellent CDR dynamics</strong> - All CDR regions show optimal flexibility and stability for antigen binding")
    elif overall >= 70:
        insights.append("✅ <strong>Good CDR performance</strong> - CDR regions maintain appropriate dynamics with minor concerns")
    elif overall >= 50:
        insights.append("⚠️ <strong>Moderate CDR quality</strong> - Some CDR stability issues detected, consider analysis of specific loops")
    else:
        insights.append("❌ <strong>Poor CDR dynamics</strong> - Significant CDR stability issues require investigation")
    
    # CDR-specific insights
    flexibility_score = scores.get('cdr_flexibility', 0)
    if flexibility_score >= 80:
        insights.append("✅ <strong>Optimal CDR flexibility</strong> - CDR regions show ideal flexibility for antigen recognition")
    elif flexibility_score < 50:
        insights.append("⚠️ <strong>CDR flexibility concerns</strong> - CDR regions may be too rigid or overly flexible for optimal binding")
    
    accessibility_score = scores.get('cdr_accessibility', 0)
    if accessibility_score >= 80:
        insights.append("✅ <strong>Good CDR accessibility</strong> - Binding sites remain accessible throughout simulation")
    elif accessibility_score < 60:
        insights.append("⚠️ <strong>CDR accessibility issues</strong> - Reduced surface accessibility may affect antigen binding")
    
    structural_score = scores.get('cdr_structural_stability', 0)
    if structural_score >= 85:
        insights.append("✅ <strong>Stable CDR structure</strong> - CDR backbone maintains structural integrity")
    elif structural_score < 60:
        insights.append("⚠️ <strong>CDR structural instability</strong> - High RMSD in CDR regions suggests conformational changes")
    
    relative_score = scores.get('relative_performance', 0)
    if relative_score >= 80:
        insights.append("✅ <strong>Balanced CDR-protein dynamics</strong> - CDR flexibility appropriately balanced with overall protein stability")
    elif relative_score < 50:
        insights.append("⚠️ <strong>Imbalanced dynamics</strong> - CDR behavior significantly different from overall protein dynamics")
    
    binding_score = scores.get('binding_dynamics', 0)
    if binding_score >= 80:
        insights.append("✅ <strong>Strong binding characteristics</strong> - Stable intermolecular interactions suggest good binding affinity")
    elif binding_score < 50:
        insights.append("⚠️ <strong>Weak binding indicators</strong> - Low hydrogen bonding suggests unstable antibody-antigen complex")
    
    # Add residue-specific recommendations
    insights.append("💡 <strong>Analysis tip:</strong> Focus on individual CDR loops (H1, H2, H3, L1, L2, L3) for detailed loop-specific behavior assessment")
    
    return "<ul>" + "".join(f"<li>{insight}</li>" for insight in insights) + "</ul>"


def get_score_color_class(score):
    """Get CSS class based on score"""
    if score >= 90:
        return "excellent"
    elif score >= 80:
        return "good"
    elif score >= 70:
        return "fair"
    elif score >= 60:
        return "poor"
    else:
        return "critical"


def get_score_grade(score):
    """Get letter grade based on score"""
    if score >= 90:
        return "A+"
    elif score >= 85:
        return "A"
    elif score >= 80:
        return "B+"
    elif score >= 75:
        return "B"
    elif score >= 70:
        return "C+"
    elif score >= 65:
        return "C"
    elif score >= 60:
        return "D"
    else:
        return "F"


def generate_cdr_html_report(xvg_files, output_file='cdr_analysis_dashboard.html'):
    """Generate CDR-focused HTML report"""
    
    plots_data = []
    parsers = []
    
    # Process each XVG file
    for xvg_file in sorted(xvg_files):
        try:
            parser = CDRAnalysisParser(xvg_file)
            plot_data = parser.to_plotly_format()
            
            if plot_data:
                plots_data.append(plot_data)
                parsers.append(parser)
                print(f"✓ Processed: {xvg_file} [{parser.analysis_type}]")
            else:
                print(f"✗ No data found in: {xvg_file}")
                
        except Exception as e:
            print(f"✗ Error processing {xvg_file}: {str(e)}")
    
    if not parsers:
        print("No plots to visualize!")
        return
    
    # Create CDR dashboard
    dashboard = CDRDashboard(parsers)
    dashboard_html = generate_cdr_dashboard_html(dashboard)
    
    # File information with analysis types
    file_info = []
    for parser in parsers:
        file_info.append({
            'filename': os.path.basename(parser.filename),
            'title': parser.title,
            'analysis_type': parser.analysis_type,
            'command': parser.metadata.get('command', 'N/A'),
            'data_points': parser.metrics['data_points']
        })
    
    # Enhanced CSS for CDR-specific styling
    cdr_css = """
        .cdr-residue-info {
            background: linear-gradient(135deg, #fff5f5 0%, #ffe8e8 100%);
            border-left: 5px solid #e74c3c;
            border-radius: 10px;
            padding: 20px;
            margin: 20px 0;
        }
        
        .cdr-residue-info h3 {
            color: #c0392b;
            margin: 0 0 15px 0;
            font-size: 1.3em;
        }
        
        .residue-details {
            display: flex;
            flex-direction: column;
            gap: 10px;
        }
        
        .residue-stat {
            font-size: 1.1em;
            color: #2c3e50;
        }
        
        .residue-list {
            color: #2c3e50;
        }
        
        .residue-numbers {
            font-family: 'Courier New', monospace;
            background: #ecf0f1;
            padding: 5px 10px;
            border-radius: 5px;
            font-size: 0.9em;
            color: #e74c3c;
            font-weight: bold;
        }
        
        .analysis-type-badge {
            display: inline-block;
            background: linear-gradient(135deg, #3498db, #2ecc71);
            color: white;
            padding: 4px 8px;
            border-radius: 12px;
            font-size: 0.8em;
            font-weight: bold;
            margin-left: 10px;
        }
        
        .cdr-specific {
            border-left: 4px solid #e74c3c !important;
        }
        
        .protein-specific {
            border-left: 4px solid #3498db !important;
        }
        
        .comparison-specific {
            border-left: 4px solid #f39c12 !important;
        }
    """
    
    # Create complete HTML with CDR focus
    html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>CDR Dynamics Analysis Dashboard</title>
    <meta charset="utf-8">
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
            margin: 0;
            padding: 0;
            background: linear-gradient(135deg, #c0392b 0%, #e74c3c 100%);
            min-height: 100vh;
        }}
        
        .main-container {{
            max-width: 1400px;
            margin: 0 auto;
            padding: 20px;
        }}
        
        .dashboard {{
            background: white;
            border-radius: 15px;
            padding: 30px;
            margin-bottom: 30px;
            box-shadow: 0 10px 30px rgba(0,0,0,0.1);
        }}
        
        .dashboard-header {{
            text-align: center;
            margin-bottom: 30px;
        }}
        
        .dashboard-header h1 {{
            color: #c0392b;
            font-size: 2.5em;
            margin: 0 0 10px 0;
        }}
        
        .dashboard-header p {{
            color: #7f8c8d;
            font-size: 1.2em;
            margin: 0;
        }}
        
        .summary-stats {{
            display: flex;
            justify-content: center;
            gap: 40px;
            margin: 30px 0;
            padding: 20px;
            background: linear-gradient(135deg, #f8f9fa 0%, #e9ecef 100%);
            border-radius: 10px;
        }}
        
        .stat-item {{
            text-align: center;
        }}
        
        .stat-number {{
            font-size: 2em;
            font-weight: bold;
            color: #e74c3c;
        }}
        
        .stat-label {{
            color: #6c757d;
            font-size: 0.9em;
            margin-top: 5px;
        }}
        
        .score-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
            gap: 20px;
            margin: 30px 0;
        }}
        
        .score-card {{
            background: white;
            border-radius: 12px;
            padding: 20px;
            text-align: center;
            box-shadow: 0 4px 15px rgba(0,0,0,0.1);
            transition: transform 0.3s, box-shadow 0.3s;
            border-left: 5px solid;
        }}
        
        .score-card:hover {{
            transform: translateY(-5px);
            box-shadow: 0 8px 25px rgba(0,0,0,0.15);
        }}
        
        .score-card.excellent {{ border-left-color: #27ae60; }}
        .score-card.good {{ border-left-color: #2ecc71; }}
        .score-card.fair {{ border-left-color: #f39c12; }}
        .score-card.poor {{ border-left-color: #e67e22; }}
        .score-card.critical {{ border-left-color: #e74c3c; }}
        
        .score-header {{
            display: flex;
            align-items: center;
            justify-content: center;
            gap: 10px;
            margin-bottom: 15px;
        }}
        
        .score-icon {{
            font-size: 1.5em;
        }}
        
        .score-title {{
            font-weight: bold;
            color: #2c3e50;
        }}
        
        .score-value {{
            font-size: 2.5em;
            font-weight: bold;
            color: #e74c3c;
            margin: 10px 0;
        }}
        
        .score-grade {{
            font-size: 1.2em;
            font-weight: bold;
            margin-bottom: 10px;
        }}
        
        .score-card.excellent .score-grade {{ color: #27ae60; }}
        .score-card.good .score-grade {{ color: #2ecc71; }}
        .score-card.fair .score-grade {{ color: #f39c12; }}
        .score-card.poor .score-grade {{ color: #e67e22; }}
        .score-card.critical .score-grade {{ color: #e74c3c; }}
        
        .score-desc {{
            color: #6c757d;
            font-size: 0.9em;
            margin-bottom: 15px;
        }}
        
        .progress-bar {{
            width: 100%;
            height: 6px;
            background-color: #ecf0f1;
            border-radius: 3px;
            overflow: hidden;
        }}
        
        .progress-fill {{
            height: 100%;
            background: linear-gradient(90deg, #e74c3c, #27ae60);
            transition: width 0.3s ease;
        }}
        
        .radar-container {{
            background: #f8f9fa;
            border-radius: 10px;
            padding: 20px;
            margin: 30px 0;
        }}
        
        .radar-container h3 {{
            text-align: center;
            color: #2c3e50;
            margin: 0 0 20px 0;
        }}
        
        .insights {{
            background: linear-gradient(135deg, #fff5f5 0%, #ffe8e8 100%);
            border-radius: 10px;
            padding: 20px;
            margin: 30px 0;
            border-left: 5px solid #e74c3c;
        }}
        
        .insights h3 {{
            color: #c0392b;
            margin: 0 0 15px 0;
        }}
        
        .insights ul {{
            margin: 0;
            padding-left: 20px;
        }}
        
        .insights li {{
            margin: 10px 0;
            line-height: 1.6;
        }}
        
        .plots-section {{
            background: white;
            border-radius: 15px;
            padding: 30px;
            box-shadow: 0 10px 30px rgba(0,0,0,0.1);
        }}
        
        .plots-header {{
            text-align: center;
            margin-bottom: 30px;
            border-bottom: 2px solid #ecf0f1;
            padding-bottom: 20px;
        }}
        
        .plots-header h2 {{
            color: #c0392b;
            font-size: 2em;
            margin: 0;
        }}
        
        .toc {{
            background: linear-gradient(135deg, #f8f9fa 0%, #e9ecef 100%);
            padding: 25px;
            border-radius: 10px;
            margin-bottom: 30px;
        }}
        
        .toc h3 {{
            color: #2c3e50;
            margin-top: 0;
            font-size: 1.3em;
        }}
        
        .toc ul {{
            list-style-type: none;
            padding-left: 0;
            margin: 0;
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
            gap: 10px;
        }}
        
        .toc li {{
            margin: 0;
        }}
        
        .toc a {{
            color: #e74c3c;
            text-decoration: none;
            padding: 12px 15px;
            display: block;
            border-radius: 5px;
            transition: all 0.3s;
            border: 1px solid transparent;
            background: white;
            position: relative;
        }}
        
        .toc a:hover {{
            background: #e74c3c;
            color: white;
            transform: translateX(5px);
        }}
        
        .plot-container {{
            background: white;
            margin: 30px 0;
            padding: 25px;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            transition: transform 0.2s;
        }}
        
        .plot-container:hover {{
            transform: translateY(-2px);
            box-shadow: 0 4px 20px rgba(0,0,0,0.15);
        }}
        
        .plot-info {{
            background: linear-gradient(135deg, #ecf0f1 0%, #e8eef2 100%);
            padding: 15px 20px;
            border-radius: 8px;
            margin-bottom: 20px;
            font-size: 14px;
            border-left: 4px solid #e74c3c;
        }}
        
        .plot-info strong {{
            color: #2c3e50;
            font-weight: 600;
        }}
        
        .plot-info code {{
            background-color: #34495e;
            color: #ecf0f1;
            padding: 2px 6px;
            border-radius: 3px;
            font-family: 'Courier New', Courier, monospace;
            font-size: 12px;
        }}
        
        h2 {{
            color: #c0392b;
            border-bottom: 3px solid #e74c3c;
            padding-bottom: 10px;
            margin-bottom: 20px;
        }}
        
        .footer {{
            text-align: center;
            margin-top: 50px;
            padding: 30px;
            color: rgba(255,255,255,0.8);
            font-size: 14px;
        }}
        
        .back-to-top {{
            position: fixed;
            bottom: 20px;
            right: 20px;
            background: linear-gradient(135deg, #e74c3c, #c0392b);
            color: white;
            padding: 15px;
            border-radius: 50%;
            text-decoration: none;
            font-size: 20px;
            opacity: 0.8;
            transition: opacity 0.3s, transform 0.3s;
            box-shadow: 0 4px 15px rgba(0,0,0,0.2);
        }}
        
        .back-to-top:hover {{
            opacity: 1;
            transform: translateY(-3px);
        }}
        
        {cdr_css}
        
        @media (max-width: 768px) {{
            .main-container {{
                padding: 10px;
            }}
            
            .summary-stats {{
                flex-direction: column;
                gap: 20px;
            }}
            
            .score-grid {{
                grid-template-columns: 1fr;
            }}
        }}
    </style>
</head>
<body>
    <div class="main-container">
        {dashboard_html}
        
        <div class="plots-section">
            <div class="plots-header">
                <h2>📊 CDR-Focused Analysis Results</h2>
            </div>
            
            <!-- Table of Contents -->
            <div class="toc">
                <h3>🗂️ CDR Analysis Overview</h3>
                <ul>
"""
    
    # Add TOC entries with analysis type badges
    for i, info in enumerate(file_info):
        analysis_class = ""
        if 'cdr' in info['analysis_type']:
            analysis_class = "cdr-specific"
        elif 'protein' in info['analysis_type']:
            analysis_class = "protein-specific"
        elif 'comparison' in info['analysis_type']:
            analysis_class = "comparison-specific"
            
        html_content += f'                    <li><a href="#plot{i}" class="{analysis_class}">{info["title"]} <span class="analysis-type-badge">{info["analysis_type"].replace("_", " ").title()}</span></a></li>\n'
    
    html_content += """                </ul>
            </div>
            
            <!-- Plots -->
"""
    
    # Add each plot with enhanced info
    for i, (plot_data, info) in enumerate(zip(plots_data, file_info)):
        analysis_class = ""
        if 'cdr' in info['analysis_type']:
            analysis_class = "cdr-specific"
        elif 'protein' in info['analysis_type']:
            analysis_class = "protein-specific"
        elif 'comparison' in info['analysis_type']:
            analysis_class = "comparison-specific"
            
        html_content += f"""
            <div class="plot-container {analysis_class}" id="plot{i}">
                <h2>{info['title']} <span class="analysis-type-badge">{info['analysis_type'].replace("_", " ").title()}</span></h2>
                <div class="plot-info">
                    <strong>📁 File:</strong> {info['filename']}<br>
                    <strong>🔬 Analysis Type:</strong> {info['analysis_type'].replace("_", " ").title()}<br>
                    <strong>📈 Data points:</strong> {info['data_points']:,}<br>
                    <strong>⚙️ GROMACS command:</strong> <code>{info['command']}</code>
                </div>
                <div id="plotDiv{i}"></div>
            </div>
"""
    
    html_content += """        </div>
    </div>
    
    <div class="footer">
        <p>🧬 CDR Dynamics Analysis Dashboard • Generated on """ + datetime.now().strftime('%Y-%m-%d %H:%M:%S') + """</p>
        <p>Specialized for 4G6K Gevokizumab CDR Analysis • Powered by Python & Plotly.js</p>
    </div>
    
    <a href="#" class="back-to-top">↑</a>
    
    <script>
"""
    
    # Add JavaScript to render plots
    for i, plot_data in enumerate(plots_data):
        html_content += f"""
        var plotData{i} = {json.dumps(plot_data)};
        Plotly.newPlot('plotDiv{i}', plotData{i}.data, plotData{i}.layout, {{responsive: true}});
"""
    
    html_content += """
        // Enhanced interactions
        document.querySelectorAll('a[href^="#"]').forEach(anchor => {
            anchor.addEventListener('click', function (e) {
                e.preventDefault();
                const target = document.querySelector(this.getAttribute('href'));
                if (target) {
                    target.scrollIntoView({ behavior: 'smooth', block: 'start' });
                }
            });
        });
        
        // Animate elements on scroll
        const observerOptions = {
            threshold: 0.1,
            rootMargin: '0px 0px -100px 0px'
        };
        
        const observer = new IntersectionObserver((entries) => {
            entries.forEach(entry => {
                if (entry.isIntersecting) {
                    const progressBars = entry.target.querySelectorAll('.progress-fill');
                    progressBars.forEach(bar => {
                        bar.style.width = bar.style.width || '0%';
                    });
                }
            });
        }, observerOptions);
        
        document.querySelectorAll('.score-card').forEach(card => {
            observer.observe(card);
        });
    </script>
</body>
</html>
"""
    
    # Write HTML file
    with open(output_file, 'w') as f:
        f.write(html_content)
    
    print(f"\n✓ CDR-focused dashboard generated: {output_file}")
    print(f"  CDR residues analyzed: {', '.join(map(str, dashboard.cdr_residues))}")
    print(f"  Dashboard includes {len(dashboard.scores)} CDR-specific quality metrics")


def main():
    """Main function for CDR analysis dashboard"""
    
    print("CDR-Focused GROMACS Analysis Dashboard")
    print("=" * 45)
    print("Specialized for antibody CDR dynamics analysis")
    print("Based on 4G6K gevokizumab structure")
    print()
    
    # Find all XVG files
    xvg_files = glob.glob("*.xvg")
    
    if not xvg_files:
        print("No XVG files found in the current directory!")
        return
    
    print(f"Found {len(xvg_files)} XVG files:")
    for f in sorted(xvg_files):
        print(f"  - {f}")
    print()
    
    # Generate CDR-focused dashboard
    print("Processing files and generating CDR dashboard...")
    generate_cdr_html_report(xvg_files)
    
    print("\n✅ CDR analysis dashboard complete!")
    print("🎯 CDR-specific features:")
    print("  • CDR structural stability scoring")
    print("  • CDR flexibility assessment")
    print("  • CDR accessibility analysis")
    print("  • Binding site dynamics evaluation")
    print("  • Comparative protein vs CDR analysis")
    print("  • Specialized CDR insights & recommendations")


if __name__ == "__main__":
    main()