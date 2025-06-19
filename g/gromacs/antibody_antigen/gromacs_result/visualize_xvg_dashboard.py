#!/usr/bin/env python3
"""
GROMACS XVG File Visualization Tool with Metrics Dashboard
Generates an interactive HTML report with key metrics dashboard and all XVG files
"""

import os
import re
import glob
import json
import math
from datetime import datetime
from pathlib import Path
from statistics import mean, stdev, median


class XVGParser:
    """Enhanced parser for GROMACS XVG files with metrics extraction"""
    
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
        self._parse()
        self._calculate_metrics()
    
    def _parse(self):
        """Parse XVG file header and data"""
        with open(self.filename, 'r') as f:
            lines = f.readlines()
        
        # Extract metadata from header
        for i, line in enumerate(lines):
            if line.startswith('#'):
                # Store metadata comments
                if 'Created by:' in line:
                    self.metadata['created_by'] = line.split('Created by:')[1].strip()
                elif 'Command line:' in line and i+1 < len(lines):
                    self.metadata['command'] = lines[i+1].strip('#').strip()
            elif line.startswith('@'):
                # Parse xmgrace commands
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
                # Parse data line
                try:
                    values = [float(x) for x in line.strip().split()]
                    if values:
                        self.data.append(values)
                except ValueError:
                    continue
        
        # If no title found, use filename
        if not self.title:
            self.title = Path(self.filename).stem.replace('_', ' ').title()
    
    def _extract_quoted(self, line):
        """Extract text between quotes"""
        match = re.search(r'"([^"]*)"', line)
        return match.group(1) if match else ""
    
    def _calculate_metrics(self):
        """Calculate key metrics for the data"""
        if not self.data:
            return
        
        # Transpose data for easier access
        columns = list(zip(*self.data))
        
        # Basic statistics for each column
        self.metrics = {
            'data_points': len(self.data),
            'time_range': [columns[0][0], columns[0][-1]] if columns else [0, 0],
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
                'stability_score': self._calculate_stability(col_data),
                'trend': self._calculate_trend(col_data)
            }
            self.metrics['columns'].append(col_metrics)
    
    def _calculate_stability(self, data):
        """Calculate stability score (0-100) based on coefficient of variation"""
        if not data or len(data) < 2:
            return 0
        
        # Use last 20% of data for stability assessment
        stable_window = max(10, len(data) // 5)
        recent_data = data[-stable_window:]
        
        if not recent_data or len(recent_data) < 2:
            return 0
        
        cv = stdev(recent_data) / abs(mean(recent_data)) if mean(recent_data) != 0 else float('inf')
        
        # Convert CV to stability score (lower CV = higher stability)
        if cv == 0:
            return 100
        elif cv < 0.01:
            return 95
        elif cv < 0.05:
            return 80
        elif cv < 0.1:
            return 60
        elif cv < 0.2:
            return 40
        else:
            return max(0, 40 - min(cv * 50, 40))
    
    def _calculate_trend(self, data):
        """Calculate trend direction and strength"""
        if len(data) < 10:
            return "insufficient_data"
        
        # Compare first and last quartiles
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
        """Convert data to Plotly.js format"""
        if not self.data:
            return None
        
        # Transpose data for easier access
        columns = list(zip(*self.data))
        
        traces = []
        x_data = list(columns[0])
        
        if len(columns) == 2:
            # Single Y column
            trace = {
                'x': x_data,
                'y': list(columns[1]),
                'type': 'scatter',
                'mode': 'lines',
                'name': self.legends[0] if self.legends else 'Data',
                'line': {'width': 2}
            }
            traces.append(trace)
        else:
            # Multiple Y columns
            for i in range(1, len(columns)):
                name = self.legends[i-1] if i-1 < len(self.legends) else f'Column {i}'
                trace = {
                    'x': x_data,
                    'y': list(columns[i]),
                    'type': 'scatter',
                    'mode': 'lines',
                    'name': name,
                    'line': {'width': 2}
                }
                traces.append(trace)
        
        layout = {
            'title': {'text': self.title, 'font': {'size': 20}},
            'xaxis': {'title': self.xlabel},
            'yaxis': {'title': self.ylabel},
            'hovermode': 'x unified',
            'showlegend': len(traces) > 1,
            'height': 500,
            'template': 'plotly_white'
        }
        
        return {'data': traces, 'layout': layout}


class MetricsDashboard:
    """Calculate overall simulation quality metrics and scores"""
    
    def __init__(self, parsers):
        self.parsers = parsers
        self.metrics = {}
        self.scores = {}
        self._calculate_dashboard_metrics()
    
    def _calculate_dashboard_metrics(self):
        """Calculate comprehensive metrics for all analyses"""
        
        # Initialize metrics
        self.metrics = {
            'simulation_time': 0,
            'total_data_points': 0,
            'analyses_count': len(self.parsers),
            'stability_scores': {},
            'key_values': {},
            'quality_indicators': {}
        }
        
        # Process each parser
        for parser in self.parsers:
            filename = Path(parser.filename).stem
            
            if parser.metrics:
                # Update total metrics
                self.metrics['total_data_points'] += parser.metrics['data_points']
                
                # Get simulation time range
                time_range = parser.metrics['time_range']
                sim_time = time_range[1] - time_range[0]
                if sim_time > self.metrics['simulation_time']:
                    self.metrics['simulation_time'] = sim_time
                
                # Store analysis-specific metrics
                for col_metric in parser.metrics['columns']:
                    key = f"{filename}_{col_metric['name']}"
                    self.metrics['stability_scores'][key] = col_metric['stability_score']
                    self.metrics['key_values'][key] = {
                        'mean': col_metric['mean'],
                        'final': col_metric['final_value'],
                        'trend': col_metric['trend']
                    }
        
        # Calculate quality scores
        self._calculate_quality_scores()
    
    def _calculate_quality_scores(self):
        """Calculate overall quality scores for the simulation"""
        
        scores = {}
        
        # 1. Structural Stability Score (based on RMSD)
        rmsd_stability = self._get_stability_score('rmsd')
        rmsd_final = self._get_final_value('rmsd')
        
        if rmsd_stability is not None and rmsd_final is not None:
            # Good RMSD: stable and < 0.3 nm
            rmsd_score = rmsd_stability * 0.7 + max(0, (0.5 - rmsd_final) / 0.5 * 30)
            scores['structural_stability'] = min(100, max(0, rmsd_score))
        else:
            scores['structural_stability'] = 0
        
        # 2. Thermodynamic Equilibrium Score
        temp_stability = self._get_stability_score('temperature')
        pressure_stability = self._get_stability_score('pressure')
        
        thermo_components = [s for s in [temp_stability, pressure_stability] if s is not None]
        scores['thermodynamic_equilibrium'] = mean(thermo_components) if thermo_components else 0
        
        # 3. Energetic Stability Score
        potential_stability = self._get_stability_score('potential')
        scores['energetic_stability'] = potential_stability if potential_stability is not None else 0
        
        # 4. Binding Affinity Score (based on hydrogen bonds and contacts)
        hbond_mean = self._get_mean_value('hbond')
        hbond_stability = self._get_stability_score('hbond')
        
        if hbond_mean is not None and hbond_stability is not None:
            # Good binding: stable hydrogen bonds > 800
            binding_score = hbond_stability * 0.6 + min(40, max(0, (hbond_mean - 600) / 10))
            scores['binding_affinity'] = min(100, max(0, binding_score))
        else:
            scores['binding_affinity'] = 0
        
        # 5. Compactness Score (based on radius of gyration)
        gyrate_stability = self._get_stability_score('gyrate')
        scores['compactness'] = gyrate_stability if gyrate_stability is not None else 0
        
        # 6. Surface Properties Score (based on SASA)
        sasa_stability = self._get_stability_score('sasa')
        scores['surface_properties'] = sasa_stability if sasa_stability is not None else 0
        
        # 7. Overall Simulation Quality Score
        all_scores = [score for score in scores.values() if score > 0]
        scores['overall_quality'] = mean(all_scores) if all_scores else 0
        
        self.scores = scores
    
    def _get_stability_score(self, analysis_type):
        """Get stability score for a specific analysis type"""
        for key, score in self.metrics['stability_scores'].items():
            if analysis_type.lower() in key.lower():
                return score
        return None
    
    def _get_final_value(self, analysis_type):
        """Get final value for a specific analysis type"""
        for key, values in self.metrics['key_values'].items():
            if analysis_type.lower() in key.lower():
                return values['final']
        return None
    
    def _get_mean_value(self, analysis_type):
        """Get mean value for a specific analysis type"""
        for key, values in self.metrics['key_values'].items():
            if analysis_type.lower() in key.lower():
                return values['mean']
        return None


def generate_dashboard_html(dashboard):
    """Generate HTML for the metrics dashboard"""
    
    scores = dashboard.scores
    metrics = dashboard.metrics
    
    # Create score cards
    score_cards = []
    
    score_definitions = {
        'overall_quality': {'title': 'Overall Quality', 'icon': '🎯', 'desc': 'Comprehensive simulation quality assessment'},
        'structural_stability': {'title': 'Structural Stability', 'icon': '🏗️', 'desc': 'Protein structure preservation (RMSD-based)'},
        'thermodynamic_equilibrium': {'title': 'Thermodynamic Equilibrium', 'icon': '⚖️', 'desc': 'Temperature and pressure stability'},
        'energetic_stability': {'title': 'Energetic Stability', 'icon': '⚡', 'desc': 'Potential energy convergence'},
        'binding_affinity': {'title': 'Binding Affinity', 'icon': '🔗', 'desc': 'Hydrogen bonding and molecular interactions'},
        'compactness': {'title': 'Compactness', 'icon': '📦', 'desc': 'Radius of gyration stability'},
        'surface_properties': {'title': 'Surface Properties', 'icon': '🌊', 'desc': 'Solvent accessible surface area stability'}
    }
    
    for key, score in scores.items():
        if key in score_definitions:
            info = score_definitions[key]
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
    
    # Create summary statistics
    summary_stats = f"""
    <div class="summary-stats">
        <div class="stat-item">
            <div class="stat-number">{metrics['analyses_count']}</div>
            <div class="stat-label">Analyses</div>
        </div>
        <div class="stat-item">
            <div class="stat-number">{metrics['simulation_time']:.1f}</div>
            <div class="stat-label">Simulation Time (ns)</div>
        </div>
        <div class="stat-item">
            <div class="stat-number">{metrics['total_data_points']:,}</div>
            <div class="stat-label">Data Points</div>
        </div>
    </div>
    """
    
    # Create radar chart data for overall assessment
    radar_data = {
        'labels': ['Structural', 'Thermodynamic', 'Energetic', 'Binding', 'Compactness', 'Surface'],
        'values': [
            scores.get('structural_stability', 0),
            scores.get('thermodynamic_equilibrium', 0),
            scores.get('energetic_stability', 0),
            scores.get('binding_affinity', 0),
            scores.get('compactness', 0),
            scores.get('surface_properties', 0)
        ]
    }
    
    dashboard_html = f"""
    <div class="dashboard">
        <div class="dashboard-header">
            <h1>🧬 Simulation Quality Dashboard</h1>
            <p>Comprehensive analysis of molecular dynamics simulation results</p>
        </div>
        
        {summary_stats}
        
        <div class="score-grid">
            {"".join(score_cards)}
        </div>
        
        <div class="radar-container">
            <h3>Quality Assessment Radar</h3>
            <div id="radarChart"></div>
        </div>
        
        <div class="insights">
            <h3>💡 Key Insights</h3>
            {generate_insights(dashboard)}
        </div>
    </div>
    
    <script>
        // Radar chart
        var radarData = {json.dumps(radar_data)};
        var radarTrace = {{
            type: 'scatterpolar',
            r: radarData.values,
            theta: radarData.labels,
            fill: 'toself',
            name: 'Quality Scores',
            line: {{ color: '#3498db' }},
            fillcolor: 'rgba(52, 152, 219, 0.3)'
        }};
        
        var radarLayout = {{
            polar: {{
                radialaxis: {{
                    visible: true,
                    range: [0, 100]
                }}
            }},
            showlegend: false,
            height: 400
        }};
        
        Plotly.newPlot('radarChart', [radarTrace], radarLayout, {{responsive: true}});
    </script>
    """
    
    return dashboard_html


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


def generate_insights(dashboard):
    """Generate textual insights based on the metrics"""
    insights = []
    scores = dashboard.scores
    
    # Overall assessment
    overall = scores.get('overall_quality', 0)
    if overall >= 85:
        insights.append("✅ <strong>Excellent simulation quality</strong> - All major indicators show stable, well-equilibrated system")
    elif overall >= 70:
        insights.append("✅ <strong>Good simulation quality</strong> - System shows acceptable stability with minor concerns")
    elif overall >= 50:
        insights.append("⚠️ <strong>Fair simulation quality</strong> - Some stability issues detected, consider extending simulation")
    else:
        insights.append("❌ <strong>Poor simulation quality</strong> - Significant stability issues require investigation")
    
    # Specific insights
    if scores.get('structural_stability', 0) < 60:
        insights.append("⚠️ <strong>Structural instability detected</strong> - High RMSD values suggest significant conformational changes")
    
    if scores.get('binding_affinity', 0) >= 80:
        insights.append("✅ <strong>Strong binding detected</strong> - Stable hydrogen bonding indicates good protein-ligand interaction")
    elif scores.get('binding_affinity', 0) < 50:
        insights.append("⚠️ <strong>Weak binding detected</strong> - Low hydrogen bond count suggests unstable complex")
    
    if scores.get('thermodynamic_equilibrium', 0) >= 85:
        insights.append("✅ <strong>System well equilibrated</strong> - Temperature and pressure show excellent stability")
    
    if scores.get('energetic_stability', 0) < 60:
        insights.append("⚠️ <strong>Energy convergence issues</strong> - Potential energy not fully stabilized")
    
    return "<ul>" + "".join(f"<li>{insight}</li>" for insight in insights) + "</ul>"


def generate_enhanced_html_report(xvg_files, output_file='xvg_dashboard_report.html'):
    """Generate comprehensive HTML report with dashboard and all visualizations"""
    
    plots_data = []
    parsers = []
    
    # Process each XVG file
    for xvg_file in sorted(xvg_files):
        try:
            parser = XVGParser(xvg_file)
            plot_data = parser.to_plotly_format()
            
            if plot_data:
                plots_data.append(plot_data)
                parsers.append(parser)
                print(f"✓ Processed: {xvg_file}")
            else:
                print(f"✗ No data found in: {xvg_file}")
                
        except Exception as e:
            print(f"✗ Error processing {xvg_file}: {str(e)}")
    
    if not parsers:
        print("No plots to visualize!")
        return
    
    # Create dashboard
    dashboard = MetricsDashboard(parsers)
    dashboard_html = generate_dashboard_html(dashboard)
    
    # File information
    file_info = []
    for parser in parsers:
        file_info.append({
            'filename': os.path.basename(parser.filename),
            'title': parser.title,
            'command': parser.metadata.get('command', 'N/A'),
            'data_points': parser.metrics['data_points']
        })
    
    # Create complete HTML
    html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>GROMACS Analysis Dashboard</title>
    <meta charset="utf-8">
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
            margin: 0;
            padding: 0;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
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
            color: #2c3e50;
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
            color: #3498db;
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
            color: #3498db;
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
            background: linear-gradient(90deg, #3498db, #2ecc71);
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
            background: linear-gradient(135deg, #e8f5e8 0%, #f0f8f0 100%);
            border-radius: 10px;
            padding: 20px;
            margin: 30px 0;
            border-left: 5px solid #27ae60;
        }}
        
        .insights h3 {{
            color: #2c3e50;
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
            color: #2c3e50;
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
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 10px;
        }}
        
        .toc li {{
            margin: 0;
        }}
        
        .toc a {{
            color: #3498db;
            text-decoration: none;
            padding: 10px 15px;
            display: block;
            border-radius: 5px;
            transition: all 0.3s;
            border: 1px solid transparent;
            background: white;
        }}
        
        .toc a:hover {{
            background: #3498db;
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
            border-left: 4px solid #3498db;
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
            color: #2c3e50;
            border-bottom: 3px solid #3498db;
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
            background: linear-gradient(135deg, #3498db, #2ecc71);
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
                <h2>📊 Detailed Analysis Plots</h2>
            </div>
            
            <!-- Table of Contents -->
            <div class="toc">
                <h3>🗂️ Analysis Overview</h3>
                <ul>
"""
    
    # Add TOC entries
    for i, info in enumerate(file_info):
        html_content += f'                    <li><a href="#plot{i}">{info["title"]} <small>({info["filename"]})</small></a></li>\n'
    
    html_content += """                </ul>
            </div>
            
            <!-- Plots -->
"""
    
    # Add each plot
    for i, (plot_data, info) in enumerate(zip(plots_data, file_info)):
        html_content += f"""
            <div class="plot-container" id="plot{i}">
                <h2>{info['title']}</h2>
                <div class="plot-info">
                    <strong>📁 File:</strong> {info['filename']}<br>
                    <strong>📈 Data points:</strong> {info['data_points']:,}<br>
                    <strong>⚙️ GROMACS command:</strong> <code>{info['command']}</code>
                </div>
                <div id="plotDiv{i}"></div>
            </div>
"""
    
    html_content += """        </div>
    </div>
    
    <div class="footer">
        <p>🧬 Advanced GROMACS Analysis Dashboard • Generated on """ + datetime.now().strftime('%Y-%m-%d %H:%M:%S') + """</p>
        <p>Powered by Python, Plotly.js, and Scientific Computing</p>
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
        // Smooth scrolling for navigation
        document.querySelectorAll('a[href^="#"]').forEach(anchor => {
            anchor.addEventListener('click', function (e) {
                e.preventDefault();
                const target = document.querySelector(this.getAttribute('href'));
                if (target) {
                    target.scrollIntoView({ behavior: 'smooth', block: 'start' });
                }
            });
        });
        
        // Animate progress bars on scroll
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
    
    print(f"\n✓ Enhanced dashboard report generated: {output_file}")
    print(f"  Dashboard includes {len(dashboard.scores)} quality metrics and comprehensive analysis")


def main():
    """Main function to process all XVG files and create dashboard"""
    
    print("GROMACS XVG Dashboard Visualization Tool")
    print("=" * 45)
    
    # Find all XVG files
    xvg_files = glob.glob("*.xvg")
    
    if not xvg_files:
        print("No XVG files found in the current directory!")
        return
    
    print(f"Found {len(xvg_files)} XVG files:")
    for f in sorted(xvg_files):
        print(f"  - {f}")
    print()
    
    # Generate enhanced HTML report with dashboard
    print("Processing files and generating dashboard...")
    generate_enhanced_html_report(xvg_files)
    
    print("\n✅ Dashboard visualization complete!")
    print("📊 Features included:")
    print("  • Quality scoring system")
    print("  • Stability assessment")
    print("  • Binding affinity analysis")
    print("  • Thermodynamic evaluation")
    print("  • Interactive radar chart")
    print("  • Automated insights")


if __name__ == "__main__":
    main()