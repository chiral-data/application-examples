#!/usr/bin/env python3
"""
GROMACS XVG File Visualization Tool
Generates an interactive HTML report with all XVG files in the current directory
"""

import os
import re
import glob
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from datetime import datetime
from pathlib import Path


class XVGParser:
    """Parser for GROMACS XVG files"""
    
    def __init__(self, filename):
        self.filename = filename
        self.title = ""
        self.xlabel = "X"
        self.ylabel = "Y"
        self.subtitle = ""
        self.legends = []
        self.data = None
        self.metadata = {}
        self._parse()
    
    def _parse(self):
        """Parse XVG file header and data"""
        with open(self.filename, 'r') as f:
            lines = f.readlines()
        
        # Extract metadata from header
        data_start = 0
        for i, line in enumerate(lines):
            if line.startswith('#'):
                # Store metadata comments
                if 'Created by:' in line:
                    self.metadata['created_by'] = line.split('Created by:')[1].strip()
                elif 'Command line:' in line:
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
            else:
                # Found first data line
                data_start = i
                break
        
        # Load numerical data
        self.data = np.loadtxt(self.filename, comments=['@', '#', '&'])
        
        # If no title found, use filename
        if not self.title:
            self.title = Path(self.filename).stem.replace('_', ' ').title()
    
    def _extract_quoted(self, line):
        """Extract text between quotes"""
        match = re.search(r'"([^"]*)"', line)
        return match.group(1) if match else ""
    
    def get_traces(self):
        """Generate Plotly traces for the data"""
        traces = []
        
        if self.data.ndim == 1:
            # Single data point
            return traces
        
        x = self.data[:, 0]
        
        if self.data.shape[1] == 2:
            # Single Y column
            trace = go.Scatter(
                x=x,
                y=self.data[:, 1],
                mode='lines',
                name=self.legends[0] if self.legends else 'Data',
                line=dict(width=2)
            )
            traces.append(trace)
        else:
            # Multiple Y columns
            for i in range(1, self.data.shape[1]):
                name = self.legends[i-1] if i-1 < len(self.legends) else f'Column {i}'
                trace = go.Scatter(
                    x=x,
                    y=self.data[:, i],
                    mode='lines',
                    name=name,
                    line=dict(width=2)
                )
                traces.append(trace)
        
        return traces


class XVGVisualizer:
    """Create interactive visualizations for XVG files"""
    
    def __init__(self):
        self.figures = []
        self.file_info = []
    
    def add_xvg_file(self, filename):
        """Add an XVG file to the visualization"""
        try:
            parser = XVGParser(filename)
            traces = parser.get_traces()
            
            if not traces:
                print(f"Warning: No data found in {filename}")
                return
            
            # Create figure
            fig = go.Figure()
            
            # Add all traces
            for trace in traces:
                fig.add_trace(trace)
            
            # Update layout
            fig.update_layout(
                title={
                    'text': parser.title,
                    'font': {'size': 20}
                },
                xaxis_title=parser.xlabel,
                yaxis_title=parser.ylabel,
                font=dict(size=14),
                hovermode='x unified',
                showlegend=True if len(traces) > 1 else False,
                height=500,
                template='plotly_white'
            )
            
            # Add subtitle as annotation if exists
            if parser.subtitle:
                fig.add_annotation(
                    xref="paper", yref="paper",
                    x=0.5, y=1.05,
                    text=parser.subtitle,
                    showarrow=False,
                    font=dict(size=12),
                    xanchor='center'
                )
            
            # Store figure and metadata
            self.figures.append(fig)
            self.file_info.append({
                'filename': os.path.basename(filename),
                'title': parser.title,
                'command': parser.metadata.get('command', 'N/A'),
                'data_points': len(parser.data)
            })
            
            print(f"✓ Processed: {filename}")
            
        except Exception as e:
            print(f"✗ Error processing {filename}: {str(e)}")
    
    def generate_html_report(self, output_file='xvg_visualization_report.html'):
        """Generate comprehensive HTML report with all visualizations"""
        
        if not self.figures:
            print("No figures to visualize!")
            return
        
        # Create HTML content
        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>GROMACS XVG Visualization Report</title>
    <meta charset="utf-8">
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        body {{
            font-family: Arial, sans-serif;
            margin: 0;
            padding: 20px;
            background-color: #f5f5f5;
        }}
        .header {{
            background-color: #2c3e50;
            color: white;
            padding: 30px;
            text-align: center;
            margin: -20px -20px 30px -20px;
            box-shadow: 0 2px 5px rgba(0,0,0,0.1);
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
        }}
        .plot-container {{
            background-color: white;
            margin: 20px 0;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 5px rgba(0,0,0,0.1);
        }}
        .plot-info {{
            background-color: #ecf0f1;
            padding: 15px;
            border-radius: 5px;
            margin-bottom: 15px;
            font-size: 14px;
        }}
        .plot-info strong {{
            color: #2c3e50;
        }}
        .toc {{
            background-color: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 5px rgba(0,0,0,0.1);
            margin-bottom: 30px;
        }}
        .toc h2 {{
            color: #2c3e50;
            margin-top: 0;
        }}
        .toc ul {{
            list-style-type: none;
            padding-left: 0;
        }}
        .toc li {{
            margin: 10px 0;
        }}
        .toc a {{
            color: #3498db;
            text-decoration: none;
            padding: 5px 10px;
            display: block;
            border-radius: 3px;
            transition: background-color 0.3s;
        }}
        .toc a:hover {{
            background-color: #ecf0f1;
        }}
        .footer {{
            text-align: center;
            margin-top: 50px;
            padding: 20px;
            color: #7f8c8d;
            font-size: 14px;
        }}
        h2 {{
            color: #2c3e50;
            border-bottom: 2px solid #3498db;
            padding-bottom: 10px;
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>GROMACS XVG Visualization Report</h1>
        <p>Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        <p>Total files analyzed: {len(self.figures)}</p>
    </div>
    
    <div class="container">
        <!-- Table of Contents -->
        <div class="toc">
            <h2>Table of Contents</h2>
            <ul>
"""
        
        # Add TOC entries
        for i, info in enumerate(self.file_info):
            html_content += f'                <li><a href="#plot{i}">{info["title"]} ({info["filename"]})</a></li>\n'
        
        html_content += """            </ul>
        </div>
        
        <!-- Plots -->
"""
        
        # Add each plot
        for i, (fig, info) in enumerate(zip(self.figures, self.file_info)):
            html_content += f"""
        <div class="plot-container" id="plot{i}">
            <h2>{info['title']}</h2>
            <div class="plot-info">
                <strong>File:</strong> {info['filename']}<br>
                <strong>Data points:</strong> {info['data_points']}<br>
                <strong>GROMACS command:</strong> <code>{info['command']}</code>
            </div>
            <div id="plotDiv{i}"></div>
        </div>
"""
        
        html_content += """    </div>
    
    <div class="footer">
        <p>This report was automatically generated from GROMACS XVG files using Python and Plotly.</p>
    </div>
    
    <script>
"""
        
        # Add JavaScript to render plots
        for i, fig in enumerate(self.figures):
            # Convert figure to JSON
            fig_json = fig.to_json()
            html_content += f"""
        var figure{i} = {fig_json};
        Plotly.newPlot('plotDiv{i}', figure{i}.data, figure{i}.layout, {{responsive: true}});
"""
        
        html_content += """    </script>
</body>
</html>
"""
        
        # Write HTML file
        with open(output_file, 'w') as f:
            f.write(html_content)
        
        print(f"\n✓ HTML report generated: {output_file}")
        print(f"  Open this file in a web browser to view the interactive plots.")


def main():
    """Main function to process all XVG files in current directory"""
    
    print("GROMACS XVG Visualization Tool")
    print("=" * 40)
    
    # Find all XVG files
    xvg_files = glob.glob("*.xvg")
    
    if not xvg_files:
        print("No XVG files found in the current directory!")
        return
    
    print(f"Found {len(xvg_files)} XVG files:")
    for f in sorted(xvg_files):
        print(f"  - {f}")
    print()
    
    # Create visualizer
    visualizer = XVGVisualizer()
    
    # Process each file
    print("Processing files...")
    for xvg_file in sorted(xvg_files):
        visualizer.add_xvg_file(xvg_file)
    
    # Generate HTML report
    print("\nGenerating HTML report...")
    visualizer.generate_html_report()
    
    print("\n✅ Visualization complete!")


if __name__ == "__main__":
    main()