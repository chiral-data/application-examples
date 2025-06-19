#!/usr/bin/env python3
"""
GROMACS XVG File Visualization Tool (No External Dependencies)
Generates an interactive HTML report with all XVG files using JavaScript
"""

import os
import re
import glob
import json
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
        self.data = []
        self.metadata = {}
        self._parse()
    
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


def generate_html_report(xvg_files, output_file='xvg_visualization_report.html'):
    """Generate comprehensive HTML report with all visualizations"""
    
    plots_data = []
    file_info = []
    
    # Process each XVG file
    for xvg_file in sorted(xvg_files):
        try:
            parser = XVGParser(xvg_file)
            plot_data = parser.to_plotly_format()
            
            if plot_data:
                plots_data.append(plot_data)
                file_info.append({
                    'filename': os.path.basename(xvg_file),
                    'title': parser.title,
                    'command': parser.metadata.get('command', 'N/A'),
                    'data_points': len(parser.data)
                })
                print(f"✓ Processed: {xvg_file}")
            else:
                print(f"✗ No data found in: {xvg_file}")
                
        except Exception as e:
            print(f"✗ Error processing {xvg_file}: {str(e)}")
    
    if not plots_data:
        print("No plots to visualize!")
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
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
            margin: 0;
            padding: 20px;
            background-color: #f5f5f5;
        }}
        .header {{
            background: linear-gradient(135deg, #2c3e50 0%, #3498db 100%);
            color: white;
            padding: 40px;
            text-align: center;
            margin: -20px -20px 30px -20px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        .header h1 {{
            margin: 0 0 10px 0;
            font-size: 2.5em;
        }}
        .header p {{
            margin: 5px 0;
            opacity: 0.9;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
        }}
        .plot-container {{
            background-color: white;
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
        .toc {{
            background-color: white;
            padding: 25px;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            margin-bottom: 30px;
        }}
        .toc h2 {{
            color: #2c3e50;
            margin-top: 0;
            font-size: 1.8em;
        }}
        .toc ul {{
            list-style-type: none;
            padding-left: 0;
            margin: 0;
        }}
        .toc li {{
            margin: 8px 0;
        }}
        .toc a {{
            color: #3498db;
            text-decoration: none;
            padding: 10px 15px;
            display: block;
            border-radius: 5px;
            transition: all 0.3s;
            border: 1px solid transparent;
        }}
        .toc a:hover {{
            background-color: #3498db;
            color: white;
            transform: translateX(5px);
        }}
        .footer {{
            text-align: center;
            margin-top: 50px;
            padding: 30px;
            color: #7f8c8d;
            font-size: 14px;
            border-top: 1px solid #e0e0e0;
        }}
        h2 {{
            color: #2c3e50;
            border-bottom: 3px solid #3498db;
            padding-bottom: 10px;
            margin-bottom: 20px;
        }}
        .back-to-top {{
            position: fixed;
            bottom: 20px;
            right: 20px;
            background-color: #3498db;
            color: white;
            padding: 10px 15px;
            border-radius: 5px;
            text-decoration: none;
            opacity: 0.8;
            transition: opacity 0.3s;
        }}
        .back-to-top:hover {{
            opacity: 1;
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>GROMACS XVG Visualization Report</h1>
        <p>Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        <p>Total files analyzed: <strong>{len(plots_data)}</strong></p>
    </div>
    
    <div class="container">
        <!-- Table of Contents -->
        <div class="toc">
            <h2>📊 Table of Contents</h2>
            <ul>
"""
    
    # Add TOC entries
    for i, info in enumerate(file_info):
        html_content += f'                <li><a href="#plot{i}">{info["title"]} <small>({info["filename"]})</small></a></li>\n'
    
    html_content += """            </ul>
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
    
    html_content += """    </div>
    
    <div class="footer">
        <p>This report was automatically generated from GROMACS XVG files using Python and Plotly.js</p>
        <p>© 2025 GROMACS Visualization Tool</p>
    </div>
    
    <a href="#" class="back-to-top">↑ Top</a>
    
    <script>
        // Plot data
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
    </script>
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
    
    print("GROMACS XVG Visualization Tool (No Dependencies Version)")
    print("=" * 55)
    
    # Find all XVG files
    xvg_files = glob.glob("*.xvg")
    
    if not xvg_files:
        print("No XVG files found in the current directory!")
        return
    
    print(f"Found {len(xvg_files)} XVG files:")
    for f in sorted(xvg_files):
        print(f"  - {f}")
    print()
    
    # Generate HTML report
    print("Processing files and generating HTML report...")
    generate_html_report(xvg_files)
    
    print("\n✅ Visualization complete!")


if __name__ == "__main__":
    main()