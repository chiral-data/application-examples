# GROMACS XVG Visualization Scripts

This directory contains three Python scripts for visualizing GROMACS XVG files in interactive HTML reports with varying levels of functionality.

## Scripts Available

### 1. `visualize_xvg_simple.py` (Recommended)
- **No external dependencies required** - uses only Python standard library
- Creates interactive plots using Plotly.js (loaded from CDN)
- Generates a comprehensive HTML report with all XVG files

### 2. `visualize_xvg.py` (Full-featured)
- Requires numpy and plotly Python packages
- More advanced parsing capabilities
- Install dependencies with: `pip install -r requirements.txt`

### 3. `visualize_xvg_dashboard.py` (Advanced Dashboard) ⭐ **RECOMMENDED**
- **No external dependencies required** - uses only Python standard library
- **Complete metrics dashboard** with quality scoring system
- **Automated insights** and stability assessment
- **Interactive radar chart** showing overall simulation quality
- **Professional scoring system** (A+ to F grades)
- **Binding affinity analysis** based on hydrogen bonds
- **Thermodynamic equilibrium** assessment
- **Structural stability** evaluation

## Usage

Simply run the script in a directory containing XVG files:

```bash
# For the advanced dashboard version (RECOMMENDED)
python3 visualize_xvg_dashboard.py

# For the simple version (no dependencies)
python3 visualize_xvg_simple.py

# For the full-featured version (requires numpy and plotly)
pip install -r requirements.txt
python3 visualize_xvg.py
```

## Output

The scripts generate HTML reports with different features:

### Dashboard Version (`xvg_dashboard_report.html`)
- **📊 Comprehensive Quality Dashboard** with scoring system
- **🎯 Overall Quality Score** (0-100 with letter grades A+ to F)
- **🏗️ Structural Stability** assessment (RMSD-based)
- **⚖️ Thermodynamic Equilibrium** evaluation
- **🔗 Binding Affinity** analysis (hydrogen bonds)
- **📦 Compactness** and **🌊 Surface Properties** metrics
- **💡 Automated Insights** with actionable recommendations
- **📈 Interactive Radar Chart** for quality visualization
- All features from the basic versions below

### Basic Versions (`xvg_visualization_report.html`)
- Interactive plots for all XVG files
- Table of contents for easy navigation
- Metadata about each analysis (command used, data points, etc.)
- Responsive design that works on all devices
- Hover tooltips for data exploration
- Zoom and pan capabilities

## Features

- Automatically detects and processes all `.xvg` files in the current directory
- Extracts plot titles, axis labels, and legends from XVG headers
- Handles multi-column data files
- Creates professional-looking visualizations
- No manual configuration needed

## Browser Compatibility

The generated HTML report works in all modern web browsers:
- Chrome/Edge
- Firefox
- Safari
- Opera

## Example XVG Files Supported

- RMSD (Root Mean Square Deviation)
- RMSF (Root Mean Square Fluctuation)
- Radius of Gyration
- SASA (Solvent Accessible Surface Area)
- Temperature
- Pressure
- Potential Energy
- Density
- Hydrogen Bonds
- And any other GROMACS analysis output in XVG format