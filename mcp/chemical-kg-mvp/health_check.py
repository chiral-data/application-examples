#!/usr/bin/env python3
"""
Health check script for Chemical Knowledge Graph MVP
Run this before demos to ensure all components are working.
"""

import os
import sys
import json
import traceback
from pathlib import Path

def check_component(name, check_func):
    """Check a component and report status"""
    try:
        result = check_func()
        print(f"✓ {name}: {result}")
        return True, result
    except Exception as e:
        print(f"✗ {name}: FAILED - {e}")
        return False, str(e)

def check_docker_container():
    """Check if Docker container is running"""
    import subprocess
    result = subprocess.run(['docker', 'ps', '--filter', 'name=chem-kg-mvp', '--format', '{{.Status}}'], 
                          capture_output=True, text=True)
    if result.returncode == 0 and 'Up' in result.stdout:
        return "Container is running"
    else:
        return "Container not running or not found"

def check_demo_files():
    """Check if demo files exist"""
    required_files = [
        "demo_data.json",
        "demo_chemical_paper.pdf",
        "demo_structures/"
    ]
    
    missing = [f for f in required_files if not os.path.exists(f)]
    
    if missing:
        return f"Missing: {missing}"
    else:
        return "All demo files present"

def check_streamlit_access():
    """Check if Streamlit is accessible"""
    import subprocess
    result = subprocess.run(['curl', '-f', 'http://localhost:8501/_stcore/health'], 
                          capture_output=True, text=True)
    if result.returncode == 0:
        return "Streamlit is accessible on port 8501"
    else:
        return "Streamlit not accessible (may not be started yet)"

def check_decimer_in_container():
    """Check DECIMER status in container"""
    import subprocess
    result = subprocess.run(['docker', 'exec', 'chem-kg-mvp', 'python3', '-c', 
                           'import DECIMER; print("DECIMER available")'], 
                          capture_output=True, text=True)
    if result.returncode == 0:
        return "DECIMER is available in container"
    else:
        return "DECIMER not available in container"

def run_health_check():
    """Run complete health check"""
    print("=" * 60)
    print("Chemical Knowledge Graph MVP - Health Check")
    print("=" * 60)
    
    checks = [
        ("Docker Container Status", check_docker_container),
        ("Demo Files", check_demo_files),
        ("DECIMER in Container", check_decimer_in_container),
        ("Streamlit Access", check_streamlit_access),
    ]
    
    results = {}
    total_checks = len(checks)
    passed_checks = 0
    
    for name, check_func in checks:
        success, result = check_component(name, check_func)
        results[name] = (success, result)
        if success:
            passed_checks += 1
    
    print("\n" + "=" * 60)
    print(f"Health Check Summary: {passed_checks}/{total_checks} checks passed")
    print("=" * 60)
    
    # Provide recommendations
    print("\nStatus Overview:")
    
    if results["Docker Container Status"][0]:
        print("+ Docker container is running")
    else:
        print("- Docker container needs to be started:")
        print("  docker run -d --name chem-kg-mvp --gpus all -p 8501:8501 chem-kg-mvp")
    
    if results["Demo Files"][0]:
        print("+ Demo files are ready")
    else:
        print("- Demo files missing - run demo_data.py in container")
    
    if results["DECIMER in Container"][0]:
        print("+ DECIMER is available for structure processing")
    else:
        print("! DECIMER not available - demo mode will be needed")
    
    if results["Streamlit Access"][0]:
        print("+ Application is accessible at http://localhost:8501")
    else:
        print("! Application not yet accessible (container may still be starting)")
    
    # Demo readiness score
    demo_ready = passed_checks >= 2  # At least container + demo files
    
    print(f"\nDemo Readiness: {'READY' if demo_ready else 'NOT READY'}")
    
    if demo_ready:
        print("\nDemo Instructions:")
        print("1. Open browser to http://localhost:8501")
        print("2. If DECIMER fails, click 'Activate Demo Mode'")
        print("3. Use demo_chemical_paper.pdf for testing uploads")
        print("4. Show structure extraction, CSV download, and Q&A features")
    else:
        print("\nFix the issues above before running the demo")
    
    return demo_ready

if __name__ == "__main__":
    success = run_health_check()
    sys.exit(0 if success else 1)