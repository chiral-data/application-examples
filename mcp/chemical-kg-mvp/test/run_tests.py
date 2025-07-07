#!/usr/bin/env python3
"""
Test runner for Chemical Knowledge Graph MVP
"""

import sys
import os
import subprocess
import argparse

# Add src directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

def run_pytest(test_path=None, verbose=False, capture=False, markers=None, timeout=None):
    """Run pytest with specified options"""
    
    cmd = ['python', '-m', 'pytest']
    
    if test_path:
        cmd.append(test_path)
    else:
        cmd.append(os.path.dirname(__file__))
    
    if verbose:
        cmd.extend(['-v', '-s'])
    
    if not capture:
        cmd.append('--tb=short')
    
    if markers:
        cmd.extend(['-m', markers])
    
    if timeout:
        cmd.extend(['--timeout', str(timeout)])
    
    # Add coverage if available
    try:
        import pytest_cov
        cmd.extend(['--cov=src', '--cov-report=term-missing'])
    except ImportError:
        pass
    
    print(f"Running: {' '.join(cmd)}")
    print("=" * 60)
    
    try:
        result = subprocess.run(cmd, cwd=os.path.dirname(os.path.dirname(__file__)))
        return result.returncode
    except FileNotFoundError:
        print("pytest not found. Installing required packages...")
        subprocess.run([sys.executable, '-m', 'pip', 'install', 'pytest'])
        result = subprocess.run(cmd, cwd=os.path.dirname(os.path.dirname(__file__)))
        return result.returncode

def run_unit_tests():
    """Run only unit tests"""
    return run_pytest('test/unit', verbose=True)

def run_integration_tests():
    """Run only integration tests"""
    return run_pytest('test/integration', verbose=True)

def run_specific_test(test_file):
    """Run a specific test file"""
    return run_pytest(test_file, verbose=True)

def check_dependencies():
    """Check if required dependencies are available"""
    print("Checking dependencies...")
    
    required_packages = [
        'pytest',
        'unittest.mock',
        'PIL',
        'sys',
        'os',
        'tempfile',
        'numpy',
        'cv2'
    ]
    
    optional_packages = [
        'DECIMER',
        'decimer_segmentation',
        'tensorflow'
    ]
    
    missing = []
    for package in required_packages:
        try:
            if package == 'PIL':
                import PIL
            elif package == 'unittest.mock':
                from unittest.mock import Mock
            elif package == 'cv2':
                import cv2
            else:
                __import__(package)
            print(f"✓ {package}")
        except ImportError:
            print(f"✗ {package}")
            missing.append(package)
    
    print("\nOptional DECIMER packages:")
    for package in optional_packages:
        try:
            if package == 'DECIMER':
                from DECIMER import predict_SMILES
            elif package == 'decimer_segmentation':
                from decimer_segmentation import segment_chemical_structures_from_file
            else:
                __import__(package)
            print(f"✓ {package}")
        except ImportError:
            print(f"? {package} (optional - tests will use mocks)")
    
    if missing:
        print(f"\nMissing required packages: {missing}")
        print("Install with: pip install pytest Pillow opencv-python numpy")
        return False
    
    print("All required dependencies available!")
    return True

def main():
    parser = argparse.ArgumentParser(description='Run tests for Chemical Knowledge Graph MVP')
    parser.add_argument('--unit', action='store_true', help='Run only unit tests')
    parser.add_argument('--integration', action='store_true', help='Run only integration tests')
    parser.add_argument('--file', type=str, help='Run specific test file')
    parser.add_argument('--check-deps', action='store_true', help='Check dependencies')
    parser.add_argument('--all', action='store_true', help='Run all tests (default)')
    
    # Hanging issue specific options
    parser.add_argument('--timeout', action='store_true', help='Run only timeout-related tests')
    parser.add_argument('--memory', action='store_true', help='Run only memory monitoring tests')
    parser.add_argument('--gpu', action='store_true', help='Run only GPU/CPU fallback tests')
    parser.add_argument('--perf', action='store_true', help='Run only performance benchmarks')
    parser.add_argument('--quick', action='store_true', help='Run only quick tests (< 30s each)')
    parser.add_argument('--stress', action='store_true', help='Run stress tests')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose output')
    
    args = parser.parse_args()
    
    if args.check_deps:
        return 0 if check_dependencies() else 1
    
    if not check_dependencies():
        return 1
    
    # Handle specific test categories
    if args.timeout:
        print("Running timeout-related tests...")
        return run_pytest('test', verbose=args.verbose, markers='timeout')
    elif args.memory:
        print("Running memory monitoring tests...")
        return run_pytest('test/unit/test_resource_monitoring.py', verbose=args.verbose)
    elif args.gpu:
        print("Running GPU/CPU fallback tests...")
        return run_pytest('test/unit/test_gpu_cpu_fallback.py', verbose=args.verbose)
    elif args.perf:
        print("Running performance benchmarks...")
        return run_pytest('test/integration/test_performance_benchmarks.py', verbose=args.verbose)
    elif args.quick:
        print("Running quick tests...")
        return run_pytest('test', verbose=args.verbose, markers='not slow', timeout=30)
    elif args.stress:
        print("Running stress tests...")
        return run_pytest('test', verbose=args.verbose, markers='stress')
    elif args.unit:
        return run_unit_tests()
    elif args.integration:
        return run_integration_tests()
    elif args.file:
        return run_specific_test(args.file)
    else:
        # Run all tests
        print("Running all tests...")
        unit_result = run_unit_tests()
        integration_result = run_integration_tests()
        
        if unit_result == 0 and integration_result == 0:
            print("\n" + "=" * 60)
            print("ALL TESTS PASSED!")
            return 0
        else:
            print("\n" + "=" * 60)
            print("SOME TESTS FAILED!")
            return 1

if __name__ == '__main__':
    sys.exit(main())