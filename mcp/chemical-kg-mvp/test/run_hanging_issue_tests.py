#!/usr/bin/env python3
"""
Comprehensive test runner specifically for PDF processing hanging issues.

This script runs all tests designed to identify and tackle the DECIMER
segmentation hanging problem.

Usage:
    python test/run_hanging_issue_tests.py [options]

Options:
    --quick     Run only quick tests (< 30s each)
    --full      Run all tests including stress tests
    --timeout   Run only timeout-related tests
    --memory    Run only memory monitoring tests
    --gpu       Run only GPU/CPU fallback tests
    --perf      Run only performance benchmarks
    --verbose   Enable verbose output
    --parallel  Run tests in parallel where possible

Example:
    python test/run_hanging_issue_tests.py --timeout --verbose
"""

import os
import sys
import argparse
import subprocess
import time
import signal
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))


class TestRunner:
    """Comprehensive test runner for hanging issue investigation"""
    
    def __init__(self, verbose=False, timeout_per_test=300):
        self.verbose = verbose
        self.timeout_per_test = timeout_per_test
        self.results = {}
        
    def run_test_module(self, module_path, test_filter=None):
        """Run a specific test module with timeout protection"""
        print(f"\n{'='*60}")
        print(f"Running: {module_path}")
        print(f"{'='*60}")
        
        cmd = ["python", "-m", "pytest", module_path, "-v"]
        
        if test_filter:
            cmd.extend(["-k", test_filter])
            
        if self.verbose:
            cmd.append("-s")
            
        # Add timeout and other options
        cmd.extend([
            "--tb=short",
            "--disable-warnings",
            f"--timeout={self.timeout_per_test}"
        ])
        
        start_time = time.time()
        
        try:
            # Run with subprocess to enable timeout
            result = subprocess.run(
                cmd,
                cwd=Path(__file__).parent.parent,
                capture_output=not self.verbose,
                text=True,
                timeout=self.timeout_per_test + 30  # Extra buffer for pytest overhead
            )
            
            duration = time.time() - start_time
            
            self.results[module_path] = {
                'status': 'PASSED' if result.returncode == 0 else 'FAILED',
                'duration': duration,
                'returncode': result.returncode,
                'stdout': result.stdout if not self.verbose else '',
                'stderr': result.stderr if not self.verbose else ''
            }
            
            if result.returncode == 0:
                print(f"✅ PASSED in {duration:.1f}s")
            else:
                print(f"❌ FAILED in {duration:.1f}s (exit code: {result.returncode})")
                if not self.verbose and result.stderr:
                    print(f"Error output: {result.stderr[-500:]}")  # Last 500 chars
                    
        except subprocess.TimeoutExpired:
            duration = time.time() - start_time
            print(f"⏰ TIMEOUT after {duration:.1f}s")
            self.results[module_path] = {
                'status': 'TIMEOUT',
                'duration': duration,
                'returncode': -1,
                'stdout': '',
                'stderr': 'Test timed out'
            }
            
        except Exception as e:
            duration = time.time() - start_time
            print(f"💥 ERROR: {e}")
            self.results[module_path] = {
                'status': 'ERROR',
                'duration': duration,
                'returncode': -2,
                'stdout': '',
                'stderr': str(e)
            }
    
    def print_summary(self):
        """Print test run summary"""
        print(f"\n{'='*80}")
        print("TEST RUN SUMMARY")
        print(f"{'='*80}")
        
        total_tests = len(self.results)
        passed = sum(1 for r in self.results.values() if r['status'] == 'PASSED')
        failed = sum(1 for r in self.results.values() if r['status'] == 'FAILED')
        timeout = sum(1 for r in self.results.values() if r['status'] == 'TIMEOUT')
        error = sum(1 for r in self.results.values() if r['status'] == 'ERROR')
        
        total_time = sum(r['duration'] for r in self.results.values())
        
        print(f"Total Tests: {total_tests}")
        print(f"Passed: {passed} ✅")
        print(f"Failed: {failed} ❌")
        print(f"Timeout: {timeout} ⏰")
        print(f"Error: {error} 💥")
        print(f"Total Time: {total_time:.1f}s")
        
        print(f"\nDetailed Results:")
        print(f"{'-'*80}")
        
        for module, result in self.results.items():
            status_emoji = {
                'PASSED': '✅',
                'FAILED': '❌', 
                'TIMEOUT': '⏰',
                'ERROR': '💥'
            }
            
            print(f"{status_emoji[result['status']]} {module:<50} {result['duration']:>6.1f}s")
            
        # Identify hanging issues
        hanging_tests = [
            module for module, result in self.results.items()
            if result['status'] == 'TIMEOUT' or result['duration'] > 60
        ]
        
        if hanging_tests:
            print(f"\n⚠️  POTENTIAL HANGING ISSUES DETECTED:")
            for test in hanging_tests:
                print(f"   - {test}")
                
        print(f"\n{'='*80}")
        
        return passed == total_tests


def main():
    parser = argparse.ArgumentParser(description="Run hanging issue investigation tests")
    parser.add_argument("--quick", action="store_true", help="Run only quick tests")
    parser.add_argument("--full", action="store_true", help="Run all tests including stress tests")
    parser.add_argument("--timeout", action="store_true", help="Run only timeout-related tests")
    parser.add_argument("--memory", action="store_true", help="Run only memory monitoring tests")
    parser.add_argument("--gpu", action="store_true", help="Run only GPU/CPU fallback tests")
    parser.add_argument("--perf", action="store_true", help="Run only performance benchmarks")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output")
    parser.add_argument("--parallel", action="store_true", help="Run tests in parallel")
    parser.add_argument("--timeout-per-test", type=int, default=300, help="Timeout per test in seconds")
    
    args = parser.parse_args()
    
    # Determine which tests to run
    test_modules = []
    
    if args.timeout or args.full or not any([args.quick, args.memory, args.gpu, args.perf]):
        test_modules.extend([
            "test/integration/test_pdf_processing_timeouts.py",
        ])
    
    if args.memory or args.full or not any([args.quick, args.timeout, args.gpu, args.perf]):
        test_modules.extend([
            "test/unit/test_resource_monitoring.py",
        ])
        
    if args.gpu or args.full or not any([args.quick, args.timeout, args.memory, args.perf]):
        test_modules.extend([
            "test/unit/test_gpu_cpu_fallback.py",
        ])
        
    if args.perf or args.full or not any([args.quick, args.timeout, args.memory, args.gpu]):
        test_modules.extend([
            "test/integration/test_performance_benchmarks.py",
        ])
        
    if not args.quick:
        test_modules.extend([
            "test/integration/test_decimer_robustness.py",
        ])
    
    # Remove duplicates while preserving order
    test_modules = list(dict.fromkeys(test_modules))
    
    print("🔬 PDF Processing Hanging Issue Investigation")
    print(f"Running {len(test_modules)} test modules with {args.timeout_per_test}s timeout per test")
    
    if args.quick:
        print("⚡ Quick mode: Running essential tests only")
    elif args.full:
        print("🔍 Full mode: Running all tests including stress tests")
    
    runner = TestRunner(verbose=args.verbose, timeout_per_test=args.timeout_per_test)
    
    try:
        for module in test_modules:
            runner.run_test_module(module)
            
        success = runner.print_summary()
        
        if success:
            print("🎉 All tests passed! No hanging issues detected.")
            sys.exit(0)
        else:
            print("⚠️  Some tests failed or timed out. Check results above.")
            sys.exit(1)
            
    except KeyboardInterrupt:
        print("\n🛑 Test run interrupted by user")
        runner.print_summary()
        sys.exit(130)


if __name__ == "__main__":
    main()