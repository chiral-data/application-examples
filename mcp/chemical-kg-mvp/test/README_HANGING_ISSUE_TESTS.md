# PDF Processing Hanging Issue Tests

This test suite is specifically designed to identify and tackle the PDF processing hanging issue observed with DECIMER segmentation.

## Problem Description

The application occasionally hangs during PDF processing, specifically during the DECIMER chemical structure segmentation step (`segment_chemical_structures_from_file()` in `pdf_processor.py:193`). This happens when:

1. Processing certain PDF files
2. DECIMER segmentation gets stuck without timeout
3. CPU usage drops to 0% but process doesn't complete
4. No error messages or logs are produced

## Test Suite Overview

### 1. Core Timeout Tests (`test_pdf_processing_timeouts.py`)

**Purpose**: Detect hanging conditions and enforce timeouts

**Key Tests**:
- `test_decimer_segmentation_timeout()` - Tests DECIMER with hard timeout
- `test_decimer_segmentation_with_mocked_hanging()` - Simulates hanging behavior
- `test_process_monitoring_and_termination()` - Tests process termination
- `test_memory_leak_detection()` - Detects memory leaks during repeated processing

**Expected Behavior**: All tests should complete within specified timeouts

### 2. DECIMER Robustness Tests (`test_decimer_robustness.py`)

**Purpose**: Test DECIMER with various problematic inputs

**Key Tests**:
- `test_large_structure_images()` - Tests with very large images
- `test_noisy_structure_images()` - Tests with noisy/corrupted images
- `test_concurrent_decimer_processing()` - Tests race conditions
- `test_decimer_with_specific_tensorflow_settings()` - Tests TF configurations

**Expected Behavior**: DECIMER should handle all inputs gracefully without hanging

### 3. Resource Monitoring Tests (`test_resource_monitoring.py`)

**Purpose**: Monitor memory usage and detect leaks

**Key Tests**:
- `test_repeated_processing_memory_leak()` - Detects memory leaks
- `test_large_pdf_memory_scaling()` - Tests memory scaling
- `test_memory_pressure_handling()` - Tests under memory pressure
- `test_concurrent_processing_memory_isolation()` - Tests memory isolation

**Expected Behavior**: Memory usage should be bounded and not leak

### 4. GPU/CPU Fallback Tests (`test_gpu_cpu_fallback.py`)

**Purpose**: Test hardware configuration handling

**Key Tests**:
- `test_gpu_detection_success()` - Tests GPU detection
- `test_no_gpu_available()` - Tests CPU-only fallback
- `test_gpu_memory_growth_error()` - Tests GPU error handling
- `test_tensorflow_import_error()` - Tests missing dependencies

**Expected Behavior**: Should work on both GPU and CPU without hanging

### 5. Performance Benchmarks (`test_performance_benchmarks.py`)

**Purpose**: Establish performance baselines and detect regressions

**Key Tests**:
- `test_simple_pdf_processing_benchmark()` - Baseline performance
- `test_decimer_segmentation_benchmark()` - DECIMER performance
- `test_full_pipeline_benchmark()` - End-to-end performance
- `test_stress_test_performance()` - Stress testing

**Expected Behavior**: All operations should complete within performance thresholds

## Running the Tests

### Quick Test (Essential tests only)
```bash
cd /home/ubuntu/chiral/application-examples/mcp/chemical-kg-mvp
python test/run_hanging_issue_tests.py --quick --verbose
```

### Full Test Suite
```bash
python test/run_hanging_issue_tests.py --full --verbose
```

### Specific Test Categories
```bash
# Timeout-related tests only
python test/run_hanging_issue_tests.py --timeout

# Memory monitoring tests only
python test/run_hanging_issue_tests.py --memory

# GPU/CPU fallback tests only
python test/run_hanging_issue_tests.py --gpu

# Performance benchmarks only
python test/run_hanging_issue_tests.py --perf
```

### Individual Test Files
```bash
# Run specific test file
python -m pytest test/integration/test_pdf_processing_timeouts.py -v

# Run specific test with timeout
python -m pytest test/integration/test_pdf_processing_timeouts.py::TestPDFProcessingTimeouts::test_decimer_segmentation_timeout -v --timeout=60
```

## Test Configuration

### Timeouts
- Default test timeout: 300 seconds (5 minutes)
- Individual test timeouts can be configured
- Hard timeouts using multiprocessing for hanging detection

### Resource Limits
- Memory usage monitoring with configurable thresholds
- CPU usage monitoring
- File handle leak detection

### Environment Variables
Tests respect these environment variables:
- `CUDA_VISIBLE_DEVICES` - Control GPU visibility
- `TF_FORCE_GPU_ALLOW_GROWTH` - TensorFlow GPU memory growth
- `TF_CPP_MIN_LOG_LEVEL` - TensorFlow logging level

## Expected Results

### Passing Tests
All tests should pass under normal conditions. Key expectations:

1. **No Timeouts**: No test should hang indefinitely
2. **Memory Bounds**: Memory usage should stay within reasonable limits
3. **Resource Cleanup**: All resources should be properly cleaned up
4. **Error Handling**: Errors should be handled gracefully with fallbacks

### Failing Tests Indicate Issues

**Timeout Failures**: 
- Indicates hanging in DECIMER segmentation
- Suggests need for timeout implementation in production code

**Memory Failures**:
- Memory leaks in repeated processing
- Excessive memory usage with large files

**Performance Failures**:
- Performance regressions
- Scaling issues with file size

## Debugging Hanging Issues

### 1. Enable Verbose Logging
```bash
python test/run_hanging_issue_tests.py --timeout --verbose
```

### 2. Run with Process Monitoring
```bash
# In separate terminal, monitor processes
watch -n 1 'ps aux | grep python'
htop
```

### 3. Memory Monitoring
```bash
# Monitor memory usage during tests
python -c "
import psutil
import time
while True:
    mem = psutil.virtual_memory()
    print(f'Memory: {mem.percent}% used')
    time.sleep(1)
"
```

### 4. GPU Monitoring (if applicable)
```bash
watch -n 1 nvidia-smi
```

## Root Cause Analysis

Based on test results, the hanging issue is likely caused by:

1. **DECIMER Segmentation**: The `segment_chemical_structures_from_file()` function lacks timeout
2. **TensorFlow/GPU Issues**: GPU memory allocation or TensorFlow configuration problems
3. **PDF Content**: Certain PDF structures or image content causes DECIMER to hang
4. **Memory Pressure**: Insufficient memory causing silent failures

## Recommended Solutions

### 1. Add Timeout to DECIMER Calls
```python
import multiprocessing as mp
import time

def segment_with_timeout(pdf_path, timeout=60):
    def worker():
        return segment_chemical_structures_from_file(pdf_path)
    
    with mp.Pool(1) as pool:
        result = pool.apply_async(worker)
        try:
            return result.get(timeout=timeout)
        except mp.TimeoutError:
            return []  # Fallback to empty results
```

### 2. Resource Monitoring
```python
import psutil
import threading

class ResourceMonitor:
    def __init__(self, memory_limit_mb=1000):
        self.memory_limit = memory_limit_mb
        self.monitoring = False
        
    def check_resources(self):
        memory_mb = psutil.Process().memory_info().rss / 1024 / 1024
        if memory_mb > self.memory_limit:
            raise MemoryError(f"Memory limit exceeded: {memory_mb:.1f}MB")
```

### 3. Graceful Fallbacks
```python
def extract_chemical_structures_with_fallback(self, pdf_path):
    try:
        # Try DECIMER with timeout
        return self.extract_with_decimer_timeout(pdf_path, timeout=60)
    except TimeoutError:
        print("DECIMER timed out, falling back to regular extraction")
        return self.extract_text_and_images()
    except Exception as e:
        print(f"DECIMER failed: {e}, falling back")
        return self.extract_text_and_images()
```

## Maintenance

### Adding New Tests
1. Create test file in appropriate directory (`unit/` or `integration/`)
2. Follow naming convention: `test_*_hanging_issue.py`
3. Add timeout decorators and resource monitoring
4. Update test runner configuration

### Updating Thresholds
Performance and resource thresholds may need adjustment based on:
- Hardware capabilities
- Expected file sizes
- Performance requirements

Edit thresholds in individual test files or centralize in a config file.

## CI/CD Integration

For continuous integration, add these tests to your pipeline:

```yaml
# .github/workflows/hanging-issue-tests.yml
name: Hanging Issue Tests
on: [push, pull_request]
jobs:
  test-hanging-issues:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      - name: Setup Python
        uses: actions/setup-python@v2
        with:
          python-version: '3.8'
      - name: Install dependencies
        run: pip install -r requirements.txt
      - name: Run hanging issue tests
        run: python test/run_hanging_issue_tests.py --quick
        timeout-minutes: 30
```

This comprehensive test suite should help identify and resolve the PDF processing hanging issue.