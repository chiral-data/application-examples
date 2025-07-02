"""
Performance benchmarking tests for PDF processing and DECIMER operations.
These tests help identify performance regressions and hanging conditions.
"""

import pytest
import os
import time
import tempfile
import shutil
import threading
import multiprocessing as mp
import statistics
from unittest.mock import patch, MagicMock
import numpy as np
from PIL import Image
import fitz  # PyMuPDF

from pdf_processor import PDFProcessor
from decimer_client import DECIMERClient


class PerformanceTimer:
    """Context manager for timing operations"""
    
    def __init__(self, name="Operation"):
        self.name = name
        self.start_time = None
        self.end_time = None
        self.duration = None
        
    def __enter__(self):
        self.start_time = time.time()
        return self
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.end_time = time.time()
        self.duration = self.end_time - self.start_time
        print(f"{self.name} took {self.duration:.3f} seconds")


class BenchmarkDataGenerator:
    """Generate test data for benchmarking"""
    
    @staticmethod
    def create_benchmark_pdf(complexity="medium", pages=5, structures_per_page=3):
        """Create PDFs with varying complexity for benchmarking"""
        temp_dir = tempfile.mkdtemp()
        pdf_path = os.path.join(temp_dir, f"benchmark_{complexity}_{pages}p.pdf")
        
        doc = fitz.open()
        
        for page_num in range(pages):
            page = doc.new_page()
            
            if complexity == "simple":
                # Simple text-only content
                for line in range(10):
                    page.insert_text((50, 50 + line * 20), 
                                    f"Simple text line {line} on page {page_num + 1}")
                    
            elif complexity == "medium":
                # Mixed text and simple structures
                page.insert_text((50, 50), f"Chemical Synthesis Report - Page {page_num + 1}")
                page.insert_text((50, 100), "Abstract: This paper describes the synthesis of organic compounds.")
                
                for i in range(structures_per_page):
                    page.insert_text((50, 150 + i * 30), f"Compound {i+1}: C{6+i}H{12+i}O{1+i}")
                    
            elif complexity == "complex":
                # Dense content with many potential chemical structures
                page.insert_text((50, 50), f"Complex Chemical Analysis - Page {page_num + 1}")
                
                # Add lots of text that might contain chemical formulas
                chemical_text = [
                    "The reaction proceeded via intermediate C6H5CH2CH2NH2",
                    "Benzene (C6H6) was used as solvent",
                    "Product analysis showed C8H10N4O2 formation",
                    "Spectroscopic data: 1H NMR (400 MHz, CDCl3)",
                    "13C NMR (100 MHz, CDCl3): δ 150.2, 134.5, 128.7",
                    "IR (KBr): 3300, 2950, 1650, 1550 cm-1",
                    "MS (EI): m/z 194 [M]+, 165, 137, 109, 77",
                    "Elemental analysis for C12H10N2O:",
                    "Calc: C 72.72, H 5.09, N 14.13",
                    "Found: C 72.68, H 5.12, N 14.09"
                ]
                
                for i, text in enumerate(chemical_text):
                    page.insert_text((50, 100 + i * 25), text)
                    
                # Add mock images/structures
                for i in range(structures_per_page * 2):  # More structures in complex PDFs
                    page.insert_text((300, 100 + i * 20), f"[Structure {i+1}]")
        
        doc.save(pdf_path)
        doc.close()
        
        return pdf_path, temp_dir
    
    @staticmethod
    def create_large_benchmark_pdf(size_mb=10):
        """Create a large PDF for memory/performance testing"""
        temp_dir = tempfile.mkdtemp()
        pdf_path = os.path.join(temp_dir, f"large_benchmark_{size_mb}mb.pdf")
        
        doc = fitz.open()
        
        # Calculate pages needed for target size
        target_size = size_mb * 1024 * 1024
        current_size = 0
        page_num = 0
        
        while current_size < target_size:
            page = doc.new_page()
            
            # Add dense text content
            for line in range(50):
                long_text = "Chemical formula analysis: " + "C6H6 " * 50
                page.insert_text((50, 20 + line * 15), long_text)
            
            # Check current size
            temp_path = pdf_path + ".tmp"
            doc.save(temp_path)
            current_size = os.path.getsize(temp_path)
            os.remove(temp_path)
            
            page_num += 1
            if page_num > 100:  # Safety limit
                break
        
        doc.save(pdf_path)
        doc.close()
        
        return pdf_path, temp_dir


@pytest.fixture
def benchmark_data():
    return BenchmarkDataGenerator()


class TestPDFProcessingPerformance:
    """Benchmark PDF processing operations"""
    
    def test_simple_pdf_processing_benchmark(self, benchmark_data):
        """Benchmark simple PDF processing"""
        pdf_path, temp_dir = benchmark_data.create_benchmark_pdf("simple", pages=3)
        
        try:
            # Warm-up run
            processor = PDFProcessor(pdf_path)
            processor.extract_text_and_images()
            processor.cleanup_temp_files()
            
            # Benchmark runs
            times = []
            for run in range(5):
                with PerformanceTimer(f"Simple PDF processing run {run+1}") as timer:
                    processor = PDFProcessor(pdf_path)
                    result = processor.extract_text_and_images()
                    processor.cleanup_temp_files()
                
                times.append(timer.duration)
                assert len(result) == 3, "Should process all pages"
            
            # Performance assertions
            avg_time = statistics.mean(times)
            max_time = max(times)
            
            assert avg_time < 5.0, f"Simple PDF processing too slow: {avg_time:.3f}s average"
            assert max_time < 10.0, f"Simple PDF processing too slow: {max_time:.3f}s max"
            
            print(f"Simple PDF benchmark: avg={avg_time:.3f}s, max={max_time:.3f}s, std={statistics.stdev(times):.3f}s")
            
        finally:
            shutil.rmtree(temp_dir)
    
    def test_medium_complexity_pdf_benchmark(self, benchmark_data):
        """Benchmark medium complexity PDF processing"""
        pdf_path, temp_dir = benchmark_data.create_benchmark_pdf("medium", pages=5, structures_per_page=3)
        
        try:
            times = []
            for run in range(3):
                with PerformanceTimer(f"Medium PDF processing run {run+1}") as timer:
                    processor = PDFProcessor(pdf_path)
                    result = processor.extract_text_and_images()
                    processor.cleanup_temp_files()
                
                times.append(timer.duration)
                assert len(result) == 5, "Should process all pages"
            
            avg_time = statistics.mean(times)
            assert avg_time < 15.0, f"Medium PDF processing too slow: {avg_time:.3f}s average"
            
            print(f"Medium PDF benchmark: avg={avg_time:.3f}s")
            
        finally:
            shutil.rmtree(temp_dir)
    
    def test_complex_pdf_processing_benchmark(self, benchmark_data):
        """Benchmark complex PDF processing"""
        pdf_path, temp_dir = benchmark_data.create_benchmark_pdf("complex", pages=3, structures_per_page=5)
        
        try:
            times = []
            for run in range(3):
                with PerformanceTimer(f"Complex PDF processing run {run+1}") as timer:
                    processor = PDFProcessor(pdf_path)
                    result = processor.extract_text_and_images()
                    processor.cleanup_temp_files()
                
                times.append(timer.duration)
                assert len(result) == 3, "Should process all pages"
            
            avg_time = statistics.mean(times)
            assert avg_time < 30.0, f"Complex PDF processing too slow: {avg_time:.3f}s average"
            
            print(f"Complex PDF benchmark: avg={avg_time:.3f}s")
            
        finally:
            shutil.rmtree(temp_dir)
    
    def test_large_pdf_memory_performance(self, benchmark_data):
        """Benchmark large PDF processing with memory monitoring"""
        pdf_path, temp_dir = benchmark_data.create_large_benchmark_pdf(size_mb=5)
        
        try:
            import psutil
            initial_memory = psutil.Process().memory_info().rss / 1024 / 1024
            
            with PerformanceTimer("Large PDF processing") as timer:
                processor = PDFProcessor(pdf_path)
                result = processor.extract_text_and_images()
                processor.cleanup_temp_files()
            
            final_memory = psutil.Process().memory_info().rss / 1024 / 1024
            memory_increase = final_memory - initial_memory
            
            # Performance assertions
            assert timer.duration < 60.0, f"Large PDF processing too slow: {timer.duration:.3f}s"
            assert memory_increase < 500, f"Excessive memory usage: {memory_increase:.2f}MB"
            assert len(result) > 0, "Should process pages"
            
            print(f"Large PDF: {timer.duration:.3f}s, memory increase: {memory_increase:.2f}MB")
            
        finally:
            shutil.rmtree(temp_dir)
    
    def test_concurrent_pdf_processing_benchmark(self, benchmark_data):
        """Benchmark concurrent PDF processing"""
        # Create multiple test PDFs
        pdf_data = []
        for i in range(4):
            pdf_path, temp_dir = benchmark_data.create_benchmark_pdf("medium", pages=2)
            pdf_data.append((pdf_path, temp_dir))
        
        try:
            def process_pdf_worker(pdf_path):
                processor = PDFProcessor(pdf_path)
                result = processor.extract_text_and_images()
                processor.cleanup_temp_files()
                return len(result)
            
            # Sequential processing
            with PerformanceTimer("Sequential processing") as seq_timer:
                seq_results = []
                for pdf_path, _ in pdf_data:
                    result = process_pdf_worker(pdf_path)
                    seq_results.append(result)
            
            # Concurrent processing
            with PerformanceTimer("Concurrent processing") as conc_timer:
                with mp.Pool(4) as pool:
                    async_results = [pool.apply_async(process_pdf_worker, (pdf_path,)) 
                                   for pdf_path, _ in pdf_data]
                    conc_results = [r.get(timeout=30) for r in async_results]
            
            # Verify results
            assert seq_results == conc_results, "Concurrent results should match sequential"
            
            # Performance comparison
            speedup = seq_timer.duration / conc_timer.duration
            print(f"Speedup: {speedup:.2f}x (sequential: {seq_timer.duration:.3f}s, concurrent: {conc_timer.duration:.3f}s)")
            
            # Should have some speedup but not require it for test to pass
            assert conc_timer.duration < seq_timer.duration * 1.5, "Concurrent processing shouldn't be much slower"
            
        finally:
            for _, temp_dir in pdf_data:
                shutil.rmtree(temp_dir)


class TestDECIMERPerformance:
    """Benchmark DECIMER operations"""
    
    def test_decimer_initialization_benchmark(self):
        """Benchmark DECIMER client initialization"""
        times = []
        
        for run in range(3):
            with PerformanceTimer(f"DECIMER initialization run {run+1}") as timer:
                client = DECIMERClient()
                # Access properties to ensure initialization
                _ = client.decimer_available
                _ = client.segmentation_available
            
            times.append(timer.duration)
        
        avg_time = statistics.mean(times)
        assert avg_time < 10.0, f"DECIMER initialization too slow: {avg_time:.3f}s average"
        
        print(f"DECIMER init benchmark: avg={avg_time:.3f}s")
    
    def test_decimer_segmentation_benchmark(self, benchmark_data):
        """Benchmark DECIMER segmentation performance"""
        pdf_path, temp_dir = benchmark_data.create_benchmark_pdf("medium", pages=2, structures_per_page=2)
        
        try:
            client = DECIMERClient()
            
            # Warm-up run
            client.segment_structures_from_pdf(pdf_path)
            
            # Benchmark runs
            times = []
            for run in range(3):
                with PerformanceTimer(f"DECIMER segmentation run {run+1}") as timer:
                    segments = client.segment_structures_from_pdf(pdf_path)
                
                times.append(timer.duration)
                assert isinstance(segments, list), "Should return list of segments"
            
            avg_time = statistics.mean(times)
            assert avg_time < 45.0, f"DECIMER segmentation too slow: {avg_time:.3f}s average"
            
            print(f"DECIMER segmentation benchmark: avg={avg_time:.3f}s")
            
        finally:
            shutil.rmtree(temp_dir)
    
    def test_decimer_batch_processing_benchmark(self):
        """Benchmark DECIMER batch processing"""
        temp_dir = tempfile.mkdtemp()
        
        try:
            # Create test images
            image_paths = []
            for i in range(5):
                img_array = np.random.randint(0, 255, (200, 200, 3), dtype=np.uint8)
                img = Image.fromarray(img_array)
                img_path = os.path.join(temp_dir, f"test_img_{i}.png")
                img.save(img_path)
                image_paths.append(img_path)
            
            client = DECIMERClient()
            
            # Benchmark batch processing
            with PerformanceTimer("DECIMER batch processing") as timer:
                results = client.batch_process(image_paths)
            
            assert len(results) == 5, "Should process all images"
            assert timer.duration < 60.0, f"Batch processing too slow: {timer.duration:.3f}s"
            
            print(f"Batch processing: {timer.duration:.3f}s for {len(image_paths)} images")
            
        finally:
            shutil.rmtree(temp_dir)
    
    def test_decimer_memory_usage_over_time(self, benchmark_data):
        """Test DECIMER memory usage over extended operation"""
        pdf_path, temp_dir = benchmark_data.create_benchmark_pdf("simple", pages=1)
        
        try:
            import psutil
            client = DECIMERClient()
            
            memory_samples = []
            initial_memory = psutil.Process().memory_info().rss / 1024 / 1024
            memory_samples.append(initial_memory)
            
            # Process multiple times to check for memory leaks
            for i in range(10):
                segments = client.segment_structures_from_pdf(pdf_path)
                
                current_memory = psutil.Process().memory_info().rss / 1024 / 1024
                memory_samples.append(current_memory)
                
                # Memory shouldn't grow significantly
                memory_increase = current_memory - initial_memory
                assert memory_increase < 200, f"Memory leak detected: {memory_increase:.2f}MB after {i+1} iterations"
            
            # Clear cache and force cleanup
            client.clear_cache()
            import gc
            gc.collect()
            
            final_memory = psutil.Process().memory_info().rss / 1024 / 1024
            total_increase = final_memory - initial_memory
            
            print(f"Memory usage over time: initial={initial_memory:.1f}MB, final={final_memory:.1f}MB, increase={total_increase:.1f}MB")
            
        finally:
            shutil.rmtree(temp_dir)
    
    def test_decimer_cpu_vs_gpu_performance(self, benchmark_data):
        """Compare DECIMER performance on CPU vs GPU"""
        pdf_path, temp_dir = benchmark_data.create_benchmark_pdf("simple", pages=1)
        
        try:
            # Test CPU performance
            with patch.dict(os.environ, {'CUDA_VISIBLE_DEVICES': ''}):
                with patch('tensorflow.config.list_physical_devices', return_value=[]):
                    with PerformanceTimer("DECIMER CPU mode") as cpu_timer:
                        cpu_client = DECIMERClient()
                        cpu_segments = cpu_client.segment_structures_from_pdf(pdf_path)
            
            # Test GPU performance (mocked)
            mock_gpu = MagicMock()
            with patch('tensorflow.config.list_physical_devices', return_value=[mock_gpu]):
                with patch('tensorflow.config.experimental.set_memory_growth'):
                    with PerformanceTimer("DECIMER GPU mode") as gpu_timer:
                        gpu_client = DECIMERClient()
                        gpu_segments = gpu_client.segment_structures_from_pdf(pdf_path)
            
            # Both should work
            assert isinstance(cpu_segments, list), "CPU mode should work"
            assert isinstance(gpu_segments, list), "GPU mode should work"
            
            # Performance comparison
            print(f"CPU: {cpu_timer.duration:.3f}s, GPU: {gpu_timer.duration:.3f}s")
            
            # Both should complete in reasonable time
            assert cpu_timer.duration < 60.0, "CPU mode too slow"
            assert gpu_timer.duration < 60.0, "GPU mode too slow"
            
        finally:
            shutil.rmtree(temp_dir)


class TestEndToEndPerformance:
    """End-to-end performance benchmarks"""
    
    def test_full_pipeline_benchmark(self, benchmark_data):
        """Benchmark complete PDF to SMILES pipeline"""
        pdf_path, temp_dir = benchmark_data.create_benchmark_pdf("medium", pages=2, structures_per_page=2)
        
        try:
            client = DECIMERClient()
            
            with PerformanceTimer("Full pipeline") as timer:
                # Complete workflow
                result = client.process_pdf_complete(pdf_path, output_dir=temp_dir)
            
            assert timer.duration < 120.0, f"Full pipeline too slow: {timer.duration:.3f}s"
            assert "total_structures" in result, "Should return structure count"
            
            print(f"Full pipeline: {timer.duration:.3f}s, structures: {result.get('total_structures', 0)}")
            
        finally:
            shutil.rmtree(temp_dir)
    
    def test_stress_test_performance(self, benchmark_data):
        """Stress test with multiple large PDFs"""
        pdf_files = []
        
        try:
            # Create multiple test PDFs
            for i in range(3):
                pdf_path, temp_dir = benchmark_data.create_benchmark_pdf("complex", pages=3)
                pdf_files.append((pdf_path, temp_dir))
            
            def process_single_pdf(pdf_data):
                pdf_path, _ = pdf_data
                processor = PDFProcessor(pdf_path)
                result = processor.extract_text_and_images()
                processor.cleanup_temp_files()
                return len(result)
            
            with PerformanceTimer("Stress test") as timer:
                # Process all PDFs sequentially
                results = []
                for pdf_data in pdf_files:
                    result = process_single_pdf(pdf_data)
                    results.append(result)
            
            assert timer.duration < 180.0, f"Stress test too slow: {timer.duration:.3f}s"
            assert all(r == 3 for r in results), "All PDFs should be processed correctly"
            
            print(f"Stress test: {timer.duration:.3f}s for {len(pdf_files)} PDFs")
            
        finally:
            for _, temp_dir in pdf_files:
                shutil.rmtree(temp_dir)
    
    def test_performance_regression_detection(self, benchmark_data):
        """Test for performance regressions"""
        pdf_path, temp_dir = benchmark_data.create_benchmark_pdf("medium", pages=2)
        
        try:
            # Baseline measurements
            baseline_times = []
            for run in range(5):
                with PerformanceTimer(f"Baseline run {run+1}") as timer:
                    processor = PDFProcessor(pdf_path)
                    result = processor.extract_text_and_images()
                    processor.cleanup_temp_files()
                
                baseline_times.append(timer.duration)
            
            baseline_avg = statistics.mean(baseline_times)
            baseline_std = statistics.stdev(baseline_times)
            
            # Performance thresholds (can be adjusted based on expected performance)
            expected_max_time = 20.0  # Maximum acceptable time
            regression_threshold = baseline_avg + 3 * baseline_std  # 3 sigma rule
            
            assert baseline_avg < expected_max_time, f"Performance below expectations: {baseline_avg:.3f}s"
            
            print(f"Performance baseline: {baseline_avg:.3f}s ± {baseline_std:.3f}s")
            print(f"Regression threshold: {regression_threshold:.3f}s")
            
            # This could be extended to compare against historical data
            
        finally:
            shutil.rmtree(temp_dir)


class TestPerformanceUnderLoad:
    """Test performance under various load conditions"""
    
    def test_memory_pressure_performance(self, benchmark_data):
        """Test performance under memory pressure"""
        pdf_path, temp_dir = benchmark_data.create_benchmark_pdf("simple", pages=1)
        
        try:
            # Create memory pressure
            memory_hogs = []
            
            def create_memory_pressure():
                for i in range(5):
                    # Allocate 100MB chunks
                    memory_hogs.append(np.zeros((100 * 1024 * 1024 // 8,), dtype=np.float64))
                    time.sleep(0.1)
            
            pressure_thread = threading.Thread(target=create_memory_pressure, daemon=True)
            pressure_thread.start()
            
            time.sleep(1)  # Let pressure build
            
            with PerformanceTimer("Processing under memory pressure") as timer:
                processor = PDFProcessor(pdf_path)
                result = processor.extract_text_and_images()
                processor.cleanup_temp_files()
            
            # Should still complete but may be slower
            assert timer.duration < 60.0, f"Performance under memory pressure too slow: {timer.duration:.3f}s"
            assert len(result) > 0, "Should still process successfully"
            
            # Clean up memory pressure
            memory_hogs.clear()
            import gc
            gc.collect()
            
        finally:
            shutil.rmtree(temp_dir)
    
    def test_concurrent_load_performance(self, benchmark_data):
        """Test performance under concurrent load"""
        # Create test PDFs
        pdf_data = []
        for i in range(6):
            pdf_path, temp_dir = benchmark_data.create_benchmark_pdf("simple", pages=1)
            pdf_data.append((pdf_path, temp_dir))
        
        try:
            def process_pdf_worker(pdf_data):
                pdf_path, _ = pdf_data
                processor = PDFProcessor(pdf_path)
                result = processor.extract_text_and_images()
                processor.cleanup_temp_files()
                return len(result), time.time()
            
            with PerformanceTimer("Concurrent load test") as timer:
                with mp.Pool(6) as pool:
                    async_results = [pool.apply_async(process_pdf_worker, (pdf,)) for pdf in pdf_data]
                    results = [r.get(timeout=45) for r in async_results]
            
            # All should complete successfully
            assert len(results) == 6, "All processes should complete"
            assert all(r[0] == 1 for r in results), "All should process 1 page"
            
            # Should complete in reasonable time even under load
            assert timer.duration < 90.0, f"Concurrent load test too slow: {timer.duration:.3f}s"
            
            print(f"Concurrent load test: {timer.duration:.3f}s for {len(pdf_data)} PDFs")
            
        finally:
            for _, temp_dir in pdf_data:
                shutil.rmtree(temp_dir)


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s", "--tb=short"])