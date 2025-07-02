"""
Memory usage monitoring and resource management tests.
These tests help identify memory leaks and resource issues that could cause hanging.
"""

import pytest
import os
import time
import threading
import psutil
import gc
import tempfile
import shutil
from unittest.mock import patch, MagicMock
import numpy as np
from PIL import Image

from pdf_processor import PDFProcessor
from decimer_client import DECIMERClient


class ResourceMonitor:
    """Monitor system resources during test execution"""
    
    def __init__(self):
        self.process = psutil.Process()
        self.monitoring = False
        self.samples = []
        self.monitor_thread = None
        
    def start_monitoring(self, interval=0.1):
        """Start monitoring resources"""
        self.monitoring = True
        self.samples = []
        self.monitor_thread = threading.Thread(target=self._monitor_loop, args=(interval,), daemon=True)
        self.monitor_thread.start()
        
    def stop_monitoring(self):
        """Stop monitoring and return statistics"""
        self.monitoring = False
        if self.monitor_thread:
            self.monitor_thread.join(timeout=1)
        
        if not self.samples:
            return {}
            
        memory_samples = [s['memory_mb'] for s in self.samples]
        cpu_samples = [s['cpu_percent'] for s in self.samples]
        
        return {
            'memory_initial': memory_samples[0],
            'memory_peak': max(memory_samples),
            'memory_final': memory_samples[-1],
            'memory_increase': max(memory_samples) - memory_samples[0],
            'cpu_peak': max(cpu_samples),
            'cpu_average': sum(cpu_samples) / len(cpu_samples),
            'samples_count': len(self.samples),
            'duration': self.samples[-1]['timestamp'] - self.samples[0]['timestamp']
        }
        
    def _monitor_loop(self, interval):
        """Internal monitoring loop"""
        while self.monitoring:
            try:
                memory_info = self.process.memory_info()
                cpu_percent = self.process.cpu_percent()
                
                sample = {
                    'timestamp': time.time(),
                    'memory_mb': memory_info.rss / 1024 / 1024,
                    'memory_vms': memory_info.vms / 1024 / 1024,
                    'cpu_percent': cpu_percent,
                    'num_threads': self.process.num_threads(),
                    'open_files': len(self.process.open_files())
                }
                
                self.samples.append(sample)
                time.sleep(interval)
                
            except Exception as e:
                print(f"Error monitoring resources: {e}")
                break


@pytest.fixture
def resource_monitor():
    """Provide resource monitoring fixture"""
    return ResourceMonitor()


@pytest.fixture
def create_memory_test_pdf():
    """Create PDFs designed to test memory usage"""
    
    def _create_pdf(size_mb=1, num_pages=1):
        temp_dir = tempfile.mkdtemp()
        pdf_path = os.path.join(temp_dir, f"memory_test_{size_mb}mb.pdf")
        
        # Create a PDF with specific memory requirements
        import fitz
        doc = fitz.open()
        
        for page_num in range(num_pages):
            page = doc.new_page()
            
            # Add text content
            text_content = "Chemical structure analysis " * 100
            page.insert_text((50, 50), text_content)
            
            # Add large image if needed
            if size_mb > 1:
                # Create large image data
                img_size = int((size_mb * 1024 * 1024 / num_pages) ** 0.5)
                img_array = np.random.randint(0, 255, (img_size, img_size, 3), dtype=np.uint8)
                img = Image.fromarray(img_array)
                
                # Save temporary image
                temp_img_path = os.path.join(temp_dir, f"large_img_{page_num}.png")
                img.save(temp_img_path)
        
        doc.save(pdf_path)
        doc.close()
        
        return pdf_path, temp_dir
    
    return _create_pdf


class TestMemoryUsageMonitoring:
    """Test memory usage patterns and potential leaks"""
    
    def test_basic_pdf_processing_memory_usage(self, create_memory_test_pdf, resource_monitor):
        """Test memory usage of basic PDF processing"""
        pdf_path, temp_dir = create_memory_test_pdf(size_mb=2, num_pages=3)
        
        try:
            resource_monitor.start_monitoring()
            
            processor = PDFProcessor(pdf_path)
            result = processor.extract_text_and_images()
            processor.cleanup_temp_files()
            
            # Force garbage collection
            gc.collect()
            
            stats = resource_monitor.stop_monitoring()
            
            # Assertions about memory usage
            assert stats['memory_increase'] < 100, f"Excessive memory usage: {stats['memory_increase']:.2f}MB"
            assert len(result) == 3, "Should process all pages"
            assert stats['duration'] > 0, "Should have monitoring data"
            
            print(f"Memory usage stats: {stats}")
            
        finally:
            shutil.rmtree(temp_dir)
    
    def test_decimer_client_memory_usage(self, create_memory_test_pdf, resource_monitor):
        """Test memory usage of DECIMER client operations"""
        pdf_path, temp_dir = create_memory_test_pdf(size_mb=1, num_pages=2)
        
        try:
            resource_monitor.start_monitoring()
            
            client = DECIMERClient()
            
            # Test segmentation
            segments = client.segment_structures_from_pdf(pdf_path)
            
            # Clear cache
            client.clear_cache()
            
            # Force garbage collection
            gc.collect()
            
            stats = resource_monitor.stop_monitoring()
            
            # Memory should not grow excessively
            assert stats['memory_increase'] < 200, f"Excessive DECIMER memory usage: {stats['memory_increase']:.2f}MB"
            assert isinstance(segments, list), "Should return list of segments"
            
            print(f"DECIMER memory stats: {stats}")
            
        finally:
            shutil.rmtree(temp_dir)
    
    def test_repeated_processing_memory_leak(self, create_memory_test_pdf, resource_monitor):
        """Test for memory leaks in repeated processing"""
        pdf_path, temp_dir = create_memory_test_pdf(size_mb=1, num_pages=1)
        
        try:
            resource_monitor.start_monitoring()
            
            initial_memory = psutil.Process().memory_info().rss / 1024 / 1024
            memory_samples = [initial_memory]
            
            # Process the same PDF multiple times
            for i in range(10):
                processor = PDFProcessor(pdf_path)
                result = processor.extract_text_and_images()
                processor.cleanup_temp_files()
                
                # Force cleanup
                del processor
                gc.collect()
                
                current_memory = psutil.Process().memory_info().rss / 1024 / 1024
                memory_samples.append(current_memory)
                
                # Check for gradual memory increase (leak)
                if i > 5:  # After warmup
                    memory_increase = current_memory - initial_memory
                    assert memory_increase < 50, f"Potential memory leak detected: {memory_increase:.2f}MB after {i+1} iterations"
            
            stats = resource_monitor.stop_monitoring()
            
            # Final memory should not be significantly higher
            final_increase = memory_samples[-1] - memory_samples[0]
            assert final_increase < 30, f"Memory leak detected: {final_increase:.2f}MB total increase"
            
            print(f"Memory leak test stats: {stats}")
            print(f"Memory progression: {[f'{m:.1f}' for m in memory_samples[:5]]}...")
            
        finally:
            shutil.rmtree(temp_dir)
    
    def test_large_pdf_memory_scaling(self, create_memory_test_pdf, resource_monitor):
        """Test memory usage scaling with PDF size"""
        memory_results = {}
        
        # Test different PDF sizes
        for size_mb in [1, 5, 10]:
            pdf_path, temp_dir = create_memory_test_pdf(size_mb=size_mb, num_pages=2)
            
            try:
                resource_monitor.start_monitoring()
                
                processor = PDFProcessor(pdf_path)
                result = processor.extract_text_and_images()
                processor.cleanup_temp_files()
                
                gc.collect()
                
                stats = resource_monitor.stop_monitoring()
                memory_results[size_mb] = stats['memory_increase']
                
                # Memory usage should scale reasonably with input size
                assert stats['memory_increase'] < size_mb * 50, f"Memory usage too high for {size_mb}MB PDF"
                
                print(f"PDF {size_mb}MB -> Memory increase: {stats['memory_increase']:.2f}MB")
                
            finally:
                shutil.rmtree(temp_dir)
        
        # Check scaling relationship
        assert memory_results[10] > memory_results[1], "Memory usage should increase with PDF size"
        
        # But not linearly (should be sub-linear due to processing optimizations)
        scaling_factor = memory_results[10] / memory_results[1]
        assert scaling_factor < 15, f"Memory scaling too aggressive: {scaling_factor}x"
    
    def test_concurrent_processing_memory_isolation(self, create_memory_test_pdf, resource_monitor):
        """Test memory isolation in concurrent processing"""
        import multiprocessing as mp
        
        pdf_path, temp_dir = create_memory_test_pdf(size_mb=2, num_pages=1)
        
        def process_pdf_worker(pdf_path):
            """Worker function for multiprocessing"""
            processor = PDFProcessor(pdf_path)
            result = processor.extract_text_and_images()
            processor.cleanup_temp_files()
            return len(result)
        
        try:
            resource_monitor.start_monitoring()
            
            # Process PDFs concurrently
            with mp.Pool(3) as pool:
                results = [pool.apply_async(process_pdf_worker, (pdf_path,)) for _ in range(3)]
                outputs = [r.get(timeout=30) for r in results]
            
            gc.collect()
            
            stats = resource_monitor.stop_monitoring()
            
            # All processes should complete successfully
            assert all(output == 1 for output in outputs), "All processes should process 1 page"
            
            # Memory usage should be reasonable even with concurrent processing
            assert stats['memory_increase'] < 150, f"Excessive concurrent memory usage: {stats['memory_increase']:.2f}MB"
            
            print(f"Concurrent processing memory stats: {stats}")
            
        finally:
            shutil.rmtree(temp_dir)
    
    def test_tensorflow_memory_management(self, resource_monitor):
        """Test TensorFlow memory management in DECIMER"""
        
        # Test GPU memory growth setting
        with patch('tensorflow.config.list_physical_devices') as mock_devices:
            mock_gpu = MagicMock()
            mock_devices.return_value = [mock_gpu]
            
            with patch('tensorflow.config.experimental.set_memory_growth') as mock_memory_growth:
                resource_monitor.start_monitoring()
                
                # Initialize DECIMER client
                client = DECIMERClient()
                
                stats = resource_monitor.stop_monitoring()
                
                # Should have attempted to set memory growth
                mock_memory_growth.assert_called_with(mock_gpu, True)
                
                # Memory usage should be reasonable
                assert stats['memory_increase'] < 100, "TensorFlow initialization uses too much memory"
    
    def test_image_processing_memory_cleanup(self, resource_monitor):
        """Test memory cleanup in image processing operations"""
        temp_dir = tempfile.mkdtemp()
        
        try:
            resource_monitor.start_monitoring()
            
            # Create multiple test images
            image_paths = []
            for i in range(10):
                img_array = np.random.randint(0, 255, (500, 500, 3), dtype=np.uint8)
                img = Image.fromarray(img_array)
                img_path = os.path.join(temp_dir, f"test_img_{i}.png")
                img.save(img_path)
                image_paths.append(img_path)
            
            client = DECIMERClient()
            
            # Process images in batch
            results = client.batch_process(image_paths)
            
            # Clear cache and force cleanup
            client.clear_cache()
            del client
            gc.collect()
            
            stats = resource_monitor.stop_monitoring()
            
            # Memory should not grow excessively
            assert stats['memory_increase'] < 200, f"Image processing memory leak: {stats['memory_increase']:.2f}MB"
            assert len(results) == 10, "Should process all images"
            
            print(f"Image processing memory stats: {stats}")
            
        finally:
            shutil.rmtree(temp_dir)
    
    def test_memory_pressure_handling(self, create_memory_test_pdf, resource_monitor):
        """Test behavior under memory pressure"""
        pdf_path, temp_dir = create_memory_test_pdf(size_mb=1, num_pages=1)
        
        try:
            # Create memory pressure
            memory_hogs = []
            
            def create_memory_pressure():
                # Allocate memory chunks to create pressure
                for i in range(3):
                    # Allocate 200MB chunks
                    memory_hogs.append(np.zeros((200 * 1024 * 1024 // 8,), dtype=np.float64))
                    time.sleep(0.1)
            
            pressure_thread = threading.Thread(target=create_memory_pressure, daemon=True)
            pressure_thread.start()
            
            time.sleep(0.5)  # Let pressure build
            
            resource_monitor.start_monitoring()
            
            # Try to process PDF under memory pressure
            processor = PDFProcessor(pdf_path)
            result = processor.extract_text_and_images()
            processor.cleanup_temp_files()
            
            stats = resource_monitor.stop_monitoring()
            
            # Should still work under memory pressure
            assert len(result) > 0, "Should process PDF even under memory pressure"
            
            # Clean up memory pressure
            memory_hogs.clear()
            gc.collect()
            
            print(f"Memory pressure test stats: {stats}")
            
        finally:
            shutil.rmtree(temp_dir)


class TestResourceCleanup:
    """Test proper cleanup of resources"""
    
    def test_pdf_processor_cleanup(self, create_memory_test_pdf):
        """Test PDFProcessor resource cleanup"""
        pdf_path, temp_dir = create_memory_test_pdf(size_mb=1, num_pages=2)
        
        try:
            processor = PDFProcessor(pdf_path)
            
            # Process and create temp files
            result = processor.extract_text_and_images()
            
            # Check that temp files were created
            temp_files_before = os.listdir(processor.temp_dir) if os.path.exists(processor.temp_dir) else []
            
            # Cleanup
            processor.cleanup_temp_files()
            
            # Check cleanup
            temp_files_after = os.listdir(processor.temp_dir) if os.path.exists(processor.temp_dir) else []
            
            # Should have fewer temp files after cleanup
            assert len(temp_files_after) <= len(temp_files_before), "Cleanup should remove temp files"
            
            # Test destructor cleanup
            del processor
            gc.collect()
            
        finally:
            shutil.rmtree(temp_dir)
    
    def test_decimer_client_cache_cleanup(self):
        """Test DECIMER client cache cleanup"""
        client = DECIMERClient()
        
        # Add some mock cache entries
        client.cache = {
            'key1': 'CCO',
            'key2': 'C6H6',
            'key3': 'CC(=O)O'
        }
        
        assert len(client.cache) == 3, "Cache should have test entries"
        
        # Clear cache
        client.clear_cache()
        
        assert len(client.cache) == 0, "Cache should be empty after clearing"
    
    def test_file_handle_cleanup(self, create_memory_test_pdf):
        """Test that file handles are properly closed"""
        pdf_path, temp_dir = create_memory_test_pdf(size_mb=1, num_pages=1)
        
        try:
            initial_open_files = len(psutil.Process().open_files())
            
            # Process PDF
            processor = PDFProcessor(pdf_path)
            result = processor.extract_text_and_images()
            
            # Cleanup
            processor.cleanup_temp_files()
            del processor
            gc.collect()
            
            final_open_files = len(psutil.Process().open_files())
            
            # Should not have significantly more open files
            file_increase = final_open_files - initial_open_files
            assert file_increase <= 2, f"File handle leak detected: {file_increase} additional files"
            
        finally:
            shutil.rmtree(temp_dir)


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])