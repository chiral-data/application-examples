"""
Comprehensive tests for PDF processing timeout and hanging issues.
These tests specifically target the DECIMER segmentation hanging problem.
"""

import pytest
import os
import tempfile
import threading
import time
import signal
import psutil
import subprocess
import numpy as np
from unittest.mock import patch, MagicMock, call
from PIL import Image
import fitz  # PyMuPDF
import multiprocessing as mp

from pdf_processor import PDFProcessor
from decimer_client import DECIMERClient


class TimeoutError(Exception):
    """Custom timeout exception"""
    pass


class ProcessTimeoutHandler:
    """Context manager for handling process timeouts"""
    
    def __init__(self, timeout_seconds):
        self.timeout_seconds = timeout_seconds
        self.timer = None
        
    def timeout_handler(self):
        raise TimeoutError(f"Process timed out after {self.timeout_seconds} seconds")
        
    def __enter__(self):
        self.timer = threading.Timer(self.timeout_seconds, self.timeout_handler)
        self.timer.start()
        return self
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.timer:
            self.timer.cancel()


@pytest.fixture
def create_test_pdf():
    """Create various test PDF files with different complexities"""
    def _create_pdf(pdf_type="simple", num_pages=1, structures_per_page=1):
        temp_dir = tempfile.mkdtemp()
        pdf_path = os.path.join(temp_dir, f"test_{pdf_type}.pdf")
        
        # Create PDF with fitz
        doc = fitz.open()
        
        for page_num in range(num_pages):
            page = doc.new_page()
            
            if pdf_type == "simple":
                # Simple text-only PDF
                page.insert_text((50, 50), f"Simple text page {page_num + 1}")
                
            elif pdf_type == "complex_images":
                # PDF with many complex images
                page.insert_text((50, 50), f"Complex page {page_num + 1}")
                
                # Add multiple mock chemical structure images
                for i in range(structures_per_page):
                    # Create a complex image that might cause DECIMER to hang
                    img_array = np.random.randint(0, 255, (200, 200, 3), dtype=np.uint8)
                    img = Image.fromarray(img_array)
                    
                    img_bytes = img.tobytes()
                    img_rect = fitz.Rect(50 + i * 100, 100, 150 + i * 100, 200)
                    
            elif pdf_type == "corrupted_images":
                # PDF with corrupted or problematic images
                page.insert_text((50, 50), "Page with corrupted images")
                
            elif pdf_type == "large_file":
                # Large PDF that might cause memory issues
                for i in range(50):  # Add lots of content
                    page.insert_text((50, 50 + i * 15), f"Line {i} with chemical formula C6H6")
                    
        doc.save(pdf_path)
        doc.close()
        
        return pdf_path, temp_dir
    
    return _create_pdf


@pytest.fixture
def monitor_system_resources():
    """Monitor system resources during test execution"""
    class ResourceMonitor:
        def __init__(self):
            self.initial_memory = psutil.virtual_memory().percent
            self.initial_cpu = psutil.cpu_percent()
            self.peak_memory = self.initial_memory
            self.peak_cpu = self.initial_cpu
            self.monitoring = False
            
        def start_monitoring(self):
            self.monitoring = True
            threading.Thread(target=self._monitor, daemon=True).start()
            
        def stop_monitoring(self):
            self.monitoring = False
            
        def _monitor(self):
            while self.monitoring:
                self.peak_memory = max(self.peak_memory, psutil.virtual_memory().percent)
                self.peak_cpu = max(self.peak_cpu, psutil.cpu_percent())
                time.sleep(0.1)
                
        def get_stats(self):
            return {
                'initial_memory': self.initial_memory,
                'peak_memory': self.peak_memory,
                'memory_increase': self.peak_memory - self.initial_memory,
                'initial_cpu': self.initial_cpu,
                'peak_cpu': self.peak_cpu
            }
    
    return ResourceMonitor()


class TestPDFProcessingTimeouts:
    """Test suite for PDF processing timeout and hanging issues"""
    
    def test_simple_pdf_processing_with_timeout(self, create_test_pdf, monitor_system_resources):
        """Test that simple PDF processing completes within reasonable time"""
        pdf_path, temp_dir = create_test_pdf("simple", num_pages=1)
        
        try:
            monitor_system_resources.start_monitoring()
            
            with ProcessTimeoutHandler(30):  # 30 second timeout
                processor = PDFProcessor(pdf_path)
                result = processor.extract_text_and_images()
                
                assert len(result) > 0
                assert result[0]['text'] is not None
                
        finally:
            monitor_system_resources.stop_monitoring()
            import shutil
            shutil.rmtree(temp_dir)
            
        stats = monitor_system_resources.get_stats()
        assert stats['memory_increase'] < 50  # Should not use more than 50% additional memory
    
    def test_decimer_segmentation_timeout(self, create_test_pdf):
        """Test DECIMER segmentation with timeout to prevent hanging"""
        pdf_path, temp_dir = create_test_pdf("complex_images", num_pages=1, structures_per_page=3)
        
        try:
            processor = PDFProcessor(pdf_path)
            
            # Test with timeout using multiprocessing
            def run_decimer_segmentation():
                return processor.extract_chemical_structures_with_decimer()
            
            # Use multiprocessing to enforce hard timeout
            with mp.Pool(1) as pool:
                result = pool.apply_async(run_decimer_segmentation)
                try:
                    # 60 second timeout for DECIMER processing
                    structures = result.get(timeout=60)
                    assert isinstance(structures, list)
                except mp.TimeoutError:
                    pytest.fail("DECIMER segmentation timed out - this indicates the hanging issue")
                    
        finally:
            import shutil
            shutil.rmtree(temp_dir)
    
    def test_decimer_segmentation_with_mocked_hanging(self, create_test_pdf):
        """Test DECIMER segmentation behavior when the underlying function hangs"""
        pdf_path, temp_dir = create_test_pdf("simple")
        
        def mock_hanging_segment_function(*args, **kwargs):
            """Mock function that simulates hanging"""
            time.sleep(120)  # Simulate hanging for 2 minutes
            return []
        
        try:
            with patch('pdf_processor.segment_chemical_structures_from_file', side_effect=mock_hanging_segment_function):
                processor = PDFProcessor(pdf_path)
                
                # This should timeout and fall back to regular extraction
                with ProcessTimeoutHandler(10):  # 10 second timeout
                    try:
                        result = processor.extract_chemical_structures_with_decimer()
                        pytest.fail("Should have timed out")
                    except TimeoutError:
                        # This is expected - the function should timeout
                        pass
                        
        finally:
            import shutil
            shutil.rmtree(temp_dir)
    
    def test_memory_leak_detection(self, create_test_pdf, monitor_system_resources):
        """Test for memory leaks during repeated PDF processing"""
        pdf_path, temp_dir = create_test_pdf("simple")
        
        try:
            monitor_system_resources.start_monitoring()
            initial_memory = psutil.Process().memory_info().rss
            
            # Process the same PDF multiple times
            for i in range(10):
                processor = PDFProcessor(pdf_path)
                result = processor.extract_text_and_images()
                processor.cleanup_temp_files()
                
                # Manually trigger garbage collection
                import gc
                gc.collect()
                
                current_memory = psutil.Process().memory_info().rss
                memory_increase = (current_memory - initial_memory) / 1024 / 1024  # MB
                
                # Memory should not increase significantly
                assert memory_increase < 100, f"Memory leak detected: {memory_increase:.2f}MB increase"
                
        finally:
            monitor_system_resources.stop_monitoring()
            import shutil
            shutil.rmtree(temp_dir)
    
    def test_gpu_cpu_fallback_behavior(self, create_test_pdf):
        """Test GPU/CPU fallback behavior and timeout handling"""
        pdf_path, temp_dir = create_test_pdf("simple")
        
        try:
            # Test with GPU unavailable
            with patch('tensorflow.config.list_physical_devices', return_value=[]):
                client = DECIMERClient()
                
                # Should still work on CPU
                with ProcessTimeoutHandler(30):
                    segments = client.segment_structures_from_pdf(pdf_path)
                    assert isinstance(segments, list)
                    
            # Test with TensorFlow import error
            with patch('decimer_client.tf', side_effect=ImportError("TensorFlow not available")):
                client = DECIMERClient()
                assert not client.decimer_available
                
        finally:
            import shutil
            shutil.rmtree(temp_dir)
    
    def test_large_pdf_processing_timeout(self, create_test_pdf, monitor_system_resources):
        """Test processing of large PDFs with timeout"""
        pdf_path, temp_dir = create_test_pdf("large_file", num_pages=5)
        
        try:
            monitor_system_resources.start_monitoring()
            
            processor = PDFProcessor(pdf_path)
            
            # Test with reasonable timeout for large files
            with ProcessTimeoutHandler(120):  # 2 minute timeout
                result = processor.extract_text_and_images()
                
                assert len(result) == 5  # Should process all 5 pages
                assert all(page['text'] for page in result)
                
        finally:
            monitor_system_resources.stop_monitoring()
            import shutil
            shutil.rmtree(temp_dir)
            
        stats = monitor_system_resources.get_stats()
        # Large files should not cause excessive memory usage
        assert stats['memory_increase'] < 80
    
    def test_concurrent_pdf_processing(self, create_test_pdf):
        """Test concurrent PDF processing to identify race conditions"""
        pdf_paths = []
        temp_dirs = []
        
        # Create multiple test PDFs
        for i in range(3):
            pdf_path, temp_dir = create_test_pdf("simple", num_pages=1)
            pdf_paths.append(pdf_path)
            temp_dirs.append(temp_dir)
        
        try:
            def process_pdf(pdf_path):
                processor = PDFProcessor(pdf_path)
                return processor.extract_text_and_images()
            
            # Process PDFs concurrently
            results = []
            with mp.Pool(3) as pool:
                async_results = [pool.apply_async(process_pdf, (pdf_path,)) for pdf_path in pdf_paths]
                
                for async_result in async_results:
                    try:
                        result = async_result.get(timeout=30)
                        results.append(result)
                    except mp.TimeoutError:
                        pytest.fail("Concurrent PDF processing timed out")
            
            assert len(results) == 3
            assert all(len(result) > 0 for result in results)
            
        finally:
            import shutil
            for temp_dir in temp_dirs:
                shutil.rmtree(temp_dir)
    
    def test_decimer_error_handling_and_fallback(self, create_test_pdf):
        """Test error handling and fallback when DECIMER fails"""
        pdf_path, temp_dir = create_test_pdf("simple")
        
        try:
            processor = PDFProcessor(pdf_path)
            
            # Test with DECIMER segmentation throwing exception
            with patch('pdf_processor.segment_chemical_structures_from_file', side_effect=RuntimeError("DECIMER error")):
                result = processor.extract_chemical_structures_with_decimer()
                
                # Should fall back to regular extraction
                assert isinstance(result, list)
                assert len(result) > 0  # Should have extracted something
            
            # Test with DECIMER not available
            with patch('pdf_processor.DECIMER_SEGMENTATION_AVAILABLE', False):
                result = processor.extract_chemical_structures_with_decimer()
                
                # Should fall back to regular extraction
                assert isinstance(result, list)
                
        finally:
            import shutil
            shutil.rmtree(temp_dir)
    
    def test_process_monitoring_and_termination(self, create_test_pdf):
        """Test process monitoring and forced termination of hanging processes"""
        pdf_path, temp_dir = create_test_pdf("simple")
        
        def hanging_process_target():
            """Function that simulates a hanging process"""
            processor = PDFProcessor(pdf_path)
            
            # Mock hanging DECIMER call
            with patch('pdf_processor.segment_chemical_structures_from_file') as mock_segment:
                mock_segment.side_effect = lambda *args, **kwargs: time.sleep(300)  # Hang for 5 minutes
                processor.extract_chemical_structures_with_decimer()
        
        try:
            # Start hanging process
            process = mp.Process(target=hanging_process_target)
            process.start()
            
            # Wait for a short time then terminate
            time.sleep(2)
            
            if process.is_alive():
                process.terminate()
                process.join(timeout=5)
                
                if process.is_alive():
                    process.kill()
                    process.join()
            
            assert not process.is_alive(), "Process should be terminated"
            
        finally:
            import shutil
            shutil.rmtree(temp_dir)


class TestDECIMERSpecificIssues:
    """Test suite specifically for DECIMER-related hanging issues"""
    
    def test_decimer_tensorflow_gpu_memory_growth(self):
        """Test TensorFlow GPU memory growth configuration"""
        with patch('tensorflow.config.list_physical_devices') as mock_list_devices:
            mock_gpu = MagicMock()
            mock_list_devices.return_value = [mock_gpu]
            
            with patch('tensorflow.config.experimental.set_memory_growth') as mock_memory_growth:
                client = DECIMERClient()
                
                # Should have attempted to set memory growth
                mock_memory_growth.assert_called_with(mock_gpu, True)
    
    def test_decimer_batch_processing_timeout(self, create_test_pdf):
        """Test batch processing with timeout to prevent hanging"""
        pdf_path, temp_dir = create_test_pdf("simple")
        
        try:
            # Create test images
            test_images = []
            for i in range(3):
                img_path = os.path.join(temp_dir, f"test_img_{i}.png")
                img = Image.new('RGB', (100, 100), color='white')
                img.save(img_path)
                test_images.append(img_path)
            
            client = DECIMERClient()
            
            # Mock image_to_smiles to simulate hanging
            def mock_hanging_image_to_smiles(image_path):
                if "test_img_1" in image_path:  # Make one image hang
                    time.sleep(60)
                return "CCO"
            
            with patch.object(client, 'image_to_smiles', side_effect=mock_hanging_image_to_smiles):
                with ProcessTimeoutHandler(10):  # Short timeout
                    try:
                        client.batch_process(test_images)
                        pytest.fail("Should have timed out")
                    except TimeoutError:
                        # Expected timeout
                        pass
                        
        finally:
            import shutil
            shutil.rmtree(temp_dir)
    
    def test_decimer_segmentation_memory_monitoring(self, create_test_pdf, monitor_system_resources):
        """Monitor memory usage during DECIMER segmentation"""
        pdf_path, temp_dir = create_test_pdf("complex_images", structures_per_page=5)
        
        try:
            monitor_system_resources.start_monitoring()
            
            client = DECIMERClient()
            
            # Monitor memory during segmentation
            initial_memory = psutil.Process().memory_info().rss / 1024 / 1024  # MB
            
            with ProcessTimeoutHandler(60):
                segments = client.segment_structures_from_pdf(pdf_path)
                
            final_memory = psutil.Process().memory_info().rss / 1024 / 1024  # MB
            memory_increase = final_memory - initial_memory
            
            # Memory increase should be reasonable
            assert memory_increase < 1000, f"Excessive memory usage: {memory_increase:.2f}MB"
            
        finally:
            monitor_system_resources.stop_monitoring()
            import shutil
            shutil.rmtree(temp_dir)


if __name__ == "__main__":
    # Run specific tests for debugging
    pytest.main([__file__, "-v", "-s"])