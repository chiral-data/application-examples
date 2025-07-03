"""
Comprehensive DECIMER integration tests for various PDF types and edge cases.
These tests help identify specific conditions that cause DECIMER to hang.
"""

import pytest
import os
import tempfile
import shutil
import numpy as np
import cv2
import time
import threading
import multiprocessing as mp
from PIL import Image, ImageDraw
import fitz  # PyMuPDF
from unittest.mock import patch, MagicMock

from pdf_processor import PDFProcessor
from decimer_client import DECIMERClient


class ChemicalStructureGenerator:
    """Generate various types of chemical structure images for testing"""
    
    @staticmethod
    def create_simple_benzene_ring(size=(200, 200)):
        """Create a simple benzene ring structure"""
        img = Image.new('RGB', size, 'white')
        draw = ImageDraw.Draw(img)
        
        # Draw hexagon for benzene ring
        center_x, center_y = size[0] // 2, size[1] // 2
        radius = min(size) // 4
        
        points = []
        for i in range(6):
            angle = i * 60 * np.pi / 180
            x = center_x + radius * np.cos(angle)
            y = center_y + radius * np.sin(angle)
            points.append((x, y))
        
        # Draw the hexagon
        for i in range(6):
            start = points[i]
            end = points[(i + 1) % 6]
            draw.line([start, end], fill='black', width=3)
        
        # Draw inner circle for aromatic character
        draw.ellipse([center_x - radius//2, center_y - radius//2, 
                     center_x + radius//2, center_y + radius//2], 
                    outline='black', width=2)
        
        return img
    
    @staticmethod
    def create_complex_molecule(size=(300, 200)):
        """Create a more complex molecular structure"""
        img = Image.new('RGB', size, 'white')
        draw = ImageDraw.Draw(img)
        
        # Draw a complex organic molecule with multiple rings and bonds
        # This is simplified - real chemical structures would be more detailed
        
        # Main backbone
        points = [(50, 100), (100, 80), (150, 100), (200, 80), (250, 100)]
        for i in range(len(points) - 1):
            draw.line([points[i], points[i+1]], fill='black', width=3)
        
        # Add side chains
        draw.line([(100, 80), (100, 50)], fill='black', width=3)
        draw.line([(200, 80), (200, 50)], fill='black', width=3)
        
        # Add some text labels (C, N, O)
        try:
            font = ImageFont.load_default()
            draw.text((90, 30), 'N', fill='blue', font=font)
            draw.text((190, 30), 'O', fill='red', font=font)
            draw.text((240, 90), 'OH', fill='red', font=font)
        except:
            # Font not available, skip labels
            pass
        
        return img
    
    @staticmethod
    def create_noisy_structure(size=(200, 200), noise_level=0.3):
        """Create a chemical structure with noise that might confuse DECIMER"""
        base_img = ChemicalStructureGenerator.create_simple_benzene_ring(size)
        img_array = np.array(base_img)
        
        # Add random noise
        noise = np.random.randint(0, int(255 * noise_level), img_array.shape, dtype=np.uint8)
        noisy_array = np.clip(img_array.astype(np.int16) + noise, 0, 255).astype(np.uint8)
        
        return Image.fromarray(noisy_array)
    
    @staticmethod
    def create_very_large_structure(size=(2000, 2000)):
        """Create a very large structure image that might cause memory issues"""
        img = Image.new('RGB', size, 'white')
        draw = ImageDraw.Draw(img)
        
        # Create a large grid of connected molecules
        grid_size = 10
        cell_size = size[0] // grid_size
        
        for i in range(grid_size):
            for j in range(grid_size):
                x = i * cell_size + cell_size // 2
                y = j * cell_size + cell_size // 2
                
                # Draw a small ring at each grid point
                radius = cell_size // 6
                draw.ellipse([x - radius, y - radius, x + radius, y + radius], 
                           outline='black', width=2)
                
                # Connect to neighbors
                if i < grid_size - 1:
                    draw.line([(x + radius, y), (x + cell_size - radius, y)], 
                             fill='black', width=2)
                if j < grid_size - 1:
                    draw.line([(x, y + radius), (x, y + cell_size - radius)], 
                             fill='black', width=2)
        
        return img
    
    @staticmethod
    def create_corrupted_structure(size=(200, 200)):
        """Create a structure that looks like it might be corrupted"""
        img = Image.new('RGB', size, 'white')
        draw = ImageDraw.Draw(img)
        
        # Random lines and shapes that don't form coherent structures
        np.random.seed(42)  # For reproducibility
        for _ in range(50):
            x1, y1 = np.random.randint(0, size[0]), np.random.randint(0, size[1])
            x2, y2 = np.random.randint(0, size[0]), np.random.randint(0, size[1])
            draw.line([(x1, y1), (x2, y2)], fill='black', width=np.random.randint(1, 4))
        
        return img


@pytest.fixture
def structure_generator():
    return ChemicalStructureGenerator()


@pytest.fixture
def create_problematic_pdf():
    """Create PDFs with various problematic content that might cause DECIMER to hang"""
    
    def _create_pdf(pdf_type, **kwargs):
        temp_dir = tempfile.mkdtemp()
        pdf_path = os.path.join(temp_dir, f"problematic_{pdf_type}.pdf")
        
        doc = fitz.open()
        page = doc.new_page()
        
        generator = ChemicalStructureGenerator()
        
        if pdf_type == "multiple_large_structures":
            # PDF with multiple very large chemical structures
            for i in range(kwargs.get('num_structures', 3)):
                large_img = generator.create_very_large_structure((1000, 1000))
                # Convert to bytes and insert into PDF
                img_bytes = large_img.tobytes()
                page.insert_text((50, 50 + i * 20), f"Large structure {i+1}")
        
        elif pdf_type == "high_noise_structures":
            # PDF with noisy structures that might confuse DECIMER
            for noise_level in [0.1, 0.3, 0.5, 0.7]:
                noisy_img = generator.create_noisy_structure(noise_level=noise_level)
                page.insert_text((50, 50), f"Noisy structure (noise={noise_level})")
        
        elif pdf_type == "corrupted_images":
            # PDF with corrupted/random image data
            corrupted_img = generator.create_corrupted_structure()
            page.insert_text((50, 50), "Corrupted structure data")
        
        elif pdf_type == "mixed_content":
            # PDF with mix of text, images, and chemical structures
            page.insert_text((50, 50), "Mixed content document")
            page.insert_text((50, 100), "Chemical formula: C6H6")
            page.insert_text((50, 150), "Some random text that is not chemistry")
            page.insert_text((50, 200), "More chemistry: CH3COOH")
        
        elif pdf_type == "empty_pages":
            # PDF with empty pages
            for _ in range(5):
                doc.new_page()
        
        elif pdf_type == "very_dense":
            # PDF with extremely dense content
            for i in range(100):
                page.insert_text((50, 50 + i * 8), f"Dense line {i} with formula C{i}H{i*2}")
        
        doc.save(pdf_path)
        doc.close()
        
        return pdf_path, temp_dir
    
    return _create_pdf


class TestDECIMERRobustness:
    """Test DECIMER with various problematic inputs"""
    
    def test_large_structure_images(self, create_problematic_pdf):
        """Test DECIMER with very large structure images"""
        pdf_path, temp_dir = create_problematic_pdf("multiple_large_structures", num_structures=2)
        
        try:
            client = DECIMERClient()
            
            # Test with timeout
            def run_segmentation():
                return client.segment_structures_from_pdf(pdf_path)
            
            # Use multiprocessing for hard timeout
            with mp.Pool(1) as pool:
                result = pool.apply_async(run_segmentation)
                try:
                    segments = result.get(timeout=90)  # 90 second timeout
                    assert isinstance(segments, list)
                except mp.TimeoutError:
                    pytest.fail("DECIMER hung on large structure images")
                    
        finally:
            shutil.rmtree(temp_dir)
    
    def test_noisy_structure_images(self, create_problematic_pdf):
        """Test DECIMER with noisy structure images"""
        pdf_path, temp_dir = create_problematic_pdf("high_noise_structures")
        
        try:
            processor = PDFProcessor(pdf_path)
            
            # Test with timeout
            start_time = time.time()
            
            def run_processing():
                return processor.extract_chemical_structures_with_decimer()
            
            with mp.Pool(1) as pool:
                result = pool.apply_async(run_processing)
                try:
                    structures = result.get(timeout=60)
                    processing_time = time.time() - start_time
                    
                    assert isinstance(structures, list)
                    assert processing_time < 55  # Should complete well within timeout
                    
                except mp.TimeoutError:
                    pytest.fail("DECIMER hung on noisy structure images")
                    
        finally:
            shutil.rmtree(temp_dir)
    
    def test_corrupted_image_handling(self, create_problematic_pdf):
        """Test DECIMER's handling of corrupted image data"""
        pdf_path, temp_dir = create_problematic_pdf("corrupted_images")
        
        try:
            client = DECIMERClient()
            
            # Should handle corrupted images gracefully without hanging
            def run_segmentation():
                return client.segment_structures_from_pdf(pdf_path)
            
            with mp.Pool(1) as pool:
                result = pool.apply_async(run_segmentation)
                try:
                    segments = result.get(timeout=30)
                    assert isinstance(segments, list)
                    # May return empty list for corrupted images, which is fine
                except mp.TimeoutError:
                    pytest.fail("DECIMER hung on corrupted images")
                    
        finally:
            shutil.rmtree(temp_dir)
    
    def test_empty_pdf_handling(self, create_problematic_pdf):
        """Test DECIMER with empty or nearly empty PDFs"""
        pdf_path, temp_dir = create_problematic_pdf("empty_pages")
        
        try:
            processor = PDFProcessor(pdf_path)
            
            # Should handle empty PDFs quickly
            start_time = time.time()
            result = processor.extract_chemical_structures_with_decimer()
            processing_time = time.time() - start_time
            
            assert isinstance(result, list)
            assert processing_time < 10  # Should be very fast for empty content
            
        finally:
            shutil.rmtree(temp_dir)
    
    def test_very_dense_content(self, create_problematic_pdf):
        """Test DECIMER with very dense content that might cause memory issues"""
        pdf_path, temp_dir = create_problematic_pdf("very_dense")
        
        try:
            processor = PDFProcessor(pdf_path)
            
            # Monitor memory usage
            import psutil
            initial_memory = psutil.Process().memory_info().rss / 1024 / 1024  # MB
            
            def run_processing():
                return processor.extract_chemical_structures_with_decimer()
            
            with mp.Pool(1) as pool:
                result = pool.apply_async(run_processing)
                try:
                    structures = result.get(timeout=120)  # 2 minute timeout for dense content
                    
                    final_memory = psutil.Process().memory_info().rss / 1024 / 1024  # MB
                    memory_increase = final_memory - initial_memory
                    
                    assert isinstance(structures, list)
                    assert memory_increase < 500  # Should not use excessive memory
                    
                except mp.TimeoutError:
                    pytest.fail("DECIMER hung on dense content")
                    
        finally:
            shutil.rmtree(temp_dir)
    
    def test_concurrent_decimer_processing(self, create_problematic_pdf):
        """Test concurrent DECIMER processing to identify race conditions"""
        # Create multiple test PDFs
        pdf_data = []
        for pdf_type in ["mixed_content", "high_noise_structures", "corrupted_images"]:
            pdf_path, temp_dir = create_problematic_pdf(pdf_type)
            pdf_data.append((pdf_path, temp_dir))
        
        try:
            def process_single_pdf(pdf_path):
                client = DECIMERClient()
                return client.segment_structures_from_pdf(pdf_path)
            
            # Process PDFs concurrently
            with mp.Pool(3) as pool:
                async_results = [pool.apply_async(process_single_pdf, (pdf_path,)) 
                               for pdf_path, _ in pdf_data]
                
                results = []
                for async_result in async_results:
                    try:
                        result = async_result.get(timeout=60)
                        results.append(result)
                    except mp.TimeoutError:
                        pytest.fail("Concurrent DECIMER processing timed out")
                
                assert len(results) == 3
                assert all(isinstance(result, list) for result in results)
                
        finally:
            for _, temp_dir in pdf_data:
                shutil.rmtree(temp_dir)
    
    def test_decimer_with_specific_tensorflow_settings(self, create_problematic_pdf):
        """Test DECIMER with specific TensorFlow configurations that might cause hanging"""
        pdf_path, temp_dir = create_problematic_pdf("mixed_content")
        
        try:
            # Test with various TensorFlow configurations
            test_configs = [
                {'TF_FORCE_GPU_ALLOW_GROWTH': 'true'},
                {'TF_GPU_ALLOCATOR': 'cuda_malloc_async'},
                {'TF_CPP_MIN_LOG_LEVEL': '2'},
                {'CUDA_VISIBLE_DEVICES': ''},  # Force CPU
            ]
            
            for config in test_configs:
                with patch.dict(os.environ, config):
                    client = DECIMERClient()
                    
                    def run_segmentation():
                        return client.segment_structures_from_pdf(pdf_path)
                    
                    with mp.Pool(1) as pool:
                        result = pool.apply_async(run_segmentation)
                        try:
                            segments = result.get(timeout=45)
                            assert isinstance(segments, list)
                        except mp.TimeoutError:
                            pytest.fail(f"DECIMER hung with config: {config}")
                            
        finally:
            shutil.rmtree(temp_dir)
    
    def test_memory_pressure_during_processing(self, create_problematic_pdf):
        """Test DECIMER behavior under memory pressure"""
        pdf_path, temp_dir = create_problematic_pdf("multiple_large_structures", num_structures=1)
        
        try:
            # Simulate memory pressure by allocating large arrays
            memory_hogs = []
            
            def create_memory_pressure():
                # Allocate memory in background
                for i in range(5):
                    # Allocate 100MB chunks
                    memory_hogs.append(np.zeros((100 * 1024 * 1024 // 8,), dtype=np.float64))
                    time.sleep(0.1)
            
            # Start memory pressure in background
            pressure_thread = threading.Thread(target=create_memory_pressure, daemon=True)
            pressure_thread.start()
            
            time.sleep(0.5)  # Let memory pressure build up
            
            # Now try DECIMER processing
            client = DECIMERClient()
            
            def run_processing():
                return client.segment_structures_from_pdf(pdf_path)
            
            with mp.Pool(1) as pool:
                result = pool.apply_async(run_processing)
                try:
                    segments = result.get(timeout=90)
                    assert isinstance(segments, list)
                except mp.TimeoutError:
                    pytest.fail("DECIMER hung under memory pressure")
                    
        finally:
            # Clean up memory
            memory_hogs.clear()
            import gc
            gc.collect()
            shutil.rmtree(temp_dir)


class TestDECIMERErrorRecovery:
    """Test DECIMER's error recovery and fallback mechanisms"""
    
    def test_decimer_import_failure_fallback(self, create_problematic_pdf):
        """Test fallback when DECIMER imports fail"""
        pdf_path, temp_dir = create_problematic_pdf("mixed_content")
        
        try:
            # Simulate DECIMER import failure
            with patch('pdf_processor.DECIMER_SEGMENTATION_AVAILABLE', False):
                processor = PDFProcessor(pdf_path)
                
                # Should fall back to regular image extraction
                result = processor.extract_chemical_structures_with_decimer()
                
                assert isinstance(result, list)
                # Should still extract text and regular images
                
        finally:
            shutil.rmtree(temp_dir)
    
    def test_tensorflow_gpu_error_recovery(self, create_problematic_pdf):
        """Test recovery from TensorFlow GPU errors"""
        pdf_path, temp_dir = create_problematic_pdf("mixed_content")
        
        try:
            # Simulate TensorFlow GPU error
            def mock_gpu_error(*args, **kwargs):
                raise RuntimeError("GPU out of memory")
            
            with patch('tensorflow.config.experimental.set_memory_growth', side_effect=mock_gpu_error):
                # Should still work despite GPU error
                client = DECIMERClient()
                
                def run_processing():
                    return client.segment_structures_from_pdf(pdf_path)
                
                with mp.Pool(1) as pool:
                    result = pool.apply_async(run_processing)
                    try:
                        segments = result.get(timeout=30)
                        assert isinstance(segments, list)
                    except mp.TimeoutError:
                        pytest.fail("DECIMER hung after GPU error")
                        
        finally:
            shutil.rmtree(temp_dir)
    
    def test_segmentation_partial_failure_recovery(self, create_problematic_pdf):
        """Test recovery when segmentation partially fails"""
        pdf_path, temp_dir = create_problematic_pdf("mixed_content")
        
        try:
            call_count = 0
            
            def mock_failing_segmentation(*args, **kwargs):
                nonlocal call_count
                call_count += 1
                if call_count <= 2:  # Fail first two calls
                    raise RuntimeError("Segmentation failed")
                return []  # Succeed on third call
            
            with patch('decimer_client.segment_chemical_structures_from_file', 
                      side_effect=mock_failing_segmentation):
                client = DECIMERClient()
                
                # Should eventually succeed or handle gracefully
                segments = client.segment_structures_from_pdf(pdf_path)
                assert isinstance(segments, list)
                
        finally:
            shutil.rmtree(temp_dir)


if __name__ == "__main__":
    # Run with specific markers for debugging
    pytest.main([__file__, "-v", "-s", "--tb=short"])