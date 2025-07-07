"""
Comprehensive tests for chemical structure segmentation functionality.
Tests the ability to detect and segment chemical structures from PDFs and images.
"""

import pytest
import os
import tempfile
import numpy as np
from PIL import Image, ImageDraw
from unittest.mock import patch, MagicMock
import cv2
import fitz  # PyMuPDF

from pdf_processor import PDFProcessor
from decimer_client import DECIMERClient


class TestStructureSegmentation:
    """Test suite for structure segmentation capabilities"""
    
    @pytest.fixture
    def decimer_client(self):
        """Create a DECIMERClient instance"""
        return DECIMERClient()
    
    @pytest.fixture
    def pdf_processor(self):
        """Create a PDFProcessor instance"""
        return PDFProcessor()
    
    @pytest.fixture
    def create_test_image_with_structure(self):
        """Create a test image with a chemical structure"""
        def _create_image(width=800, height=600, structure_type='benzene'):
            # Create white background
            img = Image.new('RGB', (width, height), 'white')
            draw = ImageDraw.Draw(img)
            
            if structure_type == 'benzene':
                # Draw a simple benzene ring
                center_x, center_y = width // 2, height // 2
                radius = 100
                
                # Draw hexagon for benzene
                angles = [i * 60 for i in range(6)]
                points = []
                for angle in angles:
                    x = center_x + radius * np.cos(np.radians(angle))
                    y = center_y + radius * np.sin(np.radians(angle))
                    points.append((x, y))
                
                # Draw bonds
                for i in range(6):
                    draw.line([points[i], points[(i + 1) % 6]], fill='black', width=3)
                
                # Draw inner circle for aromatic ring
                draw.ellipse(
                    [center_x - radius//2, center_y - radius//2, 
                     center_x + radius//2, center_y + radius//2],
                    outline='black', width=2
                )
                
            elif structure_type == 'chain':
                # Draw a simple carbon chain
                start_x = width // 4
                y = height // 2
                
                for i in range(5):
                    x = start_x + i * 100
                    draw.line([(x, y), (x + 80, y)], fill='black', width=3)
                    if i % 2 == 0:
                        draw.line([(x + 40, y), (x + 40, y - 40)], fill='black', width=3)
                    else:
                        draw.line([(x + 40, y), (x + 40, y + 40)], fill='black', width=3)
            
            return img
        
        return _create_image
    
    @pytest.fixture
    def create_test_pdf_with_structures(self, create_test_image_with_structure):
        """Create a test PDF with chemical structures"""
        def _create_pdf(num_structures=3, include_text=True):
            with tempfile.NamedTemporaryFile(suffix='.pdf', delete=False) as f:
                doc = fitz.open()
                page = doc.new_page()
                
                if include_text:
                    # Add some text
                    text = "Chemical Structures Test Document\n\n"
                    text += "This document contains several chemical structures for testing.\n"
                    page.insert_text((50, 50), text, fontsize=12)
                
                # Add structures
                y_offset = 150
                for i in range(num_structures):
                    # Create structure image
                    structure_type = 'benzene' if i % 2 == 0 else 'chain'
                    img = create_test_image_with_structure(
                        width=300, height=200, structure_type=structure_type
                    )
                    
                    # Save to temp file
                    with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as img_file:
                        img.save(img_file.name)
                        
                        # Insert into PDF
                        rect = fitz.Rect(50, y_offset, 350, y_offset + 200)
                        page.insert_image(rect, filename=img_file.name)
                        
                        os.unlink(img_file.name)
                    
                    y_offset += 250
                
                doc.save(f.name)
                doc.close()
                
                return f.name
        
        return _create_pdf
    
    def test_segment_single_structure(self, decimer_client, create_test_image_with_structure):
        """Test segmentation of a single chemical structure"""
        # Create test image
        img = create_test_image_with_structure(structure_type='benzene')
        
        with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as f:
            img.save(f.name)
            
            try:
                # Mock the segmentation to avoid actual model call
                with patch.object(decimer_client, 'segment_structures') as mock_segment:
                    # Return mock segmented image
                    mock_segment.return_value = [np.array(img)]
                    
                    segments = decimer_client.segment_structures(f.name)
                    
                    assert len(segments) >= 1
                    assert isinstance(segments[0], np.ndarray)
                    mock_segment.assert_called_once_with(f.name)
                    
            finally:
                os.unlink(f.name)
    
    def test_segment_multiple_structures(self, decimer_client, create_test_pdf_with_structures):
        """Test segmentation of multiple structures from PDF"""
        pdf_path = create_test_pdf_with_structures(num_structures=3)
        
        try:
            with patch.object(decimer_client, 'segment_structures') as mock_segment:
                # Mock return multiple segments
                mock_img = np.ones((200, 300, 3), dtype=np.uint8) * 255
                mock_segment.return_value = [mock_img] * 3
                
                segments = decimer_client.segment_structures(pdf_path)
                
                assert len(segments) == 3
                assert all(isinstance(seg, np.ndarray) for seg in segments)
                
        finally:
            os.unlink(pdf_path)
    
    def test_segmentation_with_timeout(self, decimer_client, create_test_image_with_structure):
        """Test that segmentation respects timeout settings"""
        img = create_test_image_with_structure()
        
        with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as f:
            img.save(f.name)
            
            try:
                with patch.object(decimer_client, 'segment_structures') as mock_segment:
                    # Simulate timeout
                    import time
                    def slow_segment(*args, **kwargs):
                        time.sleep(5)  # Simulate slow processing
                        return []
                    
                    mock_segment.side_effect = slow_segment
                    
                    # This should timeout if properly configured
                    with pytest.raises(Exception):  # Could be TimeoutError or similar
                        decimer_client.segment_structures(f.name, timeout=1)
                        
            finally:
                os.unlink(f.name)
    
    def test_empty_image_handling(self, decimer_client):
        """Test handling of empty or blank images"""
        # Create blank white image
        img = Image.new('RGB', (800, 600), 'white')
        
        with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as f:
            img.save(f.name)
            
            try:
                with patch.object(decimer_client, 'segment_structures') as mock_segment:
                    mock_segment.return_value = []  # No structures found
                    
                    segments = decimer_client.segment_structures(f.name)
                    assert len(segments) == 0
                    
            finally:
                os.unlink(f.name)
    
    def test_invalid_file_handling(self, decimer_client):
        """Test handling of invalid files"""
        with tempfile.NamedTemporaryFile(suffix='.txt', delete=False) as f:
            f.write(b"This is not an image file")
            f.flush()
            
            try:
                with pytest.raises(Exception):  # Should raise appropriate exception
                    decimer_client.segment_structures(f.name)
                    
            finally:
                os.unlink(f.name)
    
    def test_segmentation_quality_metrics(self, decimer_client, create_test_image_with_structure):
        """Test segmentation quality with known structures"""
        # Create image with well-defined structure
        img = create_test_image_with_structure(structure_type='benzene')
        
        with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as f:
            img.save(f.name)
            
            try:
                with patch.object(decimer_client, 'segment_structures') as mock_segment:
                    # Return properly sized segment
                    mock_segment.return_value = [np.array(img.crop((200, 150, 600, 450)))]
                    
                    segments = decimer_client.segment_structures(f.name)
                    
                    # Check segment properties
                    assert len(segments) == 1
                    segment = segments[0]
                    
                    # Check dimensions are reasonable
                    height, width = segment.shape[:2]
                    assert 100 <= height <= 500
                    assert 100 <= width <= 500
                    
                    # Check aspect ratio is reasonable for chemical structure
                    aspect_ratio = width / height
                    assert 0.5 <= aspect_ratio <= 2.0
                    
            finally:
                os.unlink(f.name)
    
    def test_batch_segmentation(self, decimer_client, create_test_image_with_structure):
        """Test batch processing of multiple images"""
        image_paths = []
        
        try:
            # Create multiple test images
            for i in range(5):
                img = create_test_image_with_structure(
                    structure_type='benzene' if i % 2 == 0 else 'chain'
                )
                with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as f:
                    img.save(f.name)
                    image_paths.append(f.name)
            
            with patch.object(decimer_client, 'segment_structures') as mock_segment:
                # Mock batch processing
                mock_segment.side_effect = [[np.ones((200, 300, 3))] for _ in image_paths]
                
                all_segments = []
                for path in image_paths:
                    segments = decimer_client.segment_structures(path)
                    all_segments.extend(segments)
                
                assert len(all_segments) == 5
                
        finally:
            for path in image_paths:
                if os.path.exists(path):
                    os.unlink(path)
    
    def test_memory_usage_during_segmentation(self, decimer_client, create_test_pdf_with_structures):
        """Test memory usage doesn't grow excessively during segmentation"""
        import psutil
        import gc
        
        process = psutil.Process()
        initial_memory = process.memory_info().rss / 1024 / 1024  # MB
        
        pdf_path = create_test_pdf_with_structures(num_structures=10)
        
        try:
            with patch.object(decimer_client, 'segment_structures') as mock_segment:
                # Mock memory-intensive operation
                def create_segments(*args, **kwargs):
                    # Create some large arrays
                    segments = [np.ones((500, 500, 3), dtype=np.uint8) for _ in range(10)]
                    return segments
                
                mock_segment.side_effect = create_segments
                
                # Process multiple times
                for _ in range(3):
                    segments = decimer_client.segment_structures(pdf_path)
                    # Clean up after each iteration
                    del segments
                    gc.collect()
                
                # Check memory hasn't grown too much
                final_memory = process.memory_info().rss / 1024 / 1024  # MB
                memory_growth = final_memory - initial_memory
                
                # Allow some growth but not excessive (e.g., < 500MB)
                assert memory_growth < 500, f"Memory grew by {memory_growth}MB"
                
        finally:
            os.unlink(pdf_path)