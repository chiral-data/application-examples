import pytest
import os
from unittest.mock import Mock, patch, MagicMock
import fitz
from PIL import Image

from pdf_processor import PDFProcessor


class TestPDFProcessor:
    """Test suite for PDFProcessor class"""
    
    def test_init_success(self, mock_pdf_path):
        """Test successful initialization"""
        with patch('fitz.open') as mock_fitz:
            mock_doc = Mock()
            mock_doc.is_encrypted = False
            mock_fitz.return_value = mock_doc
            
            processor = PDFProcessor(mock_pdf_path)
            
            assert processor.pdf_path == mock_pdf_path
            assert processor.doc == mock_doc
            assert os.path.exists(processor.temp_dir)
    
    def test_init_encrypted_pdf(self, mock_pdf_path):
        """Test initialization with encrypted PDF"""
        with patch('fitz.open') as mock_fitz:
            mock_doc = Mock()
            mock_doc.is_encrypted = True
            mock_fitz.return_value = mock_doc
            
            with pytest.raises(ValueError, match="PDF is encrypted"):
                PDFProcessor(mock_pdf_path)
    
    def test_init_invalid_pdf(self):
        """Test initialization with invalid PDF path"""
        with patch('fitz.open', side_effect=Exception("Invalid PDF")):
            with pytest.raises(ValueError, match="Failed to open PDF file"):
                PDFProcessor("invalid_path.pdf")
    
    def test_extract_text_and_images_success(self, mock_pdf_path, mock_pdf_document):
        """Test successful text and image extraction"""
        with patch('fitz.open', return_value=mock_pdf_document):
            with patch('fitz.Pixmap') as mock_pixmap:
                with patch('PIL.Image.open') as mock_image_open:
                    # Setup mocks
                    mock_pix = Mock()
                    mock_pix.n = 4
                    mock_pix.alpha = 1
                    mock_pix.tobytes.return_value = b'fake_image_data_with_sufficient_length_for_testing'
                    mock_pixmap.return_value = mock_pix
                    
                    mock_img = Mock()
                    mock_img.width = 100
                    mock_img.height = 100
                    mock_img.save = Mock()
                    mock_image_open.return_value = mock_img
                    
                    # Mock page.get_image_bbox
                    mock_pdf_document[0].get_image_bbox.return_value = fitz.Rect(10, 10, 110, 110)
                    mock_pdf_document[1].get_image_bbox.return_value = fitz.Rect(20, 20, 140, 140)
                    
                    processor = PDFProcessor(mock_pdf_path)
                    results = processor.extract_text_and_images()
                    
                    assert len(results) == 2
                    assert results[0]['page'] == 0
                    assert results[1]['page'] == 1
                    assert 'text' in results[0]
                    assert 'images' in results[0]
                    assert 'page_size' in results[0]
    
    def test_extract_text_alternative(self, mock_pdf_path):
        """Test alternative text extraction method"""
        with patch('fitz.open') as mock_fitz:
            mock_doc = Mock()
            mock_doc.is_encrypted = False
            mock_fitz.return_value = mock_doc
            
            processor = PDFProcessor(mock_pdf_path)
            
            # Mock page with text dict
            mock_page = Mock()
            mock_page.get_text.return_value = {
                'blocks': [
                    {
                        'lines': [
                            {
                                'spans': [
                                    {'text': 'Hello '},
                                    {'text': 'World'}
                                ]
                            }
                        ]
                    }
                ]
            }
            
            result = processor._extract_text_alternative(mock_page)
            assert 'Hello World' in result
    
    def test_extract_images_from_page_small_images(self, mock_pdf_path):
        """Test that small images are filtered out"""
        with patch('fitz.open') as mock_fitz:
            mock_doc = Mock()
            mock_doc.is_encrypted = False
            mock_fitz.return_value = mock_doc
            
            processor = PDFProcessor(mock_pdf_path)
            
            mock_page = Mock()
            # Mock small image (width=30, height=30)
            mock_page.get_images.return_value = [(1, 0, 30, 30, 8, 'DeviceRGB', '', 'Im1', 'DCTDecode')]
            
            images = processor._extract_images_from_page(mock_page, 0)
            assert len(images) == 0  # Small image should be filtered out
    
    def test_cleanup_temp_files(self, mock_pdf_path, temp_dir):
        """Test cleanup of temporary files"""
        with patch('fitz.open') as mock_fitz:
            mock_doc = Mock()
            mock_doc.is_encrypted = False
            mock_fitz.return_value = mock_doc
            
            processor = PDFProcessor(mock_pdf_path)
            processor.temp_dir = temp_dir
            
            # Create some test files
            test_files = [
                os.path.join(temp_dir, 'page_0_img_1.png'),
                os.path.join(temp_dir, 'temp_img_test.png'),
                os.path.join(temp_dir, 'keep_this.txt')
            ]
            
            for file_path in test_files:
                with open(file_path, 'w') as f:
                    f.write('test')
            
            processor.cleanup_temp_files()
            
            # Check that image files are removed but other files remain
            assert not os.path.exists(test_files[0])
            assert not os.path.exists(test_files[1])
            assert os.path.exists(test_files[2])
    
    def test_get_document_info(self, mock_pdf_path, mock_pdf_document):
        """Test document metadata extraction"""
        with patch('fitz.open', return_value=mock_pdf_document):
            processor = PDFProcessor(mock_pdf_path)
            info = processor.get_document_info()
            
            assert info['title'] == 'Test Chemical Paper'
            assert info['author'] == 'Test Author'
            assert info['page_count'] == 2
            assert info['is_pdf'] is True
            assert info['is_encrypted'] is False
    
    def test_destructor(self, mock_pdf_path):
        """Test that document is properly closed on destruction"""
        mock_doc = Mock()
        mock_doc.is_encrypted = False
        mock_doc.close = Mock()
        
        with patch('fitz.open', return_value=mock_doc):
            processor = PDFProcessor(mock_pdf_path)
            
        # Call destructor manually
        processor.__del__()
        mock_doc.close.assert_called_once()
    
    def test_extract_text_and_images_error_handling(self, mock_pdf_path):
        """Test error handling during extraction"""
        with patch('fitz.open') as mock_fitz:
            mock_doc = Mock()
            mock_doc.is_encrypted = False
            mock_doc.__iter__ = Mock(side_effect=Exception("Processing error"))
            mock_fitz.return_value = mock_doc
            
            processor = PDFProcessor(mock_pdf_path)
            
            with pytest.raises(RuntimeError, match="Failed to process PDF"):
                processor.extract_text_and_images()
    
    @patch('builtins.print')  # Mock print to suppress warning messages
    def test_extract_images_error_handling(self, mock_print, mock_pdf_path):
        """Test error handling during image extraction"""
        with patch('fitz.open') as mock_fitz:
            mock_doc = Mock()
            mock_doc.is_encrypted = False
            mock_fitz.return_value = mock_doc
            
            processor = PDFProcessor(mock_pdf_path)
            
            mock_page = Mock()
            mock_page.get_images.side_effect = Exception("Image extraction error")
            
            images = processor._extract_images_from_page(mock_page, 0)
            assert len(images) == 0
            mock_print.assert_called()
    
    def test_image_validation_and_filtering(self, mock_pdf_path):
        """Test image size and format validation"""
        with patch('fitz.open') as mock_fitz:
            mock_doc = Mock()
            mock_doc.is_encrypted = False
            mock_fitz.return_value = mock_doc
            
            processor = PDFProcessor(mock_pdf_path)
            
            mock_page = Mock()
            mock_page.get_images.return_value = [
                (1, 0, 100, 100, 8, 'DeviceRGB', '', 'Im1', 'DCTDecode'),  # Valid image
                (2, 0, 200, 200, 8, 'DeviceCMYK', '', 'Im2', 'DCTDecode')   # Invalid color space
            ]
            
            with patch('fitz.Pixmap') as mock_pixmap:
                # First image: valid RGB
                mock_pix1 = Mock()
                mock_pix1.n = 4  # RGB + alpha
                mock_pix1.alpha = 1
                mock_pix1.tobytes.return_value = b'valid_image_data_with_sufficient_length'
                
                # Second image: invalid color space
                mock_pix2 = Mock()
                mock_pix2.n = 5  # CMYK + alpha
                mock_pix2.alpha = 1
                
                mock_pixmap.side_effect = [mock_pix1, mock_pix2]
                
                with patch('PIL.Image.open') as mock_image_open:
                    mock_img = Mock()
                    mock_img.width = 100
                    mock_img.height = 100
                    mock_img.save = Mock()
                    mock_image_open.return_value = mock_img
                    
                    mock_page.get_image_bbox.return_value = fitz.Rect(0, 0, 100, 100)
                    
                    images = processor._extract_images_from_page(mock_page, 0)
                    
                    # Only the valid RGB image should be processed
                    assert len(images) == 1
                    assert images[0]['width'] == 100
                    assert images[0]['height'] == 100