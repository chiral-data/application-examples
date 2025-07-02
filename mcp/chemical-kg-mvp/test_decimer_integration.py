#!/usr/bin/env python3
"""
Quick test script to verify DECIMER integration after Docker build
"""

import sys
import os

def test_imports():
    """Test if all DECIMER packages can be imported"""
    print("Testing DECIMER package imports...")
    
    # Test basic imports
    try:
        import numpy as np
        print("✓ numpy")
    except ImportError as e:
        print(f"✗ numpy: {e}")
        return False
    
    try:
        import cv2
        print("✓ opencv-python")
    except ImportError as e:
        print(f"✗ opencv-python: {e}")
        return False
    
    try:
        import tensorflow as tf
        print(f"✓ tensorflow {tf.__version__}")
    except ImportError as e:
        print(f"✗ tensorflow: {e}")
        return False
    
    # Test DECIMER imports
    try:
        from DECIMER import predict_SMILES
        print("✓ DECIMER (predict_SMILES)")
    except ImportError as e:
        print(f"✗ DECIMER: {e}")
        return False
    
    try:
        from decimer_segmentation import segment_chemical_structures_from_file
        print("✓ decimer-segmentation")
    except ImportError as e:
        print(f"✗ decimer-segmentation: {e}")
        return False
    
    return True

def test_decimer_functionality():
    """Test basic DECIMER functionality with a mock image"""
    print("\nTesting DECIMER functionality...")
    
    try:
        from PIL import Image
        import tempfile
        import os
        
        # Create a simple test image
        with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
            img = Image.new('RGB', (200, 150), color='white')
            img.save(tmp.name, 'PNG')
            test_image_path = tmp.name
        
        print(f"Created test image: {test_image_path}")
        
        # Test our DECIMERClient
        sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))
        from decimer_client import DECIMERClient
        
        client = DECIMERClient()
        print(f"DECIMER available: {client.decimer_available}")
        print(f"Segmentation available: {client.segmentation_available}")
        
        if client.decimer_available:
            # Test SMILES prediction (will likely fail with empty image but should not crash)
            try:
                smiles = client.image_to_smiles(test_image_path)
                print(f"SMILES prediction result: {smiles}")
                print("✓ DECIMER prediction completed without errors")
            except Exception as e:
                print(f"DECIMER prediction error (expected for blank image): {e}")
        
        # Cleanup
        os.unlink(test_image_path)
        
        return True
        
    except Exception as e:
        print(f"✗ DECIMER functionality test failed: {e}")
        return False

def test_pdf_processing():
    """Test PDF processing capabilities"""
    print("\nTesting PDF processing...")
    
    try:
        sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))
        from pdf_processor import PDFProcessor, DECIMER_SEGMENTATION_AVAILABLE
        
        print(f"DECIMER segmentation available: {DECIMER_SEGMENTATION_AVAILABLE}")
        
        # Test with the sample PDF if it exists
        test_pdf_path = os.path.join(os.path.dirname(__file__), 'data', 'test.pdf')
        if os.path.exists(test_pdf_path):
            print(f"Testing with {test_pdf_path}")
            try:
                processor = PDFProcessor(test_pdf_path)
                doc_info = processor.get_document_info()
                print(f"✓ PDF processing works - {doc_info['page_count']} pages")
                processor.cleanup_temp_files()
                return True
            except Exception as e:
                print(f"PDF processing error: {e}")
                return False
        else:
            print(f"No test PDF found at {test_pdf_path}")
            print("✓ PDF processor class imported successfully")
            return True
            
    except Exception as e:
        print(f"✗ PDF processing test failed: {e}")
        return False

def main():
    print("=" * 60)
    print("DECIMER Integration Test")
    print("=" * 60)
    
    success = True
    
    # Test imports
    if not test_imports():
        success = False
    
    # Test functionality
    if not test_decimer_functionality():
        success = False
    
    # Test PDF processing
    if not test_pdf_processing():
        success = False
    
    print("\n" + "=" * 60)
    if success:
        print("✓ ALL TESTS PASSED - DECIMER integration is working!")
        return 0
    else:
        print("✗ SOME TESTS FAILED - Check the errors above")
        return 1

if __name__ == '__main__':
    sys.exit(main())