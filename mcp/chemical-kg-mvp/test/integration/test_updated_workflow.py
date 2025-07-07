#!/usr/bin/env python3
"""
Test the updated image-based DECIMER workflow
"""

import sys
import os

sys.path.insert(0, 'src')

def test_updated_decimer_client():
    """Test the updated DECIMERClient with image-based approach"""
    print("Testing Updated DECIMER Client")
    print("=" * 40)
    
    from decimer_client import DECIMERClient
    
    client = DECIMERClient()
    
    print(f"DECIMER available: {client.decimer_available}")
    print(f"Segmentation available: {client.segmentation_available}")
    
    if not (client.decimer_available and client.segmentation_available):
        print("ERROR: Required packages not available")
        return False
    
    # Test the new process_pdf_complete method
    pdf_path = "data/test.pdf"
    
    if not os.path.exists(pdf_path):
        print(f"ERROR: Test PDF not found: {pdf_path}")
        return False
    
    print(f"\nTesting complete PDF processing: {pdf_path}")
    
    try:
        results = client.process_pdf_complete(pdf_path, "test_workflow_output")
        
        print(f"\\nResults:")
        print(f"  Total structures: {results['total_structures']}")
        print(f"  Processing time: {results['processing_time']:.2f}s")
        
        valid_count = 0
        for structure in results['structures']:
            if structure.get('valid', False):
                valid_count += 1
                print(f"  Structure {structure['index']}: {structure['smiles']}")
        
        print(f"\\nValid SMILES generated: {valid_count}/{results['total_structures']}")
        
        return True
        
    except Exception as e:
        print(f"ERROR: {e}")
        return False

def test_pdf_processor():
    """Test the updated PDFProcessor"""
    print("\\nTesting Updated PDF Processor")
    print("=" * 40)
    
    from pdf_processor import PDFProcessor
    
    pdf_path = "data/test.pdf"
    
    try:
        processor = PDFProcessor(pdf_path)
        
        print("Testing image-based DECIMER extraction...")
        structures = processor.extract_chemical_structures_with_decimer()
        
        print(f"Found {len(structures)} structures:")
        for structure in structures[:3]:  # Show first 3
            print(f"  {structure['path']} - Page {structure.get('page', 'N/A')} - {structure['width']}x{structure['height']}")
        
        return True
        
    except Exception as e:
        print(f"ERROR: {e}")
        return False

if __name__ == "__main__":
    print("Testing Updated DECIMER Workflow")
    print("=" * 50)
    
    # Test updated DECIMER client
    client_ok = test_updated_decimer_client()
    
    # Test updated PDF processor
    processor_ok = test_pdf_processor()
    
    print("\\n" + "=" * 50)
    if client_ok and processor_ok:
        print("SUCCESS: All tests passed!")
    else:
        print("FAILED: Some tests failed")
    
    print("Test completed!")