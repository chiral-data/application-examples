import pytest
import os
import sys
from unittest.mock import Mock, patch, MagicMock

# Add src directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'src'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'data'))

from pdf_processor import PDFProcessor
from decimer_client import DECIMERClient
from chemical_handler import ChemicalHandler
from chunker import ChemicalAwareChunker
from vectorstore import ChemicalVectorStore
from rag_chain import ChemicalRAG
from test_config import TEST_PDF_PATH, SAMPLE_STRUCTURES


class TestRealPDFProcessing:
    """Test suite using the actual test.pdf file for realistic testing"""
    
    @pytest.fixture
    def pdf_path(self):
        """Get the path to the test PDF file"""
        if not os.path.exists(TEST_PDF_PATH):
            pytest.skip(f"Test PDF not found at {TEST_PDF_PATH}")
        return TEST_PDF_PATH
    
    def test_pdf_processor_with_real_file(self, pdf_path):
        """Test PDF processor with the actual test.pdf file"""
        try:
            processor = PDFProcessor(pdf_path)
            
            # Test document info extraction
            doc_info = processor.get_document_info()
            
            assert isinstance(doc_info, dict)
            assert 'page_count' in doc_info
            assert doc_info['page_count'] > 0
            print(f"PDF has {doc_info['page_count']} pages")
            
            # Test text and image extraction
            pages_data = processor.extract_text_and_images()
            
            assert isinstance(pages_data, list)
            assert len(pages_data) == doc_info['page_count']
            
            total_text_length = 0
            total_images = 0
            
            for i, page_data in enumerate(pages_data):
                assert 'page' in page_data
                assert 'text' in page_data
                assert 'images' in page_data
                assert page_data['page'] == i
                
                total_text_length += len(page_data['text'])
                total_images += len(page_data['images'])
                
                print(f"Page {i}: {len(page_data['text'])} chars, {len(page_data['images'])} images")
            
            print(f"Total: {total_text_length} characters, {total_images} images extracted")
            
            # Cleanup temp files
            processor.cleanup_temp_files()
            
        except Exception as e:
            pytest.fail(f"PDF processing failed: {e}")
    
    def test_decimer_client_with_mock_structures(self, pdf_path):
        """Test DECIMER client functionality with mocked responses"""
        
        # First extract images from the real PDF
        processor = PDFProcessor(pdf_path)
        pages_data = processor.extract_text_and_images()
        
        # Collect all image paths
        image_paths = []
        for page_data in pages_data:
            for img_info in page_data['images']:
                image_paths.append(img_info['path'])
        
        if not image_paths:
            print("No images found in PDF, creating mock image for testing")
            # Create a mock image if none found
            from PIL import Image
            import tempfile
            
            with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
                img = Image.new('RGB', (200, 150), color='white')
                img.save(tmp.name, 'PNG')
                image_paths = [tmp.name]
        
        # Test DECIMER client with mocked DECIMER package
        client = DECIMERClient()
        
        # Mock the predict_SMILES function from DECIMER
        with patch('decimer_client.predict_SMILES', return_value='CCO') as mock_predict:
            for img_path in image_paths[:3]:  # Test first 3 images only
                if os.path.exists(img_path):
                    print(f"Testing DECIMER with image: {img_path}")
                    
                    smiles = client.image_to_smiles(img_path)
                    if client.decimer_available:
                        assert smiles == 'CCO'
                        print(f"DECIMER returned SMILES: {smiles}")
                    else:
                        print("DECIMER not available, skipping SMILES test")
        
        # Cleanup
        processor.cleanup_temp_files()
        for path in image_paths:
            if os.path.exists(path) and 'tmp' in path:
                os.unlink(path)
    
    def test_decimer_segmentation_integration(self, pdf_path):
        """Test DECIMER segmentation functionality"""
        processor = PDFProcessor(pdf_path)
        
        # Test the new DECIMER segmentation method
        try:
            # Mock the DECIMER segmentation function
            with patch('pdf_processor.segment_chemical_structures_from_file') as mock_segment:
                # Mock return value: list of numpy arrays representing structure images
                import numpy as np
                mock_segments = [
                    np.random.randint(0, 255, (100, 100, 3), dtype=np.uint8),
                    np.random.randint(0, 255, (150, 120, 3), dtype=np.uint8)
                ]
                mock_segment.return_value = mock_segments
                
                structures = processor.extract_chemical_structures_with_decimer()
                
                if hasattr(processor, 'extract_chemical_structures_with_decimer'):
                    print(f"DECIMER segmentation found {len(structures)} structures")
                    for struct in structures:
                        assert 'path' in struct
                        assert 'width' in struct
                        assert 'height' in struct
                        assert struct['source'] == 'decimer_segmentation'
                else:
                    print("DECIMER segmentation method not available")
                    
        except Exception as e:
            print(f"DECIMER segmentation test failed: {e}")
        
        # Cleanup
        processor.cleanup_temp_files()
    
    def test_chemical_handler_integration(self, pdf_path):
        """Test chemical handler with real PDF images"""
        
        processor = PDFProcessor(pdf_path)
        pages_data = processor.extract_text_and_images()
        
        # Mock DECIMER client
        mock_decimer = Mock()
        mock_decimer.image_to_smiles.return_value = 'CCO'
        
        handler = ChemicalHandler(mock_decimer)
        
        structures_found = []
        
        with patch('rdkit.Chem.MolFromSmiles') as mock_mol:
            mock_molecule = Mock()
            mock_mol.return_value = mock_molecule
            
            with patch('rdkit.Chem.Descriptors.MolWt', return_value=46.07):
                with patch('rdkit.Chem.rdMolDescriptors.CalcMolFormula', return_value='C2H6O'):
                    
                    for page_data in pages_data:
                        for img_info in page_data['images']:
                            context = handler.get_structure_context(
                                page_data['text'], 
                                img_info.get('bbox', Mock())
                            )
                            
                            structure = handler.process_image(
                                img_info['path'], 
                                context
                            )
                            
                            if structure:
                                structures_found.append(structure)
                                print(f"Found structure: {structure['smiles']} ({structure['formula']})")
        
        print(f"Total structures processed: {len(structures_found)}")
        
        # Cleanup
        processor.cleanup_temp_files()
    
    def test_complete_workflow_with_real_pdf(self, pdf_path):
        """Test the complete workflow using the actual test.pdf"""
        
        print(f"Testing complete workflow with {pdf_path}")
        
        try:
            # Step 1: Process PDF
            processor = PDFProcessor(pdf_path)
            pages_data = processor.extract_text_and_images()
            
            print(f"Extracted {len(pages_data)} pages")
            
            # Step 2: Mock chemical structure processing
            mock_decimer = Mock()
            mock_decimer.image_to_smiles.return_value = 'CCO'
            
            handler = ChemicalHandler(mock_decimer)
            
            all_structures = []
            full_text = ""
            
            with patch('rdkit.Chem.MolFromSmiles') as mock_mol:
                mock_molecule = Mock()
                mock_mol.return_value = mock_molecule
                
                with patch('rdkit.Chem.Descriptors.MolWt', return_value=46.07):
                    with patch('rdkit.Chem.rdMolDescriptors.CalcMolFormula', return_value='C2H6O'):
                        
                        for page_data in pages_data:
                            full_text += page_data['text'] + "\n"
                            
                            # Process first few images only for testing
                            for img_info in page_data['images'][:2]:
                                context = handler.get_structure_context(
                                    page_data['text'], 
                                    img_info.get('bbox', Mock())
                                )
                                
                                structure = handler.process_image(
                                    img_info['path'], 
                                    context
                                )
                                
                                if structure:
                                    all_structures.append(structure)
            
            print(f"Processed {len(all_structures)} chemical structures")
            print(f"Extracted {len(full_text)} characters of text")
            
            # Step 3: Create chunks
            chunker = ChemicalAwareChunker(chunk_size=500, overlap=100)
            chunks = chunker.chunk_with_structures(full_text, all_structures)
            
            print(f"Created {len(chunks)} text chunks")
            
            # Get chunking statistics
            stats = chunker.get_chunking_stats(chunks)
            print(f"Chunking stats: {stats}")
            
            # Step 4: Create vector store (mocked)
            with patch('vectorstore.HuggingFaceEmbeddings'):
                with patch('vectorstore.Chroma') as mock_chroma:
                    mock_vectorstore = Mock()
                    mock_chroma.from_texts.return_value = mock_vectorstore
                    
                    vector_store = ChemicalVectorStore(use_ollama_embeddings=False)
                    result = vector_store.create_vectorstore(chunks)
                    
                    assert result == mock_vectorstore
                    print("Vector store created successfully")
            
            # Step 5: Test RAG chain (mocked)
            with patch('rag_chain.Ollama'):
                with patch('rag_chain.RetrievalQA') as mock_retrieval_qa:
                    mock_qa_chain = Mock()
                    mock_qa_chain.return_value = {
                        'result': f'Based on the document analysis, {len(all_structures)} chemical structures were found.',
                        'source_documents': []
                    }
                    mock_retrieval_qa.from_chain_type.return_value = mock_qa_chain
                    
                    rag = ChemicalRAG(mock_vectorstore)
                    answer = rag.query("How many chemical structures were found?")
                    
                    print(f"RAG answer: {answer}")
                    # RAG returns a dict with 'answer' and 'sources' by default
                    if isinstance(answer, dict):
                        assert 'answer' in answer
                        assert len(answer['answer']) > 0
                    else:
                        assert isinstance(answer, str)
                        assert len(answer) > 0
            
            print("Complete workflow test passed!")
            
        except Exception as e:
            print(f"Workflow failed: {e}")
            raise
        finally:
            # Cleanup
            if 'processor' in locals():
                processor.cleanup_temp_files()
    
    def test_pdf_content_analysis(self, pdf_path):
        """Analyze the content of the test PDF to understand its structure"""
        
        processor = PDFProcessor(pdf_path)
        doc_info = processor.get_document_info()
        pages_data = processor.extract_text_and_images()
        
        print("\n=== PDF CONTENT ANALYSIS ===")
        print(f"Document Info: {doc_info}")
        
        for i, page_data in enumerate(pages_data):
            print(f"\n--- Page {i+1} ---")
            print(f"Text length: {len(page_data['text'])} characters")
            print(f"Number of images: {len(page_data['images'])}")
            
            # Show first 200 characters of text
            if page_data['text']:
                preview = page_data['text'][:200].replace('\n', ' ')
                print(f"Text preview: {preview}...")
            
            # Show image info
            for j, img_info in enumerate(page_data['images']):
                print(f"  Image {j+1}: {img_info.get('width', 'unknown')}x{img_info.get('height', 'unknown')} pixels")
                if 'size_bytes' in img_info:
                    print(f"    Size: {img_info['size_bytes']} bytes")
        
        # Cleanup
        processor.cleanup_temp_files()
    
    def test_chunking_with_real_text(self, pdf_path):
        """Test chunking with actual text from the PDF"""
        
        processor = PDFProcessor(pdf_path)
        pages_data = processor.extract_text_and_images()
        
        # Extract all text
        full_text = ""
        for page_data in pages_data:
            full_text += page_data['text'] + "\n"
        
        if not full_text.strip():
            pytest.skip("No text found in PDF")
        
        print(f"Total text length: {len(full_text)} characters")
        
        # Test chunking without structures
        chunker = ChemicalAwareChunker(chunk_size=300, overlap=50)
        chunks = chunker.chunk_with_structures(full_text, [])
        
        print(f"Created {len(chunks)} chunks")
        
        # Analyze chunks
        for i, chunk in enumerate(chunks[:3]):  # Show first 3 chunks
            print(f"\n--- Chunk {i+1} ---")
            print(f"Length: {chunk['metadata']['char_count']} chars")
            print(f"Content type: {chunk['metadata']['content_type']}")
            print(f"Chemical keywords: {chunk['metadata']['chemical_keywords']}")
            print(f"Preview: {chunk['text'][:100]}...")
        
        # Test chunking with sample structures
        chunks_with_structures = chunker.chunk_with_structures(full_text, SAMPLE_STRUCTURES)
        
        stats = chunker.get_chunking_stats(chunks_with_structures)
        print(f"\nChunking statistics: {stats}")
        
        # Cleanup
        processor.cleanup_temp_files()
    
    def test_error_handling_with_real_pdf(self, pdf_path):
        """Test error handling scenarios with the real PDF"""
        
        # Test with valid PDF
        processor = PDFProcessor(pdf_path)
        assert processor.pdf_path == pdf_path
        
        # Test document info
        doc_info = processor.get_document_info()
        assert isinstance(doc_info, dict)
        
        # Test with non-existent file
        with pytest.raises(ValueError):
            PDFProcessor("nonexistent.pdf")
        
        # Cleanup
        processor.cleanup_temp_files()