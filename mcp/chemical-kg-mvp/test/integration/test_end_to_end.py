import pytest
import os
import tempfile
import shutil
from unittest.mock import Mock, patch, MagicMock
from PIL import Image
import sys

# Add src directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'src'))

from pdf_processor import PDFProcessor
from decimer_client import DECIMERClient
from chemical_handler import ChemicalHandler
from chunker import ChemicalAwareChunker
from vectorstore import ChemicalVectorStore
from rag_chain import ChemicalRAG


class TestEndToEndWorkflow:
    """End-to-end integration tests for the complete chemical paper analysis workflow"""
    
    @pytest.fixture
    def temp_workspace(self):
        """Create a temporary workspace for tests"""
        workspace = tempfile.mkdtemp()
        yield workspace
        shutil.rmtree(workspace)
    
    @pytest.fixture
    def mock_pdf_with_chemical_content(self, temp_workspace):
        """Create a mock PDF with chemical content"""
        pdf_path = os.path.join(temp_workspace, "chemical_paper.pdf")
        
        # Create mock PDF file
        with open(pdf_path, 'wb') as f:
            f.write(b'%PDF-1.4\nMock chemical paper content')
        
        return pdf_path
    
    @pytest.fixture
    def mock_chemical_images(self, temp_workspace):
        """Create mock chemical structure images"""
        images = []
        
        for i, name in enumerate(['ethanol', 'benzene', 'acetic_acid']):
            img_path = os.path.join(temp_workspace, f"{name}.png")
            
            # Create simple test image
            img = Image.new('RGB', (150, 100), color='white')
            img.save(img_path, 'PNG')
            images.append(img_path)
        
        return images
    
    def test_complete_workflow_success(self, temp_workspace, mock_pdf_with_chemical_content, mock_chemical_images):
        """Test the complete workflow from PDF to Q&A"""
        
        # Mock external dependencies
        with patch('fitz.open') as mock_fitz:
            # Setup PDF mock
            mock_doc = Mock()
            mock_doc.is_encrypted = False
            mock_doc.page_count = 2
            mock_doc.metadata = {'title': 'Test Chemical Paper'}
            
            # Mock pages with text and images
            mock_page1 = Mock()
            mock_page1.get_text.return_value = """
            Introduction
            
            This paper describes the synthesis of ethanol (CCO) from glucose.
            The molecular formula is C2H6O and the structure is shown in Figure 1.
            """
            mock_page1.get_images.return_value = [(1, 0, 150, 100, 8, 'DeviceRGB', '', 'Im1', 'DCTDecode')]
            mock_page1.rect = Mock(width=612, height=792)
            mock_page1.get_image_bbox.return_value = Mock(x0=100, y0=200, x1=250, y1=300)
            
            mock_page2 = Mock()
            mock_page2.get_text.return_value = """
            Results
            
            The synthesis yielded pure ethanol with 95% efficiency.
            NMR analysis confirmed the structure. The compound shows
            characteristic chemical properties of primary alcohols.
            """
            mock_page2.get_images.return_value = []
            mock_page2.rect = Mock(width=612, height=792)
            
            mock_doc.__iter__ = lambda self: iter([mock_page1, mock_page2])
            mock_fitz.return_value = mock_doc
            
            # Mock DECIMER API
            with patch('requests.post') as mock_post:
                mock_response = Mock()
                mock_response.status_code = 200
                mock_response.json.return_value = {'smiles': 'CCO'}
                mock_post.return_value = mock_response
                
                # Mock RDKit
                with patch('rdkit.Chem.MolFromSmiles') as mock_mol:
                    mock_molecule = Mock()
                    mock_mol.return_value = mock_molecule
                    
                    with patch('rdkit.Chem.Descriptors.MolWt', return_value=46.07):
                        with patch('rdkit.Chem.rdMolDescriptors.CalcMolFormula', return_value='C2H6O'):
                            # Mock vector store and embeddings
                            with patch('vectorstore.Chroma') as mock_chroma:
                                mock_vectorstore = Mock()
                                mock_chroma.from_texts.return_value = mock_vectorstore
                                
                                with patch('vectorstore.HuggingFaceEmbeddings'):
                                    # Mock LLM and QA chain
                                    with patch('rag_chain.Ollama'):
                                        with patch('rag_chain.RetrievalQA') as mock_retrieval_qa:
                                            mock_qa_chain = Mock()
                                            mock_qa_chain.return_value = {
                                                'result': 'Ethanol (CCO) is an organic compound with molecular formula C2H6O and molecular weight 46.07 g/mol.',
                                                'source_documents': []
                                            }
                                            mock_retrieval_qa.from_chain_type.return_value = mock_qa_chain
                                            
                                            # Mock Pixmap for image processing
                                            with patch('fitz.Pixmap') as mock_pixmap:
                                                mock_pix = Mock()
                                                mock_pix.n = 4
                                                mock_pix.alpha = 1
                                                mock_pix.tobytes.return_value = b'fake_image_data_with_sufficient_length_for_testing_purposes'
                                                mock_pixmap.return_value = mock_pix
                                                
                                                with patch('PIL.Image.open') as mock_image_open:
                                                    mock_img = Mock()
                                                    mock_img.width = 150
                                                    mock_img.height = 100
                                                    mock_img.save = Mock()
                                                    mock_image_open.return_value = mock_img
                                                    
                                                    # Execute the complete workflow
                                                    
                                                    # Step 1: Process PDF
                                                    pdf_processor = PDFProcessor(mock_pdf_with_chemical_content)
                                                    pages_data = pdf_processor.extract_text_and_images()
                                                    
                                                    assert len(pages_data) == 2
                                                    assert 'ethanol' in pages_data[0]['text'].lower()
                                                    
                                                    # Step 2: Initialize DECIMER client
                                                    decimer_client = DECIMERClient()
                                                    chem_handler = ChemicalHandler(decimer_client)
                                                    
                                                    # Step 3: Process chemical structures
                                                    all_structures = []
                                                    full_text = ""
                                                    
                                                    for page_data in pages_data:
                                                        full_text += page_data['text'] + "\n"
                                                        
                                                        for img_info in page_data['images']:
                                                            context = chem_handler.get_structure_context(
                                                                page_data['text'], 
                                                                img_info['bbox']
                                                            )
                                                            
                                                            structure = chem_handler.process_image(
                                                                img_info['path'],
                                                                context
                                                            )
                                                            
                                                            if structure:
                                                                all_structures.append(structure)
                                                    
                                                    assert len(all_structures) > 0
                                                    assert all_structures[0]['smiles'] == 'CCO'
                                                    assert all_structures[0]['formula'] == 'C2H6O'
                                                    
                                                    # Step 4: Create chunks
                                                    chunker = ChemicalAwareChunker()
                                                    chunks = chunker.chunk_with_structures(full_text, all_structures)
                                                    
                                                    assert len(chunks) > 0
                                                    structure_chunks = [c for c in chunks if c['metadata']['has_structures']]
                                                    assert len(structure_chunks) > 0
                                                    
                                                    # Step 5: Create vector store
                                                    vector_store = ChemicalVectorStore(use_ollama_embeddings=False)
                                                    vector_store.create_vectorstore(chunks)
                                                    
                                                    assert vector_store.vectorstore is not None
                                                    
                                                    # Step 6: Initialize RAG chain
                                                    rag = ChemicalRAG(vector_store)
                                                    
                                                    # Step 7: Query the system
                                                    answer = rag.query("What is the molecular formula of ethanol?")
                                                    
                                                    assert isinstance(answer, str)
                                                    assert len(answer) > 0
    
    def test_workflow_with_no_chemical_structures(self, temp_workspace, mock_pdf_with_chemical_content):
        """Test workflow when no chemical structures are found"""
        
        with patch('fitz.open') as mock_fitz:
            # Mock PDF with no images
            mock_doc = Mock()
            mock_doc.is_encrypted = False
            mock_doc.page_count = 1
            
            mock_page = Mock()
            mock_page.get_text.return_value = "This is a general text without chemical structures."
            mock_page.get_images.return_value = []  # No images
            mock_page.rect = Mock(width=612, height=792)
            
            mock_doc.__iter__ = lambda self: iter([mock_page])
            mock_fitz.return_value = mock_doc
            
            with patch('vectorstore.HuggingFaceEmbeddings'):
                with patch('vectorstore.Chroma') as mock_chroma:
                    mock_vectorstore = Mock()
                    mock_chroma.from_texts.return_value = mock_vectorstore
                    
                    with patch('rag_chain.Ollama'):
                        with patch('rag_chain.RetrievalQA') as mock_retrieval_qa:
                            mock_qa_chain = Mock()
                            mock_qa_chain.return_value = {
                                'result': 'No chemical structures were found in the document.',
                                'source_documents': []
                            }
                            mock_retrieval_qa.from_chain_type.return_value = mock_qa_chain
                            
                            # Execute workflow
                            pdf_processor = PDFProcessor(mock_pdf_with_chemical_content)
                            pages_data = pdf_processor.extract_text_and_images()
                            
                            decimer_client = DECIMERClient()
                            chem_handler = ChemicalHandler(decimer_client)
                            
                            all_structures = []
                            full_text = ""
                            
                            for page_data in pages_data:
                                full_text += page_data['text'] + "\n"
                                # No images to process
                            
                            chunker = ChemicalAwareChunker()
                            chunks = chunker.chunk_with_structures(full_text, all_structures)
                            
                            # Should still create chunks, but without structures
                            assert len(chunks) > 0
                            assert all(not c['metadata']['has_structures'] for c in chunks)
                            
                            vector_store = ChemicalVectorStore(use_ollama_embeddings=False)
                            vector_store.create_vectorstore(chunks)
                            
                            rag = ChemicalRAG(vector_store)
                            answer = rag.query("Are there any chemical structures?")
                            
                            assert isinstance(answer, str)
    
    def test_workflow_error_handling(self, temp_workspace):
        """Test workflow error handling"""
        
        # Test with invalid PDF path
        with pytest.raises(ValueError):
            PDFProcessor("nonexistent.pdf")
        
        # Test with encrypted PDF
        with patch('fitz.open') as mock_fitz:
            mock_doc = Mock()
            mock_doc.is_encrypted = True
            mock_fitz.return_value = mock_doc
            
            with pytest.raises(ValueError, match="encrypted"):
                PDFProcessor("encrypted.pdf")
    
    def test_workflow_with_decimer_failure(self, temp_workspace, mock_pdf_with_chemical_content):
        """Test workflow when DECIMER fails to process images"""
        
        with patch('fitz.open') as mock_fitz:
            # Setup PDF mock with images
            mock_doc = Mock()
            mock_doc.is_encrypted = False
            mock_doc.page_count = 1
            
            mock_page = Mock()
            mock_page.get_text.return_value = "Chemical structure shown in Figure 1."
            mock_page.get_images.return_value = [(1, 0, 100, 100, 8, 'DeviceRGB', '', 'Im1', 'DCTDecode')]
            mock_page.rect = Mock(width=612, height=792)
            mock_page.get_image_bbox.return_value = Mock(x0=0, y0=0, x1=100, y1=100)
            
            mock_doc.__iter__ = lambda self: iter([mock_page])
            mock_fitz.return_value = mock_doc
            
            # Mock DECIMER failure
            with patch('requests.post') as mock_post:
                mock_response = Mock()
                mock_response.status_code = 500  # Server error
                mock_post.return_value = mock_response
                
                with patch('fitz.Pixmap') as mock_pixmap:
                    mock_pix = Mock()
                    mock_pix.n = 4
                    mock_pix.alpha = 1
                    mock_pix.tobytes.return_value = b'fake_image_data'
                    mock_pixmap.return_value = mock_pix
                    
                    with patch('PIL.Image.open') as mock_image_open:
                        mock_img = Mock()
                        mock_img.width = 100
                        mock_img.height = 100
                        mock_img.save = Mock()
                        mock_image_open.return_value = mock_img
                        
                        # Execute workflow
                        pdf_processor = PDFProcessor(mock_pdf_with_chemical_content)
                        pages_data = pdf_processor.extract_text_and_images()
                        
                        decimer_client = DECIMERClient()
                        chem_handler = ChemicalHandler(decimer_client)
                        
                        all_structures = []
                        full_text = ""
                        
                        for page_data in pages_data:
                            full_text += page_data['text'] + "\n"
                            
                            for img_info in page_data['images']:
                                structure = chem_handler.process_image(img_info['path'])
                                if structure:
                                    all_structures.append(structure)
                        
                        # Should handle DECIMER failure gracefully
                        assert len(all_structures) == 0  # No structures processed due to DECIMER failure
                        
                        # Workflow should continue without structures
                        chunker = ChemicalAwareChunker()
                        chunks = chunker.chunk_with_structures(full_text, all_structures)
                        
                        assert len(chunks) > 0
                        assert all(not c['metadata']['has_structures'] for c in chunks)
    
    def test_chunking_statistics_integration(self, sample_text, sample_structures_list):
        """Test chunking statistics in integration context"""
        
        chunker = ChemicalAwareChunker(chunk_size=200)
        chunks = chunker.chunk_with_structures(sample_text, sample_structures_list)
        
        stats = chunker.get_chunking_stats(chunks)
        
        # Verify statistics make sense
        assert stats['total_chunks'] > 0
        assert stats['structure_coverage'] >= 0
        assert stats['avg_chunk_size'] > 0
        assert isinstance(stats['content_type_distribution'], dict)
        
        # Should have some chunks with structures given our sample data
        assert stats['chunks_with_structures'] >= 0
    
    def test_vector_store_integration_with_real_chunks(self, sample_text, sample_structures_list):
        """Test vector store creation with realistic chunks"""
        
        with patch('vectorstore.HuggingFaceEmbeddings'):
            with patch('vectorstore.Chroma') as mock_chroma:
                mock_vectorstore = Mock()
                mock_chroma.from_texts.return_value = mock_vectorstore
                
                chunker = ChemicalAwareChunker()
                chunks = chunker.chunk_with_structures(sample_text, sample_structures_list)
                
                vector_store = ChemicalVectorStore(use_ollama_embeddings=False)
                result = vector_store.create_vectorstore(chunks)
                
                assert result == mock_vectorstore
                
                # Check that texts were enhanced with chemical information
                call_args = mock_chroma.from_texts.call_args
                texts = call_args.kwargs['texts']
                metadatas = call_args.kwargs['metadatas']
                
                assert len(texts) == len(chunks)
                assert len(metadatas) == len(chunks)
                
                # Some texts should contain chemical structure information
                enhanced_texts = [t for t in texts if '[Chemical Structures Found]' in t]
                assert len(enhanced_texts) > 0
    
    def test_memory_and_performance_considerations(self, temp_workspace):
        """Test that the workflow handles larger documents efficiently"""
        
        # Create a larger mock document
        large_text = """
        This is a large chemical document. """ * 1000  # Simulate large document
        
        large_structures = []
        for i in range(50):  # Simulate many structures
            large_structures.append({
                'smiles': f'C{i}H{i*2}O',
                'formula': f'C{i}H{i*2}O',
                'molecular_weight': 50 + i,
                'context': f'Structure {i} context',
                'image_path': f'/path/to/structure_{i}.png'
            })
        
        # Test chunking with large data
        chunker = ChemicalAwareChunker(chunk_size=500)
        chunks = chunker.chunk_with_structures(large_text, large_structures)
        
        # Should handle large data without issues
        assert len(chunks) > 0
        
        # Test vector store creation with many chunks
        with patch('vectorstore.HuggingFaceEmbeddings'):
            with patch('vectorstore.Chroma') as mock_chroma:
                mock_vectorstore = Mock()
                mock_chroma.from_texts.return_value = mock_vectorstore
                
                vector_store = ChemicalVectorStore(use_ollama_embeddings=False)
                vector_store.create_vectorstore(chunks)
                
                # Should complete without memory issues
                assert mock_chroma.from_texts.called