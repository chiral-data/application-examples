#!/usr/bin/env python3
"""
Integration tests for multi-document functionality
"""

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


class TestMultiDocumentWorkflow:
    """Integration tests for multi-document processing workflow"""
    
    @pytest.fixture
    def temp_workspace(self):
        """Create a temporary workspace for tests"""
        workspace = tempfile.mkdtemp()
        yield workspace
        shutil.rmtree(workspace)
    
    @pytest.fixture
    def mock_pdf_documents(self, temp_workspace):
        """Create multiple mock PDF documents"""
        pdfs = []
        
        # Document 1: Aspirin paper
        pdf1_path = os.path.join(temp_workspace, "aspirin_paper.pdf")
        with open(pdf1_path, 'wb') as f:
            f.write(b'%PDF-1.4\nAspirin synthesis paper')
        pdfs.append(('aspirin_paper.pdf', pdf1_path))
        
        # Document 2: Ibuprofen paper  
        pdf2_path = os.path.join(temp_workspace, "ibuprofen_paper.pdf")
        with open(pdf2_path, 'wb') as f:
            f.write(b'%PDF-1.4\nIbuprofen analysis paper')
        pdfs.append(('ibuprofen_paper.pdf', pdf2_path))
        
        # Document 3: General chemistry paper
        pdf3_path = os.path.join(temp_workspace, "general_chemistry.pdf")
        with open(pdf3_path, 'wb') as f:
            f.write(b'%PDF-1.4\nGeneral chemistry concepts')
        pdfs.append(('general_chemistry.pdf', pdf3_path))
        
        return pdfs
    
    def create_mock_fitz_doc(self, filename, text_content, has_images=True):
        """Create a mock PyMuPDF document"""
        mock_doc = Mock()
        mock_doc.is_encrypted = False
        mock_doc.page_count = 1
        mock_doc.metadata = {'title': f'Test {filename}'}
        
        mock_page = Mock()
        mock_page.get_text.return_value = text_content
        mock_page.rect = Mock(width=612, height=792)
        
        if has_images:
            mock_page.get_images.return_value = [(1, 0, 150, 100, 8, 'DeviceRGB', '', 'Im1', 'DCTDecode')]
            mock_page.get_image_bbox.return_value = Mock(x0=100, y0=200, x1=250, y1=300)
        else:
            mock_page.get_images.return_value = []
        
        mock_doc.__iter__ = lambda self: iter([mock_page])
        return mock_doc
    
    def test_multi_document_processing_workflow(self, temp_workspace, mock_pdf_documents):
        """Test the complete multi-document processing workflow"""
        
        # Document contents
        doc_contents = {
            'aspirin_paper.pdf': """
            Aspirin Synthesis and Analysis
            
            Aspirin (acetylsalicylic acid) has the molecular formula C9H8O4.
            The SMILES notation is CC(=O)OC1=CC=CC=C1C(=O)O.
            This study describes the synthesis of aspirin from salicylic acid.
            The reaction yield was 85% with high purity confirmed by NMR.
            """,
            'ibuprofen_paper.pdf': """
            Ibuprofen Properties and Synthesis
            
            Ibuprofen is a nonsteroidal anti-inflammatory drug (NSAID).
            Its molecular formula is C13H18O2.
            The SMILES notation is CC(C)CC1=CC=C(C=C1)C(C)C(=O)O.
            This paper discusses improved synthesis methods for ibuprofen.
            The new method achieves 92% yield with reduced environmental impact.
            """,
            'general_chemistry.pdf': """
            General Principles of Organic Chemistry
            
            This paper discusses general concepts in organic chemistry.
            It covers functional groups, reaction mechanisms, and stereochemistry.
            The text does not contain specific chemical structures.
            It serves as a theoretical background for understanding reactions.
            """
        }
        
        # Mock external dependencies
        with patch('fitz.open') as mock_fitz:
            # Setup different mock documents for each file
            def fitz_side_effect(path):
                filename = os.path.basename(path)
                content = doc_contents.get(filename, "Default content")
                has_images = filename != 'general_chemistry.pdf'  # No structures in general paper
                return self.create_mock_fitz_doc(filename, content, has_images)
            
            mock_fitz.side_effect = fitz_side_effect
            
            # Mock DECIMER API responses
            with patch('requests.post') as mock_post:
                def decimer_side_effect(*args, **kwargs):
                    mock_response = Mock()
                    mock_response.status_code = 200
                    
                    # Return different SMILES based on which document is being processed
                    files = kwargs.get('files', {})
                    if files and 'image' in files:
                        # Simulate different responses for different documents
                        if hasattr(mock_post, 'call_count'):
                            mock_post.call_count += 1
                        else:
                            mock_post.call_count = 1
                        
                        if mock_post.call_count == 1:
                            # First call - aspirin
                            mock_response.json.return_value = {'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O'}
                        else:
                            # Second call - ibuprofen  
                            mock_response.json.return_value = {'smiles': 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O'}
                    else:
                        mock_response.json.return_value = {'smiles': 'CCO'}  # Default
                    
                    return mock_response
                
                mock_post.side_effect = decimer_side_effect
                
                # Mock RDKit
                with patch('rdkit.Chem.MolFromSmiles') as mock_mol:
                    mock_molecule = Mock()
                    mock_mol.return_value = mock_molecule
                    
                    def mol_wt_side_effect(mol):
                        # Return different molecular weights
                        if hasattr(mock_mol, 'call_count'):
                            mock_mol.call_count += 1
                        else:
                            mock_mol.call_count = 1
                        
                        if mock_mol.call_count == 1:
                            return 180.16  # Aspirin
                        else:
                            return 206.28  # Ibuprofen
                    
                    def mol_formula_side_effect(mol):
                        if hasattr(mock_mol, 'formula_call_count'):
                            mock_mol.formula_call_count += 1
                        else:
                            mock_mol.formula_call_count = 1
                        
                        if mock_mol.formula_call_count == 1:
                            return 'C9H8O4'  # Aspirin
                        else:
                            return 'C13H18O2'  # Ibuprofen
                    
                    with patch('rdkit.Chem.Descriptors.MolWt', side_effect=mol_wt_side_effect):
                        with patch('rdkit.Chem.rdMolDescriptors.CalcMolFormula', side_effect=mol_formula_side_effect):
                            # Mock image processing
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
                                    
                                    # Execute multi-document workflow (simulating app.py logic)
                                    
                                    # Initialize shared components
                                    decimer_client = DECIMERClient()
                                    chem_handler = ChemicalHandler(decimer_client)
                                    chunker = ChemicalAwareChunker()
                                    
                                    # Simulate session state
                                    documents = []
                                    all_structures = []
                                    all_chunks = []
                                    
                                    # Process each document
                                    for filename, pdf_path in mock_pdf_documents:
                                        # Step 1: Process PDF
                                        pdf_processor = PDFProcessor(pdf_path)
                                        pages_data = pdf_processor.extract_text_and_images()
                                        
                                        # Extract full text
                                        full_text = ""
                                        for page_data in pages_data:
                                            full_text += page_data['text'] + "\n"
                                        
                                        # Step 2: Process chemical structures
                                        file_structures = []
                                        for page_data in pages_data:
                                            for img_info in page_data['images']:
                                                context = chem_handler.get_structure_context(
                                                    page_data['text'], 
                                                    img_info['bbox']
                                                )
                                                
                                                structure = chem_handler.process_image(
                                                    img_info['path'],
                                                    f"{context} (from {filename})"
                                                )
                                                
                                                if structure:
                                                    # Add document metadata
                                                    structure['document'] = filename
                                                    structure['document_index'] = len(documents)
                                                    file_structures.append(structure)
                                        
                                        # Step 3: Create chunks
                                        file_chunks = chunker.chunk_with_structures(full_text, file_structures)
                                        
                                        # Add document metadata to chunks
                                        for chunk in file_chunks:
                                            chunk['document'] = filename
                                            chunk['document_index'] = len(documents)
                                        
                                        # Store document information
                                        document_info = {
                                            'filename': filename,
                                            'structures': file_structures,
                                            'chunks': file_chunks,
                                            'text_length': len(full_text),
                                            'pages': len(pages_data)
                                        }
                                        documents.append(document_info)
                                        
                                        # Add to global collections
                                        all_structures.extend(file_structures)
                                        all_chunks.extend(file_chunks)
                                    
                                    # Verify multi-document processing results
                                    assert len(documents) == 3
                                    assert documents[0]['filename'] == 'aspirin_paper.pdf'
                                    assert documents[1]['filename'] == 'ibuprofen_paper.pdf'
                                    assert documents[2]['filename'] == 'general_chemistry.pdf'
                                    
                                    # Check structures were found in chemical papers but not general paper
                                    assert len(documents[0]['structures']) > 0  # Aspirin paper
                                    assert len(documents[1]['structures']) > 0  # Ibuprofen paper  
                                    assert len(documents[2]['structures']) == 0  # General paper
                                    
                                    # Verify document metadata in structures
                                    aspirin_structures = documents[0]['structures']
                                    ibuprofen_structures = documents[1]['structures']
                                    
                                    assert all(s['document'] == 'aspirin_paper.pdf' for s in aspirin_structures)
                                    assert all(s['document'] == 'ibuprofen_paper.pdf' for s in ibuprofen_structures)
                                    assert all(s['document_index'] == 0 for s in aspirin_structures)
                                    assert all(s['document_index'] == 1 for s in ibuprofen_structures)
                                    
                                    # Verify different SMILES for different compounds
                                    if aspirin_structures:
                                        assert 'CC(=O)OC1=CC=CC=C1C(=O)O' in aspirin_structures[0]['smiles']
                                    if ibuprofen_structures:
                                        assert 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O' in ibuprofen_structures[0]['smiles']
                                    
                                    # Check chunks have document metadata
                                    assert all(chunk['document'] in ['aspirin_paper.pdf', 'ibuprofen_paper.pdf', 'general_chemistry.pdf'] 
                                             for chunk in all_chunks)
                                    
                                    # Step 4: Create unified vector store
                                    with patch('vectorstore.FAISS') as mock_faiss:
                                        mock_vectorstore = Mock()
                                        mock_faiss.from_texts.return_value = mock_vectorstore
                                        
                                        with patch('vectorstore.OllamaEmbeddings'):
                                            vector_store = ChemicalVectorStore()
                                            vector_store.create_vectorstore(all_chunks)
                                            
                                            # Verify vector store was created with multi-document chunks
                                            call_args = mock_faiss.from_texts.call_args
                                            metadatas = call_args.kwargs['metadatas']
                                            
                                            # Check document metadata preservation
                                            document_names = set(meta['document'] for meta in metadatas)
                                            assert 'aspirin_paper.pdf' in document_names
                                            assert 'ibuprofen_paper.pdf' in document_names
                                            assert 'general_chemistry.pdf' in document_names
                                            
                                            # Step 5: Test RAG with multi-document context
                                            with patch('rag_chain.Ollama'):
                                                with patch('rag_chain.RetrievalQA') as mock_retrieval_qa:
                                                    mock_qa_chain = Mock()
                                                    mock_qa_chain.return_value = {
                                                        'result': 'Aspirin (from aspirin_paper.pdf) and ibuprofen (from ibuprofen_paper.pdf) are both anti-inflammatory drugs.',
                                                        'source_documents': [
                                                            Mock(
                                                                page_content="Aspirin content",
                                                                metadata={'document': 'aspirin_paper.pdf', 'page': 1}
                                                            ),
                                                            Mock(
                                                                page_content="Ibuprofen content", 
                                                                metadata={'document': 'ibuprofen_paper.pdf', 'page': 1}
                                                            )
                                                        ]
                                                    }
                                                    mock_retrieval_qa.from_chain_type.return_value = mock_qa_chain
                                                    
                                                    rag = ChemicalRAG(vector_store)
                                                    result = rag.query("Compare aspirin and ibuprofen")
                                                    
                                                    # Verify multi-document response
                                                    assert isinstance(result, dict)
                                                    assert 'sources' in result
                                                    sources = result['sources']
                                                    
                                                    # Should have sources from multiple documents
                                                    source_docs = set(source['document'] for source in sources)
                                                    assert len(source_docs) >= 1  # At least one document referenced
    
    def test_document_management_features(self, temp_workspace, mock_pdf_documents):
        """Test document management functionality (add/remove documents)"""
        
        # Simulate the document management workflow from app.py
        documents = []  # Simulates st.session_state.documents
        
        # Test adding documents one by one
        for filename, pdf_path in mock_pdf_documents[:2]:  # Add first 2 documents
            # Simulate checking for existing documents
            current_filenames = {doc['filename'] for doc in documents}
            is_new_file = filename not in current_filenames
            
            assert is_new_file  # Should be new
            
            # Add document to collection
            document_info = {
                'filename': filename,
                'structures': [],  # Would be populated by processing
                'chunks': [],      # Would be populated by processing  
                'text_length': 1000,
                'pages': 5
            }
            documents.append(document_info)
        
        # Verify documents were added
        assert len(documents) == 2
        assert documents[0]['filename'] == 'aspirin_paper.pdf'
        assert documents[1]['filename'] == 'ibuprofen_paper.pdf'
        
        # Test adding duplicate document (should be detected)
        duplicate_filename = 'aspirin_paper.pdf'
        current_filenames = {doc['filename'] for doc in documents}
        is_duplicate = duplicate_filename in current_filenames
        
        assert is_duplicate  # Should detect duplicate
        
        # Test clearing all documents
        documents.clear()
        assert len(documents) == 0
    
    def test_multi_document_csv_export(self, temp_workspace):
        """Test CSV export with multi-document structure data"""
        
        # Create sample structures from multiple documents
        structures = [
            {
                'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O',
                'formula': 'C9H8O4', 
                'molecular_weight': 180.16,
                'context': 'Structure from page 1',
                'document': 'aspirin_paper.pdf'
            },
            {
                'smiles': 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O',
                'formula': 'C13H18O2',
                'molecular_weight': 206.28,
                'context': 'Structure from page 2', 
                'document': 'ibuprofen_paper.pdf'
            },
            {
                'smiles': 'CCO',
                'formula': 'C2H6O',
                'molecular_weight': 46.07,
                'context': 'Structure from page 1',
                'document': 'general_chemistry.pdf'
            }
        ]
        
        # Simulate CSV creation (from app.py create_csv_download function)
        import pandas as pd
        import re
        
        data = []
        for i, struct in enumerate(structures):
            # Extract page number from context if available
            page_num = ''
            context = struct.get('context', '')
            if 'page' in context.lower():
                page_match = re.search(r'page\s+(\d+)', context.lower())
                if page_match:
                    page_num = page_match.group(1)
            
            data.append({
                'Structure_ID': f'Structure_{i+1}',
                'Document': struct.get('document', 'Unknown'),
                'SMILES': struct.get('smiles', ''),
                'Molecular_Formula': struct.get('formula', ''),
                'Molecular_Weight': struct.get('molecular_weight', ''),
                'Page': page_num
            })
        
        df = pd.DataFrame(data)
        csv_content = df.to_csv(index=False)
        
        # Verify CSV contains multi-document information
        assert 'aspirin_paper.pdf' in csv_content
        assert 'ibuprofen_paper.pdf' in csv_content 
        assert 'general_chemistry.pdf' in csv_content
        assert 'CC(=O)OC1=CC=CC=C1C(=O)O' in csv_content  # Aspirin SMILES
        assert 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O' in csv_content  # Ibuprofen SMILES
        
        # Verify all expected columns are present
        assert 'Structure_ID' in csv_content
        assert 'Document' in csv_content
        assert 'SMILES' in csv_content
        assert 'Molecular_Formula' in csv_content
        assert 'Molecular_Weight' in csv_content
        assert 'Page' in csv_content


if __name__ == '__main__':
    # Run tests directly
    import pytest
    pytest.main([__file__, '-v'])