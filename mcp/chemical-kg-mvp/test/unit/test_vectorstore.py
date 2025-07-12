#!/usr/bin/env python3
"""
Test vector store functionality with extracted structures and multi-document support
"""

import sys
import os
import pytest
import tempfile
from unittest.mock import Mock, patch

sys.path.insert(0, 'src')

def test_vectorstore_basic():
    """Test basic vector store functionality"""
    print("Testing Vector Store")
    print("=" * 40)
    
    try:
        from vectorstore import ChemicalVectorStore
        print("✓ ChemicalVectorStore imported successfully")
    except Exception as e:
        print(f"✗ Failed to import ChemicalVectorStore: {e}")
        return False
    
    # Test initialization
    try:
        # Try HuggingFace embeddings first (no external dependencies)
        vector_store = ChemicalVectorStore(use_ollama_embeddings=False)
        print("✓ Vector store initialized with HuggingFace embeddings")
    except Exception as e:
        print(f"✗ Failed to initialize vector store: {e}")
        return False
    
    return True

def test_with_sample_data():
    """Test with sample chemical data"""
    print("\nTesting with Sample Data")
    print("=" * 40)
    
    try:
        from chunker import ChemicalAwareChunker
        from vectorstore import ChemicalVectorStore
        
        # Create sample data
        sample_text = """
        Aspirin (acetylsalicylic acid) is a medication used to reduce pain, fever, 
        or inflammation. It is one of the most widely used medications globally.
        
        Ibuprofen is a nonsteroidal anti-inflammatory drug (NSAID) that is used 
        for treating pain, fever, and inflammation.
        """
        
        sample_structures = [
            {
                "smiles": "CC(=O)Oc1ccccc1C(=O)O",
                "formula": "C9H8O4",
                "molecular_weight": 180.16,
                "context": "Aspirin chemical structure"
            },
            {
                "smiles": "CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O",
                "formula": "C13H18O2",
                "molecular_weight": 206.28,
                "context": "Ibuprofen chemical structure"
            }
        ]
        
        # Create chunks
        chunker = ChemicalAwareChunker()
        chunks = chunker.chunk_with_structures(sample_text, sample_structures)
        print(f"Created {len(chunks)} chunks")
        
        # Create vector store
        vector_store = ChemicalVectorStore(use_ollama_embeddings=False)
        vector_store.create_vectorstore(chunks)
        print("✓ Vector store created successfully")
        
        # Test search
        results = vector_store.similarity_search("What is aspirin?", k=2)
        print(f"✓ Search returned {len(results)} results")
        
        return True
        
    except Exception as e:
        print(f"✗ Error in sample data test: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_with_extracted_structures():
    """Test with structures extracted from PDF"""
    print("\nTesting with Extracted Structures")
    print("=" * 40)
    
    # Check if we have extracted structures
    structure_files = []
    for f in os.listdir('.'):
        if f.startswith('structure_') and f.endswith('.png'):
            structure_files.append(f)
    
    if not structure_files:
        print("No extracted structure files found")
        return False
    
    print(f"Found {len(structure_files)} structure files")
    
    try:
        from decimer_client import DECIMERClient
        from chemical_handler import ChemicalHandler
        from chunker import ChemicalAwareChunker
        from vectorstore import ChemicalVectorStore
        
        # Initialize components
        decimer_client = DECIMERClient()
        chem_handler = ChemicalHandler(decimer_client)
        
        # Process structures
        structures = []
        for i, struct_file in enumerate(structure_files[:3]):  # Test first 3
            print(f"Processing {struct_file}...")
            try:
                structure = chem_handler.process_image(struct_file, f"Structure {i+1}")
                if structure:
                    structures.append(structure)
                    print(f"  ✓ SMILES: {structure.get('smiles', 'N/A')}")
            except Exception as e:
                print(f"  ✗ Failed: {e}")
        
        if not structures:
            print("No structures processed successfully")
            return False
        
        # Create vector store
        print(f"\nCreating vector store with {len(structures)} structures...")
        
        sample_text = "Chemical structures extracted from scientific literature."
        chunker = ChemicalAwareChunker()
        chunks = chunker.chunk_with_structures(sample_text, structures)
        
        vector_store = ChemicalVectorStore(use_ollama_embeddings=False)
        vector_store.create_vectorstore(chunks)
        
        print("✓ Vector store created successfully with real structures")
        
        return True
        
    except Exception as e:
        print(f"✗ Error: {e}")
        import traceback
        traceback.print_exc()
        return False

def diagnose_type_error():
    """Diagnose the type subscription error"""
    print("\nDiagnosing Type Error")
    print("=" * 40)
    
    # Check Python version
    import sys
    print(f"Python version: {sys.version}")
    
    # Check if it's a type hint issue
    try:
        from typing import List, Dict, Any
        test_type = List[str]  # This would fail in Python < 3.9 without __future__
        print("✓ Type hints work correctly")
    except Exception as e:
        print(f"✗ Type hint error: {e}")
        print("  This might be the cause of the 'type' object is not subscriptable error")
    
    # Test specific imports that might cause issues
    try:
        import langchain
        print(f"✓ LangChain version: {langchain.__version__}")
    except:
        print("✗ LangChain not found")
    
    try:
        import chromadb
        print(f"✓ ChromaDB version: {chromadb.__version__}")
    except:
        print("✗ ChromaDB not found")

if __name__ == "__main__":
    print("Vector Store Testing Suite")
    print("=" * 50)
    
    # Run diagnostics first
    diagnose_type_error()
    
    # Test basic functionality
    if test_vectorstore_basic():
        # Test with sample data
        test_with_sample_data()
        
        # Test with real extracted structures
        test_with_extracted_structures()
    
    print("\nTest completed!")

def test_vectorstore_multi_document_metadata():
    """Test that vector store correctly handles multi-document metadata"""
    print("\nTesting Multi-Document Metadata")
    print("=" * 40)
    
    try:
        from vectorstore import ChemicalVectorStore
        from chunker import ChemicalAwareChunker
        
        # Create sample chunks from multiple documents
        sample_chunks = [
            {
                'chunk_id': 'doc1_chunk1',
                'text': 'Content from document 1 about aspirin synthesis',
                'metadata': {'page': 1, 'has_structures': True},
                'structures': [{'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O'}],
                'document': 'paper1.pdf',
                'document_index': 0
            },
            {
                'chunk_id': 'doc2_chunk1', 
                'text': 'Content from document 2 about ibuprofen',
                'metadata': {'page': 1, 'has_structures': True},
                'structures': [{'smiles': 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O'}],
                'document': 'paper2.pdf',
                'document_index': 1
            },
            {
                'chunk_id': 'doc1_chunk2',
                'text': 'More content from document 1',
                'metadata': {'page': 2, 'has_structures': False},
                'structures': [],
                'document': 'paper1.pdf', 
                'document_index': 0
            }
        ]
        
        # Mock the vectorstore creation
        with patch('vectorstore.FAISS') as mock_faiss:
            mock_vectorstore = Mock()
            mock_faiss.from_texts.return_value = mock_vectorstore
            
            with patch('vectorstore.OllamaEmbeddings') as mock_embeddings:
                mock_embeddings.return_value = Mock()
                
                vector_store = ChemicalVectorStore()
                vector_store.create_vectorstore(sample_chunks)
                
                # Verify FAISS was called with proper metadata
                call_args = mock_faiss.from_texts.call_args
                metadatas = call_args.kwargs['metadatas']
                
                # Check that document metadata is preserved
                assert len(metadatas) == 3
                assert metadatas[0]['document'] == 'paper1.pdf'
                assert metadatas[1]['document'] == 'paper2.pdf'
                assert metadatas[0]['document_index'] == 0
                assert metadatas[1]['document_index'] == 1
                
                # Check structure metadata
                assert metadatas[0]['has_structures'] == True
                assert metadatas[2]['has_structures'] == False
                
                print("✓ Multi-document metadata correctly handled")
                return True
                
    except Exception as e:
        print(f"✗ Multi-document metadata test failed: {e}")
        return False

def test_chunker_document_association():
    """Test that chunker properly associates content with source documents"""
    print("\nTesting Document Association in Chunker")
    print("=" * 40)
    
    try:
        from chunker import ChemicalAwareChunker
        
        # Simulate multi-document processing
        doc1_text = "Document 1 discusses aspirin synthesis. The molecular formula is C9H8O4."
        doc1_structures = [{
            'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O',
            'formula': 'C9H8O4',
            'document': 'paper1.pdf',
            'document_index': 0
        }]
        
        doc2_text = "Document 2 covers ibuprofen properties. Its formula is C13H18O2."
        doc2_structures = [{
            'smiles': 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O',
            'formula': 'C13H18O2', 
            'document': 'paper2.pdf',
            'document_index': 1
        }]
        
        chunker = ChemicalAwareChunker(chunk_size=100, overlap=20)
        
        # Process each document separately (as app.py does)
        doc1_chunks = chunker.chunk_with_structures(doc1_text, doc1_structures)
        doc2_chunks = chunker.chunk_with_structures(doc2_text, doc2_structures)
        
        # Manually add document metadata (as app.py does)
        for chunk in doc1_chunks:
            chunk['document'] = 'paper1.pdf'
            chunk['document_index'] = 0
            
        for chunk in doc2_chunks:
            chunk['document'] = 'paper2.pdf'
            chunk['document_index'] = 1
        
        # Verify document association
        assert all(chunk['document'] == 'paper1.pdf' for chunk in doc1_chunks)
        assert all(chunk['document'] == 'paper2.pdf' for chunk in doc2_chunks)
        assert all(chunk['document_index'] == 0 for chunk in doc1_chunks)
        assert all(chunk['document_index'] == 1 for chunk in doc2_chunks)
        
        print("✓ Document association in chunker working correctly")
        return True
        
    except Exception as e:
        print(f"✗ Document association test failed: {e}")
        return False

def test_rag_multi_document_query():
    """Test RAG chain with multi-document context"""
    print("\nTesting Multi-Document RAG Queries")
    print("=" * 40)
    
    try:
        from rag_chain import ChemicalRAG
        from vectorstore import ChemicalVectorStore
        
        # Mock vector store with multi-document data
        mock_vectorstore = Mock()
        mock_retriever = Mock()
        
        # Mock retrieved documents from different sources
        mock_docs = [
            Mock(
                page_content="Aspirin synthesis from paper1.pdf",
                metadata={'document': 'paper1.pdf', 'page': 1, 'has_structures': True}
            ),
            Mock(
                page_content="Ibuprofen properties from paper2.pdf", 
                metadata={'document': 'paper2.pdf', 'page': 2, 'has_structures': True}
            )
        ]
        
        mock_retriever.get_relevant_documents.return_value = mock_docs
        mock_vectorstore.vectorstore.as_retriever.return_value = mock_retriever
        
        with patch('rag_chain.Ollama') as mock_ollama:
            with patch('rag_chain.RetrievalQA') as mock_qa:
                mock_qa_chain = Mock()
                mock_qa_chain.return_value = {
                    'result': 'Aspirin and ibuprofen are both pain relievers found in the documents.',
                    'source_documents': mock_docs
                }
                mock_qa.from_chain_type.return_value = mock_qa_chain
                
                rag = ChemicalRAG(mock_vectorstore)
                result = rag.query("Compare aspirin and ibuprofen")
                
                # Verify multi-document sources are included
                sources = result['sources']
                assert len(sources) == 2
                assert sources[0]['document'] == 'paper1.pdf'
                assert sources[1]['document'] == 'paper2.pdf'
                
                print("✓ Multi-document RAG query working correctly")
                return True
                
    except Exception as e:
        print(f"✗ Multi-document RAG test failed: {e}")
        return False