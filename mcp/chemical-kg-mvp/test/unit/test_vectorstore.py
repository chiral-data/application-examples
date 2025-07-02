import pytest
from unittest.mock import Mock, patch, MagicMock
import os

from vectorstore import ChemicalVectorStore


class TestChemicalVectorStore:
    """Test suite for ChemicalVectorStore class"""
    
    def test_init_ollama_embeddings(self):
        """Test initialization with Ollama embeddings"""
        with patch('vectorstore.OllamaEmbeddings') as mock_ollama:
            store = ChemicalVectorStore(use_ollama_embeddings=True)
            
            mock_ollama.assert_called_once()
            assert store.vectorstore is None
            assert store.structures_db == {}
    
    def test_init_huggingface_embeddings(self):
        """Test initialization with HuggingFace embeddings"""
        with patch('vectorstore.HuggingFaceEmbeddings') as mock_hf:
            store = ChemicalVectorStore(use_ollama_embeddings=False)
            
            mock_hf.assert_called_once()
            assert store.vectorstore is None
    
    def test_init_environment_variables(self):
        """Test initialization with environment variables"""
        with patch.dict(os.environ, {
            'OLLAMA_HOST': 'test-host',
            'OLLAMA_PORT': '12345'
        }):
            with patch('vectorstore.OllamaEmbeddings') as mock_ollama:
                store = ChemicalVectorStore(use_ollama_embeddings=True)
                
                # Check that environment variables were used
                call_args = mock_ollama.call_args
                assert 'base_url' in call_args.kwargs
                assert 'test-host:12345' in call_args.kwargs['base_url']
    
    def test_create_vectorstore_basic(self, sample_structures_list):
        """Test basic vector store creation"""
        with patch('vectorstore.Chroma') as mock_chroma:
            mock_vectorstore = Mock()
            mock_chroma.from_texts.return_value = mock_vectorstore
            
            store = ChemicalVectorStore(use_ollama_embeddings=False)
            
            # Create sample chunks
            chunks = [
                {
                    'text': 'First chunk about ethanol',
                    'chunk_id': 0,
                    'structures': [sample_structures_list[0]]
                },
                {
                    'text': 'Second chunk about benzene',
                    'chunk_id': 1,
                    'structures': [sample_structures_list[1]]
                }
            ]
            
            result = store.create_vectorstore(chunks)
            
            assert result == mock_vectorstore
            assert store.vectorstore == mock_vectorstore
            assert len(store.structures_db) == 2
            assert 'CCO' in store.structures_db
            assert 'C6H6' in store.structures_db
    
    def test_enhance_chunk_text(self, sample_chemical_structure):
        """Test chunk text enhancement with chemical information"""
        store = ChemicalVectorStore(use_ollama_embeddings=False)
        
        chunk = {
            'text': 'Original text',
            'structures': [sample_chemical_structure]
        }
        
        enhanced_text = store._enhance_chunk_text(chunk)
        
        assert 'Original text' in enhanced_text
        assert '[Chemical Structures Found]' in enhanced_text
        assert 'CCO' in enhanced_text
        assert 'C2H6O' in enhanced_text
        assert '46.07' in enhanced_text
    
    def test_enhance_chunk_text_no_structures(self):
        """Test chunk text enhancement without structures"""
        store = ChemicalVectorStore(use_ollama_embeddings=False)
        
        chunk = {
            'text': 'Original text only',
            'structures': []
        }
        
        enhanced_text = store._enhance_chunk_text(chunk)
        
        assert enhanced_text == 'Original text only'
    
    def test_enhance_chunk_text_multiple_structures(self, sample_structures_list):
        """Test chunk text enhancement with multiple structures"""
        store = ChemicalVectorStore(use_ollama_embeddings=False)
        
        chunk = {
            'text': 'Text with multiple structures',
            'structures': sample_structures_list
        }
        
        enhanced_text = store._enhance_chunk_text(chunk)
        
        assert 'Structure 1:' in enhanced_text
        assert 'Structure 2:' in enhanced_text
        assert 'Structure 3:' in enhanced_text
        assert 'CCO' in enhanced_text
        assert 'C6H6' in enhanced_text
    
    def test_similarity_search_basic(self, mock_chroma_vectorstore):
        """Test basic similarity search"""
        store = ChemicalVectorStore(use_ollama_embeddings=False)
        store.vectorstore = mock_chroma_vectorstore
        
        results = store.similarity_search("ethanol molecule", k=3)
        
        mock_chroma_vectorstore.similarity_search.assert_called_once_with("ethanol molecule", k=3)
        assert len(results) == 2  # From mock
    
    def test_similarity_search_with_filter(self, mock_chroma_vectorstore):
        """Test similarity search with filter"""
        store = ChemicalVectorStore(use_ollama_embeddings=False)
        store.vectorstore = mock_chroma_vectorstore
        
        filter_dict = {'has_structures': True}
        results = store.similarity_search("ethanol", k=5, filter_dict=filter_dict)
        
        mock_chroma_vectorstore.similarity_search.assert_called_once_with(
            "ethanol", k=5, filter=filter_dict
        )
    
    def test_similarity_search_smiles_pattern(self, mock_chroma_vectorstore):
        """Test similarity search with SMILES pattern"""
        store = ChemicalVectorStore(use_ollama_embeddings=False)
        store.vectorstore = mock_chroma_vectorstore
        
        with patch.object(store, '_structure_aware_search') as mock_structure_search:
            mock_structure_search.return_value = [Mock()]
            
            # Query containing SMILES
            results = store.similarity_search("CCO molecule structure")
            
            mock_structure_search.assert_called_once()
    
    def test_structure_aware_search(self, mock_chroma_vectorstore):
        """Test structure-aware search functionality"""
        store = ChemicalVectorStore(use_ollama_embeddings=False)
        store.vectorstore = mock_chroma_vectorstore
        
        # Mock similarity_search to return different results for different calls
        def mock_similarity_search(query, k=5, filter=None):
            if filter:
                return [Mock(page_content="Filtered result")]
            else:
                return [Mock(page_content="Regular result")]
        
        mock_chroma_vectorstore.similarity_search.side_effect = mock_similarity_search
        
        results = store._structure_aware_search("CCO ethanol", k=3)
        
        assert len(results) >= 1
    
    def test_structure_aware_search_no_smiles(self, mock_chroma_vectorstore):
        """Test structure-aware search without SMILES"""
        store = ChemicalVectorStore(use_ollama_embeddings=False)
        store.vectorstore = mock_chroma_vectorstore
        
        results = store._structure_aware_search("general chemistry text", k=3)
        
        # Should fall back to regular search
        mock_chroma_vectorstore.similarity_search.assert_called()
    
    def test_create_vectorstore_metadata_generation(self):
        """Test metadata generation during vector store creation"""
        with patch('vectorstore.Chroma') as mock_chroma:
            mock_vectorstore = Mock()
            mock_chroma.from_texts.return_value = mock_vectorstore
            
            store = ChemicalVectorStore(use_ollama_embeddings=False)
            
            chunks = [
                {
                    'text': 'Test chunk',
                    'chunk_id': 0,
                    'page': 1,
                    'structures': [{'smiles': 'CCO'}]
                }
            ]
            
            store.create_vectorstore(chunks)
            
            # Check that from_texts was called with proper metadata
            call_args = mock_chroma.from_texts.call_args
            metadatas = call_args.kwargs['metadatas']
            
            assert len(metadatas) == 1
            metadata = metadatas[0]
            assert metadata['chunk_id'] == 0
            assert metadata['page'] == 1
            assert metadata['has_structures'] is True
            assert metadata['structure_count'] == 1
            assert metadata['structure_smiles'] == ['CCO']
    
    def test_create_vectorstore_persist_directory(self):
        """Test vector store creation with custom persist directory"""
        with patch('vectorstore.Chroma') as mock_chroma:
            store = ChemicalVectorStore(use_ollama_embeddings=False)
            
            chunks = [{'text': 'test', 'chunk_id': 0, 'structures': []}]
            custom_dir = "/custom/persist/dir"
            
            store.create_vectorstore(chunks, persist_directory=custom_dir)
            
            call_args = mock_chroma.from_texts.call_args
            assert call_args.kwargs['persist_directory'] == custom_dir
    
    def test_create_vectorstore_empty_chunks(self):
        """Test vector store creation with empty chunks list"""
        with patch('vectorstore.Chroma') as mock_chroma:
            mock_vectorstore = Mock()
            mock_chroma.from_texts.return_value = mock_vectorstore
            
            store = ChemicalVectorStore(use_ollama_embeddings=False)
            
            result = store.create_vectorstore([])
            
            # Should handle empty chunks gracefully
            call_args = mock_chroma.from_texts.call_args
            assert call_args.kwargs['texts'] == []
            assert call_args.kwargs['metadatas'] == []
    
    def test_enhance_chunk_text_structure_without_formula(self):
        """Test chunk enhancement with structure missing formula"""
        store = ChemicalVectorStore(use_ollama_embeddings=False)
        
        chunk = {
            'text': 'Test text',
            'structures': [{'smiles': 'CCO', 'molecular_weight': 46.07}]  # Missing formula
        }
        
        enhanced_text = store._enhance_chunk_text(chunk)
        
        assert 'CCO' in enhanced_text
        assert '46.07' in enhanced_text
        assert 'Formula=' not in enhanced_text  # Should skip missing formula
    
    def test_enhance_chunk_text_structure_without_mw(self):
        """Test chunk enhancement with structure missing molecular weight"""
        store = ChemicalVectorStore(use_ollama_embeddings=False)
        
        chunk = {
            'text': 'Test text',
            'structures': [{'smiles': 'CCO', 'formula': 'C2H6O'}]  # Missing MW
        }
        
        enhanced_text = store._enhance_chunk_text(chunk)
        
        assert 'CCO' in enhanced_text
        assert 'C2H6O' in enhanced_text
        assert 'MW=' not in enhanced_text  # Should skip missing MW
    
    def test_smiles_pattern_recognition(self):
        """Test SMILES pattern recognition in queries"""
        store = ChemicalVectorStore(use_ollama_embeddings=False)
        
        # Test various SMILES patterns
        smiles_queries = [
            "CCO ethanol",
            "C6H6 benzene ring",
            "CC(=O)O acetic acid",
            "c1ccccc1 aromatic",
            "C[C@H](N)C(=O)O amino acid"
        ]
        
        for query in smiles_queries:
            # Check if the pattern is detected (using regex from the code)
            import re
            pattern = r'[C,c][0-9A-Za-z@+\-\[\]\(\)\\=#$]+'
            has_smiles = bool(re.search(pattern, query))
            assert has_smiles, f"SMILES pattern not detected in: {query}"
    
    def test_vectorstore_with_duplicate_structures(self):
        """Test handling of duplicate structures in different chunks"""
        with patch('vectorstore.Chroma') as mock_chroma:
            mock_vectorstore = Mock()
            mock_chroma.from_texts.return_value = mock_vectorstore
            
            store = ChemicalVectorStore(use_ollama_embeddings=False)
            
            # Chunks with duplicate structures
            chunks = [
                {
                    'text': 'First mention of ethanol',
                    'chunk_id': 0,
                    'structures': [{'smiles': 'CCO', 'formula': 'C2H6O'}]
                },
                {
                    'text': 'Second mention of ethanol',
                    'chunk_id': 1,
                    'structures': [{'smiles': 'CCO', 'formula': 'C2H6O'}]  # Duplicate
                }
            ]
            
            store.create_vectorstore(chunks)
            
            # Should store the structure only once in structures_db
            assert len(store.structures_db) == 1
            assert 'CCO' in store.structures_db
    
    def test_error_handling_during_vectorstore_creation(self):
        """Test error handling during vector store creation"""
        with patch('vectorstore.Chroma') as mock_chroma:
            mock_chroma.from_texts.side_effect = Exception("Vector store creation failed")
            
            store = ChemicalVectorStore(use_ollama_embeddings=False)
            
            chunks = [{'text': 'test', 'chunk_id': 0, 'structures': []}]
            
            with pytest.raises(Exception):
                store.create_vectorstore(chunks)
    
    def test_similarity_search_error_handling(self):
        """Test error handling during similarity search"""
        store = ChemicalVectorStore(use_ollama_embeddings=False)
        
        mock_vectorstore = Mock()
        mock_vectorstore.similarity_search.side_effect = Exception("Search failed")
        store.vectorstore = mock_vectorstore
        
        with pytest.raises(Exception):
            store.similarity_search("test query")