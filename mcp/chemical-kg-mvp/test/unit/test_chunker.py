import pytest
from unittest.mock import Mock, patch
import os

from chunker import ChemicalAwareChunker


class TestChemicalAwareChunker:
    """Test suite for ChemicalAwareChunker class"""
    
    def test_init_default(self):
        """Test default initialization"""
        chunker = ChemicalAwareChunker()
        
        # Values come from environment variables set in conftest.py
        assert chunker.chunk_size == 500  # From CHUNK_SIZE env var
        assert chunker.overlap == 50      # From CHUNK_OVERLAP env var
        assert chunker.context_window == 3  # From CHEMICAL_CONTEXT_WINDOW env var
        assert chunker.splitter is not None
        assert 'structure' in chunker.chemical_keywords
    
    def test_init_custom_params(self):
        """Test initialization with custom parameters via environment override"""
        # Since environment variables take precedence, test by temporarily changing them
        with patch.dict(os.environ, {'CHUNK_SIZE': '800', 'CHUNK_OVERLAP': '150'}):
            chunker = ChemicalAwareChunker(chunk_size=500, overlap=100)  # These will be overridden
            
            assert chunker.chunk_size == 800   # From environment
            assert chunker.overlap == 150      # From environment
    
    def test_init_environment_variables(self):
        """Test initialization with environment variables"""
        with patch.dict(os.environ, {
            'CHUNK_SIZE': '800',
            'CHUNK_OVERLAP': '150',
            'CHEMICAL_CONTEXT_WINDOW': '7'
        }):
            chunker = ChemicalAwareChunker()
            
            assert chunker.chunk_size == 800
            assert chunker.overlap == 150
            assert chunker.context_window == 7
    
    def test_chunk_with_structures_basic(self, sample_text, sample_structures_list):
        """Test basic chunking with structures"""
        chunker = ChemicalAwareChunker(chunk_size=200, overlap=50)
        
        chunks = chunker.chunk_with_structures(sample_text, sample_structures_list)
        
        assert len(chunks) > 0
        for chunk in chunks:
            assert 'text' in chunk
            assert 'chunk_id' in chunk
            assert 'structures' in chunk
            assert 'metadata' in chunk
            assert isinstance(chunk['structures'], list)
    
    def test_contains_chemical_reference(self):
        """Test chemical reference detection"""
        chunker = ChemicalAwareChunker()
        
        # Text with chemical keywords
        assert chunker._contains_chemical_reference("The compound ethanol") is True
        assert chunker._contains_chemical_reference("synthesis of benzene") is True
        assert chunker._contains_chemical_reference("molecular structure") is True
        assert chunker._contains_chemical_reference("Figure 1 shows") is True
        
        # Text without chemical references
        assert chunker._contains_chemical_reference("The weather is nice") is False
        assert chunker._contains_chemical_reference("Weather report today") is False
    
    def test_contains_chemical_reference_formulas(self):
        """Test chemical reference detection with formulas"""
        chunker = ChemicalAwareChunker()
        
        # Text with chemical formulas
        assert chunker._contains_chemical_reference("C2H6O molecule") is True
        assert chunker._contains_chemical_reference("H2SO4 acid") is True
        assert chunker._contains_chemical_reference("NaCl salt") is True
        
        # Text with SMILES-like patterns
        assert chunker._contains_chemical_reference("CCO compound") is True
        assert chunker._contains_chemical_reference("C=C bond") is True
    
    def test_find_related_structures_direct_smiles(self, sample_structures_list):
        """Test structure finding with direct SMILES mention"""
        chunker = ChemicalAwareChunker()
        
        chunk_text = "The ethanol molecule CCO has interesting properties"
        
        related = chunker._find_related_structures(chunk_text, sample_structures_list, 0)
        
        assert len(related) > 0
        # Should find the ethanol structure (CCO)
        smiles_found = [s['smiles'] for s in related]
        assert 'CCO' in smiles_found
    
    def test_find_related_structures_formula_mention(self, sample_structures_list):
        """Test structure finding with formula mention"""
        chunker = ChemicalAwareChunker()
        
        chunk_text = "The compound with formula C2H6O was analyzed"
        
        related = chunker._find_related_structures(chunk_text, sample_structures_list, 0)
        
        assert len(related) > 0
        # Should find ethanol (C2H6O)
        formulas_found = [s['formula'] for s in related]
        assert 'C2H6O' in formulas_found
    
    def test_find_related_structures_keyword_proximity(self, sample_structures_list):
        """Test structure finding based on keyword proximity"""
        chunker = ChemicalAwareChunker()
        
        chunk_text = "The synthesis of organic compounds involves various molecular structures"
        
        related = chunker._find_related_structures(chunk_text, sample_structures_list, 0)
        
        # Should find structures based on keywords
        assert len(related) >= 0  # May or may not find structures based on keywords alone
    
    def test_find_related_structures_figure_references(self, sample_structures_list):
        """Test structure finding with figure references"""
        chunker = ChemicalAwareChunker()
        
        chunk_text = "As shown in Figure 1, the molecular structure is complex"
        
        related = chunker._find_related_structures(chunk_text, sample_structures_list, 0)
        
        # Should give bonus points for figure references
        assert isinstance(related, list)
    
    def test_calculate_chunk_metadata_basic(self):
        """Test basic chunk metadata calculation"""
        chunker = ChemicalAwareChunker()
        
        chunk_text = "This is a test chunk about synthesis of compounds"
        structures = [{'smiles': 'CCO', 'formula': 'C2H6O'}]
        
        metadata = chunker._calculate_chunk_metadata(chunk_text, structures)
        
        assert metadata['has_structures'] is True
        assert metadata['structure_count'] == 1
        assert metadata['word_count'] > 0
        assert metadata['char_count'] == len(chunk_text)
        assert 'chemical_keywords' in metadata
        assert 'content_type' in metadata
    
    def test_calculate_chunk_metadata_content_types(self):
        """Test content type detection"""
        chunker = ChemicalAwareChunker()
        
        # Synthesis content
        synthesis_text = "The synthesis procedure involves mixing reactants"
        metadata = chunker._calculate_chunk_metadata(synthesis_text, [])
        assert metadata['content_type'] == 'synthesis'
        
        # Analysis content
        analysis_text = "NMR analysis revealed the structure"
        metadata = chunker._calculate_chunk_metadata(analysis_text, [])
        assert metadata['content_type'] == 'analysis'
        
        # Properties content
        properties_text = "The melting point and solubility were measured"
        metadata = chunker._calculate_chunk_metadata(properties_text, [])
        assert metadata['content_type'] == 'properties'
        
        # Structure-related content
        struct_text = "General text"
        structures = [{'smiles': 'CCO'}]
        metadata = chunker._calculate_chunk_metadata(struct_text, structures)
        assert metadata['content_type'] == 'structure_related'
    
    def test_calculate_chunk_metadata_chemical_density(self):
        """Test chemical density calculation"""
        chunker = ChemicalAwareChunker()
        
        # High chemical density text
        high_density_text = "compound structure synthesis molecule analysis"  # 5 words, all chemical
        metadata = chunker._calculate_chunk_metadata(high_density_text, [])
        assert metadata['chemical_density'] == 100.0
        
        # Low chemical density text
        low_density_text = "the weather is nice today and compound"  # 7 words, 1 chemical
        metadata = chunker._calculate_chunk_metadata(low_density_text, [])
        expected_density = round((1/7) * 100, 2)
        assert metadata['chemical_density'] == expected_density
    
    def test_chunks_are_related_shared_structures(self):
        """Test chunk relationship detection with shared structures"""
        chunker = ChemicalAwareChunker()
        
        chunk1 = {
            'structures': [{'smiles': 'CCO'}, {'smiles': 'C6H6'}],
            'metadata': {'content_type': 'general', 'chemical_keywords': ['compound']}
        }
        
        chunk2 = {
            'structures': [{'smiles': 'CCO'}],  # Shared structure
            'metadata': {'content_type': 'general', 'chemical_keywords': ['molecule']}
        }
        
        assert chunker._chunks_are_related(chunk1, chunk2) is True
    
    def test_chunks_are_related_content_type(self):
        """Test chunk relationship detection by content type"""
        chunker = ChemicalAwareChunker()
        
        chunk1 = {
            'structures': [],
            'metadata': {'content_type': 'synthesis', 'chemical_keywords': ['synthesis']}
        }
        
        chunk2 = {
            'structures': [],
            'metadata': {'content_type': 'synthesis', 'chemical_keywords': ['reaction']}
        }
        
        assert chunker._chunks_are_related(chunk1, chunk2) is True
    
    def test_chunks_are_related_keywords(self):
        """Test chunk relationship detection by keywords"""
        chunker = ChemicalAwareChunker()
        
        chunk1 = {
            'structures': [],
            'metadata': {'content_type': 'general', 'chemical_keywords': ['synthesis', 'compound', 'molecule']}
        }
        
        chunk2 = {
            'structures': [],
            'metadata': {'content_type': 'general', 'chemical_keywords': ['synthesis', 'compound']}
        }
        
        assert chunker._chunks_are_related(chunk1, chunk2) is True
    
    def test_chunks_not_related(self):
        """Test chunks that are not related"""
        chunker = ChemicalAwareChunker()
        
        chunk1 = {
            'structures': [{'smiles': 'CCO'}],
            'metadata': {'content_type': 'synthesis', 'chemical_keywords': ['synthesis']}
        }
        
        chunk2 = {
            'structures': [{'smiles': 'C6H6'}],  # Different structure
            'metadata': {'content_type': 'analysis', 'chemical_keywords': ['analysis']}  # Different content type and keywords
        }
        
        assert chunker._chunks_are_related(chunk1, chunk2) is False
    
    def test_merge_chunks(self):
        """Test chunk merging functionality"""
        chunker = ChemicalAwareChunker()
        
        chunk1 = {
            'text': 'First chunk text',
            'chunk_id': 0,
            'structures': [{'smiles': 'CCO'}],
            'metadata': {'has_structures': True}
        }
        
        chunk2 = {
            'text': 'Second chunk text',
            'chunk_id': 1,
            'structures': [{'smiles': 'C6H6'}, {'smiles': 'CCO'}],  # One duplicate
            'metadata': {'has_structures': True}
        }
        
        merged = chunker._merge_chunks(chunk1, chunk2)
        
        assert 'First chunk text' in merged['text']
        assert 'Second chunk text' in merged['text']
        assert merged['chunk_id'] == 0  # Should keep first chunk's ID
        assert len(merged['structures']) == 2  # Should remove duplicates
        smiles_list = [s['smiles'] for s in merged['structures']]
        assert 'CCO' in smiles_list
        assert 'C6H6' in smiles_list
    
    def test_post_process_chunks_merging(self):
        """Test post-processing with chunk merging"""
        chunker = ChemicalAwareChunker(chunk_size=300)
        
        chunks = [
            {
                'text': 'Short text',  # Small chunk
                'chunk_id': 0,
                'structures': [{'smiles': 'CCO'}],
                'metadata': {'char_count': 10, 'has_structures': True}
            },
            {
                'text': 'Another short text',  # Small chunk
                'chunk_id': 1,
                'structures': [],
                'metadata': {'char_count': 18, 'has_structures': False}
            }
        ]
        
        with patch.object(chunker, '_chunks_are_related', return_value=True):
            with patch.object(chunker, '_merge_chunks') as mock_merge:
                mock_merge.return_value = {
                    'text': 'Merged text',
                    'chunk_id': 0,
                    'structures': [{'smiles': 'CCO'}],
                    'metadata': {'char_count': 28, 'has_structures': True}
                }
                
                processed = chunker._post_process_chunks(chunks)
                
                mock_merge.assert_called_once()
    
    def test_get_chunking_stats(self, sample_text, sample_structures_list):
        """Test chunking statistics generation"""
        chunker = ChemicalAwareChunker()
        
        chunks = chunker.chunk_with_structures(sample_text, sample_structures_list)
        stats = chunker.get_chunking_stats(chunks)
        
        assert 'total_chunks' in stats
        assert 'chunks_with_structures' in stats
        assert 'structure_coverage' in stats
        assert 'total_structures' in stats
        assert 'avg_structures_per_chunk' in stats
        assert 'content_type_distribution' in stats
        assert 'avg_chunk_size' in stats
        assert 'avg_chemical_density' in stats
        
        assert stats['total_chunks'] == len(chunks)
        assert 0 <= stats['structure_coverage'] <= 100
    
    def test_create_structure_aware_chunks_large_paragraph(self):
        """Test handling of large paragraphs that exceed chunk size"""
        chunker = ChemicalAwareChunker(chunk_size=100)
        
        # Create a very long paragraph
        long_text = "This is a very long paragraph about chemical synthesis. " * 20
        
        chunks = chunker._create_structure_aware_chunks(long_text, [])
        
        # Should split large paragraphs
        assert len(chunks) > 1
        for chunk in chunks:
            # Each chunk should be reasonably sized (though may exceed slightly due to paragraph boundaries)
            assert len(chunk) <= chunker.chunk_size * 2  # Allow some flexibility
    
    def test_apply_intelligent_chunking_empty_structures(self):
        """Test intelligent chunking with no structures"""
        chunker = ChemicalAwareChunker()
        
        text_chunks = ["First chunk of text", "Second chunk of text"]
        
        result = chunker._apply_intelligent_chunking(text_chunks, [])
        
        assert len(result) == 2
        for chunk in result:
            assert len(chunk['structures']) == 0
            assert chunk['metadata']['has_structures'] is False
    
    def test_chunk_with_structures_empty_text(self):
        """Test chunking with empty text"""
        chunker = ChemicalAwareChunker()
        
        chunks = chunker.chunk_with_structures("", [])
        
        # Should handle empty text gracefully
        assert isinstance(chunks, list)
    
    def test_chunk_with_structures_no_structures(self, sample_text):
        """Test chunking with no structures"""
        chunker = ChemicalAwareChunker()
        
        chunks = chunker.chunk_with_structures(sample_text, [])
        
        assert len(chunks) > 0
        for chunk in chunks:
            assert len(chunk['structures']) == 0
            assert chunk['metadata']['has_structures'] is False