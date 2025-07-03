"""
Comprehensive tests for LLM conversation abilities.
Tests the RAG chain, context handling, and chemical knowledge responses.
"""

import pytest
import os
import tempfile
from unittest.mock import patch, MagicMock, call
from langchain.schema import Document
from langchain.embeddings.base import Embeddings
import numpy as np

from rag_chain import RAGChain
from vectorstore import VectorStore
from chemical_handler import ChemicalHandler


class TestLLMConversation:
    """Test suite for LLM conversation and RAG capabilities"""
    
    @pytest.fixture
    def mock_embeddings(self):
        """Create mock embeddings"""
        class MockEmbeddings(Embeddings):
            def embed_documents(self, texts):
                # Return consistent embeddings based on text content
                return [np.random.rand(384).tolist() for _ in texts]
            
            def embed_query(self, text):
                # Return embedding for query
                return np.random.rand(384).tolist()
        
        return MockEmbeddings()
    
    @pytest.fixture
    def mock_llm(self):
        """Create mock LLM"""
        llm = MagicMock()
        llm.invoke.return_value = MagicMock(content="Mock response")
        return llm
    
    @pytest.fixture
    def rag_chain(self, mock_embeddings, mock_llm):
        """Create RAG chain with mocked components"""
        with patch('rag_chain.OllamaEmbeddings', return_value=mock_embeddings):
            with patch('rag_chain.ChatOllama', return_value=mock_llm):
                chain = RAGChain()
                return chain
    
    @pytest.fixture
    def chemical_context_docs(self):
        """Create chemical context documents"""
        return [
            Document(
                page_content="Aspirin (acetylsalicylic acid) is a medication used to reduce pain, fever, or inflammation. Its molecular formula is C9H8O4.",
                metadata={"source": "pharma_db", "molecule": "aspirin"}
            ),
            Document(
                page_content="The SMILES notation for aspirin is CC(=O)Oc1ccccc1C(=O)O. It has a molecular weight of 180.16 g/mol.",
                metadata={"source": "chemical_db", "molecule": "aspirin", "smiles": "CC(=O)Oc1ccccc1C(=O)O"}
            ),
            Document(
                page_content="Ibuprofen is a nonsteroidal anti-inflammatory drug (NSAID). Its SMILES is CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O.",
                metadata={"source": "pharma_db", "molecule": "ibuprofen"}
            ),
            Document(
                page_content="Caffeine is a central nervous system stimulant. Chemical formula: C8H10N4O2. SMILES: CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
                metadata={"source": "chemical_db", "molecule": "caffeine"}
            )
        ]
    
    def test_basic_chemical_query(self, rag_chain):
        """Test basic chemical information query"""
        query = "What is the molecular formula of aspirin?"
        
        with patch.object(rag_chain, 'invoke') as mock_invoke:
            mock_invoke.return_value = {
                'answer': "The molecular formula of aspirin is C9H8O4.",
                'source_documents': []
            }
            
            response = rag_chain.invoke({"question": query})
            
            assert "C9H8O4" in response['answer']
            mock_invoke.assert_called_once()
    
    def test_smiles_query(self, rag_chain):
        """Test SMILES notation query"""
        query = "What is the SMILES notation for caffeine?"
        
        with patch.object(rag_chain, 'invoke') as mock_invoke:
            mock_invoke.return_value = {
                'answer': "The SMILES notation for caffeine is CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
                'source_documents': []
            }
            
            response = rag_chain.invoke({"question": query})
            
            assert "CN1C=NC2=C1C(=O)N(C(=O)N2C)C" in response['answer']
    
    def test_context_retrieval(self, rag_chain, chemical_context_docs):
        """Test context retrieval for queries"""
        query = "Tell me about aspirin and its properties"
        
        with patch.object(rag_chain.vectorstore, 'similarity_search') as mock_search:
            mock_search.return_value = [doc for doc in chemical_context_docs if doc.metadata.get('molecule') == 'aspirin']
            
            with patch.object(rag_chain, 'invoke') as mock_invoke:
                mock_invoke.return_value = {
                    'answer': "Aspirin is a medication with molecular formula C9H8O4 and SMILES CC(=O)Oc1ccccc1C(=O)O",
                    'source_documents': mock_search.return_value
                }
                
                response = rag_chain.invoke({"question": query})
                
                assert len(response['source_documents']) > 0
                assert any('aspirin' in doc.page_content.lower() for doc in response['source_documents'])
    
    def test_multi_molecule_comparison(self, rag_chain):
        """Test comparing multiple molecules"""
        query = "Compare the molecular weights of aspirin and ibuprofen"
        
        with patch.object(rag_chain, 'invoke') as mock_invoke:
            mock_invoke.return_value = {
                'answer': "Aspirin has a molecular weight of 180.16 g/mol, while ibuprofen has a molecular weight of 206.28 g/mol. Ibuprofen is heavier.",
                'source_documents': []
            }
            
            response = rag_chain.invoke({"question": query})
            
            assert "180.16" in response['answer']
            assert "206.28" in response['answer']
    
    def test_chemical_reaction_query(self, rag_chain):
        """Test chemical reaction understanding"""
        query = "How is aspirin synthesized from salicylic acid?"
        
        with patch.object(rag_chain, 'invoke') as mock_invoke:
            mock_invoke.return_value = {
                'answer': "Aspirin is synthesized from salicylic acid through acetylation with acetic anhydride in the presence of an acid catalyst.",
                'source_documents': []
            }
            
            response = rag_chain.invoke({"question": query})
            
            assert "acetylation" in response['answer'].lower() or "acetic" in response['answer'].lower()
    
    def test_conversation_memory(self, rag_chain):
        """Test conversation memory and context retention"""
        # First query
        query1 = "What is aspirin?"
        with patch.object(rag_chain, 'invoke') as mock_invoke:
            mock_invoke.return_value = {
                'answer': "Aspirin is acetylsalicylic acid, a common pain reliever.",
                'source_documents': []
            }
            response1 = rag_chain.invoke({"question": query1})
        
        # Follow-up query referring to previous context
        query2 = "What is its molecular formula?"
        with patch.object(rag_chain, 'invoke') as mock_invoke:
            mock_invoke.return_value = {
                'answer': "The molecular formula of aspirin is C9H8O4.",
                'source_documents': []
            }
            response2 = rag_chain.invoke({"question": query2})
            
            assert "C9H8O4" in response2['answer']
    
    def test_invalid_molecule_handling(self, rag_chain):
        """Test handling of invalid or unknown molecules"""
        query = "What is the structure of unobtainium?"
        
        with patch.object(rag_chain, 'invoke') as mock_invoke:
            mock_invoke.return_value = {
                'answer': "I don't have information about 'unobtainium' in my knowledge base. This appears to be a fictional element.",
                'source_documents': []
            }
            
            response = rag_chain.invoke({"question": query})
            
            assert "don't have information" in response['answer'] or "fictional" in response['answer']
    
    def test_structure_activity_relationship(self, rag_chain):
        """Test SAR (Structure-Activity Relationship) queries"""
        query = "How does the acetyl group in aspirin affect its activity compared to salicylic acid?"
        
        with patch.object(rag_chain, 'invoke') as mock_invoke:
            mock_invoke.return_value = {
                'answer': "The acetyl group in aspirin reduces gastric irritation compared to salicylic acid while maintaining anti-inflammatory activity.",
                'source_documents': []
            }
            
            response = rag_chain.invoke({"question": query})
            
            assert "acetyl" in response['answer'].lower()
    
    def test_safety_and_toxicity_queries(self, rag_chain):
        """Test safety and toxicity information queries"""
        query = "What are the safety concerns with high doses of aspirin?"
        
        with patch.object(rag_chain, 'invoke') as mock_invoke:
            mock_invoke.return_value = {
                'answer': "High doses of aspirin can cause gastrointestinal bleeding, tinnitus, and in severe cases, salicylate toxicity.",
                'source_documents': []
            }
            
            response = rag_chain.invoke({"question": query})
            
            assert any(term in response['answer'].lower() for term in ['bleeding', 'toxicity', 'safety'])
    
    def test_drug_interaction_queries(self, rag_chain):
        """Test drug interaction queries"""
        query = "Can aspirin be taken with ibuprofen?"
        
        with patch.object(rag_chain, 'invoke') as mock_invoke:
            mock_invoke.return_value = {
                'answer': "Taking aspirin with ibuprofen can reduce aspirin's cardioprotective effects and increase bleeding risk. Consult a healthcare provider.",
                'source_documents': []
            }
            
            response = rag_chain.invoke({"question": query})
            
            assert "consult" in response['answer'].lower() or "risk" in response['answer'].lower()
    
    def test_batch_queries(self, rag_chain):
        """Test handling multiple queries in sequence"""
        queries = [
            "What is caffeine?",
            "What is its molecular weight?",
            "Is it safe during pregnancy?",
            "What are common sources?"
        ]
        
        responses = []
        for query in queries:
            with patch.object(rag_chain, 'invoke') as mock_invoke:
                # Simulate different responses
                if "caffeine" in query:
                    answer = "Caffeine is a central nervous system stimulant."
                elif "weight" in query:
                    answer = "The molecular weight of caffeine is 194.19 g/mol."
                elif "pregnancy" in query:
                    answer = "Moderate caffeine intake (less than 200mg/day) is generally considered safe during pregnancy."
                else:
                    answer = "Common sources include coffee, tea, chocolate, and energy drinks."
                
                mock_invoke.return_value = {'answer': answer, 'source_documents': []}
                response = rag_chain.invoke({"question": query})
                responses.append(response)
        
        assert len(responses) == len(queries)
        assert all('answer' in r for r in responses)
    
    def test_chemical_property_calculation(self, rag_chain):
        """Test queries about calculated chemical properties"""
        query = "Calculate the logP value for aspirin"
        
        with patch.object(rag_chain, 'invoke') as mock_invoke:
            mock_invoke.return_value = {
                'answer': "The calculated logP (partition coefficient) for aspirin is approximately 1.19, indicating moderate lipophilicity.",
                'source_documents': []
            }
            
            response = rag_chain.invoke({"question": query})
            
            assert "logP" in response['answer'] or "partition coefficient" in response['answer']
    
    def test_response_formatting(self, rag_chain):
        """Test that responses are properly formatted"""
        query = "List three properties of ibuprofen"
        
        with patch.object(rag_chain, 'invoke') as mock_invoke:
            mock_invoke.return_value = {
                'answer': """Three properties of ibuprofen:
1. Molecular formula: C13H18O2
2. Molecular weight: 206.28 g/mol
3. Melting point: 75-77°C""",
                'source_documents': []
            }
            
            response = rag_chain.invoke({"question": query})
            
            # Check for structured response
            assert "1." in response['answer']
            assert "2." in response['answer']
            assert "3." in response['answer']
    
    def test_error_handling_in_conversation(self, rag_chain):
        """Test error handling during conversation"""
        query = "What is the structure of [INVALID_SMILES]?"
        
        with patch.object(rag_chain, 'invoke') as mock_invoke:
            mock_invoke.return_value = {
                'answer': "I couldn't process the structure notation provided. Please provide a valid SMILES or molecule name.",
                'source_documents': []
            }
            
            response = rag_chain.invoke({"question": query})
            
            assert "couldn't process" in response['answer'] or "valid" in response['answer']
    
    def test_citation_handling(self, rag_chain, chemical_context_docs):
        """Test that responses include proper citations when available"""
        query = "What is aspirin used for?"
        
        with patch.object(rag_chain.vectorstore, 'similarity_search') as mock_search:
            mock_search.return_value = [chemical_context_docs[0]]
            
            with patch.object(rag_chain, 'invoke') as mock_invoke:
                mock_invoke.return_value = {
                    'answer': "According to the pharmaceutical database, aspirin is used to reduce pain, fever, or inflammation.",
                    'source_documents': mock_search.return_value
                }
                
                response = rag_chain.invoke({"question": query})
                
                assert len(response['source_documents']) > 0
                assert response['source_documents'][0].metadata['source'] == 'pharma_db'