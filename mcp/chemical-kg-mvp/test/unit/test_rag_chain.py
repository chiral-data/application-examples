import pytest
from unittest.mock import Mock, patch, MagicMock
import os

from rag_chain import ChemicalRAG


class TestChemicalRAG:
    """Test suite for ChemicalRAG class"""
    
    def test_init_default(self, mock_chroma_vectorstore):
        """Test default initialization"""
        with patch('rag_chain.Ollama') as mock_ollama:
            with patch('rag_chain.RetrievalQA') as mock_retrieval_qa:
                mock_qa_chain = Mock()
                mock_retrieval_qa.from_chain_type.return_value = mock_qa_chain
                
                rag = ChemicalRAG(mock_chroma_vectorstore)
                
                assert rag.vectorstore == mock_chroma_vectorstore
                assert rag.llm is not None
                assert rag.prompt is not None
                assert rag.qa_chain == mock_qa_chain
                
                # Check Ollama was initialized with defaults
                mock_ollama.assert_called_once()
                call_args = mock_ollama.call_args
                assert call_args.kwargs['model'] == 'mistral'
                assert call_args.kwargs['temperature'] == 0.7
    
    def test_init_custom_model(self, mock_chroma_vectorstore):
        """Test initialization with custom model"""
        with patch('rag_chain.Ollama') as mock_ollama:
            with patch('rag_chain.RetrievalQA') as mock_retrieval_qa:
                mock_retrieval_qa.from_chain_type.return_value = Mock()
                
                rag = ChemicalRAG(mock_chroma_vectorstore, model_name="llama2", streaming=False)
                
                call_args = mock_ollama.call_args
                assert call_args.kwargs['model'] == 'llama2'
                assert call_args.kwargs['callbacks'] == []  # No streaming callbacks
    
    def test_init_environment_variables(self, mock_chroma_vectorstore):
        """Test initialization with environment variables"""
        with patch.dict(os.environ, {
            'OLLAMA_HOST': 'test-host',
            'OLLAMA_PORT': '9999'
        }):
            with patch('rag_chain.Ollama') as mock_ollama:
                with patch('rag_chain.RetrievalQA') as mock_retrieval_qa:
                    mock_retrieval_qa.from_chain_type.return_value = Mock()
                    
                    rag = ChemicalRAG(mock_chroma_vectorstore)
                    
                    call_args = mock_ollama.call_args
                    assert 'test-host:9999' in call_args.kwargs['base_url']
    
    def test_prompt_template_content(self, mock_chroma_vectorstore):
        """Test that prompt template contains chemistry-specific instructions"""
        with patch('rag_chain.Ollama'):
            with patch('rag_chain.RetrievalQA') as mock_retrieval_qa:
                mock_retrieval_qa.from_chain_type.return_value = Mock()
                
                rag = ChemicalRAG(mock_chroma_vectorstore)
                
                prompt_text = rag.prompt_template
                assert 'chemist' in prompt_text.lower()
                assert 'smiles' in prompt_text.lower()
                assert 'context' in prompt_text.lower()
                assert 'chemical structures' in prompt_text.lower()
    
    def test_query_with_source_documents(self, mock_chroma_vectorstore):
        """Test query with source document return"""
        with patch('rag_chain.Ollama'):
            with patch('rag_chain.RetrievalQA') as mock_retrieval_qa:
                mock_qa_chain = Mock()
                mock_qa_chain.return_value = {
                    'result': 'Test answer about ethanol',
                    'source_documents': [
                        Mock(page_content='Source 1'),
                        Mock(page_content='Source 2')
                    ]
                }
                mock_retrieval_qa.from_chain_type.return_value = mock_qa_chain
                
                rag = ChemicalRAG(mock_chroma_vectorstore)
                result = rag.query("What is ethanol?", return_sources=True)
                
                assert isinstance(result, dict)
                assert 'answer' in result
                assert 'sources' in result
                assert result['answer'] == 'Test answer about ethanol'
                assert len(result['sources']) == 2
    
    def test_query_without_source_documents(self, mock_chroma_vectorstore):
        """Test query without source document return"""
        with patch('rag_chain.Ollama'):
            with patch('rag_chain.RetrievalQA') as mock_retrieval_qa:
                mock_qa_chain = Mock()
                mock_qa_chain.return_value = {
                    'result': 'Test answer about benzene'
                }
                mock_retrieval_qa.from_chain_type.return_value = mock_qa_chain
                
                rag = ChemicalRAG(mock_chroma_vectorstore)
                result = rag.query("What is benzene?", return_sources=False)
                
                assert isinstance(result, str)
                assert result == 'Test answer about benzene'
    
    def test_qa_chain_configuration(self, mock_chroma_vectorstore):
        """Test QA chain configuration parameters"""
        with patch('rag_chain.Ollama') as mock_ollama:
            with patch('rag_chain.RetrievalQA') as mock_retrieval_qa:
                mock_retrieval_qa.from_chain_type.return_value = Mock()
                
                rag = ChemicalRAG(mock_chroma_vectorstore)
                
                # Check RetrievalQA was configured correctly
                call_args = mock_retrieval_qa.from_chain_type.call_args
                assert call_args.kwargs['chain_type'] == 'stuff'
                assert call_args.kwargs['return_source_documents'] is True
                assert 'verbose' in call_args.kwargs['chain_type_kwargs']
                assert call_args.kwargs['chain_type_kwargs']['verbose'] is True
    
    def test_retriever_configuration(self, mock_chroma_vectorstore):
        """Test retriever configuration"""
        with patch('rag_chain.Ollama'):
            with patch('rag_chain.RetrievalQA') as mock_retrieval_qa:
                mock_retrieval_qa.from_chain_type.return_value = Mock()
                
                rag = ChemicalRAG(mock_chroma_vectorstore)
                
                # Check that vectorstore as_retriever was called with correct parameters
                mock_chroma_vectorstore.as_retriever.assert_called_once()
                call_args = mock_chroma_vectorstore.as_retriever.call_args
                assert call_args.kwargs['search_kwargs']['k'] == 5
    
    def test_ollama_parameters(self, mock_chroma_vectorstore):
        """Test Ollama LLM parameters"""
        with patch('rag_chain.Ollama') as mock_ollama:
            with patch('rag_chain.RetrievalQA') as mock_retrieval_qa:
                mock_retrieval_qa.from_chain_type.return_value = Mock()
                
                rag = ChemicalRAG(mock_chroma_vectorstore)
                
                call_args = mock_ollama.call_args
                ollama_kwargs = call_args.kwargs
                
                assert ollama_kwargs['num_predict'] == 1024
                assert ollama_kwargs['top_k'] == 40
                assert ollama_kwargs['top_p'] == 0.9
                assert ollama_kwargs['repeat_penalty'] == 1.1
    
    def test_streaming_callbacks(self, mock_chroma_vectorstore):
        """Test streaming callback configuration"""
        with patch('rag_chain.Ollama') as mock_ollama:
            with patch('rag_chain.RetrievalQA') as mock_retrieval_qa:
                with patch('rag_chain.StreamingStdOutCallbackHandler') as mock_callback:
                    mock_retrieval_qa.from_chain_type.return_value = Mock()
                    
                    # Test with streaming enabled
                    rag = ChemicalRAG(mock_chroma_vectorstore, streaming=True)
                    
                    call_args = mock_ollama.call_args
                    assert len(call_args.kwargs['callbacks']) == 1
                    mock_callback.assert_called_once()
    
    def test_no_streaming_callbacks(self, mock_chroma_vectorstore):
        """Test no callbacks when streaming is disabled"""
        with patch('rag_chain.Ollama') as mock_ollama:
            with patch('rag_chain.RetrievalQA') as mock_retrieval_qa:
                mock_retrieval_qa.from_chain_type.return_value = Mock()
                
                # Test with streaming disabled
                rag = ChemicalRAG(mock_chroma_vectorstore, streaming=False)
                
                call_args = mock_ollama.call_args
                assert call_args.kwargs['callbacks'] == []
    
    def test_prompt_input_variables(self, mock_chroma_vectorstore):
        """Test prompt template input variables"""
        with patch('rag_chain.Ollama'):
            with patch('rag_chain.RetrievalQA') as mock_retrieval_qa:
                mock_retrieval_qa.from_chain_type.return_value = Mock()
                
                rag = ChemicalRAG(mock_chroma_vectorstore)
                
                assert rag.prompt.input_variables == ['context', 'question']
    
    def test_query_error_handling(self, mock_chroma_vectorstore):
        """Test error handling during query"""
        with patch('rag_chain.Ollama'):
            with patch('rag_chain.RetrievalQA') as mock_retrieval_qa:
                mock_qa_chain = Mock()
                mock_qa_chain.side_effect = Exception("Query processing failed")
                mock_retrieval_qa.from_chain_type.return_value = mock_qa_chain
                
                rag = ChemicalRAG(mock_chroma_vectorstore)
                
                with pytest.raises(Exception):
                    rag.query("What is ethanol?")
    
    def test_query_with_missing_source_documents(self, mock_chroma_vectorstore):
        """Test query when source documents are missing from response"""
        with patch('rag_chain.Ollama'):
            with patch('rag_chain.RetrievalQA') as mock_retrieval_qa:
                mock_qa_chain = Mock()
                mock_qa_chain.return_value = {
                    'result': 'Test answer'
                    # Missing 'source_documents' key
                }
                mock_retrieval_qa.from_chain_type.return_value = mock_qa_chain
                
                rag = ChemicalRAG(mock_chroma_vectorstore)
                result = rag.query("Test question", return_sources=True)
                
                assert result['answer'] == 'Test answer'
                assert result['sources'] == []  # Should default to empty list
    
    def test_chemistry_specific_prompt_instructions(self, mock_chroma_vectorstore):
        """Test that prompt contains chemistry-specific instructions"""
        with patch('rag_chain.Ollama'):
            with patch('rag_chain.RetrievalQA') as mock_retrieval_qa:
                mock_retrieval_qa.from_chain_type.return_value = Mock()
                
                rag = ChemicalRAG(mock_chroma_vectorstore)
                
                prompt_text = rag.prompt_template.lower()
                
                # Check for chemistry-specific instructions
                assert 'base your answer only on the provided context' in prompt_text
                assert 'smiles notation' in prompt_text or 'smiles notations' in prompt_text
                assert 'chemical structures' in prompt_text
                assert 'cannot answer based on the context' in prompt_text
                assert 'key features' in prompt_text
    
    def test_different_base_urls(self, mock_chroma_vectorstore):
        """Test different base URL configurations"""
        test_cases = [
            {'OLLAMA_HOST': 'localhost', 'OLLAMA_PORT': '11434'},
            {'OLLAMA_HOST': 'remote-host', 'OLLAMA_PORT': '8080'},
            {'OLLAMA_HOST': '192.168.1.100', 'OLLAMA_PORT': '9090'}
        ]
        
        for env_vars in test_cases:
            with patch.dict(os.environ, env_vars):
                with patch('rag_chain.Ollama') as mock_ollama:
                    with patch('rag_chain.RetrievalQA') as mock_retrieval_qa:
                        mock_retrieval_qa.from_chain_type.return_value = Mock()
                        
                        rag = ChemicalRAG(mock_chroma_vectorstore)
                        
                        call_args = mock_ollama.call_args
                        expected_url = f"http://{env_vars['OLLAMA_HOST']}:{env_vars['OLLAMA_PORT']}"
                        assert call_args.kwargs['base_url'] == expected_url
    
    def test_qa_chain_prompt_injection(self, mock_chroma_vectorstore):
        """Test that custom prompt is properly injected into QA chain"""
        with patch('rag_chain.Ollama'):
            with patch('rag_chain.RetrievalQA') as mock_retrieval_qa:
                mock_retrieval_qa.from_chain_type.return_value = Mock()
                
                rag = ChemicalRAG(mock_chroma_vectorstore)
                
                call_args = mock_retrieval_qa.from_chain_type.call_args
                chain_kwargs = call_args.kwargs['chain_type_kwargs']
                
                assert 'prompt' in chain_kwargs
                assert chain_kwargs['prompt'] == rag.prompt
    
    def test_multiple_queries_same_instance(self, mock_chroma_vectorstore):
        """Test multiple queries on the same RAG instance"""
        with patch('rag_chain.Ollama'):
            with patch('rag_chain.RetrievalQA') as mock_retrieval_qa:
                mock_qa_chain = Mock()
                mock_qa_chain.side_effect = [
                    {'result': 'Answer 1', 'source_documents': []},
                    {'result': 'Answer 2', 'source_documents': []},
                    {'result': 'Answer 3', 'source_documents': []}
                ]
                mock_retrieval_qa.from_chain_type.return_value = mock_qa_chain
                
                rag = ChemicalRAG(mock_chroma_vectorstore)
                
                result1 = rag.query("Question 1", return_sources=False)
                result2 = rag.query("Question 2", return_sources=False)
                result3 = rag.query("Question 3", return_sources=False)
                
                assert result1 == 'Answer 1'
                assert result2 == 'Answer 2'
                assert result3 == 'Answer 3'
                assert mock_qa_chain.call_count == 3