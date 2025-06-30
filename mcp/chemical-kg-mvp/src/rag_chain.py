from langchain.llms import Ollama
from langchain.embeddings import OllamaEmbeddings
from langchain.chains import RetrievalQA
from langchain.prompts import PromptTemplate
from langchain.callbacks.streaming_stdout import StreamingStdOutCallbackHandler
import os

class ChemicalRAG:
    def __init__(self, vectorstore, model_name="mistral", streaming=True):
        # Configure Ollama connection
        ollama_host = os.getenv("OLLAMA_HOST", "localhost")
        ollama_port = os.getenv("OLLAMA_PORT", "11434")
        base_url = f"http://{ollama_host}:{ollama_port}"
        
        # Initialize LLM with streaming support
        callbacks = [StreamingStdOutCallbackHandler()] if streaming else []
        
        self.llm = Ollama(
            model=model_name,
            base_url=base_url,
            temperature=0.7,
            callbacks=callbacks,
            # Model-specific parameters
            num_predict=1024,  # Max tokens to generate
            top_k=40,
            top_p=0.9,
            repeat_penalty=1.1
        )
        
        self.vectorstore = vectorstore
        
        # Chemistry-specific prompt template
        self.prompt_template = """You are an expert chemist analyzing scientific papers. 
        Use the following context to answer the question accurately and comprehensively.
        If the context mentions chemical structures, include their SMILES notation and explain their significance.
        
        Context: {context}
        
        Question: {question}
        
        Instructions:
        1. Base your answer ONLY on the provided context
        2. If chemical structures are mentioned, describe their key features
        3. Include SMILES notations when available
        4. If you cannot answer based on the context, say so clearly
        5. Be concise but thorough
        
        Answer:"""
        
        self.prompt = PromptTemplate(
            template=self.prompt_template,
            input_variables=["context", "question"]
        )
        
        # Create QA chain
        self.qa_chain = RetrievalQA.from_chain_type(
            llm=self.llm,
            chain_type="stuff",
            retriever=self.vectorstore.vectorstore.as_retriever(
                search_kwargs={"k": 5}
            ),
            chain_type_kwargs={
                "prompt": self.prompt,
                "verbose": True
            },
            return_source_documents=True
        )
    
    def query(self, question, return_sources=True):
        """Process user query with source documents"""
        result = self.qa_chain({"query": question})
        
        if return_sources:
            return {
                "answer": result["result"],
                "sources": result.get("source_documents", [])
            }
        return result["result"]