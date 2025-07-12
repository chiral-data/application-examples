from langchain_community.llms import Ollama
from langchain_community.embeddings import OllamaEmbeddings
from langchain.chains import RetrievalQA
from langchain.prompts import PromptTemplate
from langchain.callbacks.streaming_stdout import StreamingStdOutCallbackHandler
import os

class ChemicalRAG:
    def __init__(self, vectorstore, model_name=None, streaming=True):
        if model_name is None:
            model_name = os.getenv("OLLAMA_MODEL", "llama3.2:latest")
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
        4. When referencing information, mention which document it comes from if available
        5. If you cannot answer based on the context, say so clearly
        6. Be concise but thorough
        
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
        """Process user query with source documents and document references"""
        result = self.qa_chain({"query": question})
        
        if return_sources:
            # Enhanced source information with document names
            sources = []
            for doc in result.get("source_documents", []):
                source_info = {
                    "content": doc.page_content,
                    "metadata": doc.metadata,
                    "document": doc.metadata.get("document", "Unknown"),
                    "page": doc.metadata.get("page", -1),
                    "has_structures": doc.metadata.get("has_structures", False)
                }
                sources.append(source_info)
            
            return {
                "answer": result["result"],
                "sources": sources
            }
        return result["result"]
    
    def _format_context_with_sources(self, docs):
        """Format retrieved documents with source attribution"""
        formatted_context = []
        for i, doc in enumerate(docs):
            document_name = doc.metadata.get("document", "Unknown Document")
            page = doc.metadata.get("page", "unknown")
            context_piece = f"[Source: {document_name}, Page: {page}]\n{doc.page_content}"
            formatted_context.append(context_piece)
        
        return "\n\n".join(formatted_context)