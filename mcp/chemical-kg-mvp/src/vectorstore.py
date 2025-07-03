# vectorstore.py - Chemical vector store using FAISS (Python 3.8 compatible)
from __future__ import annotations

from langchain_community.embeddings import OllamaEmbeddings
from langchain_community.vectorstores import FAISS
from langchain.text_splitter import RecursiveCharacterTextSplitter
import os
import pickle
import re

class ChemicalVectorStore:
    def __init__(self, use_ollama_embeddings=True):
        # Use Ollama embeddings to avoid sentence-transformers dependency
        try:
            ollama_host = os.getenv("OLLAMA_HOST", "localhost")
            ollama_port = os.getenv("OLLAMA_PORT", "11434")
            ollama_url = f"http://{ollama_host}:{ollama_port}"
            
            self.embeddings = OllamaEmbeddings(
                base_url=ollama_url,
                model=os.getenv("EMBEDDING_MODEL", "nomic-embed-text")
            )
            print(f"Using Ollama embeddings at {ollama_url}")
        except Exception as e:
            print(f"Failed to initialize Ollama embeddings: {e}")
            # Fallback to a simple hash-based embedding for development
            from langchain_community.embeddings import FakeEmbeddings
            self.embeddings = FakeEmbeddings(size=384)
            print("Using fallback fake embeddings")
        
        self.vectorstore = None
        self.structures_db = {}
        
    def create_vectorstore(self, chunks, persist_directory="./faiss_db"):
        """Create vector store from chunks using FAISS"""
        texts = []
        metadatas = []
        
        for chunk in chunks:
            # Enhance text with structure information
            enhanced_text = self._enhance_chunk_text(chunk)
            texts.append(enhanced_text)
            
            metadata = {
                'chunk_id': chunk['chunk_id'],
                'page': chunk.get('page', -1),
                'has_structures': len(chunk.get('structures', [])) > 0,
                'structure_count': len(chunk.get('structures', []))
            }
            
            # Store structure details separately
            if chunk.get('structures'):
                structure_smiles = [s['smiles'] for s in chunk['structures'] if s.get('smiles')]
                metadata['structure_smiles'] = ', '.join(structure_smiles[:3])  # First 3 SMILES
                
                # Store full structure data
                self.structures_db[chunk['chunk_id']] = chunk['structures']
            
            metadatas.append(metadata)
        
        # Create FAISS vector store
        self.vectorstore = FAISS.from_texts(
            texts=texts,
            embedding=self.embeddings,
            metadatas=metadatas
        )
        
        # Save to disk
        os.makedirs(persist_directory, exist_ok=True)
        self.vectorstore.save_local(persist_directory)
        
        # Save structures database
        with open(os.path.join(persist_directory, "structures.pkl"), "wb") as f:
            pickle.dump(self.structures_db, f)
        
        return self.vectorstore
    
    def _enhance_chunk_text(self, chunk):
        """Enhance chunk text with structure information"""
        enhanced = chunk['text']
        
        if chunk.get('structures'):
            enhanced += "\n\nChemical Structures Found:\n"
            for struct in chunk['structures']:
                if struct.get('smiles'):
                    enhanced += f"- SMILES: {struct['smiles']}\n"
                if struct.get('formula'):
                    enhanced += f"  Formula: {struct['formula']}\n"
                if struct.get('molecular_weight'):
                    enhanced += f"  MW: {struct['molecular_weight']:.2f}\n"
        
        return enhanced
    
    def similarity_search(self, query, k=4):
        """Search for similar chunks"""
        if not self.vectorstore:
            return []
        
        return self.vectorstore.similarity_search(query, k=k)
    
    def get_structures_for_chunk(self, chunk_id):
        """Get structures associated with a chunk"""
        return self.structures_db.get(chunk_id, [])
    
    def load_from_disk(self, persist_directory="./faiss_db"):
        """Load vector store from disk"""
        self.vectorstore = FAISS.load_local(persist_directory, self.embeddings)
        
        # Load structures database
        structures_path = os.path.join(persist_directory, "structures.pkl")
        if os.path.exists(structures_path):
            with open(structures_path, "rb") as f:
                self.structures_db = pickle.load(f)