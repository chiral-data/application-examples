from langchain_community.embeddings import OllamaEmbeddings, HuggingFaceEmbeddings
from langchain_community.vectorstores import Chroma
from langchain.text_splitter import RecursiveCharacterTextSplitter
import os
import re

class ChemicalVectorStore:
    def __init__(self, use_ollama_embeddings=True):
        if use_ollama_embeddings:
            # Use Ollama embeddings
            ollama_host = os.getenv("OLLAMA_HOST", "localhost")
            ollama_port = os.getenv("OLLAMA_PORT", "11434")
            base_url = f"http://{ollama_host}:{ollama_port}"
            
            self.embeddings = OllamaEmbeddings(
                model="nomic-embed-text",  # Optimized for embeddings
                base_url=base_url
            )
        else:
            # Fallback to local HuggingFace embeddings
            self.embeddings = HuggingFaceEmbeddings(
                model_name="sentence-transformers/all-MiniLM-L6-v2",
                model_kwargs={'device': 'cpu'},
                encode_kwargs={'normalize_embeddings': True}
            )
        
        self.vectorstore = None
        self.structures_db = {}
        
    def create_vectorstore(self, chunks, persist_directory="./chroma_db"):
        """Create vector store from chunks"""
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
                'structure_count': len(chunk.get('structures', [])),
                'structure_smiles': ','.join([s['smiles'] for s in chunk.get('structures', [])])
            }
            metadatas.append(metadata)
            
            # Store structures separately
            for struct in chunk.get('structures', []):
                self.structures_db[struct['smiles']] = struct
        
        # Create ChromaDB instance
        self.vectorstore = Chroma.from_texts(
            texts=texts,
            embedding=self.embeddings,
            metadatas=metadatas,
            persist_directory=persist_directory
        )
        
        return self.vectorstore
    
    def _enhance_chunk_text(self, chunk):
        """Enhance chunk text with chemical information"""
        text = chunk['text']
        
        if chunk.get('structures'):
            text += "\n\n[Chemical Structures Found]:\n"
            for i, struct in enumerate(chunk['structures']):
                text += f"Structure {i+1}: SMILES={struct['smiles']}"
                if 'formula' in struct:
                    text += f", Formula={struct['formula']}"
                if 'molecular_weight' in struct:
                    text += f", MW={struct['molecular_weight']:.2f}"
                text += "\n"
        
        return text
    
    def similarity_search(self, query, k=5, filter_dict=None):
        """Enhanced similarity search"""
        # Check if query contains SMILES
        smiles_pattern = r'[C,c][0-9A-Za-z@+\-\[\]\(\)\\=#$]+'
        
        if re.search(smiles_pattern, query):
            # Query contains SMILES, enhance search
            return self._structure_aware_search(query, k)
        
        # Regular text search
        if filter_dict:
            return self.vectorstore.similarity_search(
                query, k=k, filter=filter_dict
            )
        return self.vectorstore.similarity_search(query, k=k)
    
    def _structure_aware_search(self, query, k=5):
        """Search considering chemical structures"""
        # Extract potential SMILES from query
        smiles_pattern = r'[C,c][0-9A-Za-z@+\-\[\]\(\)\\=#$]+'
        potential_smiles = re.findall(smiles_pattern, query)
        
        # Search for chunks containing similar structures
        results = []
        
        if potential_smiles:
            # Find chunks with these SMILES
            # Since structure_smiles is now a comma-separated string, we need to search differently
            for smiles in potential_smiles:
                filter_results = self.vectorstore.similarity_search(
                    query, 
                    k=k,
                    filter=lambda metadata: smiles in metadata.get("structure_smiles", "")
                )
                results.extend(filter_results)
        
        # If not enough results, do regular search
        if len(results) < k:
            additional = self.vectorstore.similarity_search(
                query, k=k-len(results)
            )
            results.extend(additional)
        
        return results