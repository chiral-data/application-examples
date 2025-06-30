# chunker.py
from langchain.text_splitter import RecursiveCharacterTextSplitter
import re

class ChemicalAwareChunker:
    def __init__(self, chunk_size=1000, overlap=200):
        self.splitter = RecursiveCharacterTextSplitter(
            chunk_size=chunk_size,
            chunk_overlap=overlap,
            separators=["\n\n", "\n", ".", " ", ""]
        )
    
    def chunk_with_structures(self, text, structures):
        """Create chunks that maintain structure-text relationships"""
        chunks = []
        
        # Basic chunking
        text_chunks = self.splitter.split_text(text)
        
        for i, chunk in enumerate(text_chunks):
            # Check if chunk mentions any structure
            related_structures = []
            
            # Simple keyword matching (improve this)
            for struct in structures:
                if any(keyword in chunk.lower() for keyword in 
                       ['compound', 'structure', 'figure', 'scheme']):
                    related_structures.append(struct)
            
            chunks.append({
                'text': chunk,
                'chunk_id': i,
                'structures': related_structures,
                'metadata': {
                    'has_structures': len(related_structures) > 0
                }
            })
        
        return chunks