# chunker.py
from langchain.text_splitter import RecursiveCharacterTextSplitter
import re
import os
from typing import List, Dict, Any, Optional
from collections import defaultdict

class ChemicalAwareChunker:
    def __init__(self, chunk_size=1000, overlap=200):
        # Get configuration from environment
        self.chunk_size = int(os.getenv('CHUNK_SIZE', chunk_size))
        self.overlap = int(os.getenv('CHUNK_OVERLAP', overlap))
        self.context_window = int(os.getenv('CHEMICAL_CONTEXT_WINDOW', 5))
        
        # Initialize text splitter with chemistry-aware separators
        self.splitter = RecursiveCharacterTextSplitter(
            chunk_size=self.chunk_size,
            chunk_overlap=self.overlap,
            separators=[
                "\n\n\n",  # Multiple newlines (section breaks)
                "\n\n",    # Paragraph breaks
                "\n",      # Line breaks
                ". ",      # Sentence endings
                "? ",      # Question endings
                "! ",      # Exclamation endings
                "; ",      # Semicolon separators
                ", ",      # Comma separators
                " ",       # Word boundaries
                ""
            ]
        )
        
        # Chemical keywords for better association
        self.chemical_keywords = {
            'structure': ['structure', 'compound', 'molecule', 'chemical', 'formula', 'smiles'],
            'synthesis': ['synthesis', 'preparation', 'reaction', 'yield', 'procedure'],
            'analysis': ['analysis', 'characterization', 'nmr', 'ir', 'ms', 'spectr'],
            'properties': ['property', 'properties', 'melting', 'boiling', 'solubility', 'stability'],
            'reference': ['figure', 'fig', 'scheme', 'table', 'chart', 'diagram']
        }
    
    def chunk_with_structures(self, text: str, structures: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Create chunks that maintain structure-text relationships with improved association"""
        # First, create structure-aware chunks
        enhanced_chunks = self._create_structure_aware_chunks(text, structures)
        
        # Then apply intelligent chunking
        final_chunks = self._apply_intelligent_chunking(enhanced_chunks, structures)
        
        return final_chunks
    
    def _create_structure_aware_chunks(self, text: str, structures: List[Dict[str, Any]]) -> List[str]:
        """Create initial chunks while preserving chemical context"""
        # Split text into paragraphs first
        paragraphs = re.split(r'\n\s*\n', text)
        
        chunks = []
        current_chunk = ""
        
        for paragraph in paragraphs:
            # Check if paragraph contains chemical references
            contains_chem_ref = self._contains_chemical_reference(paragraph)
            
            # If adding this paragraph would exceed chunk size
            if len(current_chunk) + len(paragraph) > self.chunk_size:
                if current_chunk:
                    chunks.append(current_chunk.strip())
                    current_chunk = ""
                
                # If paragraph itself is too long, split it
                if len(paragraph) > self.chunk_size:
                    sub_chunks = self.splitter.split_text(paragraph)
                    chunks.extend(sub_chunks[:-1])
                    current_chunk = sub_chunks[-1] if sub_chunks else ""
                else:
                    current_chunk = paragraph
            else:
                current_chunk += "\n\n" + paragraph if current_chunk else paragraph
        
        if current_chunk:
            chunks.append(current_chunk.strip())
        
        return chunks
    
    def _apply_intelligent_chunking(self, text_chunks: List[str], structures: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Apply intelligent chunking with structure association"""
        chunks = []
        
        for i, chunk_text in enumerate(text_chunks):
            # Find related structures for this chunk
            related_structures = self._find_related_structures(chunk_text, structures, i)
            
            # Calculate chunk importance and context
            chunk_metadata = self._calculate_chunk_metadata(chunk_text, related_structures)
            
            # Create enhanced chunk
            chunk = {
                'text': chunk_text,
                'chunk_id': i,
                'structures': related_structures,
                'metadata': chunk_metadata
            }
            
            chunks.append(chunk)
        
        # Post-process chunks for better coherence
        return self._post_process_chunks(chunks)
    
    def _find_related_structures(self, chunk_text: str, structures: List[Dict[str, Any]], chunk_index: int) -> List[Dict[str, Any]]:
        """Find structures related to a text chunk using multiple methods"""
        related_structures = []
        chunk_lower = chunk_text.lower()
        
        for struct in structures:
            relevance_score = 0
            
            # Method 1: Direct SMILES mention
            if struct.get('smiles', '').lower() in chunk_lower:
                relevance_score += 10
            
            # Method 2: Formula mention
            if struct.get('formula', '').lower() in chunk_lower:
                relevance_score += 8
            
            # Method 3: Context from image location
            if 'context' in struct and struct['context']:
                context_words = struct['context'].lower().split()
                chunk_words = chunk_lower.split()
                common_words = set(context_words) & set(chunk_words)
                relevance_score += len(common_words) * 2
            
            # Method 4: Chemical keyword proximity
            for category, keywords in self.chemical_keywords.items():
                for keyword in keywords:
                    if keyword in chunk_lower:
                        relevance_score += 1
            
            # Method 5: Figure/scheme references
            figure_refs = re.findall(r'\b(?:fig|figure|scheme|chart)\s*\d+', chunk_lower)
            if figure_refs:
                relevance_score += 3
            
            # Method 6: Page proximity (if available)
            if 'page' in struct:
                # Assume chunks are roughly in page order
                estimated_page = chunk_index // 3  # Rough estimate
                page_distance = abs(struct['page'] - estimated_page)
                if page_distance <= 1:
                    relevance_score += 5 - page_distance
            
            # Add structure if relevance score is high enough
            if relevance_score >= 3:
                struct_copy = struct.copy()
                struct_copy['relevance_score'] = relevance_score
                related_structures.append(struct_copy)
        
        # Sort by relevance score (highest first)
        related_structures.sort(key=lambda x: x.get('relevance_score', 0), reverse=True)
        
        # Limit to top 5 most relevant structures
        return related_structures[:5]
    
    def _calculate_chunk_metadata(self, chunk_text: str, structures: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Calculate metadata for a chunk"""
        metadata = {
            'has_structures': len(structures) > 0,
            'structure_count': len(structures),
            'word_count': len(chunk_text.split()),
            'char_count': len(chunk_text),
            'chemical_keywords': [],
            'content_type': 'general'
        }
        
        # Identify chemical keywords
        chunk_lower = chunk_text.lower()
        for category, keywords in self.chemical_keywords.items():
            found_keywords = [kw for kw in keywords if kw in chunk_lower]
            if found_keywords:
                metadata['chemical_keywords'].extend(found_keywords)
        
        # Determine content type
        if any(kw in chunk_lower for kw in self.chemical_keywords['synthesis']):
            metadata['content_type'] = 'synthesis'
        elif any(kw in chunk_lower for kw in self.chemical_keywords['analysis']):
            metadata['content_type'] = 'analysis'
        elif any(kw in chunk_lower for kw in self.chemical_keywords['properties']):
            metadata['content_type'] = 'properties'
        elif structures:
            metadata['content_type'] = 'structure_related'
        
        # Calculate chemical density (chemical terms per 100 words)
        if metadata['word_count'] > 0:
            chemical_density = (len(metadata['chemical_keywords']) / metadata['word_count']) * 100
            metadata['chemical_density'] = round(chemical_density, 2)
        
        # Convert list to string for ChromaDB compatibility
        metadata['chemical_keywords'] = ','.join(metadata['chemical_keywords'])
        
        return metadata
    
    def _contains_chemical_reference(self, text: str) -> bool:
        """Check if text contains chemical references"""
        text_lower = text.lower()
        
        # Check for chemical keywords
        for keywords in self.chemical_keywords.values():
            if any(keyword in text_lower for keyword in keywords):
                return True
        
        # Check for chemical formulas (basic pattern)
        if re.search(r'\b[A-Z][a-z]?\d*(?:[A-Z][a-z]?\d*)*\b', text):
            return True
        
        # Check for SMILES-like patterns
        if re.search(r'\b[CNOPSFClBrI]\w*[()\[\]=#@+-]\w*', text):
            return True
        
        return False
    
    def _post_process_chunks(self, chunks: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Post-process chunks for better coherence"""
        # Merge very small chunks with high chemical content
        processed_chunks = []
        i = 0
        
        while i < len(chunks):
            current_chunk = chunks[i]
            
            # Check if current chunk is too small and has chemical content
            if (current_chunk['metadata']['char_count'] < self.chunk_size // 3 and 
                current_chunk['metadata']['has_structures'] and 
                i < len(chunks) - 1):
                
                next_chunk = chunks[i + 1]
                
                # Merge if next chunk is also small or related
                if (next_chunk['metadata']['char_count'] < self.chunk_size // 2 or 
                    self._chunks_are_related(current_chunk, next_chunk)):
                    
                    merged_chunk = self._merge_chunks(current_chunk, next_chunk)
                    processed_chunks.append(merged_chunk)
                    i += 2  # Skip next chunk as it's been merged
                    continue
            
            processed_chunks.append(current_chunk)
            i += 1
        
        # Renumber chunks
        for i, chunk in enumerate(processed_chunks):
            chunk['chunk_id'] = i
        
        return processed_chunks
    
    def _chunks_are_related(self, chunk1: Dict[str, Any], chunk2: Dict[str, Any]) -> bool:
        """Check if two chunks are thematically related"""
        # Check if they share structures
        struct1_smiles = {s.get('smiles', '') for s in chunk1['structures']}
        struct2_smiles = {s.get('smiles', '') for s in chunk2['structures']}
        
        if struct1_smiles & struct2_smiles:  # Common structures
            return True
        
        # Check if they have similar content types
        if (chunk1['metadata']['content_type'] == chunk2['metadata']['content_type'] and 
            chunk1['metadata']['content_type'] != 'general'):
            return True
        
        # Check keyword overlap
        keywords1_str = chunk1['metadata'].get('chemical_keywords', '')
        keywords2_str = chunk2['metadata'].get('chemical_keywords', '')
        
        # Convert back to sets for comparison
        keywords1 = set(keywords1_str.split(',')) if keywords1_str else set()
        keywords2 = set(keywords2_str.split(',')) if keywords2_str else set()
        
        if len(keywords1 & keywords2) >= 2:  # At least 2 common keywords
            return True
        
        return False
    
    def _merge_chunks(self, chunk1: Dict[str, Any], chunk2: Dict[str, Any]) -> Dict[str, Any]:
        """Merge two chunks intelligently"""
        merged_text = chunk1['text'] + "\n\n" + chunk2['text']
        
        # Combine structures, removing duplicates
        all_structures = chunk1['structures'] + chunk2['structures']
        unique_structures = []
        seen_smiles = set()
        
        for struct in all_structures:
            smiles = struct.get('smiles', '')
            if smiles not in seen_smiles:
                unique_structures.append(struct)
                seen_smiles.add(smiles)
        
        # Recalculate metadata
        merged_metadata = self._calculate_chunk_metadata(merged_text, unique_structures)
        
        return {
            'text': merged_text,
            'chunk_id': chunk1['chunk_id'],  # Keep first chunk's ID
            'structures': unique_structures,
            'metadata': merged_metadata
        }
    
    def get_chunking_stats(self, chunks: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Get statistics about the chunking process"""
        total_chunks = len(chunks)
        chunks_with_structures = sum(1 for c in chunks if c['metadata']['has_structures'])
        total_structures = sum(len(c['structures']) for c in chunks)
        
        content_types = defaultdict(int)
        for chunk in chunks:
            content_types[chunk['metadata']['content_type']] += 1
        
        avg_chunk_size = sum(c['metadata']['char_count'] for c in chunks) / max(total_chunks, 1)
        avg_chemical_density = sum(c['metadata'].get('chemical_density', 0) for c in chunks) / max(total_chunks, 1)
        
        return {
            'total_chunks': total_chunks,
            'chunks_with_structures': chunks_with_structures,
            'structure_coverage': round(chunks_with_structures / max(total_chunks, 1) * 100, 2),
            'total_structures': total_structures,
            'avg_structures_per_chunk': round(total_structures / max(total_chunks, 1), 2),
            'content_type_distribution': dict(content_types),
            'avg_chunk_size': round(avg_chunk_size),
            'avg_chemical_density': round(avg_chemical_density, 2)
        }