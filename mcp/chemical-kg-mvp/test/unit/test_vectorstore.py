#!/usr/bin/env python3
"""
Test vector store functionality with extracted structures
"""

import sys
import os

sys.path.insert(0, 'src')

def test_vectorstore_basic():
    """Test basic vector store functionality"""
    print("Testing Vector Store")
    print("=" * 40)
    
    try:
        from vectorstore import ChemicalVectorStore
        print("✓ ChemicalVectorStore imported successfully")
    except Exception as e:
        print(f"✗ Failed to import ChemicalVectorStore: {e}")
        return False
    
    # Test initialization
    try:
        # Try HuggingFace embeddings first (no external dependencies)
        vector_store = ChemicalVectorStore(use_ollama_embeddings=False)
        print("✓ Vector store initialized with HuggingFace embeddings")
    except Exception as e:
        print(f"✗ Failed to initialize vector store: {e}")
        return False
    
    return True

def test_with_sample_data():
    """Test with sample chemical data"""
    print("\nTesting with Sample Data")
    print("=" * 40)
    
    try:
        from chunker import ChemicalAwareChunker
        from vectorstore import ChemicalVectorStore
        
        # Create sample data
        sample_text = """
        Aspirin (acetylsalicylic acid) is a medication used to reduce pain, fever, 
        or inflammation. It is one of the most widely used medications globally.
        
        Ibuprofen is a nonsteroidal anti-inflammatory drug (NSAID) that is used 
        for treating pain, fever, and inflammation.
        """
        
        sample_structures = [
            {
                "smiles": "CC(=O)Oc1ccccc1C(=O)O",
                "formula": "C9H8O4",
                "molecular_weight": 180.16,
                "context": "Aspirin chemical structure"
            },
            {
                "smiles": "CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O",
                "formula": "C13H18O2",
                "molecular_weight": 206.28,
                "context": "Ibuprofen chemical structure"
            }
        ]
        
        # Create chunks
        chunker = ChemicalAwareChunker()
        chunks = chunker.chunk_with_structures(sample_text, sample_structures)
        print(f"Created {len(chunks)} chunks")
        
        # Create vector store
        vector_store = ChemicalVectorStore(use_ollama_embeddings=False)
        vector_store.create_vectorstore(chunks)
        print("✓ Vector store created successfully")
        
        # Test search
        results = vector_store.similarity_search("What is aspirin?", k=2)
        print(f"✓ Search returned {len(results)} results")
        
        return True
        
    except Exception as e:
        print(f"✗ Error in sample data test: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_with_extracted_structures():
    """Test with structures extracted from PDF"""
    print("\nTesting with Extracted Structures")
    print("=" * 40)
    
    # Check if we have extracted structures
    structure_files = []
    for f in os.listdir('.'):
        if f.startswith('structure_') and f.endswith('.png'):
            structure_files.append(f)
    
    if not structure_files:
        print("No extracted structure files found")
        return False
    
    print(f"Found {len(structure_files)} structure files")
    
    try:
        from decimer_client import DECIMERClient
        from chemical_handler import ChemicalHandler
        from chunker import ChemicalAwareChunker
        from vectorstore import ChemicalVectorStore
        
        # Initialize components
        decimer_client = DECIMERClient()
        chem_handler = ChemicalHandler(decimer_client)
        
        # Process structures
        structures = []
        for i, struct_file in enumerate(structure_files[:3]):  # Test first 3
            print(f"Processing {struct_file}...")
            try:
                structure = chem_handler.process_image(struct_file, f"Structure {i+1}")
                if structure:
                    structures.append(structure)
                    print(f"  ✓ SMILES: {structure.get('smiles', 'N/A')}")
            except Exception as e:
                print(f"  ✗ Failed: {e}")
        
        if not structures:
            print("No structures processed successfully")
            return False
        
        # Create vector store
        print(f"\nCreating vector store with {len(structures)} structures...")
        
        sample_text = "Chemical structures extracted from scientific literature."
        chunker = ChemicalAwareChunker()
        chunks = chunker.chunk_with_structures(sample_text, structures)
        
        vector_store = ChemicalVectorStore(use_ollama_embeddings=False)
        vector_store.create_vectorstore(chunks)
        
        print("✓ Vector store created successfully with real structures")
        
        return True
        
    except Exception as e:
        print(f"✗ Error: {e}")
        import traceback
        traceback.print_exc()
        return False

def diagnose_type_error():
    """Diagnose the type subscription error"""
    print("\nDiagnosing Type Error")
    print("=" * 40)
    
    # Check Python version
    import sys
    print(f"Python version: {sys.version}")
    
    # Check if it's a type hint issue
    try:
        from typing import List, Dict, Any
        test_type = List[str]  # This would fail in Python < 3.9 without __future__
        print("✓ Type hints work correctly")
    except Exception as e:
        print(f"✗ Type hint error: {e}")
        print("  This might be the cause of the 'type' object is not subscriptable error")
    
    # Test specific imports that might cause issues
    try:
        import langchain
        print(f"✓ LangChain version: {langchain.__version__}")
    except:
        print("✗ LangChain not found")
    
    try:
        import chromadb
        print(f"✓ ChromaDB version: {chromadb.__version__}")
    except:
        print("✗ ChromaDB not found")

if __name__ == "__main__":
    print("Vector Store Testing Suite")
    print("=" * 50)
    
    # Run diagnostics first
    diagnose_type_error()
    
    # Test basic functionality
    if test_vectorstore_basic():
        # Test with sample data
        test_with_sample_data()
        
        # Test with real extracted structures
        test_with_extracted_structures()
    
    print("\nTest completed!")