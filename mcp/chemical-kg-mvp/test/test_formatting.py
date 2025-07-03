#!/usr/bin/env python3
"""
Test script for the RAG response formatting system
"""
import sys
import os

# Add src directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from llm_call import format_rag_response
import json

def test_example_response():
    """Test with the example response from the user"""
    
    # Sample RAG response that matches the user's example
    sample_response = {
        'answer': 'The main chemical compounds discussed in the provided context are:\n\n1. 1,3,4-oxadiazolyl compounds:\n - Formula: C17H16N2O2 (MW: 280.33)\n - SMILES: CC1=NC2=C(C=C1C(=O)OC)C3=CC=CC(=C3N2)C4CC4\n These compounds are substituted with various R groups and exhibit inhibitory activities against PDHK2, PDHK4, and BCKDK. The exact structure of the specific compound 18 is not explicitly provided in the context.\n\n2. Pyrrolidin-1-yl-containing compound:\n - Formula: C9H14N4OS (MW: 226.30)\n - SMILES: CCC1=NN=C(C)S1.CC1=NN=C(C)O1\n This compound is part of the synthesis pathway and is not explicitly mentioned in terms of its biological activity.\n\n3. Pyrimidine-based compounds:\n - Formula: C5H10O2 (MW: 102.13)\n - SMILES: CC(C)CC(=O)O\n - Formula: C4H8NO2P (MW: 133.09)\n - SMILES: C([C@H]1COC(=O)N1)P\n These compounds appear to be intermediates in the synthesis pathway and are not explicitly mentioned in terms of their biological activity.\n\nIt is unclear what specific role these compounds play in the context of inhibiting PDHK2, PDHK4, and BCKDK.',
        'sources': [
            {
                'page_content': 'intramolecular pi-pi interaction. \nTable 3 \nInhibitory activities of substituted 1,3,4-oxadiazolyl compounds against PDHK2, \nPDHK4 and BCKDK.\ncompd \nR1 \nIC50 (μM) \nPDHK2 \nPDHK4 \nBCKDK a \n15 \nH \n0.0558 \n0.4624 \n18 \n0.0042 \n0.0241 \n>100 (3.0%) b \n19 \n0.0060 \n0.0329 \n>100 (18%) b \n20 \n0.0171 \n0.0615',
                'metadata': {'chunk_id': 25, 'page': 3, 'has_structures': True, 'structure_count': 3}
            }
        ]
    }
    
    print("Testing RAG Response Formatting")
    print("=" * 50)
    
    print("\nORIGINAL RESPONSE:")
    print("-" * 30)
    print(sample_response['answer'][:500] + "..." if len(sample_response['answer']) > 500 else sample_response['answer'])
    
    print("\nFORMATTING WITH LLAMA MODEL...")
    print("-" * 30)
    
    try:
        formatted_response = format_rag_response(sample_response)
        
        print("\nFORMATTED RESPONSE:")
        print("-" * 30)
        print(formatted_response)
        
        # Check if formatting was successful
        checks = {
            "SMILES removed": "CC1=NC2=" not in formatted_response,
            "MW removed": "MW: 280.33" not in formatted_response,
            "Formula removed": "Formula: C17H16N2O2" not in formatted_response,
            "Structure references": any(word in formatted_response.lower() for word in ["structure", "compound"]),
            "Readable format": len(formatted_response) > 50,
            "IC50 preserved": "ic50" in formatted_response.lower() or "inhibitory" in formatted_response.lower()
        }
        
        print("\nFORMATTING CHECKS:")
        print("-" * 30)
        for check, passed in checks.items():
            status = "PASS" if passed else "FAIL"
            print(f"{check}: {status}")
        
        all_passed = all(checks.values())
        print(f"\nOVERALL: {'SUCCESS' if all_passed else 'NEEDS IMPROVEMENT'}")
        
        return all_passed
        
    except Exception as e:
        print(f"\nERROR: {e}")
        print("\nTesting fallback formatting...")
        
        # Test fallback
        from llm_call import ResponseFormatter
        formatter = ResponseFormatter()
        fallback_result = formatter._basic_format(sample_response)
        
        print("FALLBACK RESULT:")
        print(fallback_result)
        
        return False

def test_simple_response():
    """Test with a simple string response"""
    
    print("\nTesting Simple String Response")
    print("=" * 50)
    
    simple_response = "This is a simple test response with SMILES: CC1=NC2=C(C=C1C(=O)OC)C3=CC=CC and Formula: C17H16N2O2 MW: 280.33"
    
    print(f"ORIGINAL: {simple_response}")
    
    try:
        formatted = format_rag_response(simple_response)
        print(f"FORMATTED: {formatted}")
        
        # Simple check
        smiles_removed = "CC1=NC2=" not in formatted
        print(f"SMILES removed: {'YES' if smiles_removed else 'NO'}")
        
        return smiles_removed
        
    except Exception as e:
        print(f"ERROR: {e}")
        return False

def test_ollama_connection():
    """Test if Ollama is accessible"""
    
    print("\nTesting Ollama Connection")
    print("=" * 50)
    
    try:
        import requests
        import os
        
        ollama_host = os.getenv("OLLAMA_HOST", "localhost")
        ollama_port = os.getenv("OLLAMA_PORT", "11434")
        base_url = f"http://{ollama_host}:{ollama_port}"
        
        print(f"Connecting to: {base_url}")
        
        # Test connection
        response = requests.get(f"{base_url}/api/tags", timeout=5)
        
        if response.status_code == 200:
            models = response.json().get('models', [])
            print("Ollama connection successful!")
            print(f"Available models: {[m.get('name', 'unknown') for m in models]}")
            
            # Check if our model is available
            model_name = os.getenv("OLLAMA_MODEL", "llama3.2:latest")
            model_available = any(model_name in m.get('name', '') for m in models)
            print(f"Model '{model_name}' available: {'YES' if model_available else 'NO'}")
            
            return model_available
        else:
            print(f"Ollama responded with status: {response.status_code}")
            return False
            
    except Exception as e:
        print(f"Connection failed: {e}")
        print("Make sure Ollama is running and accessible")
        return False

def main():
    """Run all tests"""
    
    print("RAG Response Formatting Test Suite")
    print("=" * 60)
    
    # Test Ollama connection first
    ollama_ok = test_ollama_connection()
    
    # Test complex response
    complex_ok = test_example_response()
    
    # Test simple response
    simple_ok = test_simple_response()
    
    print("\n" + "=" * 60)
    print("TEST SUMMARY")
    print("=" * 60)
    print(f"Ollama Connection: {'OK' if ollama_ok else 'FAILED'}")
    print(f"Complex Response: {'OK' if complex_ok else 'FAILED'}")
    print(f"Simple Response: {'OK' if simple_ok else 'FAILED'}")
    
    if all([ollama_ok, complex_ok, simple_ok]):
        print("\nALL TESTS PASSED! Ready to use in webapp.")
    elif complex_ok or simple_ok:
        print("\nPARTIAL SUCCESS. Formatting works but check Ollama connection.")
    else:
        print("\nTESTS FAILED. Check setup and try again.")
    
    print("\nNext steps:")
    print("   1. If tests passed, the webapp should work correctly")
    print("   2. If Ollama failed, check docker-compose services")
    print("   3. If formatting failed, check the prompt template")

if __name__ == "__main__":
    main()