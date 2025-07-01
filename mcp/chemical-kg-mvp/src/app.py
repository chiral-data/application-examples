# app.py
import streamlit as st
import os
from pdf_processor import PDFProcessor
from decimer_client import DECIMERClient
from chemical_handler import ChemicalHandler
from chunker import ChemicalAwareChunker
from vectorstore import ChemicalVectorStore
from rag_chain import ChemicalRAG

st.set_page_config(page_title="Chemical Paper Analyzer MVP", layout="wide")

# Initialize session state
if 'vectorstore' not in st.session_state:
    st.session_state.vectorstore = None
if 'structures' not in st.session_state:
    st.session_state.structures = []

st.title("Chemical Paper Analyzer - MVP")
st.markdown("Upload a chemistry paper and ask questions about structures and content")

# Sidebar for file upload
with st.sidebar:
    st.header("Upload Paper")
    uploaded_file = st.file_uploader("Choose a PDF file", type="pdf")
    
    if uploaded_file is not None:
        # Save uploaded file
        with open("temp_paper.pdf", "wb") as f:
            f.write(uploaded_file.getvalue())
        
        if st.button("Process Paper"):
            with st.spinner("Processing paper..."):
                # Initialize components
                pdf_processor = PDFProcessor("temp_paper.pdf")
                decimer_client = DECIMERClient()
                chem_handler = ChemicalHandler(decimer_client)
                chunker = ChemicalAwareChunker()
                
                # Check if user wants to use DECIMER segmentation
                use_decimer_segmentation = st.sidebar.checkbox("Use DECIMER Segmentation", value=True)
                
                # Extract content
                st.info("Extracting text and images...")
                pages_data = pdf_processor.extract_text_and_images()
                
                # Extract full text
                full_text = ""
                for page_data in pages_data:
                    full_text += page_data['text'] + "\n"
                
                # Process chemical structures
                st.info("Identifying chemical structures...")
                all_structures = []
                
                if use_decimer_segmentation:
                    # Use DECIMER to segment and process structures
                    chemical_images = pdf_processor.extract_chemical_structures_with_decimer()
                    
                    for img_info in chemical_images:
                        structure = chem_handler.process_image(
                            img_info['path'],
                            f"Chemical structure {img_info['index'] + 1}"
                        )
                        
                        if structure:
                            all_structures.append(structure)
                else:
                    # Use standard image extraction
                    for page_data in pages_data:
                        for img_info in page_data['images']:
                            context = chem_handler.get_structure_context(
                                page_data['text'], 
                                img_info['bbox']
                            )
                            
                            structure = chem_handler.process_image(
                                img_info['path'],
                                context
                            )
                            
                            if structure:
                                all_structures.append(structure)
                
                st.session_state.structures = all_structures
                
                # Create chunks
                st.info("Creating document chunks...")
                chunks = chunker.chunk_with_structures(full_text, all_structures)
                
                # Create vector store
                st.info("Building vector database...")
                vector_store = ChemicalVectorStore()
                vector_store.create_vectorstore(chunks)
                st.session_state.vectorstore = vector_store
                
                st.success("Paper processed successfully!")
                st.metric("Total Structures Found", len(all_structures))

# Main interface
col1, col2 = st.columns([2, 1])

with col1:
    st.header("Ask Questions")
    
    if st.session_state.vectorstore is not None:
        # Initialize RAG with configured model
        model_name = os.getenv("OLLAMA_MODEL", "llama3.2:latest")
        rag = ChemicalRAG(st.session_state.vectorstore, model_name=model_name)
        
        # Query input
        user_query = st.text_input("Enter your question about the paper:")
        
        if st.button("Ask") and user_query:
            with st.spinner("Searching and generating answer..."):
                answer = rag.query(user_query)
                
                st.markdown("### Answer:")
                st.write(answer)
                
                # Show relevant structures
                st.markdown("### Related Chemical Structures:")
                # Simple visualization (in production, use RDKit)
                for struct in st.session_state.structures[:3]:
                    st.code(f"SMILES: {struct['smiles']}")
                    st.text(f"Formula: {struct['formula']}")
                    st.text(f"MW: {struct['molecular_weight']:.2f}")
    else:
        st.warning("Please upload and process a paper first.")

with col2:
    st.header("Extracted Structures")
    
    if st.session_state.structures:
        for i, struct in enumerate(st.session_state.structures):
            with st.expander(f"Structure {i+1}"):
                st.code(struct['smiles'])
                st.text(f"Formula: {struct['formula']}")
                st.text(f"MW: {struct['molecular_weight']:.2f}")
                if os.path.exists(struct['image_path']):
                    st.image(struct['image_path'])

# Example queries
st.markdown("### Example Questions:")
example_queries = [
    "What are the main chemical compounds discussed in this paper?",
    "What synthesis methods are described?",
    "What are the key findings about the chemical structures?",
    "Compare the molecular weights of the compounds found.",
    "What reactions are mentioned in the paper?"
]

for query in example_queries:
    st.code(query)