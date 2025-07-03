# app.py
import streamlit as st
import os
import base64
import io
import pandas as pd
from PIL import Image
from pdf_processor import PDFProcessor
from decimer_client import DECIMERClient
from chemical_handler import ChemicalHandler
from chunker import ChemicalAwareChunker
from vectorstore import ChemicalVectorStore
from rag_chain import ChemicalRAG

# Configure page
st.set_page_config(
    page_title="Chemical Knowledge Graph - MVP", 
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Pharma color scheme CSS
st.markdown("""
<style>
    /* Main theme colors - Professional pharma palette */
    .main {
        background-color: #f8fafc;
    }
    
    /* Sidebar styling */
    .css-1d391kg {
        background-color: #1e3a5f;
    }
    
    /* Headers */
    h1 {
        color: #1e3a5f;
        font-family: 'Arial', sans-serif;
        font-weight: 600;
        border-bottom: 3px solid #0066cc;
        padding-bottom: 10px;
    }
    
    h2, h3 {
        color: #2c5282;
        font-family: 'Arial', sans-serif;
        font-weight: 500;
    }
    
    /* Buttons */
    .stButton > button {
        background-color: #0066cc;
        color: white;
        border-radius: 6px;
        border: none;
        padding: 0.5rem 1rem;
        font-weight: 500;
        transition: all 0.3s ease;
    }
    
    .stButton > button:hover {
        background-color: #0052a3;
        box-shadow: 0 4px 8px rgba(0, 102, 204, 0.2);
    }
    
    /* Download buttons */
    .download-button {
        background-color: #16a085;
        color: white;
        border-radius: 6px;
        padding: 0.5rem 1rem;
        text-decoration: none;
        font-weight: 500;
        margin: 0.25rem;
        display: inline-block;
        transition: all 0.3s ease;
    }
    
    .download-button:hover {
        background-color: #138d75;
        box-shadow: 0 4px 8px rgba(22, 160, 133, 0.2);
    }
    
    /* Cards and containers */
    .stExpander {
        border: 1px solid #e2e8f0;
        border-radius: 8px;
        background-color: white;
        box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
    }
    
    /* Metrics */
    [data-testid="metric-container"] {
        background-color: white;
        border: 1px solid #e2e8f0;
        border-radius: 8px;
        padding: 1rem;
        box-shadow: 0 2px 4px rgba(0, 0, 0, 0.05);
    }
    
    /* Success/Info messages */
    .stSuccess {
        background-color: #d4edda;
        border-color: #c3e6cb;
        color: #155724;
    }
    
    .stInfo {
        background-color: #cce7ff;
        border-color: #99d6ff;
        color: #0066cc;
    }
    
    /* Code blocks */
    .stCode {
        background-color: #f1f5f9;
        border: 1px solid #e2e8f0;
        border-radius: 6px;
    }
    
    /* File uploader */
    .stFileUploader {
        background-color: white;
        border: 2px dashed #cbd5e0;
        border-radius: 8px;
        padding: 2rem;
    }
    
    /* Sidebar content */
    .css-1d391kg .stMarkdown {
        color: white;
    }
    
    .css-1d391kg .stSelectbox label {
        color: white;
    }
</style>
""", unsafe_allow_html=True)

# Initialize session state
if 'vectorstore' not in st.session_state:
    st.session_state.vectorstore = None
if 'structures' not in st.session_state:
    st.session_state.structures = []

# 🧪 Chemical Knowledge Graph - MVP

st.markdown("""
<div style="background: linear-gradient(90deg, #0066cc 0%, #1e3a5f 100%); padding: 1.5rem; border-radius: 10px; margin-bottom: 2rem;">
    <h3 style="color: white; margin: 0; text-align: center;">
        📄 Upload chemistry papers • 🔍 Extract structures • 💬 Ask intelligent questions
    </h3>
</div>
""", unsafe_allow_html=True)

st.markdown("""
**Transform your chemical literature into an intelligent knowledge base.** 
This tool uses AI to extract chemical structures, analyze content, and answer questions about your papers.
""")

# Helper functions for downloads
def create_download_link(val, filename, link_text):
    \"\"\"Create a download link for data\"\"\"
    b64 = base64.b64encode(val).decode()
    return f'<a href="data:application/octet-stream;base64,{b64}" download="{filename}" class="download-button">{link_text}</a>'

def create_csv_download(structures):
    \"\"\"Create CSV data from structures\"\"\"
    if not structures:
        return None
    
    data = []
    for i, struct in enumerate(structures):
        data.append({
            'Structure_ID': f'Structure_{i+1}',
            'SMILES': struct.get('smiles', ''),
            'Molecular_Formula': struct.get('formula', ''),
            'Molecular_Weight': struct.get('molecular_weight', ''),
            'Context': struct.get('context', '')[:100] + '...' if len(struct.get('context', '')) > 100 else struct.get('context', '')
        })
    
    df = pd.DataFrame(data)
    return df.to_csv(index=False).encode('utf-8')

def get_image_download_link(image_path, filename):
    \"\"\"Create download link for images\"\"\"
    if os.path.exists(image_path):
        with open(image_path, \"rb\") as file:
            return create_download_link(file.read(), filename, f\"📥 Download {filename}\")
    return None

# Sidebar for file upload
with st.sidebar:
    st.markdown(\"### 📁 Upload Document\")\n    uploaded_file = st.file_uploader(\"Choose a PDF file\", type=\"pdf\", help=\"Upload a chemistry paper to analyze\")"}
    
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