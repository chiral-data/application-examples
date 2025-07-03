# app.py
# Fix for SQLite version compatibility with ChromaDB
try:
    __import__('pysqlite3')
    import sys
    sys.modules['sqlite3'] = sys.modules.pop('pysqlite3')
except ImportError:
    pass

import streamlit as st
import os
import base64
import io
import json
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
if 'demo_mode' not in st.session_state:
    st.session_state.demo_mode = False

# Demo data functions
def load_demo_data():
    """Load pre-generated demo data"""
    try:
        demo_data_path = os.path.join(os.path.dirname(__file__), '..', 'demo_data.json')
        if os.path.exists(demo_data_path):
            with open(demo_data_path, 'r') as f:
                demo_data = json.load(f)
            return demo_data
    except Exception as e:
        st.error(f"Error loading demo data: {e}")
    
    # Fallback demo data if file doesn't exist
    return {
        "structures": [
            {
                "smiles": "CC(=O)Oc1ccccc1C(=O)O",
                "formula": "C9H8O4", 
                "molecular_weight": 180.16,
                "context": "Aspirin (acetylsalicylic acid) is a medication used to reduce pain, fever, or inflammation.",
                "name": "Aspirin",
                "image_path": None
            },
            {
                "smiles": "CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O",
                "formula": "C13H18O2",
                "molecular_weight": 206.28, 
                "context": "Ibuprofen is a nonsteroidal anti-inflammatory drug (NSAID).",
                "name": "Ibuprofen",
                "image_path": None
            }
        ],
        "text": "This is demo chemical literature containing pharmaceutical compounds.",
        "source": "demo_fallback"
    }

def activate_demo_mode():
    """Activate demo mode with pre-loaded data"""
    st.session_state.demo_mode = True
    demo_data = load_demo_data()
    
    # Load demo structures
    st.session_state.structures = demo_data.get('structures', [])
    
    # Create mock vectorstore with demo text
    from chunker import ChemicalAwareChunker
    from vectorstore import ChemicalVectorStore
    
    try:
        chunker = ChemicalAwareChunker()
        chunks = chunker.chunk_with_structures(demo_data.get('text', ''), st.session_state.structures)
        
        vector_store = ChemicalVectorStore()
        vector_store.create_vectorstore(chunks)
        st.session_state.vectorstore = vector_store
        
        st.success("Demo mode activated! Sample chemical structures loaded.")
        return True
    except Exception as e:
        st.error(f"Error setting up demo mode: {e}")
        return False

# Chemical Knowledge Graph - MVP

st.markdown("""
<div style="background: linear-gradient(90deg, #0066cc 0%, #1e3a5f 100%); padding: 1.5rem; border-radius: 10px; margin-bottom: 2rem;">
    <h3 style="color: white; margin: 0; text-align: center;">
        Upload chemistry papers • Extract structures • Ask intelligent questions
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
        with open(image_path, "rb") as file:
            return create_download_link(file.read(), filename, f"Download {filename}")
    return None

# Sidebar for file upload
with st.sidebar:
    st.markdown("### Upload Document")
    uploaded_file = st.file_uploader("Choose a PDF file", type="pdf", help="Upload a chemistry paper to analyze")
    
    # Demo mode toggle
    st.markdown("---")
    st.markdown("### Demo Mode")
    st.markdown("*For quick demonstration or if DECIMER models fail*")
    
    if st.button("Activate Demo Mode", help="Load sample chemical structures for demonstration"):
        if activate_demo_mode():
            st.rerun()
    
    if st.session_state.demo_mode:
        st.success("Demo mode active")
        if st.button("Exit Demo Mode"):
            st.session_state.demo_mode = False
            st.session_state.structures = []
            st.session_state.vectorstore = None
            st.rerun()
    
    if uploaded_file is not None:
        # Save uploaded file
        with open("temp_paper.pdf", "wb") as f:
            f.write(uploaded_file.getvalue())
        
        # Processing options
        st.markdown("### Processing Options")
        use_decimer_segmentation = st.checkbox("Use DECIMER Segmentation", value=True, help="Use AI-powered structure detection")
        extract_images = st.checkbox("Extract Structure Images", value=True, help="Save structure images for download")
        
        if st.button("Process Paper", help="Start processing the uploaded paper"):
            error_occurred = False
            processing_status = st.empty()
            
            try:
                with st.spinner("Processing paper..."):
                    processing_status.info("Initializing components...")
                    
                    # Initialize components with error handling
                    try:
                        pdf_processor = PDFProcessor("temp_paper.pdf")
                        processing_status.info("PDF processor initialized ✓")
                    except Exception as e:
                        st.error(f"Failed to initialize PDF processor: {e}")
                        st.info("💡 Try uploading a different PDF file or use demo mode.")
                        error_occurred = True
                    
                    if not error_occurred:
                        try:
                            decimer_client = DECIMERClient()
                            chem_handler = ChemicalHandler(decimer_client)
                            chunker = ChemicalAwareChunker()
                            processing_status.info("DECIMER components initialized ✓")
                        except Exception as e:
                            st.warning(f"DECIMER initialization failed: {e}")
                            st.info("💡 Consider activating Demo Mode for a quick demonstration.")
                            if st.button("Activate Demo Mode Now"):
                                if activate_demo_mode():
                                    st.rerun()
                            error_occurred = True
                    
                    if not error_occurred:
                        # Extract content
                        processing_status.info("Extracting text and images...")
                        try:
                            pages_data = pdf_processor.extract_text_and_images()
                            
                            # Extract full text
                            full_text = ""
                            for page_data in pages_data:
                                full_text += page_data['text'] + "\n"
                            
                            processing_status.info(f"Extracted {len(full_text)} characters of text ✓")
                            
                        except Exception as e:
                            st.error(f"Failed to extract content from PDF: {e}")
                            st.info("💡 The PDF might be corrupted or password-protected.")
                            error_occurred = True
                    
                    if not error_occurred:
                        # Process chemical structures
                        processing_status.info("Identifying chemical structures...")
                        all_structures = []
                        
                        try:
                            if use_decimer_segmentation:
                                processing_status.info("Using DECIMER segmentation...")
                                # Use DECIMER to segment and process structures
                                try:
                                    chemical_images = pdf_processor.extract_chemical_structures_with_decimer()
                                    
                                    for img_info in chemical_images:
                                        try:
                                            structure = chem_handler.process_image(
                                                img_info['path'],
                                                f"Chemical structure {img_info['index'] + 1}"
                                            )
                                            
                                            if structure:
                                                all_structures.append(structure)
                                        except Exception as e:
                                            st.warning(f"Failed to process structure {img_info['index'] + 1}: {e}")
                                            continue
                                            
                                except Exception as e:
                                    st.warning(f"DECIMER segmentation failed: {e}")
                                    st.info("Falling back to standard image extraction...")
                                    use_decimer_segmentation = False
                            
                            if not use_decimer_segmentation:
                                processing_status.info("Using standard image extraction...")
                                # Use standard image extraction
                                for page_data in pages_data:
                                    for img_info in page_data['images']:
                                        try:
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
                                        except Exception as e:
                                            st.warning(f"Failed to process image: {e}")
                                            continue
                            
                            st.session_state.structures = all_structures
                            processing_status.info(f"Found {len(all_structures)} chemical structures ✓")
                            
                        except Exception as e:
                            st.error(f"Structure processing failed: {e}")
                            st.info("💡 Try using Demo Mode to see the interface functionality.")
                            error_occurred = True
                    
                    if not error_occurred:
                        # Create chunks
                        processing_status.info("Creating document chunks...")
                        try:
                            chunks = chunker.chunk_with_structures(full_text, all_structures)
                            processing_status.info(f"Created {len(chunks)} document chunks ✓")
                        except Exception as e:
                            st.error(f"Failed to create document chunks: {e}")
                            error_occurred = True
                    
                    if not error_occurred:
                        # Create vector store
                        processing_status.info("Building vector database...")
                        try:
                            vector_store = ChemicalVectorStore()
                            vector_store.create_vectorstore(chunks)
                            st.session_state.vectorstore = vector_store
                            processing_status.info("Vector database created ✓")
                        except Exception as e:
                            st.error(f"Failed to create vector database: {e}")
                            st.info("💡 This might be an embedding model issue. Try demo mode.")
                            error_occurred = True
                
            except Exception as e:
                st.error(f"Unexpected error during processing: {e}")
                st.info("💡 Please try Demo Mode or check the logs for more details.")
                error_occurred = True
            
            finally:
                processing_status.empty()
                
                if not error_occurred:
                
                st.success("Paper processed successfully!")
                
                # Display processing results
                col_metrics1, col_metrics2, col_metrics3 = st.columns(3)
                with col_metrics1:
                    st.metric("Total Structures Found", len(all_structures))
                with col_metrics2:
                    st.metric("Processing Method", "DECIMER" if use_decimer_segmentation else "Standard")
                with col_metrics3:
                    if all_structures:
                        avg_mw = sum(s.get('molecular_weight', 0) for s in all_structures) / len(all_structures)
                        st.metric("Avg. Molecular Weight", f"{avg_mw:.1f}")

# Main interface
col1, col2 = st.columns([2, 1])

with col1:
    st.header("Ask Questions")
    
    if st.session_state.vectorstore is not None:
        # Initialize RAG with configured model
        model_name = os.getenv("OLLAMA_MODEL", "llama3.2:latest")
        rag = ChemicalRAG(st.session_state.vectorstore, model_name=model_name)
        
        # Query input with enhanced interface
        st.markdown("**What would you like to know about the paper?**")
        user_query = st.text_area(
            "Enter your question:", 
            height=80,
            placeholder="e.g., What are the main chemical compounds? What synthesis methods are described?"
        )
        
        # Quick question buttons
        st.markdown("**Quick Questions:**")
        quick_questions = [
            "What chemical compounds are discussed?",
            "What synthesis methods are described?",
            "What are the molecular weights?",
            "Compare the structures found"
        ]
        
        cols = st.columns(2)
        for i, question in enumerate(quick_questions):
            with cols[i % 2]:
                if st.button(question, key=f"quick_{i}"):
                    user_query = question
        
        if st.button("Ask Question", type="primary") and user_query:
            with st.spinner("Analyzing paper and generating answer..."):
                answer = rag.query(user_query)
                
                st.markdown("### Answer:")
                st.markdown(f'<div style="background-color: #f8fafc; padding: 1rem; border-radius: 8px; border-left: 4px solid #0066cc;">{answer}</div>', unsafe_allow_html=True)
                
                # Show relevant structures
                if st.session_state.structures:
                    st.markdown("### Related Chemical Structures:")
                    for i, struct in enumerate(st.session_state.structures[:3]):
                        with st.expander(f"Structure {i+1} - {struct.get('formula', 'Unknown')}"):
                            col_struct1, col_struct2 = st.columns([2, 1])
                            with col_struct1:
                                st.code(f"SMILES: {struct['smiles']}")
                                st.text(f"Formula: {struct['formula']}")
                                st.text(f"MW: {struct['molecular_weight']:.2f}")
                            with col_struct2:
                                if os.path.exists(struct.get('image_path', '')):
                                    st.image(struct['image_path'], width=120)
    else:
        st.info("Upload and process a paper to start asking questions about chemical structures and content.")
        
        # Show example capabilities
        st.markdown("### What you can ask:")
        capabilities = [
            "**Structure Analysis**: Identify and compare molecular structures",
            "**Synthesis Methods**: Understand chemical synthesis procedures", 
            "**Property Queries**: Get molecular weights, formulas, and properties",
            "**Content Summary**: Summarize key findings and conclusions",
            "**Reaction Analysis**: Understand chemical reactions described"
        ]
        
        for capability in capabilities:
            st.markdown(f"• {capability}")

with col2:
    st.header("Extracted Structures")
    
    if st.session_state.structures:
        # Download all structures as CSV
        csv_data = create_csv_download(st.session_state.structures)
        if csv_data:
            st.markdown(
                create_download_link(csv_data, "chemical_structures.csv", "Download All Structures (CSV)"),
                unsafe_allow_html=True
            )
        
        st.markdown("---")
        
        for i, struct in enumerate(st.session_state.structures):
            with st.expander(f"Structure {i+1}"):
                col_a, col_b = st.columns([3, 1])
                
                with col_a:
                    st.code(struct['smiles'])
                    st.text(f"Formula: {struct['formula']}")
                    st.text(f"MW: {struct['molecular_weight']:.2f}")
                    if struct.get('context'):
                        st.text(f"Context: {struct['context'][:100]}...")
                
                with col_b:
                    if os.path.exists(struct['image_path']):
                        st.image(struct['image_path'], width=150)
                        # Download link for individual image
                        img_download = get_image_download_link(
                            struct['image_path'], 
                            f"structure_{i+1}.png"
                        )
                        if img_download:
                            st.markdown(img_download, unsafe_allow_html=True)
    else:
        st.info("No structures extracted yet. Upload and process a paper to see results.")

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