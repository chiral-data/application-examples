# app.py
# Fix for Python 3.8 type annotation compatibility
from __future__ import annotations

# Disable ChromaDB telemetry to avoid posthog issues
import os
os.environ["CHROMA_TELEMETRY"] = "false"
os.environ["ANONYMIZED_TELEMETRY"] = "false"

# Fix for SQLite version compatibility with ChromaDB
try:
    __import__('pysqlite3')
    import sys
    sys.modules['sqlite3'] = sys.modules.pop('pysqlite3')
except ImportError:
    pass

import streamlit as st
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
from llm_call import format_rag_response

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
        background-color: #0066cc !important;
        color: white !important;
        border-radius: 6px;
        padding: 0.5rem 1rem;
        text-decoration: none !important;
        font-weight: 500;
        margin: 0.25rem;
        display: inline-block;
        transition: all 0.3s ease;
        border: 2px solid #0066cc;
    }
    
    .download-button:hover {
        background-color: #0052a3 !important;
        border-color: #0052a3;
        box-shadow: 0 4px 8px rgba(0, 102, 204, 0.3);
        color: white !important;
        text-decoration: none !important;
    }
    
    .download-button:visited {
        color: white !important;
        text-decoration: none !important;
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
if 'chat_history' not in st.session_state:
    st.session_state.chat_history = []
if 'rag_chain' not in st.session_state:
    st.session_state.rag_chain = None


# Chemical Knowledge Graph - MVP

st.markdown("""
<div style="background: linear-gradient(90deg, #0066cc 0%, #1e3a5f 100%); padding: 1.5rem; border-radius: 10px; margin-bottom: 2rem;">
    <h3 style="color: white; margin: 0; text-align: center;">
        Chemistry Paper Analyzer
    </h3>
</div>
""", unsafe_allow_html=True)

st.markdown("""
**Transform your chemical literature into an intelligent knowledge base.**

This tool uses AI to extract chemical structures, analyze content, and answer questions about your papers.
""")

# Helper functions for downloads
def create_download_link(val, filename, link_text):
    """Create a download link for data"""
    b64 = base64.b64encode(val).decode()
    return f'<a href="data:application/octet-stream;base64,{b64}" download="{filename}" class="download-button">{link_text}</a>'

def create_csv_download(structures):
    """Create CSV data from structures"""
    if not structures:
        return None
    
    data = []
    for i, struct in enumerate(structures):
        # Extract page number from context if available
        page_num = ''
        context = struct.get('context', '')
        if 'page' in context.lower():
            # Try to extract page number from context like "Structure from page 2"
            import re
            page_match = re.search(r'page\s+(\d+)', context.lower())
            if page_match:
                page_num = page_match.group(1)
        
        data.append({
            'Structure_ID': f'Structure_{i+1}',
            'SMILES': struct.get('smiles', ''),
            'Molecular_Formula': struct.get('formula', ''),
            'Molecular_Weight': struct.get('molecular_weight', ''),
            'Page': page_num
        })
    
    df = pd.DataFrame(data)
    return df.to_csv(index=False).encode('utf-8')

def get_image_download_link(image_path, filename):
    """Create download link for images"""
    if os.path.exists(image_path):
        with open(image_path, "rb") as file:
            return create_download_link(file.read(), filename, f"Download {filename}")
    return None

# Sidebar for file upload
with st.sidebar:
    st.markdown("### Upload Document")
    uploaded_file = st.file_uploader("Choose a PDF file", type="pdf", help="Upload a chemistry paper to analyze")
    
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
            
            # Add progress log below spinner
            progress_container = st.container()
            with progress_container:
                progress_log = []
                log_display = st.empty()
                
                def update_progress(message):
                    progress_log.append(message)
                    # Keep last 8 messages
                    display_text = "\n".join(progress_log[-8:])
                    log_display.text(display_text)
            
            try:
                with st.spinner("Processing paper..."):
                    processing_status.info("Initializing components...")
                    update_progress("Starting PDF processing...")
                    
                    # Initialize components with error handling
                    try:
                        pdf_processor = PDFProcessor("temp_paper.pdf")
                        processing_status.info("PDF processor initialized")
                        update_progress("PDF processor ready")
                    except Exception as e:
                        st.error(f"Failed to initialize PDF processor: {e}")
                        st.info("Try uploading a different PDF file.")
                        error_occurred = True
                    
                    if not error_occurred:
                        try:
                            decimer_client = DECIMERClient()
                            chem_handler = ChemicalHandler(decimer_client)
                            chunker = ChemicalAwareChunker()
                            processing_status.info("DECIMER components initialized")
                            update_progress("DECIMER components ready")
                            update_progress(f"GPU available: {decimer_client.decimer_available}")
                        except Exception as e:
                            st.warning(f"DECIMER initialization failed: {e}")
                            st.info("Try uploading a different PDF file or check the logs for more details.")
                            error_occurred = True
                    
                    if not error_occurred:
                        # Extract content
                        processing_status.info("Extracting text and images...")
                        update_progress("Extracting content from PDF...")
                        try:
                            pages_data = pdf_processor.extract_text_and_images()
                            
                            # Extract full text
                            full_text = ""
                            for page_data in pages_data:
                                full_text += page_data['text'] + "\n"
                            
                            processing_status.info(f"Extracted {len(full_text)} characters of text")
                            update_progress(f"Extracted {len(pages_data)} pages, {len(full_text):,} characters")
                            
                        except Exception as e:
                            st.error(f"Failed to extract content from PDF: {e}")
                            st.info("The PDF might be corrupted or password-protected.")
                            error_occurred = True
                    
                    if not error_occurred:
                        # Process chemical structures
                        processing_status.info("Identifying chemical structures...")
                        all_structures = []
                        
                        try:
                            if use_decimer_segmentation:
                                processing_status.info("Using DECIMER segmentation...")
                                update_progress("Starting DECIMER AI segmentation...")
                                # Use DECIMER to segment and process structures
                                try:
                                    chemical_images = pdf_processor.extract_chemical_structures_with_decimer()
                                    update_progress(f"Found {len(chemical_images)} chemical structures")
                                    
                                    for i, img_info in enumerate(chemical_images):
                                        try:
                                            # Update progress every few structures
                                            if i % 3 == 0:
                                                update_progress(f"Processing structure {i+1}/{len(chemical_images)}...")
                                            
                                            structure = chem_handler.process_image(
                                                img_info['path'],
                                                f"Structure from page {img_info.get('page', 'unknown')}"
                                            )
                                            
                                            if structure:
                                                all_structures.append(structure)
                                                if structure.get('smiles') and (i % 3 == 0 or i == len(chemical_images)-1):
                                                    update_progress(f"Generated SMILES for structure {i+1}")
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
                            processing_status.info(f"Found {len(all_structures)} chemical structures")
                            update_progress(f"Processing complete! Found {len(all_structures)} structures")
                            
                        except Exception as e:
                            st.error(f"Structure processing failed: {e}")
                            st.info("Try uploading a different PDF file or check GPU acceleration.")
                            error_occurred = True
                    
                    if not error_occurred:
                        # Create chunks
                        processing_status.info("Creating document chunks...")
                        try:
                            chunks = chunker.chunk_with_structures(full_text, all_structures)
                            processing_status.info(f"Created {len(chunks)} document chunks")
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
                            processing_status.info("Vector database created")
                        except Exception as e:
                            st.error(f"Failed to create vector database: {e}")
                            st.info("This might be an embedding model issue. Check the logs for more details.")
                            error_occurred = True
                
            except Exception as e:
                st.error(f"Unexpected error during processing: {e}")
                st.info("Please check the logs for more details or try a different PDF.")
                error_occurred = True
            
            finally:
                processing_status.empty()
                
                if not error_occurred:
                    st.success("Paper processed successfully!")
                    
                    # Display processing results
                    st.metric("Total Structures Found", len(all_structures))
    
    # Real-time progress logging provides visibility into processing steps

# Main interface
st.header("Ask Questions")

if st.session_state.vectorstore is not None:
    # Initialize RAG with configured model (once)
    if st.session_state.rag_chain is None:
        model_name = os.getenv("OLLAMA_MODEL", "llama3.2:latest")
        st.session_state.rag_chain = ChemicalRAG(st.session_state.vectorstore, model_name=model_name)
    
    # Chat interface
    st.markdown("**Chat with your document:**")
    
    # Example questions (text only)
    st.markdown("*Example questions you can ask:*")
    st.markdown("""
    • What are the main chemical compounds discussed?
    
    • What synthesis methods are described?
    
    • What are the molecular weights of the compounds?
    
    • Compare the structures found in this paper
    
    • Explain the reaction mechanisms described
    """)
    
    # Chat history display
    chat_container = st.container()
    with chat_container:
        if st.session_state.chat_history:
            st.markdown("### Conversation:")
            for i, (question, answer) in enumerate(st.session_state.chat_history):
                # User question
                st.markdown(f'<div style="background-color: #e3f2fd; padding: 0.8rem; border-radius: 8px; margin: 0.5rem 0; border-left: 4px solid #2196f3;"><strong>You:</strong> {question}</div>', unsafe_allow_html=True)
                
                # Assistant answer
                st.markdown(f'<div style="background-color: #f8fafc; padding: 0.8rem; border-radius: 8px; margin: 0.5rem 0; border-left: 4px solid #0066cc;"><strong>Assistant:</strong> {answer}</div>', unsafe_allow_html=True)
                
                # Add some spacing
                st.markdown("---")
    
    # Chat input at the bottom
    with st.form(key="chat_form", clear_on_submit=True):
        user_input = st.text_area(
            "Ask a question about the paper:",
            height=100,
            placeholder="Type your question here...",
            key="chat_input"
        )
        submit_button = st.form_submit_button("Send", type="primary")
        
        if submit_button and user_input.strip():
            with st.spinner("Thinking..."):
                try:
                    raw_answer = st.session_state.rag_chain.query(user_input.strip())
                    # Format the answer using LLM for better readability
                    formatted_answer = format_rag_response(raw_answer)
                    # Add to chat history
                    st.session_state.chat_history.append((user_input.strip(), formatted_answer))
                    st.rerun()
                except Exception as e:
                    if "Connection refused" in str(e) or "11434" in str(e):
                        st.error("Ollama LLM service is not available. Please ensure Ollama is installed and running.")
                        st.info("The chemical structure extraction is working, but Q&A requires Ollama to be set up.")
                    else:
                        st.error(f"Error generating answer: {e}")
    
    # Clear chat button
    if st.session_state.chat_history:
        if st.button("Clear Chat History", key="clear_chat"):
            st.session_state.chat_history = []
            st.rerun()
else:
    st.info("Upload and process a paper to start chatting about chemical structures and content.")
    
    # Show capabilities without buttons
    st.markdown("### What you can ask:")
    
    st.markdown("""
    **Structure Analysis:** Identify and compare molecular structures
    - What are the main chemical compounds discussed in this paper?
    
    **Synthesis Methods:** Understand chemical synthesis procedures  
    - What synthesis methods are described?
    
    **Property Queries:** Get molecular weights, formulas, and properties
    - What are the molecular weights of the compounds found?
    
    **Content Summary:** Summarize key findings and conclusions
    - What are the key findings about the chemical structures?
    
    **Reaction Analysis:** Understand chemical reactions described
    - What reactions are mentioned in the paper?
    """)
    

# Extracted Structures section
st.markdown("---")
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
    
    # Display structures in a grid layout
    cols_per_row = 2
    for i in range(0, len(st.session_state.structures), cols_per_row):
        cols = st.columns(cols_per_row)
        for j, struct in enumerate(st.session_state.structures[i:i+cols_per_row]):
            with cols[j]:
                with st.expander(f"Structure {i+j+1}"):
                    # Structure image
                    if os.path.exists(struct.get('image_path', '')):
                        st.image(struct['image_path'], width=200)
                        # Download link for individual image
                        img_download = get_image_download_link(
                            struct['image_path'], 
                            f"structure_{i+j+1}.png"
                        )
                        if img_download:
                            st.markdown(img_download, unsafe_allow_html=True)
                    
                    # Structure details
                    st.code(struct['smiles'])
                    st.text(f"Formula: {struct['formula']}")
                    st.text(f"MW: {struct['molecular_weight']:.2f}")
else:
    st.info("No structures extracted yet. Upload and process a paper to see results.")