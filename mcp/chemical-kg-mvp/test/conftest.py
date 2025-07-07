import pytest
import os
import sys
import tempfile
import shutil
from unittest.mock import MagicMock, patch
from PIL import Image
import io

# Add src directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files"""
    temp_dir = tempfile.mkdtemp()
    yield temp_dir
    shutil.rmtree(temp_dir)

@pytest.fixture
def mock_pdf_path(temp_dir):
    """Create a mock PDF file path"""
    pdf_path = os.path.join(temp_dir, "test.pdf")
    # Create a dummy file
    with open(pdf_path, 'wb') as f:
        f.write(b'%PDF-1.4\nDummy PDF content for testing')
    return pdf_path

@pytest.fixture
def mock_image_path(temp_dir):
    """Create a mock chemical structure image"""
    img_path = os.path.join(temp_dir, "test_structure.png")
    
    # Create a simple test image
    img = Image.new('RGB', (200, 150), color='white')
    img.save(img_path, 'PNG')
    
    return img_path

@pytest.fixture
def sample_chemical_structure():
    """Sample chemical structure data"""
    return {
        'smiles': 'CCO',
        'formula': 'C2H6O',
        'molecular_weight': 46.07,
        'context': 'ethanol molecule found in figure 1',
        'image_path': '/path/to/ethanol.png',
        'page': 0
    }

@pytest.fixture
def sample_structures_list(sample_chemical_structure):
    """List of sample chemical structures"""
    return [
        sample_chemical_structure,
        {
            'smiles': 'C6H6',
            'formula': 'C6H6',
            'molecular_weight': 78.11,
            'context': 'benzene ring structure',
            'image_path': '/path/to/benzene.png',
            'page': 1
        },
        {
            'smiles': 'CC(=O)O',
            'formula': 'C2H4O2',
            'molecular_weight': 60.05,
            'context': 'acetic acid shown in scheme 2',
            'image_path': '/path/to/acetic_acid.png',
            'page': 1
        }
    ]

@pytest.fixture
def sample_text():
    """Sample scientific text with chemical content"""
    return """
    Introduction
    
    This paper describes the synthesis of several organic compounds including ethanol (CCO) 
    and benzene (C6H6). The synthesis was carried out using standard organic chemistry methods.
    
    Experimental Section
    
    Compound 1 (ethanol) was synthesized by fermentation of glucose. The molecular formula 
    is C2H6O with a molecular weight of 46.07 g/mol. The structure is shown in Figure 1.
    
    Compound 2 (benzene) has the formula C6H6 and molecular weight 78.11 g/mol. The aromatic 
    ring structure is depicted in Figure 2. This compound was purified by distillation.
    
    Results and Discussion
    
    The NMR spectra confirmed the structures of both compounds. Ethanol showed the expected 
    triplet at 1.2 ppm and quartet at 3.7 ppm in the 1H NMR spectrum.
    
    Conclusion
    
    Both compounds were successfully synthesized and characterized using spectroscopic methods.
    """

@pytest.fixture
def mock_decimer_response():
    """Mock DECIMER API response"""
    return {
        'smiles': 'CCO',
        'confidence': 0.95
    }

@pytest.fixture
def mock_ollama_response():
    """Mock Ollama LLM response"""
    return {
        'result': 'Based on the provided context, ethanol (CCO) is an important organic compound with the molecular formula C2H6O.',
        'source_documents': []
    }

@pytest.fixture
def mock_rdkit():
    """Mock RDKit functionality"""
    with patch('rdkit.Chem.MolFromSmiles') as mock_mol:
        mock_molecule = MagicMock()
        mock_molecule.GetNumAtoms.return_value = 9
        mock_mol.return_value = mock_molecule
        
        with patch('rdkit.Chem.Descriptors.MolWt') as mock_mw:
            mock_mw.return_value = 46.07
            
            with patch('rdkit.Chem.rdMolDescriptors.CalcMolFormula') as mock_formula:
                mock_formula.return_value = 'C2H6O'
                yield {
                    'mol': mock_mol,
                    'mw': mock_mw,
                    'formula': mock_formula
                }

@pytest.fixture
def mock_pdf_document():
    """Mock PyMuPDF document"""
    mock_doc = MagicMock()
    mock_doc.page_count = 2
    mock_doc.is_encrypted = False
    mock_doc.is_pdf = True
    mock_doc.needs_pass = False
    mock_doc.metadata = {
        'title': 'Test Chemical Paper',
        'author': 'Test Author',
        'subject': 'Organic Chemistry',
        'creator': 'Test Creator'
    }
    
    # Mock pages
    mock_page1 = MagicMock()
    mock_page1.get_text.return_value = "Page 1 content with ethanol structure"
    mock_page1.get_images.return_value = [(1, 0, 100, 100, 8, 'DeviceRGB', '', 'Im1', 'DCTDecode')]
    mock_page1.rect = MagicMock()
    mock_page1.rect.width = 612
    mock_page1.rect.height = 792
    
    mock_page2 = MagicMock()
    mock_page2.get_text.return_value = "Page 2 content with benzene structure"
    mock_page2.get_images.return_value = [(2, 0, 120, 120, 8, 'DeviceRGB', '', 'Im2', 'DCTDecode')]
    mock_page2.rect = MagicMock()
    mock_page2.rect.width = 612
    mock_page2.rect.height = 792
    
    mock_doc.__iter__ = lambda self: iter([mock_page1, mock_page2])
    mock_doc.__getitem__ = lambda self, idx: [mock_page1, mock_page2][idx]
    
    return mock_doc

@pytest.fixture
def mock_chroma_vectorstore():
    """Mock ChromaDB vector store"""
    mock_vectorstore = MagicMock()
    mock_vectorstore.similarity_search.return_value = [
        MagicMock(page_content="Test chunk with ethanol", metadata={'chunk_id': 0}),
        MagicMock(page_content="Test chunk with synthesis", metadata={'chunk_id': 1})
    ]
    mock_vectorstore.as_retriever.return_value = MagicMock()
    return mock_vectorstore

@pytest.fixture(autouse=True)
def mock_environment_variables():
    """Set up test environment variables"""
    test_env = {
        'CHUNK_SIZE': '500',
        'CHUNK_OVERLAP': '50',
        'CHEMICAL_CONTEXT_WINDOW': '3',
        'OLLAMA_HOST': 'localhost',
        'OLLAMA_PORT': '11434',
        'DECIMER_API_URL': 'http://test-decimer-api.com',
        'DECIMER_TIMEOUT': '10',
        'DECIMER_MAX_RETRIES': '2',
        'MAX_IMAGE_SIZE_MB': '5'
    }
    
    with patch.dict(os.environ, test_env):
        yield test_env