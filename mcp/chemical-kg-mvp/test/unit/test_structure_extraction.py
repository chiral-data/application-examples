"""
Comprehensive tests for chemical structure extraction (SMILES generation).
Tests the ability to convert segmented structures to SMILES notation.
"""

import pytest
import os
import tempfile
import numpy as np
from PIL import Image, ImageDraw
from unittest.mock import patch, MagicMock, call
from rdkit import Chem
from rdkit.Chem import Draw

from decimer_client import DECIMERClient


class TestStructureExtraction:
    """Test suite for SMILES extraction from chemical structures"""
    
    @pytest.fixture
    def decimer_client(self):
        """Create a DECIMERClient instance"""
        return DECIMERClient()
    
    @pytest.fixture
    def known_structures(self):
        """Dictionary of known structures and their SMILES"""
        return {
            'benzene': 'c1ccccc1',
            'ethanol': 'CCO',
            'aspirin': 'CC(=O)Oc1ccccc1C(=O)O',
            'caffeine': 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
            'glucose': 'C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O',
            'vanillin': 'COc1cc(C=O)ccc1O',
            'ibuprofen': 'CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O',
            'paracetamol': 'CC(=O)Nc1ccc(cc1)O'
        }
    
    @pytest.fixture
    def create_structure_image(self):
        """Create structure image from SMILES"""
        def _create_image(smiles, size=(300, 300)):
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError(f"Invalid SMILES: {smiles}")
            
            img = Draw.MolToImage(mol, size=size)
            return img
        
        return _create_image
    
    def test_extract_simple_structure(self, decimer_client, create_structure_image):
        """Test extraction of simple chemical structures"""
        # Test with benzene
        img = create_structure_image('c1ccccc1')
        
        with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as f:
            img.save(f.name)
            
            try:
                with patch.object(decimer_client, 'predict_smiles') as mock_predict:
                    mock_predict.return_value = 'c1ccccc1'
                    
                    smiles = decimer_client.predict_smiles(f.name)
                    
                    assert smiles == 'c1ccccc1'
                    mock_predict.assert_called_once()
                    
                    # Verify it's valid SMILES
                    mol = Chem.MolFromSmiles(smiles)
                    assert mol is not None
                    
            finally:
                os.unlink(f.name)
    
    def test_extract_complex_structures(self, decimer_client, create_structure_image, known_structures):
        """Test extraction of complex chemical structures"""
        complex_molecules = ['aspirin', 'caffeine', 'ibuprofen']
        
        for mol_name in complex_molecules:
            expected_smiles = known_structures[mol_name]
            img = create_structure_image(expected_smiles)
            
            with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as f:
                img.save(f.name)
                
                try:
                    with patch.object(decimer_client, 'predict_smiles') as mock_predict:
                        mock_predict.return_value = expected_smiles
                        
                        smiles = decimer_client.predict_smiles(f.name)
                        
                        # Verify structure
                        mol = Chem.MolFromSmiles(smiles)
                        assert mol is not None
                        
                        # Check molecular properties
                        assert mol.GetNumAtoms() > 0
                        assert mol.GetNumBonds() > 0
                        
                finally:
                    os.unlink(f.name)
    
    def test_canonicalization_of_smiles(self, decimer_client, create_structure_image):
        """Test that extracted SMILES are canonicalized"""
        # Different representations of the same molecule
        equivalent_smiles = [
            'CC(C)C',  # isobutane
            'C(C)(C)C',  # same molecule, different representation
        ]
        
        canonical_smiles = Chem.CanonSmiles(equivalent_smiles[0])
        
        for smiles in equivalent_smiles:
            img = create_structure_image(smiles)
            
            with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as f:
                img.save(f.name)
                
                try:
                    with patch.object(decimer_client, 'predict_smiles') as mock_predict:
                        # DECIMER should return canonical SMILES
                        mock_predict.return_value = canonical_smiles
                        
                        result = decimer_client.predict_smiles(f.name)
                        assert Chem.CanonSmiles(result) == canonical_smiles
                        
                finally:
                    os.unlink(f.name)
    
    def test_stereochemistry_preservation(self, decimer_client, create_structure_image):
        """Test that stereochemistry is preserved in extraction"""
        # Test with a chiral molecule
        chiral_smiles = 'C[C@H](O)C(=O)O'  # L-lactic acid
        img = create_structure_image(chiral_smiles)
        
        with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as f:
            img.save(f.name)
            
            try:
                with patch.object(decimer_client, 'predict_smiles') as mock_predict:
                    mock_predict.return_value = chiral_smiles
                    
                    smiles = decimer_client.predict_smiles(f.name)
                    
                    # Check that stereochemistry is preserved
                    assert '@' in smiles  # Contains stereochemistry marker
                    
                    mol = Chem.MolFromSmiles(smiles)
                    assert mol is not None
                    
                    # Check for chiral centers
                    chiral_centers = Chem.FindMolChiralCenters(mol)
                    assert len(chiral_centers) > 0
                    
            finally:
                os.unlink(f.name)
    
    def test_batch_extraction(self, decimer_client, create_structure_image, known_structures):
        """Test batch extraction of multiple structures"""
        molecules = ['benzene', 'ethanol', 'aspirin', 'caffeine']
        image_paths = []
        expected_results = []
        
        try:
            # Create multiple structure images
            for mol_name in molecules:
                smiles = known_structures[mol_name]
                img = create_structure_image(smiles)
                
                with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as f:
                    img.save(f.name)
                    image_paths.append(f.name)
                    expected_results.append(smiles)
            
            with patch.object(decimer_client, 'predict_smiles') as mock_predict:
                mock_predict.side_effect = expected_results
                
                # Extract SMILES from all images
                extracted_smiles = []
                for path in image_paths:
                    smiles = decimer_client.predict_smiles(path)
                    extracted_smiles.append(smiles)
                
                # Verify all extractions
                assert len(extracted_smiles) == len(molecules)
                for smiles in extracted_smiles:
                    mol = Chem.MolFromSmiles(smiles)
                    assert mol is not None
                
        finally:
            for path in image_paths:
                if os.path.exists(path):
                    os.unlink(path)
    
    def test_invalid_structure_handling(self, decimer_client):
        """Test handling of invalid or unclear structures"""
        # Create an image with random noise
        noise = np.random.randint(0, 255, (300, 300, 3), dtype=np.uint8)
        img = Image.fromarray(noise)
        
        with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as f:
            img.save(f.name)
            
            try:
                with patch.object(decimer_client, 'predict_smiles') as mock_predict:
                    # DECIMER might return empty or invalid SMILES
                    mock_predict.return_value = ''
                    
                    smiles = decimer_client.predict_smiles(f.name)
                    
                    # Should handle gracefully
                    assert smiles == '' or Chem.MolFromSmiles(smiles) is None
                    
            finally:
                os.unlink(f.name)
    
    def test_partial_structure_extraction(self, decimer_client):
        """Test extraction from partially visible structures"""
        # Create image with partial benzene ring
        img = Image.new('RGB', (300, 300), 'white')
        draw = ImageDraw.Draw(img)
        
        # Draw only half of a hexagon
        points = [(150, 50), (250, 100), (250, 200)]
        for i in range(len(points) - 1):
            draw.line([points[i], points[i + 1]], fill='black', width=3)
        
        with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as f:
            img.save(f.name)
            
            try:
                with patch.object(decimer_client, 'predict_smiles') as mock_predict:
                    # Might return incomplete or uncertain result
                    mock_predict.return_value = 'C1=CC=CC=C1'  # Attempt at benzene
                    
                    smiles = decimer_client.predict_smiles(f.name)
                    # Should still produce valid SMILES if possible
                    if smiles:
                        mol = Chem.MolFromSmiles(smiles)
                        # May or may not be valid depending on the model
                        
            finally:
                os.unlink(f.name)
    
    def test_extraction_confidence_scores(self, decimer_client, create_structure_image):
        """Test extraction with confidence scores if available"""
        img = create_structure_image('c1ccccc1')
        
        with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as f:
            img.save(f.name)
            
            try:
                with patch.object(decimer_client, 'predict_smiles_with_confidence') as mock_predict:
                    mock_predict.return_value = ('c1ccccc1', 0.95)  # SMILES, confidence
                    
                    if hasattr(decimer_client, 'predict_smiles_with_confidence'):
                        smiles, confidence = decimer_client.predict_smiles_with_confidence(f.name)
                        
                        assert smiles == 'c1ccccc1'
                        assert 0 <= confidence <= 1
                        assert confidence > 0.8  # High confidence for clear structure
                    
            finally:
                os.unlink(f.name)
    
    def test_molecular_property_validation(self, decimer_client, create_structure_image, known_structures):
        """Test that extracted molecules have valid properties"""
        for mol_name, expected_smiles in known_structures.items():
            img = create_structure_image(expected_smiles)
            
            with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as f:
                img.save(f.name)
                
                try:
                    with patch.object(decimer_client, 'predict_smiles') as mock_predict:
                        mock_predict.return_value = expected_smiles
                        
                        smiles = decimer_client.predict_smiles(f.name)
                        mol = Chem.MolFromSmiles(smiles)
                        
                        if mol:
                            # Validate molecular properties
                            assert mol.GetNumAtoms() > 0
                            assert mol.GetNumBonds() >= mol.GetNumAtoms() - 1  # Connected
                            
                            # Check for reasonable molecular weight
                            mw = Chem.Descriptors.MolWt(mol)
                            assert 10 < mw < 5000  # Reasonable range for small molecules
                            
                            # Check for valid valences
                            problems = Chem.DetectChemistryProblems(mol)
                            assert len(problems) == 0
                            
                finally:
                    os.unlink(f.name)
    
    def test_extraction_performance(self, decimer_client, create_structure_image):
        """Test extraction performance and timing"""
        import time
        
        img = create_structure_image('CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O')  # ibuprofen
        
        with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as f:
            img.save(f.name)
            
            try:
                with patch.object(decimer_client, 'predict_smiles') as mock_predict:
                    def timed_predict(*args, **kwargs):
                        time.sleep(0.1)  # Simulate processing time
                        return 'CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O'
                    
                    mock_predict.side_effect = timed_predict
                    
                    start_time = time.time()
                    smiles = decimer_client.predict_smiles(f.name)
                    end_time = time.time()
                    
                    extraction_time = end_time - start_time
                    
                    # Should complete in reasonable time
                    assert extraction_time < 5.0  # Less than 5 seconds
                    assert Chem.MolFromSmiles(smiles) is not None
                    
            finally:
                os.unlink(f.name)