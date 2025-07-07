import pytest
from unittest.mock import Mock, patch, MagicMock

from chemical_handler import ChemicalHandler


class TestChemicalHandler:
    """Test suite for ChemicalHandler class"""
    
    def test_init(self):
        """Test ChemicalHandler initialization"""
        mock_decimer_client = Mock()
        handler = ChemicalHandler(mock_decimer_client)
        
        assert handler.decimer == mock_decimer_client
    
    def test_process_image_success(self, mock_rdkit, mock_image_path):
        """Test successful image processing"""
        mock_decimer_client = Mock()
        mock_decimer_client.image_to_smiles.return_value = 'CCO'
        
        handler = ChemicalHandler(mock_decimer_client)
        
        with patch('rdkit.Chem.MolFromSmiles') as mock_mol_from_smiles:
            mock_mol = Mock()
            mock_mol_from_smiles.return_value = mock_mol
            
            with patch('rdkit.Chem.Descriptors.MolWt', return_value=46.07):
                with patch('rdkit.Chem.rdMolDescriptors.CalcMolFormula', return_value='C2H6O'):
                    result = handler.process_image(mock_image_path, "ethanol structure")
                    
                    assert result is not None
                    assert result['smiles'] == 'CCO'
                    assert result['molecular_weight'] == 46.07
                    assert result['formula'] == 'C2H6O'
                    assert result['context'] == "ethanol structure"
                    assert result['image_path'] == mock_image_path
    
    def test_process_image_no_smiles(self, mock_image_path):
        """Test image processing when DECIMER returns no SMILES"""
        mock_decimer_client = Mock()
        mock_decimer_client.image_to_smiles.return_value = None
        
        handler = ChemicalHandler(mock_decimer_client)
        result = handler.process_image(mock_image_path)
        
        assert result is None
    
    def test_process_image_invalid_smiles(self, mock_image_path):
        """Test image processing with invalid SMILES"""
        mock_decimer_client = Mock()
        mock_decimer_client.image_to_smiles.return_value = 'INVALID_SMILES'
        
        handler = ChemicalHandler(mock_decimer_client)
        
        with patch('rdkit.Chem.MolFromSmiles', return_value=None):
            result = handler.process_image(mock_image_path)
            
            assert result is None
    
    def test_process_image_with_rdkit_error(self, mock_image_path):
        """Test image processing when RDKit throws an error"""
        mock_decimer_client = Mock()
        mock_decimer_client.image_to_smiles.return_value = 'CCO'
        
        handler = ChemicalHandler(mock_decimer_client)
        
        with patch('rdkit.Chem.MolFromSmiles', side_effect=Exception("RDKit error")):
            result = handler.process_image(mock_image_path)
            
            assert result is None
    
    def test_get_structure_context_basic(self):
        """Test basic structure context extraction"""
        mock_decimer_client = Mock()
        handler = ChemicalHandler(mock_decimer_client)
        
        page_text = "Line 1\\nLine 2\\nLine 3\\nLine 4\\nLine 5\\nLine 6\\nLine 7"
        mock_bbox = Mock()  # Bounding box (not used in simple implementation)
        
        context = handler.get_structure_context(page_text, mock_bbox)
        
        assert isinstance(context, str)
        assert len(context) > 0
    
    def test_get_structure_context_empty_text(self):
        """Test context extraction with empty text"""
        mock_decimer_client = Mock()
        handler = ChemicalHandler(mock_decimer_client)
        
        context = handler.get_structure_context("", Mock())
        
        assert context == ""
    
    def test_process_image_empty_context(self, mock_rdkit, mock_image_path):
        """Test processing image with empty context"""
        mock_decimer_client = Mock()
        mock_decimer_client.image_to_smiles.return_value = 'CCO'
        
        handler = ChemicalHandler(mock_decimer_client)
        
        with patch('rdkit.Chem.MolFromSmiles') as mock_mol_from_smiles:
            mock_mol = Mock()
            mock_mol_from_smiles.return_value = mock_mol
            
            with patch('rdkit.Chem.Descriptors.MolWt', return_value=46.07):
                with patch('rdkit.Chem.rdMolDescriptors.CalcMolFormula', return_value='C2H6O'):
                    result = handler.process_image(mock_image_path, "")
                    
                    assert result is not None
                    assert result['context'] == ""
    
    def test_process_image_complex_molecule(self, mock_image_path):
        """Test processing a more complex molecule"""
        mock_decimer_client = Mock()
        mock_decimer_client.image_to_smiles.return_value = 'CC(C)(C)C1=CC=C(C=C1)O'  # tert-butylphenol
        
        handler = ChemicalHandler(mock_decimer_client)
        
        with patch('rdkit.Chem.MolFromSmiles') as mock_mol_from_smiles:
            mock_mol = Mock()
            mock_mol_from_smiles.return_value = mock_mol
            
            with patch('rdkit.Chem.Descriptors.MolWt', return_value=150.22):
                with patch('rdkit.Chem.rdMolDescriptors.CalcMolFormula', return_value='C10H14O'):
                    result = handler.process_image(mock_image_path, "phenol derivative")
                    
                    assert result is not None
                    assert result['smiles'] == 'CC(C)(C)C1=CC=C(C=C1)O'
                    assert result['molecular_weight'] == 150.22
                    assert result['formula'] == 'C10H14O'
    
    def test_get_structure_context_multiline(self):
        """Test context extraction with multiline text"""
        mock_decimer_client = Mock()
        handler = ChemicalHandler(mock_decimer_client)
        
        page_text = """This is the introduction section.
        
        Figure 1 shows the structure of ethanol.
        The molecular formula is C2H6O.
        This compound is commonly used as a solvent.
        
        The synthesis procedure is described below.
        Step 1: Mix the reactants.
        Step 2: Heat the mixture."""
        
        mock_bbox = Mock()
        context = handler.get_structure_context(page_text, mock_bbox)
        
        # Should return first few lines
        assert "introduction" in context.lower()
        assert isinstance(context, str)
    
    def test_process_image_special_characters_in_smiles(self, mock_image_path):
        """Test processing SMILES with special characters"""
        mock_decimer_client = Mock()
        # SMILES with brackets, charges, etc.
        mock_decimer_client.image_to_smiles.return_value = '[Na+].[Cl-]'
        
        handler = ChemicalHandler(mock_decimer_client)
        
        with patch('rdkit.Chem.MolFromSmiles') as mock_mol_from_smiles:
            mock_mol = Mock()
            mock_mol_from_smiles.return_value = mock_mol
            
            with patch('rdkit.Chem.Descriptors.MolWt', return_value=58.44):
                with patch('rdkit.Chem.rdMolDescriptors.CalcMolFormula', return_value='ClNa'):
                    result = handler.process_image(mock_image_path)
                    
                    assert result is not None
                    assert result['smiles'] == '[Na+].[Cl-]'
    
    def test_process_image_very_long_context(self, mock_rdkit, mock_image_path):
        """Test processing with very long context text"""
        mock_decimer_client = Mock()
        mock_decimer_client.image_to_smiles.return_value = 'CCO'
        
        handler = ChemicalHandler(mock_decimer_client)
        
        # Very long context
        long_context = "This is a very long context. " * 100
        
        with patch('rdkit.Chem.MolFromSmiles') as mock_mol_from_smiles:
            mock_mol = Mock()
            mock_mol_from_smiles.return_value = mock_mol
            
            with patch('rdkit.Chem.Descriptors.MolWt', return_value=46.07):
                with patch('rdkit.Chem.rdMolDescriptors.CalcMolFormula', return_value='C2H6O'):
                    result = handler.process_image(mock_image_path, long_context)
                    
                    assert result is not None
                    assert result['context'] == long_context
    
    def test_rdkit_import_error_handling(self, mock_image_path):
        """Test handling when RDKit is not available"""
        mock_decimer_client = Mock()
        mock_decimer_client.image_to_smiles.return_value = 'CCO'
        
        handler = ChemicalHandler(mock_decimer_client)
        
        # Mock import error for RDKit
        with patch('rdkit.Chem.MolFromSmiles', side_effect=ImportError("RDKit not installed")):
            result = handler.process_image(mock_image_path)
            
            # Should gracefully handle the import error
            assert result is None
    
    def test_molecular_weight_calculation_error(self, mock_image_path):
        """Test handling molecular weight calculation errors"""
        mock_decimer_client = Mock()
        mock_decimer_client.image_to_smiles.return_value = 'CCO'
        
        handler = ChemicalHandler(mock_decimer_client)
        
        with patch('rdkit.Chem.MolFromSmiles') as mock_mol_from_smiles:
            mock_mol = Mock()
            mock_mol_from_smiles.return_value = mock_mol
            
            with patch('rdkit.Chem.Descriptors.MolWt', side_effect=Exception("MW calculation error")):
                with patch('rdkit.Chem.rdMolDescriptors.CalcMolFormula', return_value='C2H6O'):
                    result = handler.process_image(mock_image_path)
                    
                    # Should handle the error gracefully
                    assert result is None
    
    def test_formula_calculation_error(self, mock_image_path):
        """Test handling formula calculation errors"""
        mock_decimer_client = Mock()
        mock_decimer_client.image_to_smiles.return_value = 'CCO'
        
        handler = ChemicalHandler(mock_decimer_client)
        
        with patch('rdkit.Chem.MolFromSmiles') as mock_mol_from_smiles:
            mock_mol = Mock()
            mock_mol_from_smiles.return_value = mock_mol
            
            with patch('rdkit.Chem.Descriptors.MolWt', return_value=46.07):
                with patch('rdkit.Chem.rdMolDescriptors.CalcMolFormula', side_effect=Exception("Formula error")):
                    result = handler.process_image(mock_image_path)
                    
                    # Should handle the error gracefully
                    assert result is None