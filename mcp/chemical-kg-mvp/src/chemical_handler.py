# chemical_handler.py
from rdkit import Chem
from rdkit.Chem import Descriptors
import json

class ChemicalHandler:
    def __init__(self, decimer_client):
        self.decimer = decimer_client
    
    def process_image(self, img_path, context_text=""):
        """Process chemical structure image"""
        smiles = self.decimer.image_to_smiles(img_path)
        
        if smiles:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                return {
                    'smiles': smiles,
                    'molecular_weight': Descriptors.MolWt(mol),
                    'formula': Chem.rdMolDescriptors.CalcMolFormula(mol),
                    'context': context_text,
                    'image_path': img_path
                }
        return None
    
    def get_structure_context(self, page_text, image_bbox):
        """Extract text context around image"""
        # Simple approach: get text near image
        # In production, use more sophisticated methods
        lines = page_text.split('\n')
        # Return nearby lines as context
        return ' '.join(lines[:5])  # Simplified