import os

# Test configuration settings
TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'data')
TEST_PDF_PATH = os.path.join(TEST_DATA_DIR, 'test.pdf')

# Sample chemical structures for testing
SAMPLE_STRUCTURES = [
    {
        'smiles': 'CCO',
        'formula': 'C2H6O',
        'molecular_weight': 46.07,
        'context': 'ethanol molecule found in figure 1',
        'image_path': '/test/path/ethanol.png',
        'page': 0
    },
    {
        'smiles': 'C6H6',
        'formula': 'C6H6',
        'molecular_weight': 78.11,
        'context': 'benzene ring structure shown in scheme 2',
        'image_path': '/test/path/benzene.png',
        'page': 1
    },
    {
        'smiles': 'CC(=O)O',
        'formula': 'C2H4O2',
        'molecular_weight': 60.05,
        'context': 'acetic acid depicted in figure 3',
        'image_path': '/test/path/acetic_acid.png',
        'page': 1
    }
]

# Sample text content for testing
SAMPLE_CHEMICAL_TEXT = """
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