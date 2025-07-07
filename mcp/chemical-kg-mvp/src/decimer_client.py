# decimer_client.py
import os
from typing import Optional, Dict, Any, List
from PIL import Image
import time

try:
    import tensorflow as tf
    # Force GPU usage if available
    gpus = tf.config.list_physical_devices('GPU')
    if gpus:
        try:
            # Enable memory growth for GPU
            for gpu in gpus:
                tf.config.experimental.set_memory_growth(gpu, True)
            print(f"DECIMER: GPU acceleration enabled with {len(gpus)} GPU(s)")
        except RuntimeError as e:
            print(f"DECIMER: GPU setup error: {e}")
    
    from DECIMER import predict_SMILES
    DECIMER_AVAILABLE = True
except ImportError:
    DECIMER_AVAILABLE = False
    print("Warning: DECIMER package not installed. Install with: pip install DECIMER")

try:
    from decimer_segmentation import segment_chemical_structures_from_file, segment_chemical_structures
    import cv2
    import pdf2image
    SEGMENTATION_AVAILABLE = True
except ImportError:
    SEGMENTATION_AVAILABLE = False
    print("Warning: decimer-segmentation package not installed. Install with: pip install decimer-segmentation")

class DECIMERClient:
    def __init__(self):
        """Initialize DECIMER client for chemical structure recognition"""
        self.decimer_available = DECIMER_AVAILABLE
        self.segmentation_available = SEGMENTATION_AVAILABLE
        
        # Cache for processed images
        self.cache = {}
        self.cache_enabled = True
        
        # Statistics tracking
        self.stats = {
            'processed_images': 0,
            'cache_hits': 0,
            'errors': 0,
            'segmented_structures': 0
        }
        
        if not self.decimer_available:
            print("DECIMER not available. Chemical structure recognition will be skipped.")
    
    def convert_pdf_to_images(self, pdf_path: str, output_dir: str = "temp_pdf_pages", dpi: int = 300) -> List[str]:
        """Convert PDF to high-resolution images"""
        try:
            print(f"Converting PDF to images: {pdf_path}")
            
            # Create output directory
            os.makedirs(output_dir, exist_ok=True)
            
            # Convert PDF to images
            images = pdf2image.convert_from_path(pdf_path, dpi=dpi)
            
            image_paths = []
            for i, image in enumerate(images):
                image_path = os.path.join(output_dir, f"page_{i+1}.png")
                image.save(image_path)
                image_paths.append(image_path)
                print(f"Saved page {i+1}: {image_path} ({image.size})")
            
            print(f"Converted PDF to {len(image_paths)} images")
            return image_paths
            
        except Exception as e:
            print(f"Error converting PDF to images: {e}")
            self.stats['errors'] += 1
            return []
    
    def segment_structures_from_pdf(self, pdf_path: str, expand: bool = True, output_dir: str = "temp_pdf_pages") -> List[Any]:
        """Extract chemical structure images from PDF using image-based approach"""
        if not self.segmentation_available:
            print("DECIMER segmentation not available")
            return []
        
        try:
            print(f"Processing PDF with image-based approach: {pdf_path}")
            
            # Step 1: Convert PDF to images
            image_paths = self.convert_pdf_to_images(pdf_path, output_dir)
            
            if not image_paths:
                print("Failed to convert PDF to images")
                return []
            
            # Step 2: Segment structures from each image
            all_segments = []
            for image_path in image_paths:
                print(f"Segmenting structures from: {image_path}")
                segments = self.segment_structures_from_image(image_path, expand=expand)
                all_segments.extend(segments)
                print(f"Found {len(segments)} structures in {image_path}")
            
            self.stats['segmented_structures'] += len(all_segments)
            print(f"Total structures found: {len(all_segments)}")
            return all_segments
            
        except Exception as e:
            print(f"Error processing PDF {pdf_path}: {e}")
            self.stats['errors'] += 1
            return []
    
    def segment_structures_from_image(self, image_path: str, expand: bool = True) -> List[Any]:
        """Extract chemical structures from a single image"""
        if not self.segmentation_available:
            print("DECIMER segmentation not available")
            return []
        
        try:
            # Read image using cv2
            image = cv2.imread(image_path)
            if image is None:
                print(f"Failed to read image: {image_path}")
                return []
            
            segments = segment_chemical_structures(image, expand=expand)
            self.stats['segmented_structures'] += len(segments)
            return segments
        except Exception as e:
            print(f"Error segmenting image {image_path}: {e}")
            self.stats['errors'] += 1
            return []
    
    def image_to_smiles(self, image_path: str) -> Optional[str]:
        """Convert chemical structure image to SMILES using DECIMER"""
        if not self.decimer_available:
            print("DECIMER not available for SMILES prediction")
            return None
        
        # Check cache first
        if self.cache_enabled:
            cache_key = self._get_cache_key(image_path)
            if cache_key in self.cache:
                self.stats['cache_hits'] += 1
                return self.cache[cache_key]
        
        try:
            # Validate image exists
            if not os.path.exists(image_path):
                print(f"Image file not found: {image_path}")
                return None
            
            # Predict SMILES using DECIMER
            print(f"Processing image: {image_path}")
            start_time = time.time()
            
            smiles = predict_SMILES(image_path)
            
            processing_time = time.time() - start_time
            print(f"DECIMER prediction completed in {processing_time:.2f} seconds")
            
            if smiles and self._validate_smiles(smiles):
                self.stats['processed_images'] += 1
                
                # Cache the result
                if self.cache_enabled:
                    self.cache[cache_key] = smiles
                
                return smiles
            else:
                print(f"Invalid SMILES returned: {smiles}")
                return None
                
        except Exception as e:
            print(f"Error processing image {image_path}: {e}")
            self.stats['errors'] += 1
            return None
    
    def batch_process(self, image_paths: List[str]) -> Dict[str, Optional[str]]:
        """Process multiple images and return SMILES for each"""
        results = {}
        
        for i, image_path in enumerate(image_paths):
            print(f"Processing image {i+1}/{len(image_paths)}")
            try:
                smiles = self.image_to_smiles(image_path)
                results[image_path] = smiles
            except Exception as e:
                print(f"Failed to process {image_path}: {e}")
                results[image_path] = None
                self.stats['errors'] += 1
        
        return results
    
    def _validate_smiles(self, smiles: str) -> bool:
        """Validate SMILES string using RDKit, allowing R-group placeholders"""
        if not smiles or len(smiles) < 2:
            return False
        
        # Try RDKit validation if available
        try:
            from rdkit import Chem
            
            # Check if SMILES contains R-group placeholders (common in medicinal chemistry)
            r_group_patterns = ['[R1]', '[R2]', '[R3]', '[R4]', '[R5]', '[R]']
            has_r_groups = any(pattern in smiles for pattern in r_group_patterns)
            
            if has_r_groups:
                # For R-group containing SMILES, do basic syntax check instead of RDKit validation
                # RDKit doesn't recognize R-group placeholders as valid atoms
                print(f"SMILES contains R-group placeholders - accepting as valid: {smiles}")
                return True
            
            # For regular SMILES, use RDKit validation
            mol = Chem.MolFromSmiles(smiles)
            return mol is not None
        except ImportError:
            # RDKit not available, trust DECIMER output
            return True
        except Exception as e:
            print(f"SMILES validation error: {e}")
            return False
    
    def _get_cache_key(self, image_path: str) -> str:
        """Generate cache key for image"""
        try:
            # Use file modification time and size as key
            stat = os.stat(image_path)
            return f"{image_path}_{stat.st_mtime}_{stat.st_size}"
        except:
            return image_path
    
    def clear_cache(self):
        """Clear the response cache"""
        self.cache.clear()
        print("Cache cleared")
    
    def get_stats(self) -> Dict[str, Any]:
        """Get usage statistics"""
        return self.stats.copy()
    
    def process_pdf_complete(self, pdf_path: str, output_dir: str = "extracted_structures") -> Dict[str, Any]:
        """Complete workflow: segment structures from PDF and convert to SMILES using image-based approach"""
        if not (self.decimer_available and self.segmentation_available):
            print("Both DECIMER and segmentation packages required for complete PDF processing")
            return {"error": "Missing required packages"}
        
        start_time = time.time()
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Segment structures from PDF using image-based approach
        segments = self.segment_structures_from_pdf(pdf_path, output_dir=f"{output_dir}_temp_pages")
        
        results = {
            "pdf_path": pdf_path,
            "total_structures": len(segments),
            "processing_time": 0,
            "structures": []
        }
        
        # Process each segmented structure
        for i, segment in enumerate(segments):
            try:
                # Save segment as image
                structure_path = os.path.join(output_dir, f"structure_{i+1}.png")
                cv2.imwrite(structure_path, segment)
                
                print(f"Processing structure {i+1}/{len(segments)}: {structure_path}")
                
                # Convert to SMILES
                smiles = self.image_to_smiles(structure_path)
                
                structure_info = {
                    "index": i + 1,
                    "image_path": structure_path,
                    "smiles": smiles,
                    "valid": smiles is not None,
                    "image_size": segment.shape if hasattr(segment, 'shape') else None
                }
                
                results["structures"].append(structure_info)
                
                if smiles:
                    print(f"  Success: {smiles}")
                else:
                    print(f"  Failed: No valid SMILES generated")
                
            except Exception as e:
                print(f"Error processing segment {i+1}: {e}")
                results["structures"].append({
                    "index": i + 1,
                    "error": str(e)
                })
        
        results["processing_time"] = time.time() - start_time
        
        # Summary
        valid_structures = sum(1 for s in results["structures"] if s.get("valid", False))
        print(f"\nProcessing complete: {valid_structures}/{len(segments)} structures successfully converted to SMILES")
        print(f"Total processing time: {results['processing_time']:.2f} seconds")
        
        return results