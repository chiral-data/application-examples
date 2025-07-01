# pdf_processor.py
import fitz  # PyMuPDF
from PIL import Image
import os
import io
from typing import List, Dict, Any, Optional

try:
    from decimer_segmentation import segment_chemical_structures_from_file
    import cv2
    DECIMER_SEGMENTATION_AVAILABLE = True
except ImportError:
    DECIMER_SEGMENTATION_AVAILABLE = False

class PDFProcessor:
    def __init__(self, pdf_path):
        self.pdf_path = pdf_path
        try:
            self.doc = fitz.open(pdf_path)
        except Exception as e:
            raise ValueError(f"Failed to open PDF file: {e}")
        
        if self.doc.is_encrypted:
            raise ValueError("PDF is encrypted and cannot be processed")
        
        self.temp_dir = "temp"
        os.makedirs(self.temp_dir, exist_ok=True)
    
    def extract_text_and_images(self):
        """Extract text and images with location info and better error handling"""
        results = []
        
        try:
            for page_num, page in enumerate(self.doc):
                # Extract text with better handling
                try:
                    text = page.get_text()
                    if not text.strip():
                        # Try OCR or alternative text extraction if needed
                        text = self._extract_text_alternative(page)
                except Exception as e:
                    print(f"Warning: Failed to extract text from page {page_num}: {e}")
                    text = ""
                
                # Extract images with improved handling
                images = self._extract_images_from_page(page, page_num)
                
                results.append({
                    'page': page_num,
                    'text': text,
                    'images': images,
                    'page_size': {'width': page.rect.width, 'height': page.rect.height}
                })
                
        except Exception as e:
            raise RuntimeError(f"Failed to process PDF: {e}")
        
        return results
    
    def _extract_text_alternative(self, page):
        """Alternative text extraction method"""
        try:
            # Try different text extraction methods
            text_dict = page.get_text("dict")
            text = ""
            for block in text_dict.get("blocks", []):
                if "lines" in block:
                    for line in block["lines"]:
                        for span in line.get("spans", []):
                            text += span.get("text", "") + " "
                        text += "\n"
            return text
        except:
            return ""
    
    def _extract_images_from_page(self, page, page_num):
        """Extract images from a page with better error handling"""
        images = []
        
        try:
            image_list = page.get_images(full=True)
            
            for img_index, img in enumerate(image_list):
                try:
                    # Get image reference and metadata
                    xref = img[0]
                    
                    # Check if image is large enough to be meaningful
                    if len(img) > 2 and img[2] < 50 and img[3] < 50:
                        continue  # Skip very small images
                    
                    # Extract image
                    pix = fitz.Pixmap(self.doc, xref)
                    
                    # Only process RGB, RGBA, or grayscale images
                    if pix.n - pix.alpha < 4:
                        img_data = pix.tobytes("png")
                        
                        # Validate image data
                        if len(img_data) < 1000:  # Skip very small images
                            pix = None
                            continue
                        
                        try:
                            img_pil = Image.open(io.BytesIO(img_data))
                            
                            # Check image dimensions
                            if img_pil.width < 50 or img_pil.height < 50:
                                pix = None
                                continue
                            
                            # Save with unique filename
                            img_filename = f"page_{page_num}_img_{img_index}_{xref}.png"
                            img_path = os.path.join(self.temp_dir, img_filename)
                            img_pil.save(img_path, "PNG")
                            
                            # Get image position on page
                            try:
                                bbox = page.get_image_bbox(img)
                            except:
                                bbox = fitz.Rect(0, 0, img_pil.width, img_pil.height)
                            
                            images.append({
                                'path': img_path,
                                'page': page_num,
                                'bbox': bbox,
                                'width': img_pil.width,
                                'height': img_pil.height,
                                'size_bytes': len(img_data),
                                'format': 'PNG'
                            })
                            
                        except Exception as e:
                            print(f"Warning: Failed to process image {img_index} on page {page_num}: {e}")
                    
                    # Clean up pixmap
                    pix = None
                    
                except Exception as e:
                    print(f"Warning: Failed to extract image {img_index} from page {page_num}: {e}")
                    continue
                    
        except Exception as e:
            print(f"Warning: Failed to extract images from page {page_num}: {e}")
        
        return images
    
    def cleanup_temp_files(self):
        """Clean up temporary image files"""
        try:
            if os.path.exists(self.temp_dir):
                for filename in os.listdir(self.temp_dir):
                    if filename.startswith(('page_', 'temp_img_')):
                        file_path = os.path.join(self.temp_dir, filename)
                        try:
                            os.remove(file_path)
                        except Exception as e:
                            print(f"Warning: Failed to remove temp file {file_path}: {e}")
        except Exception as e:
            print(f"Warning: Failed to cleanup temp files: {e}")
    
    def get_document_info(self):
        """Get document metadata"""
        try:
            metadata = self.doc.metadata
            return {
                'title': metadata.get('title', ''),
                'author': metadata.get('author', ''),
                'subject': metadata.get('subject', ''),
                'creator': metadata.get('creator', ''),
                'producer': metadata.get('producer', ''),
                'creation_date': metadata.get('creationDate', ''),
                'modification_date': metadata.get('modDate', ''),
                'page_count': self.doc.page_count,
                'is_pdf': self.doc.is_pdf,
                'needs_pass': self.doc.needs_pass,
                'is_encrypted': self.doc.is_encrypted
            }
        except Exception as e:
            print(f"Warning: Failed to get document info: {e}")
            return {'page_count': self.doc.page_count if self.doc else 0}
    
    def extract_chemical_structures_with_decimer(self, expand: bool = True) -> List[Dict[str, Any]]:
        """Extract chemical structures using DECIMER segmentation"""
        if not DECIMER_SEGMENTATION_AVAILABLE:
            print("DECIMER segmentation not available. Using standard image extraction.")
            return self.extract_text_and_images()
        
        try:
            print(f"Using DECIMER to segment chemical structures from {self.pdf_path}")
            
            # Get chemical structure segments
            segments = segment_chemical_structures_from_file(self.pdf_path, expand=expand)
            print(f"DECIMER found {len(segments)} chemical structures")
            
            # Save segments as images
            structures = []
            for i, segment in enumerate(segments):
                try:
                    # Save segment image
                    img_filename = f"decimer_structure_{i}.png"
                    img_path = os.path.join(self.temp_dir, img_filename)
                    cv2.imwrite(img_path, segment)
                    
                    # Get dimensions
                    height, width = segment.shape[:2]
                    
                    structures.append({
                        'path': img_path,
                        'index': i,
                        'width': width,
                        'height': height,
                        'source': 'decimer_segmentation'
                    })
                    
                except Exception as e:
                    print(f"Error saving segment {i}: {e}")
            
            return structures
            
        except Exception as e:
            print(f"DECIMER segmentation failed: {e}. Falling back to standard extraction.")
            return self.extract_text_and_images()
    
    def __del__(self):
        """Cleanup when object is destroyed"""
        try:
            if hasattr(self, 'doc') and self.doc:
                self.doc.close()
        except:
            pass