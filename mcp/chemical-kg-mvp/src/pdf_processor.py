# pdf_processor.py
import fitz  # PyMuPDF
from PIL import Image
import os
import io
from typing import List, Dict, Any, Optional

try:
    from decimer_segmentation import segment_chemical_structures_from_file, segment_chemical_structures
    import cv2
    import pdf2image
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
        """Extract chemical structures using DECIMER segmentation with image-based approach"""
        if not DECIMER_SEGMENTATION_AVAILABLE:
            print("DECIMER segmentation not available. Using standard image extraction.")
            # Return images from standard extraction for processing
            pages_data = self.extract_text_and_images()
            images = []
            for page_data in pages_data:
                images.extend(page_data.get('images', []))
            return images
        
        try:
            print(f"Using DECIMER image-based segmentation for {self.pdf_path}")
            
            # Step 1: Convert PDF to high-resolution images
            try:
                import pdf2image
                print("Converting PDF to images...")
                
                # Create temp directory for page images
                pages_dir = os.path.join(self.temp_dir, "pdf_pages")
                os.makedirs(pages_dir, exist_ok=True)
                
                # Convert PDF to images
                page_images = pdf2image.convert_from_path(self.pdf_path, dpi=300)
                print(f"Converted PDF to {len(page_images)} page images")
                
                # Save page images
                page_image_paths = []
                for i, page_image in enumerate(page_images):
                    page_path = os.path.join(pages_dir, f"page_{i+1}.png")
                    page_image.save(page_path)
                    page_image_paths.append(page_path)
                    print(f"Saved page {i+1}: {page_path}")
                
            except Exception as e:
                print(f"PDF to image conversion failed: {e}")
                raise
            
            # Step 2: Segment structures from each page image
            all_structures = []
            structure_count = 0
            
            for page_idx, page_image_path in enumerate(page_image_paths):
                try:
                    print(f"Segmenting structures from page {page_idx + 1}...")
                    
                    # Read page image
                    page_image = cv2.imread(page_image_path)
                    if page_image is None:
                        print(f"Failed to read page image: {page_image_path}")
                        continue
                    
                    # Segment structures from this page
                    from decimer_segmentation import segment_chemical_structures
                    segments = segment_chemical_structures(page_image, expand=expand)
                    
                    print(f"Found {len(segments)} structures on page {page_idx + 1}")
                    
                    # Save each segment
                    for seg_idx, segment in enumerate(segments):
                        try:
                            # Save segment image
                            img_filename = f"structure_page_{page_idx+1}_{seg_idx+1}.png"
                            img_path = os.path.join(self.temp_dir, img_filename)
                            cv2.imwrite(img_path, segment)
                            
                            # Get dimensions
                            height, width = segment.shape[:2]
                            
                            all_structures.append({
                                'path': img_path,
                                'index': structure_count,
                                'page': page_idx + 1,
                                'segment_on_page': seg_idx + 1,
                                'width': width,
                                'height': height,
                                'source': 'decimer_image_segmentation'
                            })
                            
                            structure_count += 1
                            
                        except Exception as e:
                            print(f"Error saving segment {seg_idx} from page {page_idx + 1}: {e}")
                    
                except Exception as e:
                    print(f"Error processing page {page_idx + 1}: {e}")
                    continue
            
            print(f"DECIMER image-based segmentation completed: {len(all_structures)} structures found")
            return all_structures
            
        except Exception as e:
            print(f"DECIMER image-based segmentation failed: {e}. Falling back to standard extraction.")
            # Return images from standard extraction for processing
            pages_data = self.extract_text_and_images()
            images = []
            for page_data in pages_data:
                images.extend(page_data.get('images', []))
            return images
    
    def __del__(self):
        """Cleanup when object is destroyed"""
        try:
            if hasattr(self, 'doc') and self.doc:
                self.doc.close()
        except:
            pass