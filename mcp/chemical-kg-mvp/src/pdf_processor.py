# pdf_processor.py
import fitz  # PyMuPDF
from PIL import Image
import os

class PDFProcessor:
    def __init__(self, pdf_path):
        self.pdf_path = pdf_path
        self.doc = fitz.open(pdf_path)
    
    def extract_text_and_images(self):
        """Extract text and images with location info"""
        results = []
        
        for page_num, page in enumerate(self.doc):
            # Extract text
            text = page.get_text()
            
            # Extract images
            image_list = page.get_images()
            images = []
            
            for img_index, img in enumerate(image_list):
                xref = img[0]
                pix = fitz.Pixmap(self.doc, xref)
                
                if pix.n - pix.alpha < 4:  # GRAY or RGB
                    img_data = pix.tobytes("png")
                    img = Image.open(io.BytesIO(img_data))
                    
                    # Save temporarily
                    img_path = f"temp_img_{page_num}_{img_index}.png"
                    img.save(img_path)
                    
                    images.append({
                        'path': img_path,
                        'page': page_num,
                        'bbox': page.get_image_bbox(img[7])
                    })
                
                pix = None
            
            results.append({
                'page': page_num,
                'text': text,
                'images': images
            })
        
        return results