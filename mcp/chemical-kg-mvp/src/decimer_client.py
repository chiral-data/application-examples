# decimer_client.py
import requests
import base64
from PIL import Image
import io

class DECIMERClient:
    def __init__(self):
        # Use DECIMER web API for MVP (faster than local setup)
        self.api_url = "https://www.decimer.ai/process_image"
    
    def image_to_smiles(self, image_path):
        """Convert chemical structure image to SMILES"""
        try:
            with open(image_path, 'rb') as img_file:
                img_base64 = base64.b64encode(img_file.read()).decode()
            
            response = requests.post(
                self.api_url,
                json={'image': img_base64}
            )
            
            if response.status_code == 200:
                return response.json().get('smiles', None)
            return None
        except Exception as e:
            print(f"Error: {e}")
            return None