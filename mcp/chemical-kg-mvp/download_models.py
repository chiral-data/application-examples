#!/usr/bin/env python3
"""Download DECIMER models during Docker build to cache them in the image."""

# Fix for SQLite version compatibility with ChromaDB
try:
    __import__('pysqlite3')
    import sys
    sys.modules['sqlite3'] = sys.modules.pop('pysqlite3')
except ImportError:
    pass

import os
import sys

# Set environment variable to download models to a specific location
os.environ['DATA_DIR'] = '/root/.data'

try:
    # Import DECIMER to trigger model download
    import DECIMER
    from decimer_segmentation import segment_chemical_structures_from_file
    
    print("DECIMER models downloaded successfully!")
    print(f"Models stored in: {os.environ['DATA_DIR']}")
    
    # List downloaded files to verify
    model_dir = os.path.join(os.environ['DATA_DIR'], 'DECIMER-V2')
    if os.path.exists(model_dir):
        print(f"Contents of {model_dir}:")
        for item in os.listdir(model_dir):
            print(f"  - {item}")
    
except Exception as e:
    print(f"Error downloading DECIMER models: {e}")
    sys.exit(1)