import pytest
import os
import json
import time
from unittest.mock import Mock, patch, MagicMock
import requests
from PIL import Image

from decimer_client import DECIMERClient


class TestDECIMERClient:
    """Test suite for DECIMERClient class"""
    
    def test_init_default(self):
        """Test default initialization"""
        client = DECIMERClient()
        
        assert client.use_local is False
        assert client.timeout == 30
        assert client.max_retries == 3
        assert 'decimer.ai' in client.api_url
        assert client.cache_enabled is True
        assert client.stats['api_calls'] == 0
    
    def test_init_custom_params(self):
        """Test initialization with custom parameters"""
        client = DECIMERClient(
            use_local=True,
            api_url="http://custom-api.com",
            timeout=60,
            max_retries=5
        )
        
        assert client.use_local is True
        assert client.timeout == 60
        assert client.max_retries == 5
        assert client.api_url == "http://custom-api.com"
    
    def test_init_environment_variables(self):
        """Test initialization with environment variables"""
        with patch.dict(os.environ, {
            'DECIMER_API_URL': 'http://env-api.com',
            'DECIMER_TIMEOUT': '45',
            'DECIMER_MAX_RETRIES': '4'
        }):
            client = DECIMERClient()
            assert client.api_url == "http://env-api.com"
    
    @patch('requests.post')
    def test_api_image_to_smiles_success(self, mock_post, mock_image_path):
        """Test successful API call for SMILES conversion"""
        # Mock successful response
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {'smiles': 'CCO'}
        mock_post.return_value = mock_response
        
        client = DECIMERClient()
        
        with patch.object(client, '_prepare_image', return_value='base64_image_data'):
            with patch.object(client, '_validate_smiles', return_value=True):
                smiles = client._api_image_to_smiles(mock_image_path)
                
                assert smiles == 'CCO'
                assert client.stats['api_calls'] == 1
                mock_post.assert_called_once()
    
    @patch('requests.post')
    def test_api_image_to_smiles_retry_on_rate_limit(self, mock_post, mock_image_path):
        """Test retry logic on rate limiting"""
        # First call returns 429 (rate limited), second succeeds
        mock_response_429 = Mock()
        mock_response_429.status_code = 429
        
        mock_response_200 = Mock()
        mock_response_200.status_code = 200
        mock_response_200.json.return_value = {'smiles': 'CCO'}
        
        mock_post.side_effect = [mock_response_429, mock_response_200]
        
        client = DECIMERClient(max_retries=2)
        
        with patch.object(client, '_prepare_image', return_value='base64_image_data'):
            with patch.object(client, '_validate_smiles', return_value=True):
                with patch('time.sleep'):  # Mock sleep to speed up test
                    smiles = client._api_image_to_smiles(mock_image_path)
                    
                    assert smiles == 'CCO'
                    assert mock_post.call_count == 2
    
    @patch('requests.post')
    def test_api_image_to_smiles_timeout(self, mock_post, mock_image_path):
        """Test timeout handling"""
        mock_post.side_effect = requests.exceptions.Timeout()
        
        client = DECIMERClient(max_retries=2)
        
        with patch.object(client, '_prepare_image', return_value='base64_image_data'):
            with patch('time.sleep'):
                smiles = client._api_image_to_smiles(mock_image_path)
                
                assert smiles is None
                assert client.stats['errors'] == 1
                assert mock_post.call_count == 2
    
    def test_prepare_image_success(self, mock_image_path):
        """Test successful image preparation"""
        client = DECIMERClient()
        
        # Mock os.path.exists and os.path.getsize
        with patch('os.path.exists', return_value=True):
            with patch('os.path.getsize', return_value=1024):  # 1KB file
                with patch('PIL.Image.open') as mock_image_open:
                    mock_img = Mock()
                    mock_img.mode = 'RGB'
                    mock_img.width = 100
                    mock_img.height = 100
                    mock_img.save = Mock()
                    mock_image_open.return_value.__enter__ = Mock(return_value=mock_img)
                    mock_image_open.return_value.__exit__ = Mock(return_value=None)
                    
                    with patch('base64.b64encode', return_value=b'encoded_data'):
                        result = client._prepare_image(mock_image_path)
                        
                        assert result == 'encoded_data'
    
    def test_prepare_image_file_not_found(self):
        """Test image preparation with missing file"""
        client = DECIMERClient()
        
        with patch('os.path.exists', return_value=False):
            result = client._prepare_image("nonexistent.png")
            assert result is None
    
    def test_prepare_image_too_large(self, mock_image_path):
        """Test image preparation with oversized file"""
        client = DECIMERClient()
        
        with patch('os.path.exists', return_value=True):
            # Mock file size larger than limit (default 10MB)
            with patch('os.path.getsize', return_value=20 * 1024 * 1024):
                result = client._prepare_image(mock_image_path)
                assert result is None
    
    def test_prepare_image_resize_large_image(self, mock_image_path):
        """Test that large images are resized"""
        client = DECIMERClient()
        
        with patch('os.path.exists', return_value=True):
            with patch('os.path.getsize', return_value=1024):
                with patch('PIL.Image.open') as mock_image_open:
                    mock_img = Mock()
                    mock_img.mode = 'RGB'
                    mock_img.width = 3000  # Large width
                    mock_img.height = 2000  # Large height
                    mock_img.thumbnail = Mock()
                    mock_img.save = Mock()
                    mock_image_open.return_value.__enter__ = Mock(return_value=mock_img)
                    mock_image_open.return_value.__exit__ = Mock(return_value=None)
                    
                    with patch('base64.b64encode', return_value=b'encoded_data'):
                        result = client._prepare_image(mock_image_path)
                        
                        mock_img.thumbnail.assert_called_once()
                        assert result == 'encoded_data'
    
    def test_validate_smiles_valid(self):
        """Test SMILES validation with valid SMILES"""
        client = DECIMERClient()
        
        # Test with RDKit available
        with patch('rdkit.Chem.MolFromSmiles') as mock_mol:
            mock_mol.return_value = Mock()  # Non-None molecule
            
            assert client._validate_smiles('CCO') is True
            assert client._validate_smiles('C6H6') is True
    
    def test_validate_smiles_invalid(self):
        """Test SMILES validation with invalid SMILES"""
        client = DECIMERClient()
        
        assert client._validate_smiles('') is False
        assert client._validate_smiles('X') is False
        assert client._validate_smiles('invalid_characters_!@#') is False
        
        # Test with RDKit returning None (invalid molecule)
        with patch('rdkit.Chem.MolFromSmiles', return_value=None):
            assert client._validate_smiles('INVALID') is False
    
    def test_validate_smiles_no_rdkit(self):
        """Test SMILES validation without RDKit"""
        client = DECIMERClient()
        
        with patch('rdkit.Chem.MolFromSmiles', side_effect=ImportError):
            # Should fall back to basic validation
            assert client._validate_smiles('CCO') is True
            assert client._validate_smiles('invalid_!@#') is False
    
    def test_caching_functionality(self, mock_image_path):
        """Test response caching"""
        client = DECIMERClient()
        
        with patch.object(client, '_api_image_to_smiles', return_value='CCO') as mock_api:
            with patch.object(client, '_get_cache_key', return_value='test_key'):
                # First call should hit API
                result1 = client.image_to_smiles(mock_image_path)
                assert result1 == 'CCO'
                assert mock_api.call_count == 1
                
                # Second call should use cache
                result2 = client.image_to_smiles(mock_image_path)
                assert result2 == 'CCO'
                assert mock_api.call_count == 1  # Still only one API call
                assert client.stats['cache_hits'] == 1
    
    def test_local_decimer_setup_success(self):
        """Test successful local DECIMER setup"""
        with patch('DECIMER', create=True) as mock_decimer:
            client = DECIMERClient(use_local=True)
            assert client.local_decimer == mock_decimer
    
    def test_local_decimer_setup_failure(self):
        """Test local DECIMER setup failure"""
        with patch('builtins.__import__', side_effect=ImportError):
            client = DECIMERClient(use_local=True)
            assert client.local_decimer is None
            assert client.use_local is False
    
    @patch('builtins.print')
    def test_local_image_to_smiles(self, mock_print, mock_image_path):
        """Test local DECIMER processing"""
        mock_decimer = Mock()
        mock_decimer.predict_SMILES.return_value = 'CCO'
        
        client = DECIMERClient()
        client.local_decimer = mock_decimer
        
        with patch('PIL.Image.open') as mock_image_open:
            mock_img = Mock()
            mock_img.mode = 'RGB'
            mock_image_open.return_value = mock_img
            
            with patch.object(client, '_validate_smiles', return_value=True):
                result = client._local_image_to_smiles(mock_image_path)
                
                assert result == 'CCO'
                mock_decimer.predict_SMILES.assert_called_once_with(mock_img)
    
    def test_get_cache_key(self, mock_image_path):
        """Test cache key generation"""
        client = DECIMERClient()
        
        with patch('os.stat') as mock_stat:
            mock_stat.return_value = Mock(st_mtime=1234567890, st_size=1024)
            
            key = client._get_cache_key(mock_image_path)
            expected = f"{mock_image_path}_1234567890_1024"
            assert key == expected
    
    def test_clear_cache(self):
        """Test cache clearing"""
        client = DECIMERClient()
        client.cache['test_key'] = 'test_value'
        
        client.clear_cache()
        assert len(client.cache) == 0
    
    def test_get_stats(self):
        """Test statistics retrieval"""
        client = DECIMERClient()
        client.stats['api_calls'] = 5
        client.stats['cache_hits'] = 2
        
        stats = client.get_stats()
        assert stats['api_calls'] == 5
        assert stats['cache_hits'] == 2
        assert isinstance(stats, dict)
    
    def test_batch_process(self, temp_dir):
        """Test batch processing multiple images"""
        # Create test image files
        image_paths = []
        for i in range(3):
            img_path = os.path.join(temp_dir, f"test_{i}.png")
            img = Image.new('RGB', (100, 100), color='white')
            img.save(img_path)
            image_paths.append(img_path)
        
        client = DECIMERClient()
        
        with patch.object(client, 'image_to_smiles') as mock_image_to_smiles:
            mock_image_to_smiles.side_effect = ['CCO', 'C6H6', None]
            
            results = client.batch_process(image_paths)
            
            assert len(results) == 3
            assert results[image_paths[0]] == 'CCO'
            assert results[image_paths[1]] == 'C6H6'
            assert results[image_paths[2]] is None
    
    def test_image_to_smiles_with_local_fallback(self, mock_image_path):
        """Test image_to_smiles with local DECIMER and API fallback"""
        client = DECIMERClient(use_local=True)
        client.local_decimer = Mock()
        
        # Local DECIMER fails, should fall back to API
        with patch.object(client, '_local_image_to_smiles', return_value=None):
            with patch.object(client, '_api_image_to_smiles', return_value='CCO'):
                with patch.object(client, '_get_cache_key', return_value='test_key'):
                    result = client.image_to_smiles(mock_image_path)
                    
                    assert result == 'CCO'
    
    @patch('builtins.print')
    def test_error_handling_during_processing(self, mock_print, mock_image_path):
        """Test error handling during image processing"""
        client = DECIMERClient()
        
        with patch.object(client, '_prepare_image', side_effect=Exception("Processing error")):
            result = client._api_image_to_smiles(mock_image_path)
            assert result is None
            assert client.stats['errors'] == 1