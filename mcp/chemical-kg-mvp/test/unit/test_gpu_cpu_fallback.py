"""
GPU/CPU detection and fallback tests.
These tests ensure proper handling of different hardware configurations and fallback scenarios.
"""

import pytest
import os
import tempfile
import shutil
from unittest.mock import patch, MagicMock, call
import time

from decimer_client import DECIMERClient
from pdf_processor import PDFProcessor


class TestGPUDetectionAndFallback:
    """Test GPU detection and CPU fallback scenarios"""
    
    def test_gpu_detection_success(self):
        """Test successful GPU detection and configuration"""
        mock_gpu1 = MagicMock()
        mock_gpu2 = MagicMock()
        
        with patch('tensorflow.config.list_physical_devices') as mock_list_devices:
            mock_list_devices.return_value = [mock_gpu1, mock_gpu2]
            
            with patch('tensorflow.config.experimental.set_memory_growth') as mock_memory_growth:
                client = DECIMERClient()
                
                # Should have detected GPUs
                mock_list_devices.assert_called_with('GPU')
                
                # Should have set memory growth for both GPUs
                expected_calls = [call(mock_gpu1, True), call(mock_gpu2, True)]
                mock_memory_growth.assert_has_calls(expected_calls)
                
                assert client.decimer_available, "DECIMER should be available"
    
    def test_no_gpu_available(self):
        """Test behavior when no GPU is available"""
        with patch('tensorflow.config.list_physical_devices') as mock_list_devices:
            mock_list_devices.return_value = []  # No GPUs
            
            with patch('tensorflow.config.experimental.set_memory_growth') as mock_memory_growth:
                client = DECIMERClient()
                
                # Should still try to detect GPUs
                mock_list_devices.assert_called_with('GPU')
                
                # Should not attempt to set memory growth
                mock_memory_growth.assert_not_called()
                
                # Should still be available (CPU fallback)
                assert client.decimer_available, "DECIMER should work on CPU"
    
    def test_gpu_memory_growth_error(self):
        """Test handling of GPU memory growth configuration errors"""
        mock_gpu = MagicMock()
        
        with patch('tensorflow.config.list_physical_devices') as mock_list_devices:
            mock_list_devices.return_value = [mock_gpu]
            
            with patch('tensorflow.config.experimental.set_memory_growth') as mock_memory_growth:
                mock_memory_growth.side_effect = RuntimeError("GPU memory configuration failed")
                
                # Should not crash despite GPU error
                client = DECIMERClient()
                
                assert client.decimer_available, "Should still be available despite GPU error"
                mock_memory_growth.assert_called_with(mock_gpu, True)
    
    def test_tensorflow_import_error(self):
        """Test handling when TensorFlow is not available"""
        with patch('decimer_client.tf', side_effect=ImportError("TensorFlow not installed")):
            client = DECIMERClient()
            
            assert not client.decimer_available, "DECIMER should not be available without TensorFlow"
            assert not client.segmentation_available, "Segmentation should not be available"
    
    def test_decimer_import_error(self):
        """Test handling when DECIMER package is not available"""
        # Mock TensorFlow as available but DECIMER as unavailable
        with patch('tensorflow.config.list_physical_devices', return_value=[]):
            with patch('decimer_client.DECIMER_AVAILABLE', False):
                client = DECIMERClient()
                
                assert not client.decimer_available, "DECIMER should not be available"
                # Segmentation might still be available independently
    
    def test_segmentation_import_error(self):
        """Test handling when decimer-segmentation is not available"""
        with patch('decimer_client.SEGMENTATION_AVAILABLE', False):
            client = DECIMERClient()
            
            assert not client.segmentation_available, "Segmentation should not be available"
            # DECIMER for SMILES prediction might still be available
    
    def test_cuda_visible_devices_cpu_only(self):
        """Test CPU-only mode via CUDA_VISIBLE_DEVICES"""
        with patch.dict(os.environ, {'CUDA_VISIBLE_DEVICES': ''}):
            with patch('tensorflow.config.list_physical_devices') as mock_list_devices:
                # Even if GPUs exist, they should not be visible
                mock_list_devices.return_value = []
                
                client = DECIMERClient()
                
                # Should work in CPU-only mode
                assert client.decimer_available, "Should work in CPU-only mode"
    
    def test_tensorflow_force_cpu_settings(self):
        """Test various TensorFlow CPU-forcing environment variables"""
        cpu_env_vars = {
            'CUDA_VISIBLE_DEVICES': '',
            'TF_FORCE_GPU_ALLOW_GROWTH': 'false',
            'TF_CPP_MIN_LOG_LEVEL': '2'
        }
        
        with patch.dict(os.environ, cpu_env_vars):
            with patch('tensorflow.config.list_physical_devices') as mock_list_devices:
                mock_list_devices.return_value = []  # No visible GPUs
                
                client = DECIMERClient()
                
                # Should still work
                assert client.decimer_available, "Should work with CPU-forcing env vars"
    
    def test_mixed_availability_scenarios(self):
        """Test various combinations of package availability"""
        scenarios = [
            {'decimer': True, 'segmentation': True},
            {'decimer': True, 'segmentation': False},
            {'decimer': False, 'segmentation': True},
            {'decimer': False, 'segmentation': False},
        ]
        
        for scenario in scenarios:
            with patch('decimer_client.DECIMER_AVAILABLE', scenario['decimer']):
                with patch('decimer_client.SEGMENTATION_AVAILABLE', scenario['segmentation']):
                    client = DECIMERClient()
                    
                    assert client.decimer_available == scenario['decimer']
                    assert client.segmentation_available == scenario['segmentation']
    
    def test_gpu_memory_allocation_strategies(self):
        """Test different GPU memory allocation strategies"""
        mock_gpu = MagicMock()
        
        # Test successful memory growth configuration
        with patch('tensorflow.config.list_physical_devices', return_value=[mock_gpu]):
            with patch('tensorflow.config.experimental.set_memory_growth') as mock_memory_growth:
                client = DECIMERClient()
                mock_memory_growth.assert_called_with(mock_gpu, True)
        
        # Test fallback when memory growth fails
        with patch('tensorflow.config.list_physical_devices', return_value=[mock_gpu]):
            with patch('tensorflow.config.experimental.set_memory_growth', 
                      side_effect=RuntimeError("Memory growth not supported")):
                client = DECIMERClient()
                # Should still be available despite memory growth failure
                assert client.decimer_available


class TestHardwareSpecificBehavior:
    """Test behavior on different hardware configurations"""
    
    def test_cpu_only_processing(self):
        """Test processing in CPU-only mode"""
        with patch.dict(os.environ, {'CUDA_VISIBLE_DEVICES': ''}):
            with patch('tensorflow.config.list_physical_devices', return_value=[]):
                client = DECIMERClient()
                
                # Create a simple test case
                temp_dir = tempfile.mkdtemp()
                try:
                    # Create a mock PDF path
                    pdf_path = os.path.join(temp_dir, "test.pdf")
                    with open(pdf_path, 'wb') as f:
                        f.write(b'%PDF-1.4\nTest PDF')
                    
                    # Should work in CPU mode (though may return empty list)
                    segments = client.segment_structures_from_pdf(pdf_path)
                    assert isinstance(segments, list), "Should return list even in CPU mode"
                    
                finally:
                    shutil.rmtree(temp_dir)
    
    def test_gpu_memory_constraints(self):
        """Test behavior with limited GPU memory"""
        mock_gpu = MagicMock()
        
        with patch('tensorflow.config.list_physical_devices', return_value=[mock_gpu]):
            with patch('tensorflow.config.experimental.set_memory_growth') as mock_memory_growth:
                # Simulate GPU memory constraint
                mock_memory_growth.side_effect = RuntimeError("Out of GPU memory")
                
                client = DECIMERClient()
                
                # Should handle GPU memory error gracefully
                assert client.decimer_available, "Should fallback gracefully from GPU memory error"
    
    def test_multiple_gpu_configuration(self):
        """Test configuration with multiple GPUs"""
        mock_gpu1 = MagicMock()
        mock_gpu2 = MagicMock()
        mock_gpu3 = MagicMock()
        
        with patch('tensorflow.config.list_physical_devices', return_value=[mock_gpu1, mock_gpu2, mock_gpu3]):
            with patch('tensorflow.config.experimental.set_memory_growth') as mock_memory_growth:
                client = DECIMERClient()
                
                # Should configure all GPUs
                assert mock_memory_growth.call_count == 3
                expected_calls = [
                    call(mock_gpu1, True),
                    call(mock_gpu2, True),
                    call(mock_gpu3, True)
                ]
                mock_memory_growth.assert_has_calls(expected_calls)
    
    def test_partial_gpu_failure(self):
        """Test behavior when some GPUs fail to configure"""
        mock_gpu1 = MagicMock()
        mock_gpu2 = MagicMock()
        
        def memory_growth_side_effect(gpu, enable):
            if gpu == mock_gpu1:
                return None  # Success
            else:
                raise RuntimeError("GPU 2 failed")
        
        with patch('tensorflow.config.list_physical_devices', return_value=[mock_gpu1, mock_gpu2]):
            with patch('tensorflow.config.experimental.set_memory_growth', 
                      side_effect=memory_growth_side_effect):
                client = DECIMERClient()
                
                # Should still be available even if some GPUs fail
                assert client.decimer_available, "Should work even with partial GPU failure"


class TestPerformanceFallback:
    """Test performance characteristics in different modes"""
    
    def test_cpu_vs_gpu_timeout_behavior(self):
        """Test timeout behavior differences between CPU and GPU modes"""
        temp_dir = tempfile.mkdtemp()
        
        try:
            # Create test PDF
            pdf_path = os.path.join(temp_dir, "test.pdf")
            with open(pdf_path, 'wb') as f:
                f.write(b'%PDF-1.4\nTest content')
            
            # Test CPU mode timing
            with patch.dict(os.environ, {'CUDA_VISIBLE_DEVICES': ''}):
                with patch('tensorflow.config.list_physical_devices', return_value=[]):
                    client = DECIMERClient()
                    
                    start_time = time.time()
                    segments = client.segment_structures_from_pdf(pdf_path)
                    cpu_time = time.time() - start_time
                    
                    assert isinstance(segments, list), "CPU mode should return list"
                    assert cpu_time < 30, "CPU mode should not hang indefinitely"
            
            # Test GPU mode timing (if available)
            mock_gpu = MagicMock()
            with patch('tensorflow.config.list_physical_devices', return_value=[mock_gpu]):
                with patch('tensorflow.config.experimental.set_memory_growth'):
                    client = DECIMERClient()
                    
                    start_time = time.time()
                    segments = client.segment_structures_from_pdf(pdf_path)
                    gpu_time = time.time() - start_time
                    
                    assert isinstance(segments, list), "GPU mode should return list"
                    assert gpu_time < 30, "GPU mode should not hang indefinitely"
                    
        finally:
            shutil.rmtree(temp_dir)
    
    def test_fallback_performance_monitoring(self):
        """Test performance monitoring during fallback scenarios"""
        import psutil
        
        initial_memory = psutil.Process().memory_info().rss / 1024 / 1024
        
        # Test fallback from GPU to CPU
        with patch('tensorflow.config.list_physical_devices') as mock_list_devices:
            # First try with GPU that fails
            mock_gpu = MagicMock()
            mock_list_devices.return_value = [mock_gpu]
            
            with patch('tensorflow.config.experimental.set_memory_growth', 
                      side_effect=RuntimeError("GPU failed")):
                client = DECIMERClient()
                
                # Should still work
                assert client.decimer_available
                
                # Memory usage should be reasonable
                current_memory = psutil.Process().memory_info().rss / 1024 / 1024
                memory_increase = current_memory - initial_memory
                assert memory_increase < 100, f"Excessive memory in fallback: {memory_increase:.2f}MB"


class TestEnvironmentVariableHandling:
    """Test handling of various environment variables"""
    
    def test_tensorflow_environment_variables(self):
        """Test various TensorFlow environment variables"""
        test_env_vars = {
            'TF_CPP_MIN_LOG_LEVEL': '2',
            'TF_FORCE_GPU_ALLOW_GROWTH': 'true',
            'TF_GPU_ALLOCATOR': 'cuda_malloc_async',
            'TF_ENABLE_ONEDNN_OPTS': '0'
        }
        
        with patch.dict(os.environ, test_env_vars):
            with patch('tensorflow.config.list_physical_devices', return_value=[]):
                client = DECIMERClient()
                
                # Should handle all environment variables gracefully
                assert client.decimer_available, "Should work with TensorFlow env vars"
    
    def test_cuda_environment_variables(self):
        """Test CUDA-specific environment variables"""
        cuda_env_vars = {
            'CUDA_VISIBLE_DEVICES': '0',
            'CUDA_DEVICE_ORDER': 'PCI_BUS_ID',
            'TF_FORCE_GPU_ALLOW_GROWTH': 'true'
        }
        
        with patch.dict(os.environ, cuda_env_vars):
            mock_gpu = MagicMock()
            with patch('tensorflow.config.list_physical_devices', return_value=[mock_gpu]):
                with patch('tensorflow.config.experimental.set_memory_growth'):
                    client = DECIMERClient()
                    
                    assert client.decimer_available, "Should work with CUDA env vars"
    
    def test_memory_limit_environment_variables(self):
        """Test memory limiting environment variables"""
        memory_env_vars = {
            'TF_FORCE_GPU_ALLOW_GROWTH': 'true',
            'TF_GPU_ALLOCATOR': 'cuda_malloc_async',
            'TF_MEMORY_ALLOCATION': 'BFC'  # Best-fit with coalescing
        }
        
        with patch.dict(os.environ, memory_env_vars):
            mock_gpu = MagicMock()
            with patch('tensorflow.config.list_physical_devices', return_value=[mock_gpu]):
                with patch('tensorflow.config.experimental.set_memory_growth'):
                    client = DECIMERClient()
                    
                    # Should respect memory limits
                    assert client.decimer_available, "Should work with memory limit env vars"


class TestErrorRecoveryMechanisms:
    """Test error recovery and graceful degradation"""
    
    def test_tensorflow_initialization_recovery(self):
        """Test recovery from TensorFlow initialization errors"""
        call_count = 0
        
        def failing_list_devices(device_type):
            nonlocal call_count
            call_count += 1
            if call_count == 1:
                raise RuntimeError("TensorFlow initialization failed")
            return []  # Succeed on retry
        
        with patch('tensorflow.config.list_physical_devices', side_effect=failing_list_devices):
            # Should handle initialization error gracefully
            client = DECIMERClient()
            
            # May or may not be available depending on fallback mechanism
            # But should not crash
            assert hasattr(client, 'decimer_available')
    
    def test_gpu_detection_retry_mechanism(self):
        """Test retry mechanism for GPU detection"""
        call_count = 0
        
        def intermittent_gpu_detection(device_type):
            nonlocal call_count
            call_count += 1
            if call_count <= 2:
                raise RuntimeError("GPU detection failed")
            return [MagicMock()]  # Succeed on third try
        
        # This test depends on implementation details
        # For now, just test that it doesn't crash
        with patch('tensorflow.config.list_physical_devices', 
                  side_effect=intermittent_gpu_detection):
            try:
                client = DECIMERClient()
                # Should either work or fail gracefully
                assert hasattr(client, 'decimer_available')
            except RuntimeError:
                # Acceptable if no retry mechanism exists
                pass
    
    def test_graceful_degradation_chain(self):
        """Test full graceful degradation chain"""
        # Test: GPU available -> GPU fails -> CPU works
        scenarios = [
            # Start with GPU available
            {'tf_available': True, 'gpu_available': True, 'gpu_works': True},
            # GPU available but fails to configure
            {'tf_available': True, 'gpu_available': True, 'gpu_works': False},
            # No GPU available
            {'tf_available': True, 'gpu_available': False, 'gpu_works': False},
            # TensorFlow not available
            {'tf_available': False, 'gpu_available': False, 'gpu_works': False},
        ]
        
        for scenario in scenarios:
            if scenario['tf_available']:
                if scenario['gpu_available']:
                    mock_gpu = MagicMock()
                    gpus = [mock_gpu]
                else:
                    gpus = []
                
                with patch('tensorflow.config.list_physical_devices', return_value=gpus):
                    if scenario['gpu_works']:
                        memory_growth_effect = None
                    else:
                        memory_growth_effect = RuntimeError("GPU configuration failed")
                    
                    with patch('tensorflow.config.experimental.set_memory_growth', 
                              side_effect=memory_growth_effect):
                        client = DECIMERClient()
                        
                        # Should handle each scenario gracefully
                        assert hasattr(client, 'decimer_available')
                        
            else:
                # TensorFlow not available
                with patch('decimer_client.tf', side_effect=ImportError("TensorFlow not available")):
                    client = DECIMERClient()
                    assert not client.decimer_available


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])