#!/usr/bin/env python3
"""
Quick Validation: Real Data Integration Test
=============================================
Tests the complete pipeline: NIfTI → Simulation → Output
"""
import sys
import tempfile
import numpy as np
from pathlib import Path

# Add project root to path
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

# Test 1: NIfTI loader
def test_nifti_loader():
    print("Test 1: NIfTI loader...")
    try:
        import nibabel as nib
        from src.dicom_loader import load_nifti, resample_to_grid
        
        # Create synthetic NIfTI
        with tempfile.NamedTemporaryFile(suffix='.nii.gz', delete=False) as f:
            data = np.random.rand(100, 100, 50).astype(np.float32)
            affine = np.eye(4)
            affine[0,0] = affine[1,1] = affine[2,2] = 1.0
            img = nib.Nifti1Image(data, affine)
            nib.save(img, f.name)
            temp_path = f.name
        
        # Load
        loaded_data, loaded_affine, meta = load_nifti(temp_path)
        assert loaded_data.shape == (100, 100, 50)
        print(f"  PASS: load_nifti: shape={loaded_data.shape}")
        
        # Resample
        resampled, target_affine = resample_to_grid(loaded_data, loaded_affine, (64, 64, 64))
        assert resampled.shape == (64, 64, 64)
        print(f"  PASS: resample_to_grid: {resampled.shape}")
        
        # Cleanup
        Path(temp_path).unlink()
        return True
    except Exception as e:
        print(f"  FAIL: {e}")
        return False


# Test 2: BraTS patient loader
def test_brats_loader():
    print("Test 2: BraTS patient loader...")
    try:
        from src.dicom_loader import load_brats_patient, prepare_simulation_inputs
        
        # Create fake BraTS directory structure
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            patient_dir = tmpdir / "Brats2021_00001"
            patient_dir.mkdir()
            
            # Create dummy modalities
            import nibabel as nib
            for modality in ['t1', 't1ce', 't2', 'flair', 'seg']:
                data = np.random.rand(100, 100, 50).astype(np.float32)
                if modality == 'seg':
                    data = np.random.randint(0, 4, (100, 100, 50)).astype(np.float32)
                affine = np.eye(4)
                affine[0,0] = affine[1,1] = affine[2,2] = 1.0
                img = nib.Nifti1Image(data, affine)
                nib.save(img, patient_dir / f"Brats2021_00001_{modality}.nii.gz")
            
            # Load
            patient_data = load_brats_patient(patient_dir, target_shape=(64, 64, 64))
            print(f"  PASS: load_brats_patient: keys={list(patient_data.keys())}")
            
            # Prepare simulation
            sim_inputs = prepare_simulation_inputs(patient_data, target_shape=(64, 64, 64))
            print(f"  PASS: prepare_simulation_inputs: u0={sim_inputs['u0'].shape}, tensor={sim_inputs['tensor_field'].shape}")
            
        return True
    except Exception as e:
        print(f"  FAIL: {e}")
        import traceback
        traceback.print_exc()
        return False


# Test 3: Inverse estimation with NIfTI
def test_inverse_estimation():
    print("Test 3: Inverse estimation with NIfTI...")
    try:
        # Create synthetic NIfTI segmentations
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            import nibabel as nib
            
            # Create two segmentations with known volume change
            for i, (vol_factor, name) in enumerate([(1.0, 't0'), (1.5, 't1')]):
                data = np.zeros((50, 50, 30), dtype=np.float32)
                # Sphere with radius proportional to volume
                center = (25, 25, 15)
                radius = int(5 * (vol_factor ** (1/3)))
                z, y, x = np.mgrid[0:50, 0:50, 0:30]
                dist = np.sqrt((x-center[0])**2 + (y-center[1])**2 + (z-center[2])**2)
                data[dist <= radius] = 1.0
                
                affine = np.eye(4)
                affine[0,0] = affine[1,1] = affine[2,2] = 1.0
                img = nib.Nifti1Image(data.astype(np.float32), affine)
                nib.save(img, f"/tmp/seg_{name}.nii.gz")
            
            # Test volume loading - import directly from file
            import importlib.util
            spec = importlib.util.spec_from_file_location("inv_est", PROJECT_ROOT / "src" / "51_inverse_parameter_estimation.py")
            inv_est = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(inv_est)
            
            V0, V1 = inv_est._load_volumes_from_brats(Path("/tmp/seg_t0.nii.gz"), Path("/tmp/seg_t1.nii.gz"))
            print(f"  PASS: Volume loading: V0={V0:.1f}, V1={V1:.1f} mm3")
            print(f"  Volume ratio: {V1/V0:.3f} (expected ~1.5)")
            
            return True
    except Exception as e:
        print(f"  FAIL: {e}")
        import traceback
        traceback.print_exc()
        return False


# Test 4: End-to-end 42_anisotropic_pde with real mask
def test_anisotropic_pde_real():
    print("Test 4: 42_anisotropic_pde with real mask...")
    try:
        # Create synthetic real mask and tensor
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            
            # Create mask NPZ
            mask = np.zeros((64, 64), dtype=np.float32)
            y, x = np.mgrid[0:64, 0:64]
            mask[(y-32)**2 + (x-32)**2 < 5**2] = 1.0
            np.savez(tmpdir / "real_mask.npz", u0=mask)
            
            # Create tensor NPZ
            D_xx = np.ones((64, 64), dtype=np.float32) * 0.0013
            D_xy = np.zeros((64, 64), dtype=np.float32)
            D_yy = np.ones((64, 64), dtype=np.float32) * 0.0013
            np.savez(tmpdir / "real_tensor.npz", D_xx=D_xx, D_xy=D_xy, D_yy=D_yy)
            
            # Run 42_anisotropic_pde with real data
            import subprocess
            result = subprocess.run([
                sys.executable, "src/42_anisotropic_pde.py",
                "--real-mask", str(tmpdir / "real_mask.npz"),
                "--real-tensor", str(tmpdir / "real_tensor.npz"),
                "--use-real-only",
                "--output-dir", str(tmpdir / "out"),
                "--grid-size", "64"
            ], capture_output=True, text=True, timeout=60, cwd=PROJECT_ROOT)
            
            if result.returncode == 0:
                print(f"  PASS: Real data simulation completed")
                print(f"  Output: {result.stdout[-500:]}")
                return True
            else:
                print(f"  FAIL: {result.stderr[:500]}")
                return False
    except Exception as e:
        print(f"  FAIL: {e}")
        return False


# Main
if __name__ == "__main__":
    print("=" * 60)
    print("REAL DATA INTEGRATION VALIDATION")
    print("=" * 60)
    
    tests = [
        test_nifti_loader,
        test_brats_loader,
        test_inverse_estimation,
        test_anisotropic_pde_real,
    ]
    
    passed = 0
    for test in tests:
        if test():
            passed += 1
        print()
    
    print("=" * 60)
    print(f"RESULTS: {passed}/{len(tests)} tests passed")
    print("=" * 60)
    
    if passed == len(tests):
        print("PASS: All tests passed!")
        exit(0)
    else:
        print("FAIL: Some tests failed")
        exit(1)