#!/usr/bin/env python3
"""
Transform Harvard-Oxford atlas from FSL's MNI152 space to MNI152NLin6Asym space at 2mm resolution.

This script:
1. Gets the MNI152NLin6Asym template from templateflow
2. Uses FSL's flirt to calculate the transformation
3. Transforms all Harvard-Oxford atlas files to the new space

Requirements:
- FSL installed and FSLDIR environment variable set
- templateflow Python package: pip install templateflow
"""

import os
import subprocess
import shutil
from pathlib import Path
from templateflow import api as tflow

mask_dir = '/project/3018051.01/ruggero/derivatives/masks'

def run_command(cmd, description=""):
    """Run a shell command and handle errors."""
    if description:
        print(f"{description}...")
    
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"Error: {result.stderr}")
        raise RuntimeError(f"Command failed: {cmd}")
    
    return result.stdout


def apply_transform(input_mask, output_mask, binary=True):
    """
    Apply the FSL->NLin6Asym transformation to a mask file.
    
    Simple wrapper that uses the transformation matrix created by main().
    
    Parameters:
    -----------
    input_mask : str
        Path to input mask file in FSL MNI152 space
    output_mask : str
        Path for output mask file in MNI152NLin6Asym space
    binary : bool
        If True, uses nearest neighbor interpolation (for binary/label masks)
        If False, uses trilinear interpolation (for probabilistic masks)
    
    Example:
    --------
    # After running the main script:
    apply_transform('my_mask.nii.gz', 'my_mask_MNI152NlLin6Asym.nii.gz')
    apply_transform('prob_mask.nii.gz', 'prob_mask_MNI152NlLin6Asym.gz', binary=False)
    """
    transform_matrix = os.path.join(mask_dir, "MNI152_to_MNI152NLin6Asym_2mm.mat")
    reference_template = os.path.join(mask_dir, "MNI152NLin6Asym_T1_2mm_brain.nii.gz")
    
    # Check if required files exist
    if not os.path.exists(transform_matrix):
        raise FileNotFoundError(
            f"{transform_matrix} not found. Please run the main script first."
        )
    if not os.path.exists(reference_template):
        raise FileNotFoundError(
            f"{reference_template} not found. Please run the main script first."
        )
    if not os.path.exists(input_mask):
        raise FileNotFoundError(f"Input mask {input_mask} not found.")
    
    # Choose interpolation method
    interp = 'nearestneighbour' if binary else 'trilinear'
    
    # Apply transformation
    flirt_cmd = (
        f"flirt -in {input_mask} "
        f"-ref {reference_template} "
        f"-applyxfm -init {transform_matrix} "
        f"-interp {interp} "
        f"-out {output_mask}"
    )
    
    run_command(flirt_cmd, f"Transforming {input_mask}")
    print(f"  Output saved: {output_mask}")
    
    
def main():
    # Check if FSL is installed
    fsl_dir = os.environ.get('FSLDIR')
    if not fsl_dir:
        raise EnvironmentError("FSLDIR not set. Please install FSL and set FSLDIR environment variable.")
    
    print("=" * 60)
    print("Harvard-Oxford Atlas Transformation to MNI152NLin6Asym (2mm)")
    print("=" * 60)
    print()
    
    # Define paths
    fsl_template = os.path.join(fsl_dir, "data/standard/MNI152_T1_2mm_brain.nii.gz")
    ho_dir = os.path.join(fsl_dir, "data/atlases/HarvardOxford")
    output_dir = Path(os.path.join(mask_dir, "harvardoxford_nlin6asym_2mm"))
    output_dir.mkdir(exist_ok=True)
    
    # Step 1: Get MNI152NLin6Asym template from templateflow
    print("Step 1: Getting MNI152NLin6Asym template from templateflow...")
    try:
        nlin6_template_path = tflow.get(
            'MNI152NLin6Asym', 
            resolution=2, 
            desc='brain', 
            suffix='T1w', 
            extension='nii.gz'
        )
        print(f"  Template found: {nlin6_template_path}")
        
        # Copy to working directory for easy reference
        nlin6_template_local = "MNI152NLin6Asym_T1_2mm_brain.nii.gz"
        shutil.copy(nlin6_template_path, nlin6_template_local)
        print(f"  Template copied to: {nlin6_template_local}")
    except Exception as e:
        print(f"Error getting template: {e}")
        raise
    
    print()
    
    # Step 2: Calculate transformation matrix
    print("Step 2: Calculating transformation between MNI spaces...")
    transform_matrix = "FSL_to_NLin6Asym_2mm.mat"
    
    flirt_cmd = (
        f"flirt -in {fsl_template} "
        f"-ref {nlin6_template_local} "
        f"-omat {transform_matrix} "
        f"-dof 12 "
        f"-searchrx -30 30 "
        f"-searchry -30 30 "
        f"-searchrz -30 30"
    )
    
    run_command(flirt_cmd, "  Running FLIRT registration")
    print(f"  Transformation matrix saved: {transform_matrix}")
    print()
    
    # Step 3: Transform Harvard-Oxford atlases
    print("Step 3: Transforming Harvard-Oxford atlases...")
    
    # Track successful transformations
    transformed_files = []
    
    # Transform cortical maxprob atlases
    print("  Transforming cortical atlases...")
    for thr in [0, 25, 50]:
        input_file = os.path.join(ho_dir, f"HarvardOxford-cort-maxprob-thr{thr}-2mm.nii.gz")
        output_file = output_dir / f"HarvardOxford-cort-maxprob-thr{thr}-2mm.nii.gz"
        
        if os.path.exists(input_file):
            flirt_cmd = (
                f"flirt -in {input_file} "
                f"-ref {nlin6_template_local} "
                f"-applyxfm -init {transform_matrix} "
                f"-interp nearestneighbour "
                f"-out {output_file}"
            )
            run_command(flirt_cmd)
            print(f"    ✓ Cortical maxprob threshold {thr}")
            transformed_files.append(str(output_file))
        else:
            print(f"    ✗ Cortical maxprob threshold {thr} (not found)")
    
    # Transform subcortical maxprob atlases
    print("  Transforming subcortical atlases...")
    for thr in [0, 25, 50]:
        input_file = os.path.join(ho_dir, f"HarvardOxford-sub-maxprob-thr{thr}-2mm.nii.gz")
        output_file = output_dir / f"HarvardOxford-sub-maxprob-thr{thr}-2mm.nii.gz"
        
        if os.path.exists(input_file):
            flirt_cmd = (
                f"flirt -in {input_file} "
                f"-ref {nlin6_template_local} "
                f"-applyxfm -init {transform_matrix} "
                f"-interp nearestneighbour "
                f"-out {output_file}"
            )
            run_command(flirt_cmd)
            print(f"    ✓ Subcortical maxprob threshold {thr}")
            transformed_files.append(str(output_file))
        else:
            print(f"    ✗ Subcortical maxprob threshold {thr} (not found)")
    
    # Transform cortical probabilistic maps
    print("  Transforming cortical probabilistic maps...")
    cortical_prob_files = list(Path(ho_dir).glob("HarvardOxford-cort-prob-*.nii.gz"))
    for idx, probmap in enumerate(cortical_prob_files, 1):
        output_file = output_dir / probmap.name
        
        flirt_cmd = (
            f"flirt -in {probmap} "
            f"-ref {nlin6_template_local} "
            f"-applyxfm -init {transform_matrix} "
            f"-interp trilinear "
            f"-out {output_file}"
        )
        run_command(flirt_cmd)
        transformed_files.append(str(output_file))
        
        # Print progress
        if idx % 10 == 0:
            print(f"    ✓ Transformed {idx}/{len(cortical_prob_files)} cortical maps")
    
    print(f"    ✓ All {len(cortical_prob_files)} cortical probabilistic maps transformed")
    
    # Transform subcortical probabilistic maps
    print("  Transforming subcortical probabilistic maps...")
    subcortical_prob_files = list(Path(ho_dir).glob("HarvardOxford-sub-prob-*.nii.gz"))
    for idx, probmap in enumerate(subcortical_prob_files, 1):
        output_file = output_dir / probmap.name
        
        flirt_cmd = (
            f"flirt -in {probmap} "
            f"-ref {nlin6_template_local} "
            f"-applyxfm -init {transform_matrix} "
            f"-interp trilinear "
            f"-out {output_file}"
        )
        run_command(flirt_cmd)
        transformed_files.append(str(output_file))
    
    print(f"    ✓ All {len(subcortical_prob_files)} subcortical probabilistic maps transformed")
    print()
    
    # Summary
    print("=" * 60)
    print("TRANSFORMATION COMPLETE!")
    print("=" * 60)
    print(f"Total files transformed: {len(transformed_files)}")
    print(f"Output directory: {output_dir.absolute()}")
    print(f"Transformation matrix: {transform_matrix}")
    print()
    print("Key output files:")
    print(f"  - Cortical: HarvardOxford-cort-maxprob-thr{{0,25,50}}-2mm.nii.gz")
    print(f"  - Subcortical: HarvardOxford-sub-maxprob-thr{{0,25,50}}-2mm.nii.gz")
    print(f"  - Individual probabilistic maps for each region")
    print()
    print("To verify alignment in FSLeyes:")
    print(f"  fsleyes {nlin6_template_local} \\")
    print(f"           your_fmriprep_bold.nii.gz \\")
    print(f"           {output_dir}/HarvardOxford-cort-maxprob-thr25-2mm.nii.gz -cm random")
    print()

if __name__ == "__main__":
    try:
        # main()
        
        # apply transformation to existing masks in anatomical directory
        masks = ['IPLD.nii.gz', 'OFG_L.nii.gz', 'OFG_R.nii.gz', 'OFG.nii.gz', 'VOT_L.nii.gz']

        for mask in masks:
            mask = os.path.join(mask_dir, 'anatomical', mask)
            output = mask.replace('.nii.gz', '_MNI152NlLin6Asym.nii.gz')
            apply_transform(mask, output)
            
    except Exception as e:
        print(f"\nError: {e}")
        exit(1)