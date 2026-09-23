from config.pythonConsoleAutoImport import *
import numpy as np
import os

#Use script to get npy and origin and spacing from roi data in dragonfly

def extract_roi_to_numpy(roi_id, save_path):
    """Extract ROI data to NumPy file and print spacing and origin"""
    os.makedirs(save_path, exist_ok=True)
    
    roi = orsObj(roi_id)

    # Get numpy array
    numpy_array = roi.getNDArray()
    
    # Get spacing 
    x_spacing = roi.getXSpacing()
    y_spacing = roi.getYSpacing()
    z_spacing = roi.getZSpacing()
    spacing = (x_spacing, y_spacing, z_spacing)
    
    # Get position in space
    origin = roi.getOrigin()
    
    # Save just the numpy array
    npy_path = os.path.join(save_path, 'points.npy')
    np.save(npy_path, numpy_array)
    
    # Print spacing and origin to console
    print(f"Spacing: {spacing}")
    print(f"Origin: {origin}")
    print(f"Data saved to {npy_path}")
    
    return npy_path

# Run in dragonfly python with:  exec(open(r'C:\Users\Admin\Desktop\michael\UTM\extract_roi.py').read())
if __name__ == "__main__":
    roi_id = 'FF3E66C03D004BF1800669798C13A0B6CxvVolume_ROI' #change
    save_path = r'C:\Users\Admin\Desktop\michael\UTM' #change
    
    extract_roi_to_numpy(roi_id, save_path) 