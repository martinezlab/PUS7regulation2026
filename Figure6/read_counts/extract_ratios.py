import os
import glob
import numpy as np
import pandas as pd
import tifffile as tiff
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches
from skimage import measure
from skimage.measure import regionprops_table
from cellpose import models
import gc

# ==========================================
# USER SETTINGS
# ==========================================
BASE_DIR = r"PUS7_localization"
DIR_IMAGE = "rep2_PUS7localization"
IMAGE_DIR = os.path.join(BASE_DIR, DIR_IMAGE)

# Output Directories
DIR_INDIVIDUAL = os.path.join(BASE_DIR, "rep2_individual_results")
DIR_MERGED = os.path.join(BASE_DIR, "rep2_merged_results")

os.makedirs(DIR_INDIVIDUAL, exist_ok=True)
os.makedirs(DIR_MERGED, exist_ok=True)

# Channel Definitions (0, 1, 2)
CHANNEL_NUCLEI_IDX = 0    # DAPI
CHANNEL_SIGNAL_IDX = 1    # FITC/PUS7
CHANNEL_CYTO_IDX   = 2    # TRANS/Whole Cell

# Cellpose Settings
DIAMETER_NUCLEI = 60
DIAMETER_CYTO = 100

# ==========================================
# INITIALIZE MODELS
# ==========================================
print("Initializing Cellpose models on GPU...")
model_nuc  = models.Cellpose(gpu=True, model_type='nuclei')
model_cyto = models.Cellpose(gpu=True, model_type='cyto')

# ==========================================
# FILE DISCOVERY
# ==========================================
search_pattern = os.path.join(IMAGE_DIR, "*_dapi.tiff")
dapi_files = glob.glob(search_pattern)
unique_image_names = [os.path.basename(f).replace("_dapi.tiff", "") for f in dapi_files]

print(f"Found {len(unique_image_names)} images to process.")

# ==========================================
# MAIN PROCESSING LOOP
# ==========================================
for i, IMAGE_NAME in enumerate(unique_image_names):
    
    # --- CHECKPOINT: RESUME CAPABILITY ---
    expected_csv = os.path.join(DIR_INDIVIDUAL, f"{IMAGE_NAME}_quantification.csv")
    if os.path.exists(expected_csv):
        print(f"[{i+1}/{len(unique_image_names)}] Skipping {IMAGE_NAME} (Result exists).")
        continue
    
    print(f"\n[{i+1}/{len(unique_image_names)}] Processing: {IMAGE_NAME}")

    try:
        # 1. CONSTRUCT FILENAMES
        FILENAME_DAPI  = f"{IMAGE_NAME}_dapi.tiff"
        FILENAME_FITC  = f"{IMAGE_NAME}_fitc.tiff"
        FILENAME_TRANS = f"{IMAGE_NAME}_trans.tiff"

        FILE_DAPI  = os.path.join(IMAGE_DIR, FILENAME_DAPI)
        FILE_FITC  = os.path.join(IMAGE_DIR, FILENAME_FITC)
        FILE_TRANS = os.path.join(IMAGE_DIR, FILENAME_TRANS)

        # 2. LOAD IMAGES
        if not (os.path.exists(FILE_DAPI) and os.path.exists(FILE_FITC) and os.path.exists(FILE_TRANS)):
            print(f"  -> Skipping {IMAGE_NAME}: Missing one or more channels.")
            continue

        img_nuc  = tiff.imread(FILE_DAPI)
        img_sig  = tiff.imread(FILE_FITC)
        img_cyto = tiff.imread(FILE_TRANS)

        # Stack and Move Axis to (Y, X, Channels)
        image_stack = np.stack([img_nuc, img_sig, img_cyto], axis=0)
        image = np.moveaxis(image_stack, 0, -1)

        # Extract Analysis Channels
        img_nuc_analysis  = image[:, :, CHANNEL_NUCLEI_IDX]
        img_sig_analysis  = image[:, :, CHANNEL_SIGNAL_IDX]
        img_cyto_analysis = image[:, :, CHANNEL_CYTO_IDX]

        # 3. SEGMENTATION
        # Nuclei
        mask_nuc, _, _, _ = model_nuc.eval(
            img_nuc_analysis, 
            diameter=DIAMETER_NUCLEI
        )
        
        # Cytoplasm
        mask_cyto, _, _, _ = model_cyto.eval(
            img_cyto_analysis, 
            diameter=DIAMETER_CYTO
        )

        # 4. RELABEL NUCLEI (Link Nucleus to Cytoplasm)
        nuc_labels_relabeled = np.zeros_like(mask_nuc)
        props = measure.regionprops(mask_nuc)

        for region in props:
            coords = region.coords
            cyto_vals_under_nuc = mask_cyto[coords[:, 0], coords[:, 1]]
            cyto_vals_under_nuc = cyto_vals_under_nuc[cyto_vals_under_nuc > 0]
            
            if len(cyto_vals_under_nuc) > 0:
                unique, counts = np.unique(cyto_vals_under_nuc, return_counts=True)
                dominant_cyto_label = unique[np.argmax(counts)]
                nuc_labels_relabeled[coords[:, 0], coords[:, 1]] = dominant_cyto_label

        # 5. MEASURE INTENSITY (Raw)
        df_cyto = pd.DataFrame(regionprops_table(
            mask_cyto, 
            intensity_image=img_sig_analysis, 
            properties=('label', 'mean_intensity', 'area')
        ))

        df_nuc = pd.DataFrame(regionprops_table(
            nuc_labels_relabeled, 
            intensity_image=img_sig_analysis, 
            properties=('label', 'mean_intensity', 'area')
        ))

        # Calculate Total Intensity
        df_cyto['total_int'] = df_cyto['mean_intensity'] * df_cyto['area']
        df_nuc['total_int']  = df_nuc['mean_intensity'] * df_nuc['area']

        # 6. CALCULATE RATIOS
        merged = pd.merge(df_nuc, df_cyto, on='label', suffixes=('_nuc', '_cyto'))
        
        # Filter: Cytoplasm must be larger than Nucleus (Geometric Sanity Check)
        valid = merged[merged['area_cyto'] > merged['area_nuc']].copy()

        # Calculate Cytosol Only
        valid['total_int_cytosol_only'] = valid['total_int_cyto'] - valid['total_int_nuc']

        # Ratios
        valid['ratio_nuc_vs_whole'] = valid['total_int_nuc'] / valid['total_int_cyto']
        valid['ratio_nuc_vs_cytosol'] = valid['total_int_nuc'] / valid['total_int_cytosol_only']
        valid['ratio_cytosol_vs_nuc'] = valid['total_int_cytosol_only'] / valid['total_int_nuc']
        
        # Log2 Ratios for plotting
        valid['log2_ratio_nuc'] = np.log2(valid['ratio_nuc_vs_cytosol'])
        valid['log2_ratio_cyto'] = np.log2(valid['ratio_cytosol_vs_nuc'])

        # extract image name
        parts = IMAGE_NAME.split('_')
        if len(parts) >= 2:
            group_id = f"{parts[0]}_{parts[1]}"
        else:
            group_id = "Unknown"

        valid['Image_Name'] = IMAGE_NAME
        valid['Group_ID'] = group_id

        # 7. SAVE INDIVIDUAL CSV
        valid.to_csv(os.path.join(DIR_INDIVIDUAL, f"{IMAGE_NAME}_quantification.csv"), index=False)

        # 8. VISUALIZATION - OVERLAY (Saved, not shown)
        fig, ax = plt.subplots(figsize=(10, 10))
        ax.imshow(img_sig_analysis, cmap='gray')
        ax.contour(mask_cyto, levels=np.unique(mask_cyto), colors='cyan', linewidths=0.5, alpha=0.5)
        ax.contour(nuc_labels_relabeled, levels=np.unique(nuc_labels_relabeled), colors='red', linewidths=0.5, alpha=0.8)
        
        ax.set_title(f"PUS7 Signal Overlay\n{IMAGE_NAME}")
        ax.axis('off')
        
        # Custom Legend
        red_patch = mpatches.Patch(color='red', label='Nuclei')
        cyan_patch = mpatches.Patch(color='cyan', label='Cytoplasm')
        plt.legend(handles=[red_patch, cyan_patch], loc='upper right', frameon=True, facecolor='black', labelcolor='white')
        
        plt.tight_layout()
        plt.savefig(os.path.join(DIR_INDIVIDUAL, f"{IMAGE_NAME}_overlay.png"), dpi=300, bbox_inches='tight')
        plt.close()

        # 9. VISUALIZATION - VIOLIN PLOTS (Saved, not shown)
        
        # N/C Plot
        plt.figure(figsize=(7, 6))
        sns.violinplot(data=valid, y='log2_ratio_nuc', color='skyblue', inner='quartile')
        plt.axhline(0, color='red', linestyle='--', alpha=0.7)
        plt.ylabel("Log2 (Nucleus / Cytoplasm)")
        plt.title(f"N/C Distribution\n{IMAGE_NAME}")
        plt.savefig(os.path.join(DIR_INDIVIDUAL, f"{IMAGE_NAME}_NC_violin.png"), dpi=300, bbox_inches='tight')
        plt.close()

        print(f"  -> Finished {IMAGE_NAME}. Found {len(valid)} cells.")

        # --- FORCE MEMORY CLEANUP ---
        del img_nuc, img_sig, img_cyto, image_stack, image
        del mask_nuc, mask_cyto, nuc_labels_relabeled
        del df_nuc, df_cyto, merged, valid
        gc.collect()

    except Exception as e:
        print(f"  -> ERROR processing {IMAGE_NAME}: {e}")
        gc.collect()

# ==========================================
# FINAL AGGREGATION
# ==========================================
print("\nProcessing complete. Aggregating results...")
# Find all CSVs produced by the loop
all_csv_files = glob.glob(os.path.join(DIR_INDIVIDUAL, "*_quantification.csv"))

if all_csv_files:
    master_df = pd.concat((pd.read_csv(f) for f in all_csv_files), ignore_index=True)
    
    unique_groups = master_df['Group_ID'].unique()
    print(f"Found {len(unique_groups)} unique conditions: {unique_groups}")
    
    for group in unique_groups:
        group_data = master_df[master_df['Group_ID'] == group]
        save_path = os.path.join(DIR_MERGED, f"{group}_merged.csv")
        group_data.to_csv(save_path, index=False)
        print(f"  -> Saved {group} ({len(group_data)} cells) to: {save_path}")

    # Save a global master file just in case
    master_df.to_csv(os.path.join(DIR_MERGED, "all_data.csv"), index=False)

else:
    print("No CSV files found in output directory.")