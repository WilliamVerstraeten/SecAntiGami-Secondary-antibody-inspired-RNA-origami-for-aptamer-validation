#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Start in python environment from terminal.
# This code opens an interactive window for each microscopy image to be analyzed.
# The user can chose the GUVS that should be analyzed by hand.
# Those GUVs are automatically analyzed and appended to the output list.
# This code was set up with the aid of ChatGPT (chatgpt.com).

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Patch
from skimage import io, morphology, measure
from skimage.filters import gaussian, threshold_multiotsu
from skimage.morphology import disk
from scipy import ndimage
import gc

# =========================
# PARAMETERS
# =========================
parentfolder = "XXX"
newpath = os.path.join(parentfolder, "processed")
os.makedirs(newpath, exist_ok=True)
file_output = os.path.join(newpath, "output_Bio_12.csv")

min_object_size = 1000
ring_width = 2
save_qc = True
offsets = [0, 3, 6]

# =========================
# HELPER FUNCTIONS
# =========================
def select_guvs(Lipid_image, labels, regions):
    """
    Opens and interactve Window for choosing the GUVs to be analyzed.
    Green = Chosen, Red = ignored
    """
    selected_labels = set()
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.imshow(Lipid_image, cmap='gray')
    ax.set_title("Click GUVs to select (green=selected, red=ignored)")
    plt.axis('off')

    def redraw():
        ax.clear()
        ax.imshow(Lipid_image, cmap='gray')
        for r in regions:
            mask = labels == r.label
            color = 'green' if r.label in selected_labels else 'red'
            ax.contour(mask, levels=[0.5], colors=color, linewidths=1.5)
            y, x = r.centroid
            ax.text(x, y, str(r.label), color='white', fontsize=9, ha='center')
        fig.canvas.draw_idle()

    def on_click(event):
        if event.inaxes != ax:
            return
        y, x = int(event.ydata), int(event.xdata)
        lbl = labels[y, x]
        if lbl == 0:
            return
        if lbl in selected_labels:
            selected_labels.remove(lbl)
        else:
            selected_labels.add(lbl)
        redraw()

    fig.canvas.mpl_connect('button_press_event', on_click)
    redraw()
    plt.show(block=True)

    print("Please select GUVs and close Window for progress.")
    return selected_labels

# =========================
# MAIN LOOP
# =========================
image_counter = 0
header_written = False

for img in sorted(os.listdir(parentfolder)):
    if not img.lower().endswith(".tif"):
        continue

    image_counter += 1
    img_path = os.path.join(parentfolder, img)
    print(f"\nProcessing {img}")

    # -------------------------
    # IMAGE LOADING
    # -------------------------
    temp = io.imread(img_path)

    if temp.ndim < 3 or temp.shape[2] < 2:
        print(f"{img} does not have neough channels. Skipping.")
        continue

    RNA_image = temp[0, :, :] # fit to data (Biotin data: (C, X, Y); GFP data: (X, Y, C))
    Lipid_image = temp[1, :, :] # fit to data

    background = np.median(RNA_image)
    background_mean = np.mean(RNA_image)

    # -------------------------
    # GUV DETECTION
    # -------------------------
    Lipid_norm = Lipid_image / np.max(Lipid_image)
    smoothed = gaussian(Lipid_norm, sigma=2, preserve_range=True)
    thresholds = threshold_multiotsu(smoothed, classes=3)
    binary = smoothed > thresholds[0]
    binary = morphology.remove_small_objects(binary, min_size=min_object_size)
    binary = ndimage.binary_fill_holes(binary)

    labels = measure.label(binary)
    regions = measure.regionprops(labels)

    if len(regions) == 0:
        print("No GUVs found – Skipping image.")
        del temp, RNA_image, Lipid_image, binary, labels, regions
        gc.collect()
        continue

    # -------------------------
    # INTERACTIVE SELECTION
    # -------------------------
    selected_labels = select_guvs(Lipid_image, labels, regions)

    # -------------------------
    # SAVE FINAL SELECTED GUVS
    # -------------------------
    qc_selection_dir = os.path.join(newpath, "selected_areas")
    os.makedirs(qc_selection_dir, exist_ok=True)

    fig_sel, ax_sel = plt.subplots(figsize=(6,6))
    ax_sel.imshow(Lipid_image, cmap='gray')
    plt.axis('off')

    for region in regions:
        mask = labels == region.label
        color = 'green' if region.label in selected_labels else 'red'
        ax_sel.contour(mask, levels=[0.5], colors=color, linewidths=1.5)
        y, x = region.centroid
        ax_sel.text(x, y, str(region.label), color='white', fontsize=9, ha='center')

    ax_sel.set_title(f"Selected GUVs – Image {image_counter}")
    plt.savefig(os.path.join(qc_selection_dir, f"{image_counter}_{img}_selected_areas.png"),
                dpi=300, bbox_inches='tight')
    plt.close(fig_sel)

    if len(selected_labels) == 0:
        print("No GUVs chosen – skipping image.")
        del temp, RNA_image, Lipid_image, binary, labels, regions
        gc.collect()
        continue

    # -------------------------
    # ANALYSIS
    # -------------------------
    rows = []
    offset_list = offsets + [-o for o in offsets if o != 0]
    colors = cm.viridis(np.linspace(0, 1, len(offset_list)))
    offset_to_color = dict(zip(offset_list, colors))

    fig_qc, ax_qc = plt.subplots(figsize=(6, 6))
    ax_qc.imshow(Lipid_image, cmap='gray')
    plt.axis('off')

    for region in regions:
        if region.label not in selected_labels:
            continue
        base_mask = labels == region.label

        for offset in offset_list:
            mask = base_mask.copy()
            if offset > 0:
                mask = morphology.binary_dilation(mask, disk(offset))
            elif offset < 0:
                mask = morphology.binary_erosion(mask, disk(abs(offset)))

            if np.sum(mask) == 0:
                continue

            outer = morphology.binary_dilation(mask, disk(ring_width))
            inner = morphology.binary_erosion(mask, disk(ring_width))
            ring = outer ^ inner
            if np.sum(ring) == 0:
                continue

            intensities = RNA_image[ring]

            rows.append({
                "image_nr": image_counter,
                "img_name": img,
                "guv_label": region.label,
                "offset_px": offset,
                "ring_width_px": ring_width,
                "intensity_mean": np.mean(intensities),
                "intensity_median": np.median(intensities),
                "background_median": background,
                "background_mean": background_mean,
                "area_ring_px": np.sum(ring)
            })

            if save_qc:
                overlay = np.zeros((*ring.shape, 4))
                overlay[..., :3] = offset_to_color[offset][:3]
                overlay[..., 3] = ring * 0.4
                ax_qc.imshow(overlay)

    # -------------------------
    # SAVE QC IMAGE
    # -------------------------
    if save_qc:
        legend = [Patch(facecolor=offset_to_color[o][:3], label=f"{o} px", alpha=0.4)
                  for o in offset_list]
        ax_qc.legend(handles=legend, fontsize=8)
        ax_qc.set_title(f"Selected GUVs – Image {image_counter}")

        qc_dir = os.path.join(newpath, "synthetic_rings")
        os.makedirs(qc_dir, exist_ok=True)
        plt.savefig(os.path.join(qc_dir, f"{image_counter}_{img}_selected.png"),
                    dpi=300, bbox_inches='tight')
        plt.close(fig_qc)

    # -------------------------
    # SAVE CSV
    # -------------------------
    df = pd.DataFrame(rows)
    if not header_written:
        df.to_csv(file_output, index=False, mode='w')
        header_written = True
    else:
        df.to_csv(file_output, index=False, mode='a', header=False)

    # -------------------------
    # CLEANUP
    # -------------------------
    del temp, RNA_image, Lipid_image, binary, labels, regions, rows
    gc.collect()
    print(f"{img} completely processed.")
