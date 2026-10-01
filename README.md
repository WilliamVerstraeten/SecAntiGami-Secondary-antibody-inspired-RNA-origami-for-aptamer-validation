# SecAntiGami — analysis code

Analysis code for the figures in:

> **SecAntiGami: A Secondary-Antibody-Inspired RNA Origami for Qualitative Aptamer Validation**
> Linus S. Frey, William Verstraeten, Mai P. Tran, Dana Venkert, Luca Monari, Cody Geary, Kerstin Göpfrich
> *Nano Letters* (2026). DOI: [10.1021/acs.nanolett.6c03022](https://doi.org/10.1021/acs.nanolett.6c03022)

The repository holds two independent analysis pipelines: dynamic light
scattering (DLS) and confocal microscopy of giant unilamellar vesicles (GUVs).
They share no code and can be run separately.

```
.
├── DLS/          dynamic light scattering — includes the raw data
└── Confocal/     GUV image processing and quantification — raw images on request
```

---

## DLS

Reproduces four panels from the raw Zetasizer exports, which are included in
`DLS/data/`.

| Section of the notebook | Panel | Output file |
|---|---|---|
| 1 | Figure 2E — titration of SecAntiGami against Biotin-KL | `combined_titration.svg` |
| 2 | Figure 2D — Biotin-KL and negative controls | `combined_Biotin-KL.svg` |
| 3 | Figure 3E — aptamer variants V1–V3 | `combined_Variants.svg` |
| 4 | Supplementary Figure S8d — SecAntiGami-1KL and N-Biotin-KL | `4KL_1KL_NBio_combined.svg` |

```bash
cd DLS
pip install -r requirements.txt
jupyter lab DLS_analysis.ipynb
```

Run the sections in order — section 3 reuses the tables written by section 2.
Processed tables land in `output/` and figures in `output/graphs/`; both are
created automatically.

**Method.** The instrument exports two data types; this pipeline uses the
per-measurement export only, taking the machine-computed main peak from the
`Pk 1 Mean Int (d.nm)` column. Three technical repeats are averaged to give one
replicate value, with the standard deviation and standard error of the mean
computed across them. The reported diameter *increase* is the mean peak size
after RNA addition minus the mean peak size of the same sample before RNA
addition, with uncertainties propagated through the subtraction using the
`uncertainties` package. Throughout, *n* refers to independent replicates
performed on different days, not to technical repeats.

Measurements that the instrument flagged as faulty (attenuator too high, count
rate too low, and similar) and measurements not used in the analysis were
removed from the `.txt` files before they entered this workflow.

The data behind Supplementary Figure S8d were acquired on a Zetasizer Ultra
rather than the Zetasizer Nano ZS used elsewhere, with a different
spectrophotometer and different RNA purification columns. See the Supporting
Information for details.

---

## Confocal

Quantifies SecAntiGami fluorescence at the GUV membrane. Two steps, run in
order.

**1 · `Preprocess_images_final.ipynb`** converts the microscope's `.czi` files
to `.tif`, collapsing singleton dimensions. Set `input_folder` to the directory
holding the `.czi` files; `.tif` files are written to a `tiff_output`
subfolder.

**2 · `Analyze_GUVs_final.py`** measures the membrane intensity. Run it from a
terminal in an interactive Python environment — it opens a window per image and
will not work headless:

```bash
cd Confocal
python Analyze_GUVs_final.py
```

Set `parentfolder` to the directory of `.tif` files before running. For each
image the script normalises the lipid channel, smooths it, segments the GUVs by
multi-Otsu thresholding, discards objects smaller than 1000 px and fills holes.
The detected GUVs are then presented in a clickable window: green outlines are
selected, red are ignored. This manual step exists so that GUVs which are out
of focus, deformed, aggregated, or already counted in another image can be
excluded. Every selection is written to `processed/selected_areas/` as a PNG,
so the choices remain auditable after the fact.

Each selected GUV is eroded and dilated by 2 px and the two areas subtracted,
giving a 4 px ring at the membrane. Mean and median intensities of that ring in
the RNA channel are written to `processed/output_*.csv`, together with the image
median and mean used as background.

Two details worth knowing before you run it:

- The script also measures rings displaced by ±3 px and ±6 px from the membrane
  (`offsets`), as a control for how sensitive the readout is to ring placement.
  **The values reported in the paper are the `offset_px == 0` rows.** Filter on
  that column when reproducing the figures.
- Channel order differs between datasets. The biotin images are `(C, X, Y)`
  and the GFP images are `(X, Y, C)`, so the two lines assigning `RNA_image`
  and `Lipid_image` must be adjusted to match the data being analysed.

Background subtraction is the image mean subtracted from the ring mean. For the
Broccoli SecAntiGami (Supplementary Figure S11) no correction was applied, as
the background intensity was zero.

### Raw confocal images are not included

The `.czi` and `.tif` microscopy images are not deposited here — the full set
runs to tens of gigabytes, which is impractical for a Git repository. **They
are available from the corresponding authors on reasonable request**
(c.geary@zmbh.uni-heidelberg.de; k.goepfrich@zmbh.uni-heidelberg.de).

The code in `Confocal/` is complete and will run on any equivalently structured
two-channel dataset; only the input paths need setting.

---

## Requirements

DLS: `pandas`, `numpy`, `matplotlib`, `uncertainties`, `jupyterlab` — see
`DLS/requirements.txt`.

Confocal: `numpy`, `pandas`, `matplotlib`, `scikit-image`, `scipy`, `czifile`,
`tifffile`. An interactive matplotlib backend is required for the GUV selection
step.


## Citation

If you use this code, please cite the paper above.
