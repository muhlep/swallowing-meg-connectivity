# swallowing-meg-connectivity
MATLAB / FieldTrip pipeline for source-level MEG connectivity analysis of swallowing

# Source-level MEG Connectivity Analysis Pipeline

This repository contains MATLAB scripts implementing a complete source-level MEG connectivity analysis pipeline based on FieldTrip. The pipeline is designed for data recorded on CTF MEG systems and enables the estimation of frequency-specific functional and directed connectivity between cortical regions of interest (ROIs).

The code supports end-to-end processing starting from raw MEG datasets and culminates in group-level statistical analyses of large-scale cortical networks.

---

## Overview of the Pipeline

The pipeline implements the following main analysis steps:

1. **Trial Definition and Preprocessing**
   - Raw MEG data recorded on CTF systems (`.ds` format)
   - Trials defined based on manually marked event triggers
   - Analysis focuses on EMG2 (active swallowing condition) and EMG3 (control condition)
   - Optional band-pass filtering and line-noise removal

2. **Source Reconstruction**
   - Linearly Constrained Minimum Variance (LCMV) beamformer
   - Source reconstruction performed using FieldTrip
   - Fixed-orientation source estimates on a warped MNI grid

3. **ROI-Level Signal Extraction**
   - ROI-level time series are computed internally during beamforming
   - Source estimates are aggregated across predefined cortical ROIs
   - No standalone ROI time-series extraction module is provided

4. **Connectivity Analysis**
   - Weighted Phase Lag Index (wPLI) for undirected functional connectivity
   - Phase Slope Index (PSI) for directed connectivity
   - Frequency-specific analysis across standard bands:
     - theta, alpha, beta, low gamma, high gamma

5. **Statistical Analysis**
   - Group-level contrasts (EMG2 vs. EMG3)
   - Laterality indices
   - Cluster-based permutation testing (CBPT) at the ROI × ROI level

---

## Data Requirements

- **MEG system:** CTF
- **Data format:** Native CTF `.ds` datasets
- **Event markers:** Manually defined triggers (EMG2, EMG3)
- **Pipeline entry point:** Raw `.ds` files (no preprocessed input required)

---

## Software Requirements

- **MATLAB**
- **FieldTrip** (tested with version **2023-12-20**)

FieldTrip must be installed and added to the MATLAB path prior to running the scripts.

---

## Repository Structure
- preprocessing/     Trial definition and preprocessing
- beamformer/        LCMV beamformer source reconstruction and ROI aggregation
- connectivity/      wPLI and PSI connectivity analysis
- statistics/        Statistical analysis and permutation testing

Each subfolder contains a dedicated `README.md` describing its role within the pipeline.

---

## Notes

- The code is provided for transparency and reproducibility of the analysis approach.
- Scripts are organized modularly and can be adapted to related MEG connectivity studies.
- Visualization scripts used for figure generation are intentionally not included as core pipeline components.

---

## Citation

If you use this code or parts of the analysis pipeline in your work, please cite:

Oostenveld, R., Fries, P., Maris, E., & Schoffelen, J.-M. (2011).  
*FieldTrip: Open Source Software for Advanced Analysis of MEG, EEG, and Invasive Electrophysiological Data*.  
Computational Intelligence and Neuroscience, 2011, Article ID 156869.

Additional methodological citations should follow the corresponding primary publications.
