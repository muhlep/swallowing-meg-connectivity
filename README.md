# swallowing-meg-connectivity
MATLAB/FieldTrip pipeline for source-level MEG connectivity analysis of swallowing

# Source-level MEG Connectivity Analysis Pipeline

This repository contains MATLAB scripts implementing a complete source-level MEG connectivity analysis pipeline based on FieldTrip. The pipeline is designed for data recorded on CTF MEG systems and enables the estimation of frequency-specific functional and directed connectivity between cortical regions of interest (ROIs).

The code supports end-to-end processing starting from raw MEG datasets and culminates in group-level statistical analyses and visualization of large-scale cortical networks.

---

## Overview of the Pipeline

The pipeline implements the following main analysis steps:

1. **Trial Definition and Preprocessing**
   - Raw MEG data recorded on CTF systems (`.ds` format)
   - Trials defined based on manually marked event triggers
   - Analysis focuses on EMG2 (active condition) and EMG3 (control condition)

2. **Source Reconstruction**
   - Linearly Constrained Minimum Variance (LCMV) beamformer
   - Source-level time series reconstructed using FieldTrip
   - Fixed orientation source estimates

3. **ROI Time Series Extraction**
   - Aggregation of source estimates within predefined cortical ROIs
   - Generation of representative ROI-level time series

4. **Connectivity Analysis**
   - Weighted Phase Lag Index (wPLI)
   - Phase Slope Index (PSI)
   - Frequency-specific analysis across standard bands (theta, alpha, beta, low gamma, high gamma)

5. **Statistical Analysis**
   - Condition contrasts and time-window comparisons
   - Laterality indices
   - Cluster-based permutation testing

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

FieldTrip must be installed and correctly added to the MATLAB path prior to running the scripts.

---

## Repository Structure
preprocessing/     Trial definition and preprocessing
beamformer/        LCMV beamformer source reconstruction
roi_timeseries/    ROI-level source time series extraction
connectivity/      wPLI and PSI connectivity analysis
statistics/        Statistical analysis and permutation testing

Each subfolder contains a dedicated `README.md` describing its role within the pipeline.

---

## Notes

- The code is provided for transparency and reproducibility of the analysis approach.
- Scripts are organized modularly and can be adapted to related MEG connectivity studies.
- Visualization scripts used for figure generation are intentionally not included as core pipeline components.

---

## Citation

If you use this code or parts of the analysis pipeline in your work, please cite:

- Oostenveld et al., 2011. *FieldTrip: Open Source Software for Advanced Analysis of MEG, EEG, and Invasive Electrophysiological Data*. Computational Intelligence and Neuroscience.

Additional methodological citations should follow the corresponding primary publications.
