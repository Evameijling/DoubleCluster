# Siblings or strangers; Determining the distance and age of the Double Cluster h and χ Persei

This repository contains the code developed for my **Bachelor thesis** at Amsterdam University College, which investigates the **Double Cluster (NGC 869 and NGC 884)** in the Perseus constellation. The research aimed to determine the **distance and age** of these clusters and explore whether they originate from the same molecular cloud.

## Research Overview

Open clusters are valuable for studying stellar formation and evolution, as their stars are born at the same time and under similar conditions. The **Double Cluster** presents an interesting case due to the proximity of its two components, leading to the question:  
➡ *Are NGC 869 and NGC 884 truly related siblings, or just chance (visual) neighbors?*

### **Observations & Data Collection**
- Conducted own observations at the **Anton Pannekoek Observatory** using a 51cm Richey-Chretien and a 35cm Meade Schmidt-Cassegrain telescope.
- Captured images using **B, V, and R broadband filters** with a total exposure time of **2 hours and 7 minutes**.
- Used **Gaia Data Release 3** to supplement photometric data with **parallax measurements**.

### **Analysis Workflow**
1. **Star Detection** – Identified sources in astronomical images.
2. **Astrometric Calibration** – Mapped pixel coordinates to celestial coordinates.
3. **Photometric Calibration** – Measured and calibrated stellar magnitudes.
4. **Data Filtering** – Removed low-quality detections and non-cluster stars.
5. **Distance & Age Determination** – Applied **color-magnitude diagram (CMD) fitting** and **zero-age main sequence (ZAMS) models** to estimate cluster parameters.

## Key Findings
- **Distance to NGC 869**: **2565 ± 85 pc**  
- **Distance to NGC 884**: **2380 ± 70 pc**  
- **Age of both clusters**: **20–40 million years**  
- The clusters **may** have originated from the same molecular cloud, but their spatial separation suggests they are distinct clusters.

## Repository Structure
- `1_Starfinder.py` – Detects stars in images.
- `2_Pltsolved.py` – Performs astrometric calibration.
- `3_Magcalibrate.py` – Conducts photometric calibration.
- `4_Combined.py` – Merges multiple observations.
- `5_Filter.py` – Filters data to remove noise and outliers.
- `6_Stardistribution.py` – Analyzes the spatial distribution of stars.
- `7_Results.py` – Generates CMD plots and calculates distances.
- `8_Membership.py` – Evaluates cluster membership probabilities.

## Requirements
- **Python 3.x**
- `astropy`
- `numpy`
- `matplotlib`
- `pandas`
- `astroquery`
- `scipy`
- `photutils`

*This research was conducted as part of my Bachelor thesis at Amsterdam University College. For more details, contact me at evameijling@gmail.com.*
