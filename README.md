Hubble Space Telescope Image Analysis:

# Contents
1.Cosmic Ray Removal

2.Star-Finding

3.Photometry

4.Creating a Hertzsprung–Russell Diagram

This project implements a complete astronomical image analysis designed to process Hubble Space Telescope (HST) images 
-construct a Hertzsprung-Russell (HR) diagram. 
-This handles raw astronomical images from two HST filters (F336W and F555W), 
-processes them to remove artifacts, 
-identifies stellar sources, performs photometric measurements, and ultimately produces a color-magnitude diagram that reveals the stellar population characteristics.
-The HR diagram is one of the most important tools in astrophysics, showing the relationship between a star's color (temperature) and its brightness (luminosity). This allows astronomers to understand stellar evolution, age, and composition of star clusters.

# The Filters:
-*F336W*
-*F555W*

#Requirements:

This Pipeline requires the following Python libraries:
**Package purposes:**
- `numpy` - Numerical operations and array handling
- `pandas` - Data organization and catalog management
- `astropy` - FITS file handling and astronomical data structures
- `scipy` - Image processing and filtering operations
- `matplotlib` - Visualization and plotting

# Directory Structure

Before running the pipeline, organize your data as follows:
```
your-project-directory/
├── Data/
│   ├── F336W/              # Directory containing F336W filter FITS images
│   ├── F555W/              # Directory containing F555W filter FITS images
│   └── [outputs will be saved here]
├── optical.py         # Main pipeline script
└── README.md               # This file

# Task 1: Cosmic Ray Removal
**Purpose**: Cosmic rays are high-energy particles from space that create bright spots in astronomical images. By combining multiple exposures of the same field, we can identify and remove these artifacts.

# Task 2: Star-Finding
**Purpose**: Automatically identify stellar sources in the image field for subsequent analysis.

*different thresholds*
The F555W filter typically captures more sources at fainter magnitudes in this particular dataset, so a lower threshold helps match the detection counts between filters.

# Task 3: Aperture Photometry
**Purpose**: Measure the brightness (flux) of each detected star to calculate magnitudes.
for:
```
   magnitude = -2.5 × log₁₀(flux) + zeropoint
```
   where zeropoint = -21.1 (from hst workbook)

# Task 4: Hertzsprung-Russell Diagram
**Purpose**: Visualize the color-magnitude relationship of the stellar population.

**Outputs**:
- `hr_diagram.png` - Final color-magnitude diagram

## Output Files Description

| File | Description |
|------|-------------|
| `F336W_combined.fits` | Cosmic ray cleaned F336W image (FITS format) |
| `F555W_combined.fits` | Cosmic ray cleaned F555W image (FITS format) |
| `F336W_comparison.png` | Side-by-side comparison: single exposure vs. combined |
| `F555W_comparison.png` | Side-by-side comparison: single exposure vs. combined |
| `F336W_detected.png` | F336W image with all detected sources circled in red |
| `F555W_detected.png` | F555W image with all detected sources circled in red |
| `star_catalog.csv` | Basic catalog: source IDs and positions |
| `photometry_catalog.csv` | Full catalog: IDs, positions, fluxes, magnitudes |
| `hr_diagram.png` | Hertzsprung-Russell diagram |

'''
# Conclusion

This pipeline transforms raw HST images into a scientifically valuable HR diagram through systematic cosmic ray removal, source detection, and photometry. The modular design makes it adaptable for different datasets while maintaining code clarity and reproducibility. The resulting analysis reveals stellar population characteristics that provide insights into the age, composition, and evolutionary state of the observed star cluster, demonstrating the power of automated astronomical image processing.
'''
