# Cosmic Ray Removal
# Star-Finding
# Photometry
# Creating a Hertzsprung–Russell Diagram

import numpy as np
import pandas as pd
from astropy.io import fits
from pathlib import Path
from scipy import ndimage
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt

# Step 1: Cosmic Ray Removal

def median_combine_fits(directory, output_path=None, show_plot=False):
    # Median-combine all FITS files in a directory to remove cosmic rays.

    directory = Path(directory)

    fits_files = sorted(list(directory.glob("*.fits")) + list(directory.glob("*.fit")))

    # Load all images
    images = []
    header = None

    for fits_file in fits_files:
        with fits.open(fits_file) as hdul:
            # Find the first HDU with data
            data = None
            for hdu in hdul:
                if hdu.data is not None:
                    data = hdu.data
                    if header is None:
                        header = hdu.header.copy()
                    break

            if data is not None:
                images.append(data)
            else:
                print(f"Warning: No data found in {fits_file.name}")

    if len(images) == 0:
        raise ValueError(f"No valid image data found in {directory}")

    # Median combine
    combined_image = np.median(images, axis=0)
    print(f"Median-combined {len(images)} images")

    # Save The Result
    if output_path is not None:
        hdu = fits.PrimaryHDU(data=combined_image, header=header)
        hdu.writeto(output_path, overwrite=True)
        print(f"Saved to {output_path}")

    # Show the Plot for both filters
    if show_plot:

        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        # Smooth the images
        smoothed_single = gaussian_filter(images[0], sigma=1)
        smoothed_combined = gaussian_filter(combined_image, sigma=1)

        vmin, vmax = np.percentile(smoothed_single, [1, 99])
        axes[0].imshow(smoothed_single, cmap='gray', vmin=vmin, vmax=vmax, origin='lower', interpolation='bilinear')
        axes[0].set_title('Single Exposure (with cosmic rays)')
        axes[0].axis('off')

        vmin, vmax = np.percentile(smoothed_combined, [1, 99])
        axes[1].imshow(smoothed_combined, cmap='gray', vmin=vmin, vmax=vmax, origin='lower', interpolation='bilinear')
        axes[1].set_title('Median Combined (cosmic rays removed)')
        axes[1].axis('off')

        plt.tight_layout()
        plot_path = str(output_path).replace('.fits', '_comparison.png')
        plt.savefig(plot_path, dpi=150, bbox_inches='tight')
        print(f"Saved plot to {plot_path}")
        plt.close()

    return combined_image, header

# Step 2: Star-Finding

def find_stars(image, threshold_sigma=0.7, min_separation=3):
    # Finding stellar sources in an image, using threshold and peak detection.

    # To Calculate threshold
    background = np.median(image)
    threshold = background + threshold_sigma * np.std(image)

    # Masking the image
    mask = image > threshold
    labeled, num = ndimage.label(mask)

    # Get peak position in each region
    source = []
    border = 40  # pixels to exclude from edges(It is noisy)
    ny, nx = image.shape

    for i in range(1, num + 1):
        region = labeled == i
        y_coords, x_coords = np.where(region)
        values = image[region]
        max_idx = np.argmax(values)

        x = x_coords[max_idx]
        y = y_coords[max_idx]

        # Skip stars near LEFT and BOTTOM edges only (Used AI for this part)
        if x < border or y < border:
            continue

        # Check minimum separation
        if source:
            xs = np.array([s['x_center'] for s in source])
            ys = np.array([s['y_center'] for s in source])
            dists = np.sqrt((x - xs)**2 + (y - ys)**2)
            if np.all(dists >= min_separation):
                source.append({'id': len(source) + 1, 'x_center': x, 'y_center': y})
        else:
            source.append({'id': 1, 'x_center': x, 'y_center': y})

    return pd.DataFrame(source)


def match_catalogs(cat1, cat2, max_distance=3.0): # Used AI for this part
    # Match sources between two catalogs by position.

    matched = []
    used = set()

    for _, row1 in cat1.iterrows():
        best_idx = None
        best_dist = max_distance

        for idx2, row2 in cat2.iterrows():
            if idx2 in used:
                continue

            dist = np.sqrt((row1['x_center'] - row2['x_center'])**2 + (row1['y_center'] - row2['y_center'])**2)

            if dist < best_dist:
                best_dist = dist
                best_idx = idx2

        if best_idx is not None:
            matched.append({
                'id': len(matched) + 1,
                'x_center': (row1['x_center'] + cat2.loc[best_idx, 'x_center']) / 2,
                'y_center': (row1['y_center'] + cat2.loc[best_idx, 'y_center']) / 2
            })
            used.add(best_idx)

    return pd.DataFrame(matched)

# STEP 3: Photometry

def circular_aperture(image, x, y, radius):
    # Extract flux within a circular aperture.

    ny, nx = image.shape
    yy, xx = np.ogrid[:ny, :nx]
    mask = (xx - x)**2 + (yy - y)**2 <= radius**2
    return np.sum(image[mask])


def annulus_background(image, x, y, inner_radius, outer_radius):
    # Estimate background using an annulus.

    ny, nx = image.shape
    yy, xx = np.ogrid[:ny, :nx]
    distance = np.sqrt((xx - x)**2 + (yy - y)**2)
    annulus_mask = (distance >= inner_radius) & (distance <= outer_radius)

    if np.sum(annulus_mask) == 0:
        return 0

    return np.median(image[annulus_mask])


def aperture_photometry(image, x, y, aperture_radius=5, inner_bg=8, outer_bg=12):
    # Perform aperture photometry with local background subtraction.

    source_flux = circular_aperture(image, x, y, aperture_radius)
    bg_per_pixel = annulus_background(image, x, y, inner_bg, outer_bg)
    n_pixels = np.pi * aperture_radius**2
    flux = source_flux - (bg_per_pixel * n_pixels)
    return flux


def flux_to_magnitude(flux, zeropoint = -21.1):
    # Convert flux to magnitude.

    flux = np.array(flux)
    mag = -2.5 * np.log10(flux) + abs(zeropoint)
    return mag


def photometry_catalog(image_f336w, image_f555w, catalog, aperture_radius=5):
    # Perform photometry on all sources in catalog for both filters.

    flux_f336w = []
    flux_f555w = []

    for _, source in catalog.iterrows():
        x = source['x_center']
        y = source['y_center']

        f336w_flux = aperture_photometry(image_f336w, x, y, aperture_radius)
        f555w_flux = aperture_photometry(image_f555w, x, y, aperture_radius)

        flux_f336w.append(f336w_flux)
        flux_f555w.append(f555w_flux)

    catalog['flux_f336w'] = flux_f336w
    catalog['flux_f555w'] = flux_f555w
    catalog['mag_f336w'] = flux_to_magnitude(catalog['flux_f336w'])
    catalog['mag_f555w'] = flux_to_magnitude(catalog['flux_f555w'])
    catalog['aperture_radius'] = aperture_radius

    return catalog

# STEP 4: HR Diagram

def create_hr_diagram(catalog, output_path="hr_diagram.png"):
    # Create a Hertzsprung-Russell diagram.

    # Calculate color
    color = catalog['mag_f336w'] - catalog['mag_f555w']
    magnitude = catalog['mag_f336w']

    # Remove invalid values
    valid = np.isfinite(color) & np.isfinite(magnitude)
    color = color[valid]
    magnitude = magnitude[valid]

    print(f"Plotting {len(color)} stars in HR diagram")

    # Create plot
    plt.figure(figsize=(10, 8))
    plt.scatter(color, magnitude, s=7, alpha=0.6, edgecolors='none', c='black')

    plt.xlabel('Color (F336W - F555W)', fontsize=12)
    plt.ylabel('F336W Magnitude', fontsize=12)
    plt.title('Hertzsprung-Russell Diagram', fontsize=14)

    plt.gca().invert_yaxis()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Saved HR diagram to {output_path}")
    plt.close()

# Calling main pipeline

def main():
    # Run the complete HST analysis pipeline.

    print("HST IMAGE ANALYSIS PIPELINE")
    print("_"*50)

    # Configuration
    data_dir = "Data"
    f336w_dir = f"{data_dir}/F336W"
    f555w_dir = f"{data_dir}/F555W"

    # Step 1: Cosmic Ray Removal
    print("TASK 1: COSMIC RAY REMOVAL")
    print("_"*50)

    f336w_combined, _ = median_combine_fits(f336w_dir, f"{data_dir}/F336W_combined.fits", show_plot=True)
    f555w_combined, _ = median_combine_fits(f555w_dir, f"{data_dir}/F555W_combined.fits", show_plot=True)

    # Step 2: Star-Finding
    print("TASK 2: STAR-FINDING")
    print("_"*50)

    print("\nFinding stars in F336W...")
    cat_f336w = find_stars(f336w_combined)
    print(f"Found {len(cat_f336w)} stars")

    print("\nFinding stars in F555W...")
    cat_f555w = find_stars(f555w_combined, threshold_sigma = 0.5) # need lower threshold for this f555w
    print(f"Found {len(cat_f555w)} stars")

    print("\nMatching catalogs...")
    catalog = match_catalogs(cat_f336w, cat_f555w)
    print(f"Matched {len(catalog)} stars")

    # F336W visualization
    fig, ax = plt.subplots(figsize=(12, 12))
    smoothed = gaussian_filter(f336w_combined, sigma=1)
    vmin, vmax = np.percentile(smoothed, [1, 99])
    ax.imshow(smoothed, cmap='gray', vmin=vmin, vmax=vmax, origin='lower', interpolation='bilinear')
    for _, star in cat_f336w.iterrows():
        ax.add_patch(plt.Circle((star['x_center'], star['y_center']), 8, color='red', fill=False, linewidth=1, alpha=0.7))
    ax.set_title(f'Detected Stars ({len(cat_f336w)} sources)', fontsize=14, color='white', pad=10)
    ax.axis('off')
    plt.tight_layout()
    plt.savefig(f"{data_dir}/F336W_detected.png", dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()

    # F555W visualization
    fig, ax = plt.subplots(figsize=(12, 12))
    smoothed = gaussian_filter(f555w_combined, sigma=1)
    vmin, vmax = np.percentile(smoothed, [1, 99])
    ax.imshow(smoothed, cmap='gray', vmin=vmin, vmax=vmax, origin='lower', interpolation='bilinear')
    for _, star in cat_f555w.iterrows():
        ax.add_patch(plt.Circle((star['x_center'], star['y_center']), 8, color='red', fill=False, linewidth=1, alpha=0.7))
    ax.set_title(f'Detected Stars ({len(cat_f555w)} sources)', fontsize=14, color='white', pad=10)
    ax.axis('off')
    plt.tight_layout()
    plt.savefig(f"{data_dir}/F555W_detected.png", dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()

    catalog.to_csv(f"{data_dir}/star_catalog.csv", index=False)
    print(f"Saved star catalog")

    # Step 3: Photometry
    print("PHOTOMETRY")
    print("_"*50)

    catalog_phot = photometry_catalog(f336w_combined, f555w_combined, catalog)
    catalog_phot.to_csv(f"{data_dir}/photometry_catalog.csv", index=False)
    print(f"Saved photometry catalog with {len(catalog_phot)} sources")

    # Step 4: HR Diagram
    print("HR DIAGRAM")
    print("_"*50)

    create_hr_diagram(catalog_phot, f"{data_dir}/hr_diagram.png")

    print("ANALYSIS COMPLETE!")
    print("_"*50)
    print(f"\nOutput files saved to '{data_dir}/':")
    print(" - F336W_combined.fits")
    print(" - F336W_detected.png")
    print(" - F555W_combined.fits")
    print(" - F555W_detected.png")
    print(" - star_catalog.csv")
    print(" - photometry_catalog.csv")
    print(" - hr_diagram.png")


if __name__ == "__main__":
    main()
