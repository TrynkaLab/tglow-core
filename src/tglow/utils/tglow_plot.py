import numpy as np
from skimage import exposure
from matplotlib import pyplot as plt
import pandas as pd

# Take greyscale images and make a composite
def composite_images(imgs, equalize=False, aggregator=np.mean):

    if equalize:
        imgs = [exposure.equalize_hist(img) for img in imgs]

    imgs = [img / img.max() for img in imgs]

    if len(imgs) < 3:
        imgs = [np.zeros(shape=imgs[0].shape)] * (3-len(imgs)) + imgs

    imgs = np.dstack(imgs)
    return imgs


# Plot before and after of an imageset for a given index in the arrays
def plot_registration_imgs(before, after, filename):
    
    #print(f"[DEBUG] b:{before.dtype} a: {after.dtype}")
    
    fig, axes = plt.subplots(1, 2, figsize=(6, 3))
    im = axes[0].imshow((before * 255).astype(np.uint8), cmap='gray', vmin=0, vmax=1)
    axes[0].set_title("Before registration")
    im = axes[1].imshow((after * 255).astype(np.uint8), cmap='gray', vmin=0, vmax=1)
    axes[1].set_title("After registration")
    
    fig.tight_layout()
    fig.savefig(filename, dpi=600)
    plt.close()
    
    
def plot_grey_as_magma(image, filename):
    
    # Create the figure and axis
    fig, ax = plt.subplots()

    # Create the image plot
    img = ax.imshow(image, cmap='magma')

    # Add colorbar
    cbar = plt.colorbar(img)

    # Optional: Add labels and title
    ax.set_title("Beads")
    cbar.set_label("Values")

    # Save the figure
    fig.tight_layout()
    fig.savefig(filename, dpi=600)
    plt.close()


def plot_blobs(image, filename, blobs):
    
    fig, ax = plt.subplots()

    ax.set_title("Identified blobs")
    ax.imshow(image)
    for blob in blobs:
        y, x, r = blob
        c = plt.Circle((x, y), r, color="yellow", linewidth=2, fill=False)
        ax.add_patch(c)
    ax.set_axis_off()

    fig.tight_layout()
    fig.savefig(filename, dpi=600)
    plt.close()


def plot_histogram(data, filename, bins=30, title="Histogram", xlabel="Value", ylabel="Frequency"):
    """
    Create and save a histogram to a file.
    
    Parameters:
    - data: array-like, the data to plot
    - filename: str, the output filename (including path if needed)
    - bins: int or sequence, number of bins or bin edges (default: 10)
    - title: str, title of the histogram (default: "Histogram")
    - xlabel: str, label for x-axis (default: "Value")
    - ylabel: str, label for y-axis (default: "Frequency")
    """
    
    # Create the histogram
    plt.figure(figsize=(10, 6))
    plt.hist(data, bins=bins, edgecolor='black')
    
    # Set labels and title
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    
    # Add grid
    #plt.grid(True, linestyle='--', alpha=0.7)
    
    # Tight layout to prevent clipping of labels
    plt.tight_layout()
    
    # Save the plot to a file
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    
    # Close the plot to free up memory
    plt.close()
    
    


def plot_histogram_df(data, filename, bins=30, title_prefix="Histogram", xlabel="Value", ylabel="Frequency", ncols=3):
    """
    Create and save histograms for each column in a pandas DataFrame with mean and median lines.
    
    Parameters:
    - data: pandas DataFrame, columns will be plotted individually
    - filename: str, base filename (e.g., 'histograms' → 'histograms_col1.png', etc.)
    - bins: int or sequence, number of bins or bin edges (default: 30)
    - title_prefix: str, prefix for each histogram title
    - xlabel: str, label for x-axis
    - ylabel: str, label for y-axis
    - ncols: int, number of columns in subplot grid (default: 3)
    """
    import numpy as np
    import matplotlib.pyplot as plt
    
    # Select only numeric columns
    numeric_cols = data.select_dtypes(include='number').columns
    
    if len(numeric_cols) == 0:
        print("No numeric columns found in DataFrame")
        return
    
    nrows = int(np.ceil(len(numeric_cols) / ncols))
    
    # Create figure with subplots
    fig, axes = plt.subplots(nrows, ncols, figsize=(5*ncols, 4*nrows))
    axes = axes.flatten() if nrows > 1 or ncols > 1 else [axes]
    
    # Plot histogram for each numeric column
    for i, col in enumerate(numeric_cols):
        col_data = data[col].dropna()
        axes[i].hist(col_data, bins=bins, edgecolor='black', alpha=0.7)
        
        # Calculate mean and median
        mean_val = col_data.mean()
        median_val = col_data.median()
        
        # Add vertical lines for mean and median
        axes[i].axvline(mean_val, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_val:.2f}')
        axes[i].axvline(median_val, color='blue', linestyle='--', linewidth=2, label=f'Median: {median_val:.2f}')
        
        axes[i].set_title(f"{title_prefix}: {col}")
        axes[i].set_xlabel(xlabel)
        axes[i].set_ylabel(ylabel)
        axes[i].legend()
        #axes[i].grid(True, linestyle='--', alpha=0.7)
    
    # Hide unused subplots
    for i in range(len(numeric_cols), len(axes)):
        axes[i].set_visible(False)
    
    plt.tight_layout()
    
    # Save single combined figure
    plt.savefig(f"{filename}.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{filename}.pdf", format='pdf', bbox_inches='tight', dpi=300)
    plt.close()
