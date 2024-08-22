import numpy as np
from skimage import exposure
from matplotlib import pyplot as plt

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