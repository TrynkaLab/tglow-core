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
