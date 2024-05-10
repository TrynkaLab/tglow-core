import os
import re
import tifffile
import string
import argparse
import random
import numpy as np
import logging
from basicpy import BaSiC
from matplotlib import pyplot as plt
from scipy import ndimage as ndi
from tglow.io.tglow_io import AICSImageReader
from tglow.utils.tglow_utils import float_to_16bit_unint


# Plot results from basicpy fit
def plot_basic_results(basic, filename):
    fig, axes = plt.subplots(1, 3, figsize=(9, 3))
    im = axes[0].imshow(basic.flatfield)
    fig.colorbar(im, ax=axes[0])
    axes[0].set_title("Flatfield")
    im = axes[1].imshow(basic.darkfield)
    fig.colorbar(im, ax=axes[1])
    axes[1].set_title("Darkfield")
    axes[2].plot(basic.baseline)
    axes[2].set_xlabel("Frame")
    axes[2].set_ylabel("Baseline")
    fig.tight_layout()
    fig.savefig(filename, dpi=300)
    plt.close()

# Plot before and after of an imageset for a given index in the arrays
def plot_before_after(images, images_transformed, i, filename):
    fig, axes = plt.subplots(1, 2, figsize=(6, 3))
    im = axes[0].imshow(images[i])
    fig.colorbar(im, ax=axes[0])
    axes[0].set_title("Original")
    im = axes[1].imshow(images_transformed[i])
    fig.colorbar(im, ax=axes[1])
    axes[1].set_title("Corrected")
    fig.suptitle(f"frame {i}")
    fig.tight_layout()
    fig.savefig(filename, dpi=300)
    plt.close()

# Plot before and after for two pre-picked images, usefull for max projections
def plot_before_after_mp(images,images_transformed, filename):
    fig, axes = plt.subplots(1, 2, figsize=(6, 3))
    im = axes[0].imshow(images)
    fig.colorbar(im, ax=axes[0])
    axes[0].set_title("Original")
    im = axes[1].imshow(images_transformed)
    fig.colorbar(im, ax=axes[1])
    axes[1].set_title("Corrected")
    fig.suptitle(f"Max projected")
    fig.tight_layout()
    fig.savefig(filename, dpi=300)
    plt.close()


# Logging
logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

class BasicpyTrainer():
    
    def __init__(self, path, output_dir, output_prefix,  channel,  nimg, merge_n, tune, fit_darkfield, all_planes=False, max_project=False, plates=None, fields=None, planes=None):
        
        self.path=path
        self.channel=channel
        self.plates=plates
        self.fields=fields
        self.planes=planes
        
        self.nimg=nimg
        self.merge_n=merge_n
        self.all_planes=all_planes
        self.tune=tune
        self.fit_darkfield=fit_darkfield
        
        self.output_dir=output_dir
        self.output_prefix=output_prefix
        
        self.max_project=max_project
        
        if self.max_project:
            self.all_planes=True
            
    def train(self):
        
        # Init the image reader
        ac_reader = AICSImageReader(self.path, self.plates, self.fields)
        training_imgs = []
        i=0

        while i < self.merge_n:
            i=i+1
            # Select subset of images
            merge_files = random.choices(ac_reader.images, k=self.nimg)
            
            stack_dims = ac_reader.get_img(ac_reader.images[0]).dims
            
            log.info(f"Reading {len(merge_files)} randomly selected files")

            # Loop over the raw tiffs and read them as a list of 2d numpy arrays
            j = 0
            training_imgs_tmp = []
            while j < len(merge_files): 
                #print(merge_files[j])
                #training_imgs_tmp.append(tifffile.imread(merge_files[j]))
                q = merge_files[j]
                q.channel = self.channel
                
                if self.all_planes:
                    img = ac_reader.read_image(q)
                    if self.max_project:
                        training_imgs_tmp.append(np.max(img, axis=0))
                    else:
                        training_imgs_tmp.append(list(img))
                else:
                    if self.planes is not None:
                        rplane = random.sample(self.planes, k=1)[0]
                    else:
                        rplane = random.sample(range(0, stack_dims["Z"][0]), k=1)[0]
                    q.plane = rplane
                    img = ac_reader.read_image(q)
                    training_imgs_tmp.append(img)
                log.debug(f"Read stack of {img.shape}")
                j=j+1
                
            if self.merge_n > 1:
                log.info("Merging " + str(self.nimg) + " random images for training")
                training_imgs.append(np.max(np.array(training_imgs_tmp), axis=0))
            else:
                # This gave issues as it puts an array of arrays, giving the wrong shape
                #training_imgs.append(training_imgs_tmp)
                training_imgs=training_imgs_tmp
                
        # Convert into 3d numpy array
        merged = np.array(training_imgs)
        log.info(f"Read training files into array of shape {str(merged.shape)}")

        # Init basicpy object with autosegmentation (needed for 3d)
        basic = BaSiC(get_darkfield=self.fit_darkfield, smoothness_flatfield=1, autosegment=False)

        # Optimze parameters
        if self.tune:
            log.info("Tuning model")
            basic.autotune(merged)
        
        # Fit parameters
        log.info("Fitting model")
        basic.fit(merged)
        
        # Save output
        out = f"{self.output_dir}/{self.output_prefix}_ch{str(self.channel)}"
        
        if not os.path.exists(out):
            os.makedirs(out)
            
        basic.save_model(out, overwrite=True)
        
        # Apply correction
        merged_corrected = basic.transform(merged)
        
        # Convert data back to original 16bit uint
        merged_corrected = float_to_16bit_unint(merged_corrected)
        
        # Plot results
        plot_basic_results(basic, out + "/flat_and_darkfield.png")
        plot_before_after_mp(np.max(merged, axis=0), np.max(merged_corrected, axis=0), filename=out + "/all_imgs_max_proj_pre_post.png")
        plot_before_after(merged, merged_corrected, 0, filename=out + "/img0_pre_post.png")
        plot_before_after(merged, merged_corrected, 1, filename=out + "/img1_pre_post.png")
                

if __name__ == "__main__":
    
    # CLI 
    parser = argparse.ArgumentParser(description="Train a basicpy model on raw HCI images orgnaized into <plate>/<row>/<col>/<field>.ome.tiff stacks with CZYX")
    parser.add_argument('-i','--input', help='Base dir to input organized <plate>/<row>/<col>/<field>.ome.tiff')
    parser.add_argument('-o','--output', help='Output folder')
    parser.add_argument('--output_prefix', help='Output prefix name', default="basicpy")
    parser.add_argument('--no_tune', help="Do not tune the basicpy model", action='store_true', default=False)
    parser.add_argument('--fit_darkfield', help="Fit the darkfield component. For tglow not reccomended", action='store_true', default=False)
    parser.add_argument('--nimg', help="Number of random images to train on. If --merge_n is specified this is the number of images that are max projected into one. Sampled with replacement", default=None, required=True)
    parser.add_argument('--merge_n', help="Number of times to max project random selection --nimg number of images. If > 1 basicpy is run on this number of max projected images.", default=1)
    parser.add_argument('-p','--plate', help='Plate to process', nargs='+')
    parser.add_argument('-c','--channel', help="Channel number to correct", required=True)
    parser.add_argument('--fields', help='Fields to use', nargs='+', default=None)
    parser.add_argument('--planes', help='Z planes to use', nargs='+', default=None)
    #parser.add_argument('--gpu', help="Use the GPU", action='store_true', default=False)
    parser.add_argument('--all_planes', help="Instead of randomly picking one plane for a stack, use them all", action='store_true', default=False)
    parser.add_argument('--max_project', help="Calculate flatfields on max projections. Automatically activates --all_planes", action='store_true', default=False)
    args = parser.parse_args()
    
    input=args.input
        
    # Set subfolders to process
    if args.plate == None:
        plates=[ name for name in os.listdir(input) if os.path.isdir(os.path.join(input, name)) ]
    else:
        plates=args.plate
    
    # If max projecting all_planes also needs to be true so the full stack is read
    # This is also forced in the constructor, but this is so its printed well here
    # TODO: make the printing a method
    if args.max_project:
        args.all_planes=True
    
    print("-----------------------------------------------------------")
    print("Input plates:\t" + str(plates))
    print("Input fields:\t" + str(args.fields))
    print("Input planes:\t" + str(args.planes))    
    print("Input:\t\t" + input)
    print("Output:\t\t" + args.output)
    print("Output pre:\t" + args.output_prefix)
    print("channel:\t" + str(args.channel))
    print("n img:\t\t" + str(args.nimg))
    print("tune basicpy:\t" + str(not args.no_tune))
    print("merge n:\t" + str(args.merge_n))
    print("darkfied:\t" + str(args.fit_darkfield))
    print("max project:\t" + str(args.max_project))
    print("-----------------------------------------------------------")

    trainer = BasicpyTrainer(path=input,
                            output_dir=args.output,
                            output_prefix=args.output_prefix,
                            channel=args.channel,
                            nimg=int(args.nimg),
                            tune=not args.no_tune,
                            merge_n=int(args.merge_n),
                            fit_darkfield=args.fit_darkfield,
                            all_planes=args.all_planes,
                            max_project=args.max_project,
                            plates=plates,
                            fields=args.fields)
                            
    trainer.train()






    
    
    
    