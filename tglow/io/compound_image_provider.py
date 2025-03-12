import logging
import time
import random
import numpy as np

from tglow.io.image_query import ImageQuery
from tglow.io.tglow_io import AICSImageReader, BlacklistReader
from tglow.utils.tglow_utils import float_to_16bit_unint_scaled, float_to_16bit_unint

from skimage import data, filters
from tqdm import tqdm

# Logging
logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class CompoundImageProvider():
    
    def __init__(self,  path, nimg,  channel, blacklist=None, plates=None, fields=None, planes=None, pseudoreplicates=0, merge_n=1, max_project=False, all_planes=False):
        self.path=path
        self.blacklist=blacklist
        self.nimg=nimg

        # Read the placklist
        if self.blacklist is not None:
            bl_reader = BlacklistReader(self.blacklist)
            bl = bl_reader.read_blacklist()
        else:
            bl = None
            
        # Subset to specific channels
        self.channel=channel
        self.plates=plates
        self.fields=fields
        self.planes=planes
        
        # Init the image reader
        self.ac_reader = AICSImageReader(self.path, self.plates, self.fields, bl)
        self.all_imgs = [x for v in self.ac_reader.images.values() for x in v]
        
        # Grab the image dimension, assuming all images are the same shape
        self.stack_dims = self.ac_reader.get_img(self.all_imgs[0]).dims
        
        # Grab the max value in the image, assuming all image are the same dtype
        self.max_value = np.iinfo(self.ac_reader.get_img(self.all_imgs[0]).dtype).max

        self.pseudoreplicates = pseudoreplicates
        self.merge_n=merge_n
        self.all_planes=all_planes
        
        self.max_project=max_project
        
        if self.max_project:
            self.all_planes=True
        
        
    def fetch_training_images(self):
        training_imgs_tmp = []
        training_imgs = []
        i=0

        start_time = time.time()
        
        if self.pseudoreplicates > 0:
            sample_n = 1
        else:
            sample_n = self.merge_n
        
        pb = tqdm(total=self.nimg, desc='Reading', unit='image')
        while i < self.nimg:
            i=i+1
            training_imgs_tmp.extend(self.fetch_compound(sample_n))
            pb.update(1)
            #log.info(f"[{i}/{self.nimg}]")
        pb.close()

        log.info(f"Reading took { round(((time.time() - start_time)/60), 2)} minutes, read {len(training_imgs_tmp)} images")
        
        if self.pseudoreplicates > 0:
            log.info(f"Pseudoreplicating {self.pseudoreplicates} compound images from {self.nimg} read images, with each compound image consisting of {self.merge_n} images")
            log.info(f"Creating {self.pseudoreplicates} pseudoreplicates from {len(training_imgs_tmp)} images")
            
            i = 0
            while i < self.pseudoreplicates:
                merge_files = random.choices(training_imgs_tmp, k=self.merge_n)
                training_imgs.append(np.max(np.array(merge_files), axis=0, keepdims=False))
                i+=1
        else:
            training_imgs = training_imgs_tmp
            
        return training_imgs
        
    def fetch_compound_image_mean_max(self, max_project=False):
        sum_image=None
        sum_mask=None
        max_image=None
        i=0
        start_time = time.time()
        total_imgs = 0
        
        if self.pseudoreplicates > 0:
            log.warning("Pseudoreplicating not supported for fetching average image, option is ignored")

        pb = tqdm(total=self.nimg, desc='Reading', unit='image')
        while i < self.nimg:
            i=i+1
            for final_img in self.fetch_compound(self.merge_n):
                #if rescale:
                #    final_img = final_img.astype(np.float32) 
                #    final_img = final_img / self.max_value
                    
                if sum_image is None:
                    sum_image = np.zeros_like(final_img, dtype=np.float32)
                    sum_mask = np.ones_like(final_img, dtype=np.float32)
                    max_image = np.zeros_like(final_img, dtype=np.float32)

                sum_image = sum_image + final_img
                
                max_image[final_img > max_image] = final_img[final_img > max_image]
                    
                #if threshold:
                #    thresh = filters.threshold_otsu(final_img)
                #    sum_mask += (final_img > thresh)
                    
                total_imgs += 1
            pb.update(1)
            #log.info(f"[{i}/{self.nimg}]")
            
        pb.close()
        #if threshold:
        #    avg_image = sum_image / sum_mask
        #else:
        avg_image = sum_image / self.nimg
        log.info(f"Reading took { round(((time.time() - start_time)/60), 2)} minutes")
        
        if max_project:  
            return(max_image)
        else:
            return(avg_image)

    def fetch_compound(self, sample_n):
        
        # Select subset of images
        merge_files = random.choices(self.all_imgs, k=sample_n)

        # Loop over the raw tiffs and read them as a list of 2d numpy arrays
        j = 0
        
        training_imgs_tmp = []
        final_output = []
        
        while j < len(merge_files): 
            q = merge_files[j]
            q.channel = self.channel
            training_imgs_tmp.append(self.fetch_image(q))
            j += 1
                        
        #log.debug(f"Read {len(training_imgs_tmp)} stacks of {training_imgs_tmp[j-1].shape} into array")
        
        if self.merge_n > 1:
            final_output.append(np.max(np.array(training_imgs_tmp), axis=0))
            # Mean projection
            #mean_tmp = np.mean(np.array(training_imgs_tmp), axis=0)
            #mean_tmp[mean_tmp < 0] = 0
            #mean_tmp = float_to_16bit_unint(mean_tmp)
            #final_output.append(mean_tmp)
        else:
            #final_output.extend(training_imgs_tmp)
            final_output = training_imgs_tmp
                    
        return(final_output)

    def fetch_image(self, q):
            
            if self.all_planes:
                img = self.ac_reader.read_image(q)
                if self.max_project:
                    img = np.max(img, axis=0)
                else:
                    return(list(img))
            else:
                if self.planes is not None:
                    rplane = random.sample(self.planes, k=1)[0]
                else:
                    rplane = random.sample(range(0, self.stack_dims["Z"][0]), k=1)[0]
                q.plane = rplane
                img = self.ac_reader.read_image(q)
                
            return(img)
        


        
    
    
    