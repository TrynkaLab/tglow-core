import logging
import tifffile
import os
import argparse
import RedLionfishDeconv as rl
import numpy as np

from tglow.io.tglow_io import AICSImageReader, BlacklistReader, AICSImageWriter
from tglow.io.image_query import ImageQuery
from tglow.utils.tglow_utils import float_to_16bit_unint, float_to_32bit_unint, dict_to_str, float_to_16bit_unint_scaled

root_log = logging.getLogger()
root_log.setLevel(logging.INFO)

# Logging
logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

class RedLionFish:
    
    def __init__(self, args):
        
        self.input = args.input
        self.plate = args.plate
        self.well = args.well
        self.fields = args.fields
        
        self.output = args.output

        self.blacklist = args.blacklist
        self.niter = int(args.niter)
        self.mode = args.mode
        self.max_project = args.max_project
        self.uint32 = args.uint32
        self.clip_max=int(args.clip_max)

        if self.blacklist is not None:
            bl_reader = BlacklistReader(self.blacklist)
            bl = bl_reader.read_blacklist()
        else:
            bl = None
        
        self.reader = AICSImageReader(self.input, plates_filter=args.plate, fields_filter=self.fields, blacklist=bl)
        self.writer = AICSImageWriter(self.output)
        
        if self.fields is None:
            self.fields = self.reader.fields[self.plate]
                    
        if args.channels is None:
            img = self.reader.get_img(self.reader.images[self.plate][0])
            self.channels = range(0, img.dims['C'][0])
        else:
            self.channels = [int(x) for x in args.channels]
        
        self.psf_string = args.psf
        
        # Read PSFs
        self.psfs = {}
        for val in args.psf:
            keypair = val.split("=")
            log.info(f"Reading PSF: {keypair}")
            self.psfs[int(keypair[0])] = tifffile.imread(keypair[1])
 
    
    def run(self):
        
        for field in self.fields:
            self.run_decon(field)
    
    def run_decon(self, field):
        
        row, col = ImageQuery.well_id_to_index(self.well)
        iq = ImageQuery(self.plate, row, col, field)
        image = self.reader.read_image(iq)
        img = self.reader.get_img(iq)

        decons = []
        for channel in self.channels:
            
            if (channel in self.psfs):  
                log.info(f"Deconvoluting field {field} channel {channel}")
                cur_decon = self.deconvolute(image, channel)
                
                log.debug(f"dtype: {cur_decon.dtype} max: {np.max(cur_decon)}")
            else:
                log.info(f"NOT deconvoluting field {field} channel {channel}")
                cur_decon = image[channel,:,:,:]
                
            if self.max_project:
                decons.append(np.max(cur_decon, axis=0, keepdims=True))
            else:
                decons.append(cur_decon)
        
        
        if self.uint32:        
            decon = float_to_32bit_unint(np.array(decons))
        else:
            decon = float_to_16bit_unint_scaled(np.array(decons), self.clip_max)
        
        log.debug(f"Final dtype: {decon.dtype} max: {np.max(decon)}")

        
        log.info(f"Writing stack of dims {decon.shape}")
        self.writer.write_stack(decon, iq, img.channel_names, img.physical_pixel_sizes)
               
    def deconvolute(self, image, channel):
    
        psf = self.psfs[channel]
        decon = rl.doRLDeconvolutionFromNpArrays(image[channel,:,:,:],
                                        psf,
                                        niter=self.niter,
                                        method=self.mode)
        
        return(decon)
        
    
    def printParams(self):
        
        log.info(f"Input:\t\t{self.input}")
        log.info(f"Input plate:\t{str(self.plate)}")
        log.info(f"Input well:\t{str(self.well)}")
        log.info(f"Input fields:\t{str(self.fields)}")

        log.info(f"Output:\t\t{self.output}")
        log.info(f"psfs:\t\t{self.psf_string}")
        log.info(f"niter:\t\t{self.niter}")
        log.info(f"mode:\t\t{self.mode}")
        log.info(f"channels:\t\t{self.channels}")
        log.info(f"max project:\t\t{self.max_project}")
        log.info(f"Output 32bit:\t{str(self.uint32)}")    
        
     

# Main loop
if __name__ == "__main__":
            
    # CLI 
    parser = argparse.ArgumentParser(description="Merge raw Phenix tiff files organized per well into per channel tiffs, one page per z-stack")
    parser.add_argument('-i','--input', help='Base dir to raw input')
    parser.add_argument('-o','--output', help='Output prefix')
    parser.add_argument('-p','--plate', help='Subfolder in raw dir to process')
    parser.add_argument('-w', '--well', help='Well ID to merge')
    parser.add_argument('--fields', help='Fields to use. <field #1> | [<field #1> <field #2> <field #n>]', nargs='+', default=None)
    parser.add_argument('--psf', help='Point spread functions <zero index channel id>=/path/to/tiff', nargs='+')
    parser.add_argument('--niter', help='Number of RL iterations', default=10)
    parser.add_argument('--mode', help='Mode of caluclating RL. cpu|gpu', default="gpu")
    parser.add_argument('--max_project', help="Output max projection over Z", action='store_true', default=False)
    parser.add_argument('--blacklist', help='TSV file with "<plate>  <well>" on each row descrbing what to ignore', default=None)
    parser.add_argument('--clip_max', help='Clip values above this prior to 32>16 bit conversion. Values below this will be preserved, but scaled and rounded to 16 bit range', default=65535)

    parser.add_argument('--channels', help='Channels to output, zero indexed. Only channels specified in --psf are deconveluted!', nargs='+', default=None)
    parser.add_argument('--uint32', help="Write as 32 bit unsigned integer instead of clipping to 16 bit uint after applying basicpy model", action='store_true', default=False)
    args = parser.parse_args()
  
    log.info("-----------------------------------------------------------")
    if not os.path.exists(args.output):
        os.makedirs(args.output)
        log.info(f"Folder created: {args.output}")
    
    
    log.debug(f"FIELDS: {args.fields}")
    runner=RedLionFish(args)
    runner.printParams()
    log.info("-----------------------------------------------------------")
    
    runner.run()