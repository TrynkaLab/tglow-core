from cellpose import models, io
import numpy as np
import time
import logging
import argparse
from tglow.io.tglow_io import AICSImageReader
from tglow.io.image_query import ImageQuery
import os
import math

# Cellpose logger
#logger = io.logger_setup()
logging.getLogger("cellpose").setLevel(logging.INFO)

# Logging
logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class CellposeRunner():
    
    def __init__(self, path, plate, output, model, nucl_channel, other_channel, diameter, diameter_nucl, do_3d, anisotropy=None):
        
        self.path=path
        self.plate=plate
        
        self.anisotropy=anisotropy
        self.do_3d=do_3d
        self.__init_reader__()
        
        self.output=output
        self.model=model
        self.nucl_channel=nucl_channel
        self.other_channel=other_channel
        self.diameter=diameter
        self.diameter_nucl=diameter_nucl
        
        # Minimum object area in pixels
        self.min_size= math.pi * ((self.diameter/4) ** 2)
        
        if self.diameter_nucl is not None:
            self.min_size_nucl= math.pi * ((self.diameter_nucl/4) ** 2)

        
    def __init_reader__(self):
        
        self.reader = AICSImageReader(self.path,
                                      plates_filter=self.plate)
        
        if self.anisotropy is None and self.do_3d is True:
            log.info("Estimating anisotropy based on first image")
            img1 = self.reader.images[0]
            
            meta = self.reader.get_img(img1)            
            
            if meta.physical_pixel_sizes is not None:
                px = meta.physical_pixel_sizes
                log.info(f"Pixel sizes from ome {px}")
                self.anisotropy = px[0] / px[1]
                log.info(f"Estimated anisotropy {self.anisotropy}")
            else:
                log.error("Must provide --anisotropy when running in 3d or ome metadata must contain PhysicalPixelSizes")
                raise TypeError("Must provide --anisotropy when running in 3d or ome metadata must contain PhysicalPixelSizes")
        
    def run_model(self, query):
        
        
        log.debug(f"Processing query {query.to_string()}, row: {query.row},row: {query.get_row_letter()}, col: {query.col}, well: {query.get_well_id()}")
        
        
        # Read channel with cell outlines
        query.channel = self.other_channel
        data_cell = self.reader.read_image(query)
        
        # If max projection is true
        if not self.do_3d:
            data_cell = np.max(data_cell, axis=0)
            log.debug(f"cell shape: {data_cell.shape}")
            
        # Read channel with nucleus signal
        if self.nucl_channel:
            query.channel = self.nucl_channel
            data_nucl = self.reader.read_image(query)
            
            if not self.do_3d:
                data_nucl = np.max(data_nucl, axis=0)
                log.debug(f"nucl shape: {data_nucl.shape}")
        
            # Combine into stack
            img = np.stack([data_cell, data_nucl], axis=-1)
            channel_axis=3
            # define CHANNELS to run segementation on
            channels = [1,2]
            # grayscale=0, R=1, G=2, B=3
            # channels = [cytoplasm, nucleus]
            # if NUCLEUS channel does not exist, set the second channel to 0
            # channels = [0,0]
            # IF ALL YOUR IMAGES ARE THE SAME TYPE, you can give a list with 2 elements
            # channels = [0,0] # IF YOU HAVE GRAYSCALE
            # channels = [2,3] # IF YOU HAVE G=cytoplasm and B=nucleus
            # channels = [2,1] # IF YOU HAVE G=cytoplasm and R=nucleus
            # or if you have different types of channels in each image
        else:
            # In the case there is no nucleus channel
            img = data_cell
            channels = [0,0]
            channel_axis=None
            
        log.debug(f"Read images in stack of shape {img.shape}")
         
        start_time = time.time()
        masks, flows, styles, diams = model.eval(img,
                                                diameter=self.diameter,
                                                channels=channels,
                                                channel_axis=channel_axis,
                                                do_3D=self.do_3d,
                                                anisotropy=self.anisotropy,
                                                min_size=self.min_size)
        log.info("Cell mask running time %s seconds" % round(time.time() - start_time))

        if not os.path.exists(f"{self.output}/{query.plate}/{query.get_row_letter()}/{query.col}/"):
            os.makedirs(f"{self.output}/{query.plate}/{query.get_row_letter()}/{query.col}")

        # save results as png
        io.save_masks(img, masks, flows, f"{self.output}/{query.plate}/{query.get_row_letter()}/{query.col}/{query.field}_cell_mask_d{str(self.diameter)}_ch{self.other_channel}", tif=True, png=False, save_outlines=True)
        
        start_time = time.time()
        
        if (self.nucl_channel is not None):
            if (self.do_3d):
                nucl = img[:,:,:,1]
            else:
                nucl = img[:,:,1]
            
            masks, flows, styles, diams = model.eval(nucl,
                                                    diameter=self.diameter_nucl,
                                                    channels=[0,0],
                                                    do_3D=self.do_3d,
                                                    anisotropy=self.anisotropy,
                                                    min_size=self.min_size_nucl)
            log.info("Nucleus running time %s seconds" % round(time.time() - start_time))

            # save results as png
            io.save_masks(img, masks, flows, f"{self.output}/{query.plate}/{query.get_row_letter()}/{query.col}/{query.field}_nucl_mask_d{self.diameter_nucl}_ch{self.nucl_channel}", tif=True, png=False, save_outlines=True)


if __name__ == "__main__":
    
    # CLI 
    parser = argparse.ArgumentParser(description="Train a basicpy model on raw HCI images orgnaized into <plate>/<row>/<col>/<field>.ome.tiff stacks with CZYX")
    parser.add_argument('-i','--input', help='Base dir to input organized <plate>/<row>/<col>/<field>.ome.tiff', required=True)
    parser.add_argument('-p','--plate', help='Plates to process (at least one)', nargs='+', required=True)
    parser.add_argument('-w','--well', help='Wells to process (at least one)', nargs='+', required=True)
    parser.add_argument('-o','--output', help='Output folder', default="./")
    parser.add_argument('--nucl_channel', help='Channel for nucleus signal', required=False, default=None)
    parser.add_argument('--cell_channel', help='Channel for cell segmentation signal', required=True)
    parser.add_argument('--model', help='Cellpose model', default="cyto2")
    parser.add_argument('--gpu', help="Use the GPU", action='store_true', default=False)
    parser.add_argument('--anisotropy', help="Ratio between z / xy resolution", default=None)
    parser.add_argument('--diameter', help="Estimated cellsize", default=None)
    parser.add_argument('--diameter_nucl', help="Estimated nucleus size", default=None)
    parser.add_argument('--no_3d', help="Don't run in 3d mode", action='store_true', default=False)

    args = parser.parse_args()
    
    # Init cellpose model
    model = models.Cellpose(gpu=args.gpu, model_type=args.model)
    
    if not args.no_3d and args.diameter is None:
        e = "Must provide --diameter if running in 3d mode"
        log.error(e)
        raise TypeError(e)
    
    if not args.no_3d and (args.diameter_nucl is None and args.nucl_channel is not None):
        e = "Must provide --diameter_nucl if running in 3d mode and --nucl_channel is specified"
        log.error(e)
        raise TypeError(e)
    
    if args.no_3d:
        log.info("Max projecting image stacks before running cellpose")
        
    # Cellpose runner class
    runner = CellposeRunner(args.input,
                            args.plate,
                            args.output,
                            model,
                            args.nucl_channel,
                            args.cell_channel,
                            int(args.diameter) if args.diameter is not None else args.diameter,
                            int(args.diameter_nucl) if args.diameter_nucl is not None else args.diameter_nucl,
                            not args.no_3d,
                            args.anisotropy)
    
    # Loop, ideally one plate and well at the time is supplied, but can run all
    for plate in args.plate:
        for well in args.well:
            row_col = ImageQuery.well_id_to_index(well)    
            #log.debug(f"Detected row_col {row_col} for well {well} in {plate}")
            #log.debug(f"{runner.reader.index.keys()}")

            fields = runner.reader.index[plate][str(row_col[0])][str(row_col[1])].keys()
            log.info(f"Detected fields {fields} for well {well} in {plate}")
            
            for field in fields:
                q = ImageQuery(plate, row_col[0], row_col[1], field)
                log.info(f"Running for {q.to_string()}")
                runner.run_model(q)