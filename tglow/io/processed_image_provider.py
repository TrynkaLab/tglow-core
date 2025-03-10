import logging
import copy
import pandas as pd
import numpy as np
import os
import pickle

from tglow.io.image_query import ImageQuery
from tglow.io.tglow_io import AICSImageReader, BlacklistReader
from tglow.utils.tglow_utils import apply_registration, float_to_32bit_unint, float_to_16bit_unint
from basicpy import BaSiC

# Logging
logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class ProcessedImageProvider():
        
    def __init__(self, path, plate, blacklist=None, plate_merge=None, registration_dir=None, flatfields=None, scaling_factors=None, mask_channels=None, mask_dir=None, mask_pattern=None, uint32=False, verbose=True):
        
        # Main dir to index
        self.path = path
        
        # Blacklist with wells to ignore
        if blacklist is not None:
            self.blacklist = BlacklistReader(blacklist).read_blacklist()
        
        # Plates to process
        if len(plate) != 1:
            raise RuntimeError("Can only provide one index plate at the time with the ProccessedImageProvider")
        
        self.plates=plate

        # Save results as uint32 instead of uint16
        self.uint32=uint32
        
        self.verbose=verbose
        
        #---------------------------------------------------------------------
        # Channels for alignment and plates to merge
        self.plates_merge=plate_merge   
        self.registration_dir=registration_dir 
        
        if (self.plates_merge is None and self.registration_dir is not None):
            raise RuntimeError("Need to supply plates for registration when registration dir is set")
        
        if self.plates_merge is not None:
            # Reader for plates
            self.plate_reader=AICSImageReader(self.path, plates_filter=self.plates + self.plates_merge, blacklist=blacklist)
        
            if (self.registration_dir is None):
                raise RuntimeError("Need to supply registration_dir for registration")
            
        else:
            # Reader for plates
            self.plate_reader=AICSImageReader(self.path, plates_filter=self.plates, blacklist=blacklist)
            
        #---------------------------------------------------------------------
        # Build dict with pretrained flatfields models channel is the key
        
        if self.plates_merge is None:
            self.plates_merge = []
        
        if flatfields is None:
            self.flatfields=None
        else:
            self.flatfields={}
            for val in flatfields:
                keypair = val.split("=")

                for curp in self.plates + self.plates_merge:
                    if curp in keypair[0]:
                        log.info(f"Adding basicpy model: {keypair}")
                        self.flatfields[keypair[0]] = BaSiC.load_model(keypair[1])
                    
        #--------------------------------------------------------------------- 
        # Build dict with scaling factors
        if scaling_factors == None:
            self.scaling_factors=None
        else:
            self.scaling_factors={}
            for val in scaling_factors:
                keypair = val.split("=")
                
                for curp in self.plates + self.plates_merge:
                    if curp in keypair[0]:
                        log.info(f"Adding scaling factors: {keypair}")
                        self.scaling_factors[keypair[0]] = float(keypair[1])

        if len(self.plates_merge) == 0:
            self.plates_merge = None
        #---------------------------------------------------------------------
        # Optional demultiplexing
        self.mask_channels=None
        if mask_channels is not None:
            
            # Old way, not per plate
            #self.mask_channels = [int(c) for c in mask_channels]
            self.mask_channels = {}
            
            for val in mask_channels:
                keypair = val.split("=")
                log.info(f"Adding channel mask {keypair}")
                if keypair[0] not in self.mask_channels:
                    self.mask_channels[keypair[0]]=set()
                self.mask_channels[keypair[0]].add(int(keypair[1]))    
                
            if mask_dir is not None:
                self.mask_reader = AICSImageReader(mask_dir,
                                                self.plates,
                                                pattern=mask_pattern)
            else:
                raise RuntimeError("Must provide --mask_channels and --mask_dir")
        elif mask_dir is not None:
            raise RuntimeError("Must provide --mask_channels and --mask_dir")   
    
        #----------------------------------------------------------------------
        # Build a pandas dataframe with what the output will look like
        self.build_channel_index()
    
    # Build an index of what the new image channels will look like given these settings    
    def build_channel_index(self):
    
        # Peak at an image for channel dims and names
        img = self.plate_reader.get_img(self.plate_reader.images[self.plates[0]][0])      
        df = pd.DataFrame(columns=["ref_plate", "plate", "cycle", "channel", "name", "orig_channel", "orig_name"])
        
        cycle = 1
        channel_id = 0
        for channel in range(0, img.dims['C'][0]):
            df.at[channel, "ref_plate"] = self.plates[0]
            df.at[channel, "plate"] = self.plates[0]
            df.at[channel, "channel"] = channel
            df.at[channel, "name"] = img.channel_names[channel]
            df.at[channel, "cycle"] = cycle
            df.at[channel, "orig_channel"] = channel
            df.at[channel, "orig_name"] = img.channel_names[channel]

            channel_id += 1  
        
        # Add the channel names and IDs in the merged plate
        if self.plates_merge is not None:
            for plate_merge in self.plates_merge:
                cycle += 1
                img = self.plate_reader.get_img(self.plate_reader.images[plate_merge][0])
                
                for channel in range(0, img.dims['C'][0]):
                    df.at[channel_id, "ref_plate"] = self.plates[0]
                    df.at[channel_id, "plate"] = plate_merge
                    df.at[channel_id, "channel"] = channel_id
                    df.at[channel_id, "name"] = f"ch{channel_id} - {img.channel_names[channel]}"
                    df.at[channel_id, "cycle"] = cycle
                    df.at[channel_id, "orig_channel"] = channel
                    df.at[channel_id, "orig_name"] = img.channel_names[channel]
                    channel_id += 1
                    
        # If there are channel to mask
        if self.mask_channels is not None:
            if self.plates[0] in self.mask_channels:
                for mask_channel in self.mask_channels[self.plates[0]]:
                        mask_channel = int(mask_channel)
                        df.at[channel_id, "ref_plate"] = df.iloc[mask_channel]["ref_plate"]
                        df.at[channel_id, "plate"] = df.iloc[mask_channel]["plate"]
                        df.at[channel_id, "channel"] = channel_id
                        df.at[channel_id, "name"] = f"ch{channel_id} - {df.iloc[mask_channel]['name']} mask inclusive"
                        df.at[channel_id, "cycle"] = df.iloc[mask_channel]["cycle"]
                        df.at[channel_id, "orig_channel"] = df.iloc[mask_channel]["channel"]
                        df.at[channel_id, "orig_name"] = df.iloc[mask_channel]["name"]
                        channel_id += 1
                    
                        df.at[channel_id, "ref_plate"] = df.iloc[mask_channel]["ref_plate"]
                        df.at[channel_id, "plate"] = df.iloc[mask_channel]["plate"]
                        df.at[channel_id, "channel"] = channel_id
                        df.at[channel_id, "name"] = f"ch{channel_id} - {df.iloc[mask_channel]['name']} mask exclusive"
                        df.at[channel_id, "cycle"] = df.iloc[mask_channel]["cycle"]
                        df.at[channel_id, "orig_channel"] = df.iloc[mask_channel]["channel"]
                        df.at[channel_id, "orig_name"] = df.iloc[mask_channel]["name"]
                        channel_id += 1

        df.index = df["plate"] + "_ch" + str(df["orig_channel"])
        self.channel_index = df
        self.dims = {}
        
        for k in list(img.dims.order):
            self.dims[k] = img.dims[k][0]

        self.dims['C'] = len(self.channel_index)
   
            
    # Fetch a single image 
    def fetch_image(self, iq):
                
        # Pre allocate a nupy array to avoid creating copies (quite a big speedup)
        # Updated to init at float32 straight away, and only coverting back at the end
        stack = np.zeros((self.dims['C'], self.dims['Z'], self.dims['Y'], self.dims['X']), dtype=np.float32)
        
        # Read data in
        img = self.plate_reader.get_img(iq)
        stack[range(0, img.dims['C'][0]),:,:,:] = self.plate_reader.read_image(iq)
        
        # Keep track of the last channel populated
        last_channel = img.dims['C'][0]
        if self.verbose: log.info(f"Read up to channel {last_channel} into image stack of shape {stack.shape}, {stack.dtype}")
        #if self.verbose: log.debug(f"Stack dtype {stack.dtype}")

        #-------------------------------------------------
        # If plates need to be merged, and or registred do that here
        if self.plates_merge is not None:
            for m_plate in self.plates_merge:
                m_iq = copy.deepcopy(iq)
                m_iq.plate = m_plate
                m_dims = self.plate_reader.get_img(m_iq).dims
                
                cur_range = range(last_channel, last_channel + m_dims['C'][0])
                stack[cur_range,:,:,:] = self.plate_reader.read_image(m_iq)
                
                if self.registration_dir is not None:
                    if self.verbose: log.info(f"Applying registration")
                    reg = self.fetch_registration(iq)
                    #if self.verbose: log.info(f"beep")
                    #if self.verbose: log.debug(f"Before reg min/max /{np.max(stack, axis=(1,2,3))}")
                    stack[cur_range,:,:,:] = apply_registration(stack[cur_range,:,:,:], reg[m_plate])
                    #if self.verbose: log.debug(f"After reg  min/max {np.max(stack, axis=(1,2,3))}")
                    #stack[cur_range,:,:,:] = tmp
                
                last_channel = last_channel +  m_dims['C'][0]
            if self.verbose: log.info(f"Read up to channel {last_channel} into image stack of shape {stack.shape}, {stack.dtype}")

        
        #-------------------------------------------------
        # Apply flatfield corrections
        if self.flatfields is not None:
            
            for index, row in self.channel_index.iterrows():
                channel = int(row["channel"])
                bp_key = f"{row['plate']}_ch{row['orig_channel']}"
                
                if str(bp_key) in self.flatfields.keys():
                    if self.verbose: log.info(f"Applying flatfield {bp_key} to: {iq.plate}/{iq.get_well_id()}/ch{str(row['channel'])}/f{iq.field}")                    
                    basic_model = self.flatfields[str(bp_key)]

                    # modify in place
                    stack[channel,:,:,:] -= basic_model.darkfield
                    stack[channel,:,:,:] /= basic_model.flatfield
                    
                    #if self.uint32:
                    #    stack[channel,:,:,:]=float_to_32bit_unint(basic_model.transform(stack[channel,:,:,:]))
                    #else:
                    #    stack[channel,:,:,:]=float_to_16bit_unint(basic_model.transform(stack[channel,:,:,:]))
                                                    
            if self.verbose: log.info(f"Applied flatfieds to stack of shape {stack.shape}, {stack.dtype}")
        
        #-------------------------------------------------
        # Apply masks to the image to demultiplex
        if self.mask_channels is not None:
            if iq.plate in self.mask_channels:
                for mask_channel in self.mask_channels[iq.plate]:
                    
                    # Read mask and convert to binary
                    mask = self.mask_reader.read_image(iq)
                    if self.verbose: log.debug(f"Read mask of shape {mask.shape} with {np.max(mask)} objects and dtype {mask.dtype}")
                    
                    if (mask.shape[1] != stack.shape[1]):
                        raise RuntimeError(f"Mask and image must have the same dimensions, mask: {mask.shape}, image: {stack.shape}")
                    
                    # You can actually do operations with boolean directly in numpy
                    mask = mask > 0
                    #mask[mask > 0] = 1
                    
                    #if self.verbose: log.debug(f"Binarized mask of shape {mask.shape} with {np.max(mask)} max value and dtype {mask.dtype}")
                    #mask_inv = np.abs(1-mask)
                    
                    if self.verbose: log.debug(f"Masking channel {mask_channel}")
                    #stack[last_channel, :,:,:] = stack[mask_channel,:,:,:] * mask
                    stack[last_channel, :,:,:] = stack[mask_channel,:,:,:]
                    stack[last_channel, :,:,:] *= mask[0,:,:,:]
                    last_channel = last_channel + 1

                    #stack[last_channel, :,:,:] = stack[mask_channel,:,:,:] * ~mask
                    stack[last_channel, :,:,:] = stack[mask_channel,:,:,:]
                    stack[last_channel, :,:,:] *= ~mask[0,:,:,:]
                    last_channel = last_channel + 1

                if self.verbose: log.info(f"Masked channels {self.mask_channels}. Resulting stack {stack.shape}, {stack.dtype}")
    
        #-------------------------------------------------
        # Apply scaling factors
        if self.scaling_factors is not None:
            # Convert stack to 32 bit float
            #stack = stack.astype(np.float32)
            
            for index, row in self.channel_index.iterrows():
                channel = int(row["channel"])
                scale_key = f"{row['ref_plate']}_ch{row['channel']}"
                     
                if str(scale_key) in self.scaling_factors.keys():
                    factor = self.scaling_factors[scale_key]
                    if self.verbose: log.debug(f"Scaling {scale_key} by factor {factor} for {row['plate']}, ch{channel}")
                    #if self.verbose: log.debug(f"Pre-scale min/max {np.min(stack[channel,:,:,:])}/{np.max(stack[channel,:,:,:])} dtype:{stack[channel,:,:,:].dtype}")
                    # Divide
                    stack[channel,:,:,:] /= factor
                    #if self.verbose: log.debug(f"Post-scale min/max {np.min(stack[channel,:,:,:])}/{np.max(stack[channel,:,:,:])} dtype:{stack[channel,:,:,:].dtype}")
                else:
                    if self.verbose: log.warning(f"Scale key {scale_key} not found! NOT APPLYING SCALING FOR: {row['plate']}, ch{channel}")
            
            
        log.info(f"Processed stack min: {np.min(stack)} max: {np.max(stack)}")
        
        # Scale back to uint
        if self.uint32:
            #stack=float_to_32bit_unint(stack)
            stack = np.clip(stack, 0, np.iinfo(np.uint32).max, out=stack)
            stack = stack.astype(np.uint32, copy=False)
            #for channel in range(0,stack.shape[0]):
            #    if self.verbose: log.debug(f"Post clip min/max channel {channel} {np.min(stack[channel,:,:,:])}/{np.max(stack[channel,:,:,:])} dtype:{stack[channel,:,:,:].dtype}")
        else:
            #stack=float_to_16bit_unint(stack)
            
            stack = np.round(stack, out=stack)
            log.debug("rounded")
            stack = np.clip(stack, 0, np.iinfo(np.uint16).max, out=stack)
            log.debug("clipped")
            stack = stack.astype(np.uint16)
            log.debug("cast")
            log.debug(f"{stack.dtype}, max: {np.max(stack)}")
            #for channel in range(0,stack.shape[0]):
            #    if self.verbose: log.debug(f"Post clip min/max channel {channel} {np.min(stack[channel,:,:,:])}/{np.max(stack[channel,:,:,:])} dtype:{stack[channel,:,:,:].dtype}")
    
        return stack
    
    
    # Fetch the registration matrices for a well as a dict
    def fetch_registration(self, iq) -> dict:
        pickle_path = f"{self.registration_dir}/{iq.plate}/{ImageQuery.ID_TO_ROW[iq.row]}/{iq.col}/{iq.field}.pickle"
        
        if os.path.isfile(pickle_path):
            pickle_file = open(pickle_path, "rb")
            alignment_matrices = pickle.load(pickle_file)
            pickle_file.close()
        else:
            raise RuntimeError(f"Alignment matrix for {pickle_path} does not exist")

        return alignment_matrices  
  
    
    def get_wells(self):
        return self.plate_reader.get_wells()


    def write_channel_index(self, output):
        
        # Write the output
        outdir = f"{output}/{self.plates[0]}"
        if not os.path.exists(outdir):
            os.makedirs(outdir)
            log.info(f"Folder created: {outdir}")    
        
        if not os.path.isfile(f"{outdir}/channel_indices.tsv"):
            self.channel_index.to_csv(f"{outdir}/channel_indices.tsv", sep="\t", index=False)