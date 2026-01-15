"""Higher-level processed image provider.

`ProcessedImageProvider` assembles image stacks from one or more plates,
applies optional flatfield correction, registration, masking and scaling and
returns a combined CZYX (or similar) numpy array ready for downstream
processing.
"""

import logging
import copy
import pandas as pd
import numpy as np
import os
import pickle

from tglow.io.image_query import ImageQuery
from tglow.io.tglow_io import AICSImageReader, BlacklistReader
from tglow.utils.tglow_utils import apply_registration, apply_registration_cv, float_to_32bit_unint, float_to_16bit_unint
from basicpy import BaSiC

# Logging
logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def sigmoid(x, slope, bias):
    return 1 / (1 + np.exp(-slope * (x - bias)))


class ProcessedImageProvider():
        """Builds processed image stacks for a plate.

        The provider can merge channels from multiple plates (cycles), apply
        flatfield corrections, registration transforms and optional masking
        to produce a single multi-channel stack per well/field.
        """
        
    def __init__(self, path, plate, blacklist=None, plate_merge=None, registration_dir=None, flatfields=None, scaling_factors=None, mask_channels=None, mask_dir=None, mask_pattern=None, uint32=False, verbose=True, scaling_slope=0.001, scaling_bias=None):
        
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

        # Do some input checking to avoid runtime issues
        if ((scaling_slope is None) and (scaling_bias is not None)) | ((scaling_slope is not None) and (scaling_bias is None)):
            log.error("Either scaling_bias or scaling_slope is None. Either both must be set to a value, or must both be set to None")
            raise RuntimeError("Must supply both scaling_bias and scaling slope")
            
        if ((scaling_slope is not None) and (scaling_factors is None)):
            log.error("Must supply scaling_factors when using scaling_slope & scaling_bias")
            raise RuntimeError("Must supply scaling_factors when using scaling_slope & scaling_bias")            

        # These control the shape of the sigmoid curve of intensity that is 
        # optionally applied to the scaling factor
        if scaling_slope == None:
            self.scaling_slopes=None
        else:
            self.scaling_slopes={}
            for val in scaling_slope:
                keypair = val.split("=")
                for curp in self.plates + self.plates_merge:
                    if curp in keypair[0]:
                        log.info(f"Adding scaling slope: {keypair}")
                        self.scaling_slopes[keypair[0]] = float(keypair[1])
        
        if scaling_bias == None:
            self.scaling_biases=None
        else:
            self.scaling_biases={}
            for val in scaling_bias:
                keypair = val.split("=")
                
                for curp in self.plates + self.plates_merge:
                    if curp in keypair[0]:
                        log.info(f"Adding scaling bias: {keypair}")
                        self.scaling_biases[keypair[0]] = float(keypair[1])

        #---------------------------------------------------------------------
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
        """Construct a pandas DataFrame describing output channels.

        The `channel_index` maps output channel ids to the source plate/cycle and
        original channel id and name. It is used to populate and transform the
        assembled stack.
        """
    
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
        """Fetch and return a processed image stack for the provided `ImageQuery`.

        Returns a numpy array (C,Z,Y,X) with corrections and merges applied.
        """
 
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
        # Read any additional cycles
        if self.plates_merge is not None:
            for m_plate in self.plates_merge:                
                m_iq = copy.deepcopy(iq)
                m_iq.plate = m_plate
              
                # Check the fields are available
                m_fields = self.plate_reader.get_fields(m_iq)
                        
                if iq.field not in m_fields:
                    if self.verbose: log.warning(f"Field {iq.field} not found for {m_plate}. Returning None")
                    return None
                
                m_dims = self.plate_reader.get_img(m_iq).dims
                cur_range = range(last_channel, last_channel + m_dims['C'][0])
                
                stack[cur_range,:,:,:] = self.plate_reader.read_image(m_iq)
                if self.verbose: log.info(f"Read images for {m_plate}")
                
                last_channel = last_channel + m_dims['C'][0]
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

                    basic_model.darkfield = basic_model.darkfield.astype(np.float32)
                    basic_model.flatfield = basic_model.flatfield.astype(np.float32)

                    # modify in place
                    stack[channel,:,:,:] -= basic_model.darkfield[np.newaxis,:,:]
                    stack[channel,:,:,:] /= basic_model.flatfield[np.newaxis,:,:]
                                                    
            if self.verbose: log.info(f"Applied flatfieds to stack of shape {stack.shape}, {stack.dtype}")
        
        #-------------------------------------------------
        # Apply registration
        if self.registration_dir is not None and self.plates_merge is not None:
            # Reset last channel to first cycle max
            last_channel = img.dims['C'][0]

            for m_plate in self.plates_merge:
                m_iq = copy.deepcopy(iq)
                m_iq.plate = m_plate
                m_dims = self.plate_reader.get_img(m_iq).dims
                cur_range = range(last_channel, last_channel + m_dims['C'][0])
                                
                if self.verbose: log.info(f"Applying registration for {m_plate}")
                reg = self.fetch_registration(iq)
                
                # Make sure the registration has the right dtype
                # This is a 3x3 homography matrix
                # In the tglow pipeline, only the translation is set
                # [1, 0, x]
                # [0, 1, y]
                # [0, 0, 1]
                reg_mat = reg[m_plate].astype(np.float32)

                # https://stackoverflow.com/questions/78691652/applying-an-affine-transformation-using-opencv-skimage-and-scipy-returns-diffe
                # As the matrices were saved following the StackReg convention, we need to flip them again
                reg_mat = np.linalg.inv(reg_mat)
                
                # Use open cv warpPerspective instead of skimage, as it is much faster
                stack[cur_range,:,:,:] = apply_registration_cv(stack[cur_range,:,:,:], reg_mat)
                last_channel = last_channel +  m_dims['C'][0]

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
                    
                    if self.scaling_biases is None:
                        # Divide
                        stack[channel,:,:,:] /= factor
                    else:
                        scaling_bias = self.scaling_biases[scale_key]
                        scaling_slope = self.scaling_slopes[scale_key]

                        if self.verbose: log.debug(f"Weighing scale factor with sigmoid for {scale_key} with slope, bias {scaling_slope}, {scaling_bias}")
                        
                        # Determine the weight of the scaling for each pixel value, this forms a "soft threshold"
                        # which means the scaling will be softenend for low intensities. scaling_bias
                        # sets the point where the sigmoid returns 0.5.
                        # scaling_slope controls the slope or smoothnes of the transition
                        # Setting it too smooth will have a bad impact on the data. Setting the bias too low is equivalent
                        # to scaling equally over all pixels. 
                        # Optimal values are pre-caclulated outside this script.
                        scale_weight = sigmoid(stack[channel,:,:,:], scaling_slope, scaling_bias)
                        
                        # This ensures when the weight is 0, no scaling is applied and when the weight
                        # is one, the scaling is equal to factor. It also ensures if scaling is >0<1
                        # it works in the way thats intended, i.e. the values get closer to 1 when the
                        # weight goes down
                        stack[channel,:,:,:] /= ((scale_weight * (factor-1)) + 1)
                        
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
            stack = stack.astype(np.uint16, copy=False)
            log.debug("cast")
            log.debug(f"{stack.dtype}, max: {np.max(stack)}")
            #for channel in range(0,stack.shape[0]):
            #    if self.verbose: log.debug(f"Post clip min/max channel {channel} {np.min(stack[channel,:,:,:])}/{np.max(stack[channel,:,:,:])} dtype:{stack[channel,:,:,:].dtype}")
    
        return stack
    
    
    # Fetch the registration matrices for a well as a dict
    def fetch_registration(self, iq) -> dict:
        """Load registration matrices for a well from the registration directory.

        Expects a pickle containing a dict of plate->3x3 homography matrices.
        """
        pickle_path = f"{self.registration_dir}/{iq.plate}/{ImageQuery.ID_TO_ROW[iq.row]}/{iq.col}/{iq.field}.pickle"
        
        if os.path.isfile(pickle_path):
            pickle_file = open(pickle_path, "rb")
            alignment_matrices = pickle.load(pickle_file)
            pickle_file.close()
        else:
            raise RuntimeError(f"Alignment matrix for {pickle_path} does not exist")

        return alignment_matrices  
  
    
    def get_wells(self):
        """Return the set of wells indexed for the main plate."""
        return self.plate_reader.get_wells()


    def write_channel_index(self, output):
        """Write `channel_index` to `output/<plate>/channel_indices.tsv`.

        If the output folder does not exist it will be created.
        """
        # Write the output
        outdir = f"{output}/{self.plates[0]}"
        if not os.path.exists(outdir):
            os.makedirs(outdir)
            log.info(f"Folder created: {outdir}")    
        
        if not os.path.isfile(f"{outdir}/channel_indices.tsv"):
            self.channel_index.to_csv(f"{outdir}/channel_indices.tsv", sep="\t", index=False)