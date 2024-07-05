import os
import re
import tifffile
import argparse
import numpy as np
import copy
from basicpy import BaSiC
from skimage import transform
from matplotlib import pyplot as plt
from pystackreg import StackReg

import pickle
import logging
from tglow.io.tglow_io import AICSImageReader
from tglow.io.image_query import ImageQuery
from tglow.utils.tglow_utils import float_to_16bit_unint, float_to_32bit_unint, dict_to_str

# Logging
logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

#------------------------------------------------------------------
# Main runner
class MergeAndAlign:
    
    # Constructor, and parse arguments
    def __init__(self, args):
          
        self.input=args.input
        self.output=args.output
        self.write_zstack=not args.no_zstack
        self.write_max_projection_onefile=args.max_project_onefile
        self.write_max_projection=args.max_project
        self.uint32=args.uint32
        
        
        if (self.write_max_projection):
            self.write_zstack=False
            log.info("Writing max projection as individual files, not saving z-stacks")
        
        # Plates to process
        self.plates=args.plate
        
        # Channels for alignment and plates to merge
        self.plates_merge=args.plate_merge   
        self.registration_dir=args.registration_dir 
        
        if (self.plates_merge is None and self.registration_dir is not None):
            raise RuntimeError("Need to supply plates for registration when registration dir is set")
        
        if self.plates_merge is not None:
            # Reader for plates
            self.plate_reader=AICSImageReader(self.input, plates_filter=self.plates + self.plates_merge)
        
            if (self.registration_dir is None):
                raise RuntimeError("Need to supply registration_dir for registration")
            
        else:
            # Reader for plates
            self.plate_reader=AICSImageReader(self.input, plates_filter=self.plates)
                    
        # Build dict with pretrained BaSiCpy models channel is the key
        if args.basicpy_model == None:
            self.flatfields=None
        else:
            self.flatfields={}
            for val in args.basicpy_model:
                keypair = val.split("=")
                log.info(f"{keypair}")
                self.flatfields[keypair[0]] = BaSiC.load_model(keypair[1])
                # Hack arround missing baseline (this should've been fixed in the dev version)
                #self.flatfields[keypair[0]].baseline=""     
            
        if (not (self.write_max_projection | self.write_zstack)):
            raise RuntimeError("No output option specified")
            #log.info(f"No output specified, exiting")
            #return
            
            
        if args.fields is None:
            self.fields = self.plate_reader.fields[self.plates[0]]
            log.info(f"Autodetected fields {self.fields} based on first plate and image")
        else:
            self.fields=args.fields
            

    # Main loops
    def run(self, well):

        # Setup for registration
        if (self.plates_merge != None):
            
            if (len(self.plates) !=1):
                raise RuntimeError("When merging channels --plate can only have one value as it is used as the reference plate")
            
            log.info(f"Writing channel info")
            channels = self.write_updated_channel_map()
        else:
            channels = None

        log.info(f"Staring merging")

        # Loop over possible plates
        for plate in self.plates:
        
            channel_index = self.get_bp_channel_keys(plate)
        
            row, col = ImageQuery.well_id_to_index(well)
            out_dir_final= f"{self.output}/{plate}/{ImageQuery.ID_TO_ROW[str(row)]}/{col}"
            
            if not os.path.exists(out_dir_final):
                os.makedirs(out_dir_final)
                log.info(f"Folder created: {out_dir_final}")
                
                    
            for field in self.fields:
                
                # Array to store the fields for the max projection
                max_projection = []
                
                iq = ImageQuery(plate, row, col, field)
                
                # If channels is none (when not merging)
                if channels is None:
                    channels = range(0, self.plate_reader.get_img(iq).dims['C'][0])
                
                # Read data in
                stack = self.plate_reader.read_image(iq)

                # If plates need to be merged, and or registred
                # do that here
                if self.plates_merge is not None:
                    
                    for m_plate in self.plates_merge:
                        m_iq = copy.deepcopy(iq)
                        m_iq.plate = m_plate
                        m_stack = self.plate_reader.read_image(m_iq)
                        
                        if self.registration_dir is not None:
                            log.info("Applying registration")
                            reg = self.fetch_registration(iq)
                            m_stack = self.apply_registration(m_stack, reg[m_plate])
                            
                        stack = np.concatenate((stack, m_stack), axis=0)

                log.info(f"Read merged image stack of shape {stack.shape}")
                
                # Loop over channels
                for channel in channels:
                
                    merged = stack[channel,:,:,:]
                
                    log.info(f"{plate} well: {well} channel: {str(channel)} field: {field}")
                    
                    # Final output filename for channel
                    cur_out = f"{out_dir_final}/{field}_{plate}_{well}_ch{str(channel)}.tiff"      
                
                    # If flatfields are specified, apply them here
                    bp_key=channel_index[f"{plate}_ch{channel}"]
                
                    if not self.flatfields == None:
                        if str(bp_key) in self.flatfields.keys():
                            log.info(f"Applying BaSiC model for channel: {str(bp_key)}")
                            basic_model = self.flatfields[str(bp_key)]
                            merged = basic_model.transform(merged)
                            
                            # Convert to either 32 or 16 bit unsigned int
                            if self.uint32:
                                merged=float_to_32bit_unint(merged)
                            else:
                                merged=float_to_16bit_unint(merged)
                                
                    # Write out the final tiff file, saves as ZYX
                    if self.write_zstack:
                        tifffile.imwrite(cur_out, merged, shape=merged.shape, imagej=True, photometric='MINISBLACK', metadata={'axes': 'ZYX'})
                        
                    # Store max projection accross the first axis (Z) for current channel 
                    if (self.write_max_projection | self.write_max_projection_onefile):
                        mp = np.max(merged, axis=0)
                        
                        if (self.write_max_projection):
                            log.info(f"{plate} well: {well} channel: {str(channel)} field: {field}, max projection")
                        
                            # Final output filename for channel
                            cur_out = f"{out_dir_final}/{field}_{plate}_{well}_ch{str(channel)}.tiff"
                            tifffile.imwrite(cur_out, mp, shape=mp.shape, imagej=True, metadata={'axes': 'YX'})        
                                
                        else:
                            max_projection.append(mp)

                # Write the max projection, CYX
                if (self.write_max_projection_onefile):
                    log.info(f"{plate} well: {well} field: {field} max projection")
                    cur_out = f"{out_dir_final}/{plate}_{well}_{field}_max_projection.tiff"
                    
                    max_merged = np.array(max_projection)
                    tifffile.imwrite(cur_out, max_merged, shape=max_merged.shape, imagej=True, metadata={'axes': 'CYX'})

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
        
    # Apply 2d registration 
    def apply_registration(self, stack, alignment_matrix):
            
        # Apply transformation matrix
        tform = transform.AffineTransform(matrix=alignment_matrix)

        # If the stack is 2d YX
        if len(stack.shape) == 2:
            # transform image using the saved transformation matrix
            stack = transform.warp(stack, tform, preserve_range=True)
            
        # If the stack is 3d CYX or ZYX
        elif len(stack.shape) == 3:
            
            for i in range(0, stack.shape[0]):
                #log.debug(f"Aligning stack of shape {stack[i,:,:].shape}")
                stack[i,:,:] = transform.warp(stack[i,:,:], tform, preserve_range=True)
        
        # If the stack is CZYX
        elif len(stack.shape) == 4:
            for i in range(0, stack.shape[0]):
                for j in range(0, stack.shape[1]):
                    #log.debug(f"Aligning stack of shape {stack[i,j,:,:].shape}")
                    stack[i,j,:,:] = transform.warp(stack[i,j,:,:], tform, preserve_range=True)
        
        return stack
                  
    # Write a json file with the channel id's in the merged tiffs
    def write_updated_channel_map(self):
    
         # Peak at an image for channel names
        img = self.plate_reader.get_img(self.plate_reader.images[self.plates[0]][0])
        
        channel_names_old=img.channel_names
                
        if channel_names_old is None:
            log.warning("Skipping writing new channel names as none were found")
            return
        
        plate_names = [self.plates[0] for channel in channel_names_old]
        
        if self.plates_merge != None:
            for plate_merge in self.plates_merge:
                img = self.plate_reader.get_img(self.plate_reader.images[plate_merge][0])
                channel_names_old += img.channel_names
                plate_names += [plate_merge for channel in img.channel_names]
                            
        
        outdir = f"{self.output}/{self.plates[0]}"
        if not os.path.exists(outdir):
            os.makedirs(outdir)
            log.info(f"Folder created: {outdir}")    
        
        if not os.path.isfile(f"{outdir}/channel_indices.tsv"):
            channel_index = open(f"{outdir}/channel_indices.tsv", "w")
            channel_index.write(f"channel\tchannel_name\told_channel_name\tplate\n")     
            
            rm_ch = re.compile("^ch\d+ - ")
            
            for i in range(0, len(channel_names_old)):
                log.debug(f"{i}: {channel_names_old[i]}")
                new_name = rm_ch.sub("", channel_names_old[i])
                channel_name_new = f"ch{i} - {new_name}"
                channel_index.write(f"{i}\t{channel_name_new}\t{channel_names_old[i]}\t{plate_names[i]}\n")
                
            channel_index.flush()
            channel_index.close()
        else:
            log.warning("channel_indices.tsv already exists for this plate, skipping writing writing again")
        
        return range(0, len(channel_names_old))
    
    # Make a map of new and old channel ids
    def get_bp_channel_keys(self, plate):
        
        channel_index = {}
        channels = range(0, self.plate_reader.get_img(self.plate_reader.images[plate][0]).dims['C'][0])

        new_idx=0
        
        for ch in channels:
            channel_index[f"{plate}_ch{new_idx}"] = f"{plate}_ch{ch}"
            new_idx+=1
            
        if self.plates_merge != None:
            for plate_merge in self.plates_merge:
                channels = range(0, self.plate_reader.get_img(self.plate_reader.images[plate_merge][0]).dims['C'][0])
                for ch in channels:
                    channel_index[f"{plate}_ch{new_idx}"] = f"{plate_merge}_ch{ch}"
                    new_idx+=1
                    
        return channel_index
            
    
    def printParams(self):
        
        log.info(f"Input plates:\t{str(self.plates)}")
        #log.info(f"Input fields:\t{str(self.fields)}")
        #log.info(f"Input planes:\t{str(self.planes)}")   
        #log.info(f"Input channel:\t{str(self.channels)}")     
        log.info(f"Input:\t\t{self.input}")
        log.info(f"Output:\t\t{self.output}")    
        log.info(f"Basicpy models:\t{dict_to_str(self.flatfields)}")    
        log.info(f"Output 32bit:\t{str(self.uint32)}")    
        
        log.info(f"----------------------------")    
        log.info(f"Merge plates:\t{str(self.plates_merge)}")     
        #log.info(f"Merge channels:\t{str(self.channels_merge)}")  
        #log.info(f"Ref channel:\t{str(self.ref_channel)}")  
        #log.info(f"Qry channel:\t{str(self.qry_channel)}")
        #log.info(f"Ref channel eval:\t{str(self.ref_channel_eval)}")  
        #log.info(f"Qry channel eval:\t{str(self.qry_channel_eval)}")  
        #log.info(f"Merge hist match:\t{str(self.eq_merge)}")  
        #log.info(f"Merge invert:\t{str(self.inv_merge)}")  
        #log.info(f"Transform:\t{str(self.transform)}")        
        #log.info(f"Save registr:\t{str(self.save_reg)}")     
        #log.info(f"Plot:\t\t{str(self.plot)}")     




# Main loop
if __name__ == "__main__":
            
    # CLI 
    parser = argparse.ArgumentParser(description="Merge raw Phenix tiff files organized per well into per channel tiffs, one page per z-stack")
    parser.add_argument('-w', '--well', help='Well ID to merge')
    parser.add_argument('-i','--input', help='Base dir to raw input')
    parser.add_argument('-o','--output', help='Output prefix')
    parser.add_argument('-p','--plate', help='Subfolder in raw dir to process', nargs='+')
    parser.add_argument('-m','--plate_merge', help='Plates to align using pystackreg (e.g. multiple cycles of the same plate). Output will be added as additional channels.', nargs='+', default=None)
    parser.add_argument('--no_zstack', help="Do not output one tiff file per channel containing all z-stacks", action='store_true', default=False)
    parser.add_argument('--max_project', help="Output max projection as per channel YX tiffs", action='store_true', default=False)
    parser.add_argument('--max_project_onefile', help="Output max projection one CYX tiff. Specify --no_zstack to skip per channel tiffs.", action='store_true', default=False)
    parser.add_argument('--registration_dir', help="Path to registration root storing <plate>/<row>/<col>/<field>.pickle", default=None)
    parser.add_argument('--fields', help='Fields to use. <field #1> | [<field #1> <field #2> <field #n>]', nargs='+', default=None)
    #parser.add_argument('--planes', help='Z planes to use. <plane #1> | [<plane #1> <plane #2> <plane #n>]', nargs='+', default=None)
    #parser.add_argument('--channels', help='Channels to use. <channel #1> | [<channel #1> <channel #2> <channel #n>]', nargs='+', default=None)
    #parser.add_argument('--channels_merge', help='Channels to for merging from second plate. <channel #1> | [<channel #1> <channel #2> <channel #n>]', nargs='+', default=None)
    #parser.add_argument('-r','--ref_channel', help='Reference channel for registration (nucleus)', nargs=1, default=None)
    #parser.add_argument('-q','--qry_channel', help='Query channel for registration (nucleus)', nargs=1, default=None)
    parser.add_argument('--basicpy_model', help='Basicpy model dir for a channel. If merging channels ids are assigned seqeuntially for extra channels <channel #1>=/path/to/model | [<channel #1>=/path/to/model <channel #n>=/path/to/model]', nargs='+', default=None)
    parser.add_argument('--uint32', help="Write as 32 bit unsigned integer instead of clipping to 16 bit uint after applying basicpy model", action='store_true', default=False)
    args = parser.parse_args()
  
    log.info("-----------------------------------------------------------")
    if not os.path.exists(args.output):
        os.makedirs(args.output)
        log.info(f"Folder created: {args.output}")
    

    runner=MergeAndAlign(args)
    log.info(f"Well:\t{args.well}")
    runner.printParams()
    log.info("-----------------------------------------------------------")
    
    runner.run(args.well)

    #main(well, input, output, ref_channel, qry_channel, write_zstack, write_max_projection, basicpy_models, uint32)







