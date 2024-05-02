import sys
import numpy as np
import tifffile
import collections
import re
import os
import glob
import logging

from aicsimageio import AICSImage
from aicsimageio.writers import OmeTiffWriter

from parse_xml import PerkinElmerParser
from tglow_utils import build_well_index, default_to_regular

# Logging
logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

# Informal interface for plate
#class plateIo:
#    
#    def read_well(self, well, channel=None) -> np.array:
#        """Read a well as a 5d array. If channel is specified, read a specific channel"""
#        pass
#    
#    def write_well(self, well, channel=None):
#        """Write a well as a 5d tiff stack in CZYX. If channel is specified write a specific channel"""
#        pass
#    
#    def get_index(self) -> dict:
#        """Get the plate index which gives the path to files"""
#        pass
    
#class imageProvider:
#    
#    def get_image(self, plate, row, col, field, channel, plane) -> np.array:
#        """Get an image"""
#        pass

class ImageQuery:
    
    ROW_TO_ID=build_well_index(invert=False, rows_as_string=True)
    ID_TO_ROW=build_well_index(invert=True, rows_as_string=True)

    def __init__(self, plate, row, col, field, channel=None, plane=None):
        
        if type(plate) is str:
            self.plate = plate
        else:
            raise TypeError("Not a valid type for plate")
        
        if type(row) is int:
            self.row = str(row)
        elif type(row) is str:
            self.row = row
        else:
            raise TypeError("Not a valid type for row")
            
        if type(col) is int:
            self.col = str(col)
        elif type(col) is str:
            self.col = col
        else:
            raise TypeError("Not a valid type for col")
        
        if type(field) is int:
            self.field = str(field)
        elif type(field) is str:
            self.field = field
        else:
            raise TypeError("Not a valid type for field")
        
        if channel is not None:
            if type(channel) is int:
                self.channel = str(channel)
            elif type(channel) is str:
                self.channel = channel
            else:
                raise TypeError("Not a valid type for channel")
        else:
            self.channel=None
            
        if plane is not None:
            if type(plane) is int:
                self.plane = str(plane)
            elif type(plane) is str:
                self.plane = plane
            else:
                raise TypeError("Not a valid type for plane")
        else:
            self.plane=None
            
    def get_well_id(self):
        return f"{ImageQuery.ID_TO_ROW[self.row]}{self.col.zfill(2)}"
    
    def well_id_to_index(well_id) -> tuple:
        
        pat = re.compile(r'^([a-z])(\d+)', flags=re.IGNORECASE)
        match = re.match(pat, well_id)
        if match:
            row = match.group(1)
            col = int(match.group(2))
        else:
            raise Exception(f"No match found for {well_id}, does not match ^[a-Z]\d+")
        return (int(ImageQuery.ROW_TO_ID[row]), col)
        
    def to_string(self):
        
        if (self.channel is None and self.plane is None):
            return f"r{self.row}c{self.col}f{self.field}"
        elif (self.channel is None and self.plane is not None):
            return f"r{self.row}c{self.col}f{self.field}p{self.plane}"
        elif (self.channel is not None and self.plane is None):
            return f"r{self.row}c{self.col}f{self.field}ch{self.channel}"
        else:
            return f"r{self.row}c{self.col}f{self.field}ch{self.channel}p{self.plane}"

    
class IndexedImageReader:
    """Read raw image data from a dictionary of /plate/row/col/field/channel/plane/path_to_file.tiff"""
    
    def __init__(self, index, path, dtype=np.uint16, resolution=None) -> None:
        """Prefix path to the subpaths in index"""
        self.path = path
        self.index = index
        self.dtype = dtype
        
        self.__peak_metadata__()
        
        if resolution is not None:
            self.resolution = resolution
        
        log.info(f"Inititalized reader with")
        log.info(f"plate: {self.plate_order}")
        log.info(f"rows: {self.row_order}")
        log.info(f"cols: {self.col_order}")
        log.info(f"fields: {self.field_order}")
        
        log.info(f"C: {self.channel_order}")
        log.info(f"Z: {self.plane_order}")

        log.info(f"estimated resolution in ZYX: {self.resolution}")

    def __peak_metadata__(self):
        """Get the resolution for the first image in the index, assume this is the same for the rest of the files"""
        
        self.plate_order = list(self.index.keys())        
        plate = self.plate_order[0]
        
        self.row_order = list(self.index[plate].keys())        
        row = list(self.index[plate].keys())[0]
        
        self.col_order = list(self.index[plate][row].keys())        
        col = self.col_order[0]
        
        self.field_order = list(self.index[plate][row][col].keys())        
        field = self.field_order[0]
        
        self.channel_order = list(self.index[plate][row][col][field].keys())
        channel = self.channel_order[0]
        
        self.plane_order = list(self.index[plate][row][col][field][channel].keys())
        plane = self.plane_order[0]

        img = tifffile.imread(self.index[plate][row][col][field][channel][plane])
        self.resolution = (len(self.index[plate][row][col][field][channel].keys()), img.shape[0], img.shape[1])
                    
    def get_image(self, query) -> np.ndarray:
        """Get a single image"""
        if not isinstance(query, ImageQuery):
            raise TypeError("Query is not of type ImageQuery")
        
        cur_img = tifffile.imread(self.index[query.plate][query.row][query.col][query.field][query.channel][query.plane])
        return(cur_img)

    def read_stack(self, query) -> np.ndarray:
        """Read an image stack into a CZYX array"""
        if not isinstance(query, ImageQuery):
            raise TypeError("Query is not of type ImageQuery")
        
        i=0
        arrays=[] 
        
        channels=[]
        if (query.channel is None):
            channels = self.get_channel_order()
        else:
            channels = [query.channel]
        
        #log.info(f"Identified {len(channels)} channels")        
        
        planes=[]
        if (query.plane is None):
            planes = self.get_plane_order()
        else:
            planes = [query.plane]
            
        #log.info(f"Identified {len(planes)} planes")        

        # Init empty 4d numpy array
        stack = np.empty(shape=(len(channels),
                                len(planes),
                                self.resolution[1],
                                self.resolution[2]),
                                dtype=self.dtype)
        
        for channel in range(len(channels)):
            for plane in range(len(planes)):
                stack[channel][plane] = tifffile.imread(self.index[query.plate][query.row][query.col][query.field][channels[channel]][planes[plane]])
        
        log.info(f"read stack of {stack.shape}")        
        
        return stack
    
    def get_channel_order(self) -> list:
        return self.channel_order

    def get_plane_order(self) -> list:
        return self.plane_order
        

class PerkinElmerRawReader(IndexedImageReader):
    """Read raw image data from a PerkinElmer export or from a re-formatted output as long as it has an index file"""
    
    def __init__(self, index_xml, path, dtype=np.uint16, resolution=None) -> None:        
        self.pe_index=PerkinElmerParser(index_xml)
        self.path = path
        cur_index = PerkinElmerRawReader.__convert_pe_index__(self)
        super().__init__(index=cur_index, path=path, dtype=dtype, resolution=resolution)
    

    def __convert_pe_index__(self):
        plate = self.pe_index.plate["id"]
        nested_dict = lambda: collections.defaultdict(nested_dict)
        index = nested_dict()
        
        for file in self.pe_index.images.values():
            index[plate][str(file['row'])][str(file['col'])][str(file['field'])][str(file['channel'])][str(file['plane'])] = self.path + file['file']
            
        return default_to_regular(index)


class AICSImageReader():
    """Reads image data from ome tiffs in a folder structure /plate/row/col/field.ome.tiff where field.ome.tiff is a CZYX array"""
    
    def __init__(self, path, dtype=np.uint16, resolution=None) -> None:
        self.path = path

    def _deprecated_build_index__(self):
        
        wells = os.listdir(self.path)
        pattern = re.compile(r"^[A-Z]\d+$")
        wells = [s for s in wells if pattern.match(s)]
        
        plate = os.path.basename(self.path)
        
        nested_dict = lambda: collections.defaultdict(nested_dict)
        
        index = nested_dict()
        
        for well in wells:
            row, col = ImageQuery.well_id_to_index(well)
                        
            fov = glob.glob(f"{self.path}/{well}/*.tiff")
            fov.sort()
            for field in range(len(fov)):
                index[str(plate)][str(row)][str(col)][str(field)] = fov[field]
                
        self.index = default_to_regular(index)
    
    def get_img(self, query):
        img = AICSImage(f"{self.path}/{query.plate}/{ImageQuery.ID_TO_ROW[query.row]}/{query.col}/{query.field}.ome.tiff")
        return img
    
    def read_image(self, query) -> np.ndarray:
        """Get a single image"""
        if not isinstance(query, ImageQuery):
            raise TypeError("Query is not of type ImageQuery")    
        
        img = AICSImage(f"{self.path}/{query.plate}/{ImageQuery.ID_TO_ROW[query.row]}/{query.col}/{query.field}.ome.tiff")
        
        if (query.channel is None) and (query.plane is None):
            # returns 4D CZYX numpy array
            return img.get_image_data("CZYX", T=0) 
        
        if (query.channel is not None) and (query.plane is None):
            # returns 3D ZYX numpy array
            return img.get_image_data("ZYX", T=0, C=int(query.channel)) 
        
        if (query.channel is None) and (query.plane is not None):
            # returns 3D CYX numpy array
            return img.get_image_data("CYX", T=0, Z=int(query.plane))
        
        if (query.channel is not None) and (query.plane is not None):
            # returns 2D YX numpy array
            return img.get_image_data("YX", T=0, Z=int(query.plane), C=int(query.channel))

    def read_stack(self, query) -> np.ndarray:
        """Read an image stack into a CZYX array"""
        if not isinstance(query, ImageQuery):
            raise TypeError("Query is not of type ImageQuery")
        
        query.channel = None
        query.plane = None
        
        return self.read_image(query)
        
        
class AICSImageWriter():
    """Writes image data from ome tiffs in a folder structure /plate/row/col/field.ome.tiff where field.ome.tiff is a CZYX array"""
    
    def __init__(self, path) -> None:
        self.path = path
    
    def write_stack(self, stack, query, channel_names=None, physical_pixel_sizes=None):
        """Write a CZYX array into the folder structure"""

        if not isinstance(query, ImageQuery):
            raise TypeError("Query is not of type ImageQuery")
        
        if not isinstance(stack, np.ndarray):
            raise TypeError("Stack is not of type np array")
        
        if not len(stack.shape) == 4:
            raise TypeError("Stack is not 4 dimensional")
        
        outdir = f"{self.path}/{query.plate}/{ImageQuery.ID_TO_ROW[query.row]}/{query.col}"
        
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        
        OmeTiffWriter.save(stack,
         f"{outdir}/{query.field}.ome.tiff",
          dim_order="CZYX",
          channel_names=channel_names,
          physical_pixel_sizes=physical_pixel_sizes,
          )
