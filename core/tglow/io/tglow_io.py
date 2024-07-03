import sys
import numpy as np
import tifffile
import collections
import re
import os
import glob
import logging
import string
import csv

from aicsimageio import AICSImage
from aicsimageio.writers import OmeTiffWriter
from aicsimageio.readers.ome_tiff_reader import OmeTiffReader


from tglow.io.image_query import ImageQuery
from tglow.io.perkin_elmer_parser import PerkinElmerParser
from tglow.utils.tglow_utils import default_to_regular

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


class NamedImageQuery(ImageQuery):
    """Image query with an index that convers string names to indices"""
    
    def __init__(self, index, plate, row, col, field, channel=None, plane=None):
        super().__init__(index["plate"][plate],
                         index["row"][row],
                         index["col"][col],
                         index["field"][field],
                         index["channel"][channel],
                         index["plane"][plane])


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

        log.info(f"estimated dimensions in ZYX: {self.resolution}")
        

    def __peak_metadata__(self):
        """Get the resolution for the first image in the index, assume this is the same for the rest of the files"""
        
        self.plate_order = list(self.index.keys())        
        plate = self.plate_order[0]
        
        self.row_order = list(self.index[plate].keys())  
        self.row_order.sort(key=int)
        row = list(self.index[plate].keys())[0]
        
        self.col_order = list(self.index[plate][row].keys())
        self.col_order.sort(key=int)        
        col = self.col_order[0]
        
        self.field_order = list(self.index[plate][row][col].keys())
        self.field_order.sort(key=int)    
        field = self.field_order[0]
    
        self.channel_order = list(self.index[plate][row][col][field].keys())
        self.channel_order.sort(key=int)
        channel = self.channel_order[0]
        
        self.plane_order = list(self.index[plate][row][col][field][channel].keys())
        self.plane_order.sort(key=int)
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
    
    def get_filenames(self, query) -> list:
        """Get the filenames associated with an image query in the order they are read"""
        
        if not isinstance(query, ImageQuery):
            raise TypeError("Query is not of type ImageQuery")
                
        channels=[]
        if (query.channel is None):
            channels = self.get_channel_order()
        else:
            channels = [query.channel]
                
        planes=[]
        if (query.plane is None):
            planes = self.get_plane_order()
        else:
            planes = [query.plane]
        
        filenames=[]
        
        for channel in range(len(channels)):
            for plane in range(len(planes)):
                filenames.append(os.path.basename(self.index[query.plate][query.row][query.col][query.field][channels[channel]][planes[plane]]))
        
        return filenames
    
    def get_channel_order(self) -> list:
        return self.channel_order

    def get_plane_order(self) -> list:
        return self.plane_order
        

class PerkinElmerRawReader(IndexedImageReader):
    """Read raw image data from a PerkinElmer export or from a re-formatted output as long as it has an index file"""
    
    def __init__(self, index_xml, path, dtype=np.uint16, resolution=None,  plates_filter=None, fields_filter=None, blacklist=None) -> None:        
        self.pe_index=PerkinElmerParser(index_xml)
        self.path = path
        self.pixel_sizes=self.pe_index.estimate_pixel_sizes()
        log.info(f"Estimated pixels sizes to be {self.pixel_sizes} um (z, y, x)")
        cur_index = self.__convert_pe_index__()
        super().__init__(index=cur_index, path=path, dtype=dtype, resolution=resolution,  plates_filter=plates_filter, fields_filter=fields_filter)
    

    def __convert_pe_index__(self):
        plate = self.pe_index.plate["id"]
        nested_dict = lambda: collections.defaultdict(nested_dict)
        index = nested_dict()
        
        for file in self.pe_index.images.values():
            index[plate][str(file['row'])][str(file['col'])][str(file['field'])][str(file['channel'])][str(file['plane'])] = f"{self.path}/{file['file']}"
            
        return default_to_regular(index)
        

class AICSImageReader():
    """Reads image data from ome tiffs in a folder structure /plate/row/col/field.ome.tiff where field.ome.tiff is a CZYX array"""
    
    def __init__(self, path, plates_filter=None, fields_filter=None, blacklist=None, dtype=np.uint16, resolution=None) -> None:
        self.path = path
        self.plates_filter = plates_filter
        self.fields_filter = fields_filter
        
        if blacklist is None:
            self.blacklist = []
        else:
            self.blacklist = blacklist
        
        self.__build_index__()

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
    
    def __build_index__(self):
        
        plates = self.__list_directories__(self.path)
        
        if self.plates_filter is not None:
            plates = [plate for plate in plates if plate in self.plates_filter]
            
        nested_dict = lambda: collections.defaultdict(nested_dict)
        index = nested_dict()
        
        chars = string.ascii_uppercase
        nums = list(range(1, len(chars)))
        conv = dict(zip(chars,nums))

        nplates=0
        nwells=0
        nwells_blacklisted=0
        wells={}
        
        self.plates = set(plates)
        self.rows = {}
        self.cols = {}
        self.fields = {}
        self.images = {}

        for plate in plates:
            rows = self.__list_directories__(f"{self.path}/{plate}")
            
            if plate not in self.rows:
                self.rows[plate] = set()
            self.rows[plate].update(rows)
            
            nplates = nplates+1
            
            wells[plate] = set()
            for row in rows: 
                cols = self.__list_directories__(f"{self.path}/{plate}/{row}")
                
                if plate not in self.cols:
                    self.cols[plate] = set()
                    
                self.cols[plate].update(cols)
                
                for col in cols:
                    well = f"{row}{col.zfill(2)}"
                
                    # Ignore blacklisted wells in indexing
                    if f"{plate}:{well}" in self.blacklist:
                        nwells_blacklisted=nwells_blacklisted+1
                        continue
                        
                        
                    nwells = nwells+1
                    wells[plate].add(well)
                    
                    fields = glob.glob(f"{self.path}/{plate}/{row}/{col}/*.ome.tiff")
                    fields = [os.path.normpath(f) for f in fields]
                    fields = [os.path.basename(f) for f in fields]
                    fields = [f.replace(".ome.tiff", "") for f in fields]
                    
                    if self.fields_filter is not None:
                        fields = [field for field in fields if field in self.fields_filter]
                        
                    if plate not in self.fields:
                        self.fields[plate] = set()
                    self.fields[plate].update(fields)
                    
                    if plate not in self.images:
                        self.images[plate] = []
                                            
                    for field in fields:
                        iq = ImageQuery(plate, conv[row], col, field)
                        index[str(plate)][str(conv[row])][str(col)][str(field)] = iq #f"{self.path}/{plate}/{row}/{col}/{field}.ome.tiff"
                        self.images[plate].append(iq)
                        
        self.index = default_to_regular(index)
        self.wells = wells
                    
        #if self.fields_filter is not None:
        #    self.fields = [field for field in self.fields if field in self.fields_filter]
        
        log.info(f"Indexed {nwells} wells in {nplates} plates with {len(fields)} fields each, skipped {nwells_blacklisted} blacklisted wells")
        log.info(plates)
        log.info(wells)
    
    def __list_directories__(self, path):
        return [ name for name in os.listdir(path) if os.path.isdir(os.path.join(path, name)) ]
    
    def get_wells(self, plate):
        return self.wells[plate]
    
    def get_img(self, query):
        img = AICSImage(f"{self.path}/{query.plate}/{ImageQuery.ID_TO_ROW[query.row]}/{query.col}/{query.field}.ome.tiff")
        return img
    
    def read_image(self, query) -> np.ndarray:
        """Get a single image"""
        if not isinstance(query, ImageQuery):
            raise TypeError("Query is not of type ImageQuery")    
        
        #img = AICSImage(f"{self.path}/{query.plate}/{ImageQuery.ID_TO_ROW[query.row]}/{query.col}/{query.field}.ome.tiff")
        img = OmeTiffReader(f"{self.path}/{query.plate}/{ImageQuery.ID_TO_ROW[query.row]}/{query.col}/{query.field}.ome.tiff")

        if (query.channel is None) and (query.plane is None):
            # returns 4D CZYX numpy array
            #return img.get_image_data("CZYX", T=0)
            return img.get_image_dask_data("CZYX", T=0).compute()

        if (query.channel is not None) and (query.plane is None):
            # returns 3D ZYX numpy array
            #return img.get_image_data("ZYX", T=0, C=int(query.channel))
            return img.get_image_dask_data("ZYX", T=0, C=int(query.channel)).compute()

        if (query.channel is None) and (query.plane is not None):
            # returns 3D CYX numpy array
            #return img.get_image_data("CYX", T=0, Z=int(query.plane))
            return img.get_image_data("CYX", T=0, Z=int(query.plane)).compute()

        if (query.channel is not None) and (query.plane is not None):
            # returns 2D YX numpy array
            #return img.get_image_data("YX", T=0, Z=int(query.plane), C=int(query.channel))
            return img.get_image_dask_data("YX", T=0, Z=int(query.plane), C=int(query.channel)).compute()

    def read_stack(self, query) -> np.ndarray:
        """Read an image stack into a CZYX array"""
        if not isinstance(query, ImageQuery):
            raise TypeError("Query is not of type ImageQuery")
        
        query.channel = None
        query.plane = None
        
        return self.read_image(query)


class BlacklistReader():
    
    def __init__(self, path, sep="\t"):
        self.path=path
        self.sep=sep
        
    def read_blacklist(self, separator=":"):
    
        result = []
        with open(self.path, 'r') as file:
            reader = csv.reader(file, delimiter=self.sep)
            for row in reader:
                if len(row) >= 2:
                    combined_string = separator.join(row[:2])
                    result.append(combined_string)
        return result


        
class AICSImageWriter():
    """Writes image data from ome tiffs in a folder structure /plate/row/col/field.ome.tiff where field.ome.tiff is a CZYX array"""
    
    def __init__(self, path, channel_names=None, physical_pixel_sizes=None) -> None:
        self.path = path
        self.channel_names=channel_names
        self.physical_pixel_sizes=physical_pixel_sizes
        
    
    def write_stack(self, stack, query, channel_names=None, physical_pixel_sizes=None, image_names=None):
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
        
        if channel_names is None and self.channel_names is not None:
            channel_names = self.channel_names
        
        if physical_pixel_sizes is None and self.physical_pixel_sizes is not None:
            physical_pixel_sizes = self.physical_pixel_sizes
        
        if query.channel is not None and channel_names is not None:
            channel_names = [channel_names[int(query.channel)]]
        
        log.debug(f"pps: {physical_pixel_sizes}")
        log.debug(f"cnn: {channel_names}")
        log.debug(f"imn: {image_names}")
        log.debug(f"shp: {stack.shape}")


        OmeTiffWriter.save(stack,
         f"{outdir}/{query.field}.ome.tiff",
          dim_order="CZYX",
          channel_names=channel_names,
          physical_pixel_sizes=physical_pixel_sizes,
          image_names=image_names
          )
