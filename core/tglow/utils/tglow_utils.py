import numpy as np
import xml.etree.ElementTree as ET
import re
import json
import string
import collections

# Encode numpy formats as JSON
class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)

# Get channel info easliy queriable in python as json format
# adapted from Matin Prete's ParseXML.py
def get_channel_channel_info(index_xml) -> dict:
    tree = ET.parse(index_xml)
    root = tree.getroot()

    default_namespace = re.findall(r'^{(.*)}',root.tag)[0]
    NS = {
        "PE": default_namespace,
        "OME": "http://www.openmicroscopy.org/Schemas/OME/2016-06"
    }

    entries=root.findall("PE:Maps/PE:Map/PE:Entry", NS)

    # Adapted & extended from Martins script
    channel_info={}
    for e in entries:
        cn = e.findall('./PE:ChannelName',NS)
        
        if cn:
            channel = {
                "id": e.attrib['ChannelID'],
                "name": cn[0].text,
                "image_type": e.find('./PE:ImageType',NS).text,
                "acquisition_type": e.find('./PE:AcquisitionType',NS).text,
                "illumination_type": e.find('./PE:IlluminationType',NS).text,
                "channel_type": e.find('./PE:ChannelType',NS).text,
                "binning_x": int(e.find('./PE:BinningX',NS).text),
                "binning_y": int(e.find('./PE:BinningY',NS).text),
                "emmisson": int(e.find('./PE:MainEmissionWavelength',NS).text),
                "excitation": int(e.find('./PE:MainExcitationWavelength',NS).text),
                "max_intensity": int(e.find('./PE:MaxIntensity',NS).text),
                "numerical_apeture": float(e.find('./PE:ObjectiveNA',NS).text),
                "objective_magnification": float(e.find('./PE:ObjectiveMagnification',NS).text),
                "exposure_time":{"unit": e.find('./PE:ExposureTime',NS).attrib['Unit'],
                                "value": float(e.find('./PE:ExposureTime',NS).text)},
                "size": (
                    int(e.find('./PE:ImageSizeX',NS).text),
                    int(e.find('./PE:ImageSizeY',NS).text)
                )}
            x = e.find('./PE:ImageResolutionX',NS)
            y = e.find('./PE:ImageResolutionY',NS)
            channel["image_resolution"]: {
                "x": {"unit":x.attrib["Unit"], "value":float(x.text),},
                "y": {"unit":y.attrib["Unit"], "value":float(y.text),},
            }
            # more info inside :
            ff_selector = f"./PE:Maps/PE:Map/PE:Entry[@ChannelID='{channel['id']}']/PE:FlatfieldProfile"
            channel["flatfield"] = root.find(ff_selector,NS).text
            channel_info[channel['id']] = channel
            
    return channel_info

# Build dict linking letters to number
def build_well_index(upper=True, invert=False, rows_as_string=False) -> dict:
    if (upper):
        chars = string.ascii_uppercase
    else:
        chars = string.ascii_lowercase
    
    nums = list(range(1, len(chars)))
    
    if rows_as_string:
        for i in range(0, len(nums)):
            nums[i] = str(nums[i])
        #nums=list([str(nums)for num in nums])

    if (invert):
        idx = dict(zip(nums,chars))
    else:
        idx = dict(zip(chars,nums))
    
    return(idx)

# Covert a defualt dict tree to regular dict             
def default_to_regular(d):
    if isinstance(d, collections.defaultdict):
        d = {k: default_to_regular(v) for k, v in d.items()}
    return d

# Convert a numpy 32bit float to a 16 bit int, clip values to zero and 65535
def float_to_16bit_unint(matrix) -> np.array:

    # Convert to 16 bit to keep consistency with output
    # Round to nearest int            
    matrix = np.rint(matrix)
    
    # Set negatives to zero
    matrix[matrix < 0] = 0
    
    # Clip values at the 16 bit max
    matrix[matrix > np.iinfo(np.uint16).max] = np.iinfo(np.uint16).max
    
    matrix = matrix.astype(np.uint16)
    
    return matrix

# Convert a numpy 32bit float to a 16 bit int, clip values to 
def float_to_32bit_unint(matrix) -> np.array:

    # Round to nearest int            
    matrix = np.rint(matrix)
    
    # Set negatives to zero
    matrix[matrix < 0] = 0
    
    # Clip values at the 32 bit max
    matrix[matrix > np.iinfo(np.uint32).max] = np.iinfo(np.uint32).max
    
    matrix = matrix.astype(np.uint32)
    
    return matrix

# Encode numpy formats as JSON
class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)
