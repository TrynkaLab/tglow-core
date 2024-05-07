#!/usr/bin/env python 
#
# Adapted from Martin Prete's scripts
####################################################################
# Script invocation:
# python actions/parse_xml.py \
#    --input_file "/path/to/Images/Index.xml" \
#    --output_path "/path/to/output"
####################################################################
# Store Index.xml information in JSON format for easier parsing
# and re-use of the pate/well/images metadata.
####################################################################

import os
import re
import json
import logging
import argparse
from xml.etree import ElementTree as ET
from image_query import ImageQuery

# Setup logging
logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

class PerkinElmerParser(object):
    """Parse PerkinElmer Index XML into a easy to access dictionary"""
    def __init__(self, index_path):
        self.index_file = index_path
        log.info(f"[+] Reading Index file: {index_path}")
        self.xml = ET.parse(index_path)
        self.NS = {"PE": re.findall(r'^{(.*)}',self.xml.getroot().tag)[0]}
        log.info(f"[+] Using namespace '{self.NS}'")
        self.parse_plate()
        self.parse_images()
        self.parse_wells()
        self.parse_channels()
        self.parse_planes()
        
    def estimate_pixel_sizes(self):
        """
        Extracts the pixel sizes from a PE index xml in (z, y, x) and returns in in microns
        Returns none if could not be estimated
        """
        img0=None
        img1=None

        for img in self.wells[0]["images"]:
            if img["plane"] == '1':
                img0=img
                next
            if img["plane"] == '2':
                img1=img
                break

        zres = abs(img0["position"]["z"]["value"] - img1["position"]["z"]["value"])
        zunit = img1["position"]["z"]["unit"]
    
        # Convert to microns
        if zunit == "m":
            zres=zres*1e6
        else:
            log.warn(f"Could not estimate z resolution, invalid unit {zunit}")
            zres=None

        # Y resolution
        if self.channels is not None and len(self.channels) > 0:
            yres = self.channels[0]["image_resolution"]["y"]["value"]
            yunit = self.channels[0]["image_resolution"]["y"]["unit"]
            
            # Convert to microns
            if yunit == "m":
                yres=yres*1e6
            else:
                log.warn(f"Could not estimate y resolution, invalid unit {yunit}")
                zres=None
        else:
            log.warn(f"Could not estimate y resolution, no channel info")
            zres=None

        # X resolution
        if self.channels is not None and len(self.channels) > 0:
            xres = self.channels[0]["image_resolution"]["x"]["value"]
            xunit = self.channels[0]["image_resolution"]["x"]["unit"]

            # Convert to microns
            if xunit == "m":
                xres=xres*1e6
            else:
                log.warn(f"Could not estimate x resolution, invalid unit {xunit}")
                xres=None
        else:
            log.warn(f"Could not estimate x resolution, no channel info")
            xres=None

        if zres is not None and yres is not None and xres is not None:
            return [zres, yres, xres]
        else:
            return None         
    
    def parse_channels(self):
        """Get PE channels"""
        log.info(f"[+] Reading Channels metadata")
        self.channels = []
        for e in self.xml.findall("./PE:Maps/PE:Map/PE:Entry", self.NS):
            cn = e.findall('./PE:ChannelName',self.NS)
            if cn:
                channel = {
                    "id": e.attrib['ChannelID'],
                    "name": cn[0].text,
                    "image_type": e.find('./PE:ImageType',self.NS).text,
                    "acquisition_type": e.find('./PE:AcquisitionType',self.NS).text,
                    "illumination_type": e.find('./PE:IlluminationType',self.NS).text,
                    "channel_type": e.find('./PE:ChannelType',self.NS).text,
                    "size": (
                        int(e.find('./PE:ImageSizeX',self.NS).text),
                        int(e.find('./PE:ImageSizeY',self.NS).text)
                    ),
                }
                x = e.find('./PE:ImageResolutionX',self.NS)
                y = e.find('./PE:ImageResolutionY',self.NS)
                channel["image_resolution"] = {
                    "x": {"unit":x.attrib["Unit"], "value":float(x.text),},
                    "y": {"unit":y.attrib["Unit"], "value":float(y.text),},
                }
                                
                # more info inside :
                # ff_selector = f"./PE:Maps/PE:Map/PE:Entry[@ChannelID='{channel['id']}']/PE:FlatfieldProfile"
                # channel["flatfield"] = self.xml.find(ff_selector,self.NS).text
                self.channels.append(channel)
                log.info(f" ├ {channel['id']} = {channel['name']}")
        
        log.info(f" └ Channel count = {len(self.channels)}")

        if len(self.channels) == 0:
            log.warn("No channel map found! This is the case for Operetta idx (TODO)")
            self.channels=None

    def parse_planes(self):
        """Get PE planes as a sorted list"""
        log.info(f"[+] Reading PlaneIDs from Images")
        all_planes = self.xml.findall("./PE:Images/PE:Image/PE:PlaneID", self.NS)
        self.planes = sorted(list(set([p.text for p in all_planes])))
        log.info(f" └ Found planes: {self.planes}")

    def parse_images(self):
        """Get all PE images"""
        log.info(f"[+] Reading images metadata")
        self.images={}
        for image in self.xml.findall("./PE:Images/PE:Image", self.NS):
            i = {
                "id": image.find("./PE:id", self.NS).text,
                "file": image.find("./PE:URL",self.NS).text,
                "state": image.find("./PE:State",self.NS).text,
                "row": int(image.find("./PE:Row",self.NS).text),
                "col": int(image.find("./PE:Col",self.NS).text),
                "field": image.find("./PE:FieldID",self.NS).text,
                "plane": image.find("./PE:PlaneID",self.NS).text,
                "channel": image.find("./PE:ChannelID",self.NS).text,
                "timepoint": image.find("./PE:TimepointID",self.NS).text,
                "sequence": image.find("./PE:SequenceID",self.NS).text if image.find("./PE:SequenceID",self.NS) is not None else None,
                "time_offset": image.find("./PE:MeasurementTimeOffset",self.NS).text,
                "time_abs": image.find("./PE:AbsTime",self.NS).text
            }
            x = image.find("./PE:PositionX",self.NS)
            y = image.find("./PE:PositionY",self.NS)
            z = image.find("./PE:PositionZ",self.NS)
            z_abs = image.find("./PE:PositionZ",self.NS)
            i["position"] = {
                "x": {"unit":x.attrib["Unit"], "value":float(x.text)},
                "y": {"unit":y.attrib["Unit"], "value":float(y.text)},
                "z": {"unit":y.attrib["Unit"], "value":float(z.text)},
                "z_abs": {"unit":z_abs.attrib["Unit"], "value":float(z_abs.text)}
            }
            self.images[i['id']] = i
        log.info(f" └ Images: {len(self.images)}")
    
    def parse_plate(self):
        """Get all PE Plate"""
        log.info(f"[+] Reading Plate metadata")
        plate = self.xml.find("./PE:Plates/PE:Plate", self.NS)
        self.plate = {
                "id": plate.find("./PE:PlateID", self.NS).text,
                "measurement": plate.find("./PE:MeasurementID",self.NS).text,
                "time": plate.find("./PE:MeasurementStartTime",self.NS).text,
                "name": plate.find("./PE:Name",self.NS).text,
                "type": plate.find("./PE:PlateTypeName",self.NS).text,
                "rows": int(plate.find("./PE:PlateRows",self.NS).text),
                "cols": int(plate.find("./PE:PlateColumns",self.NS).text)
        }
        log.info(f" ├ ID   = '{self.plate['id']}'")
        log.info(f" ├ Type = '{self.plate['type']}'")
        log.info(f" ├ Time = {self.plate['time']}")
        log.info(f" ├ Rows = {self.plate['rows']}")
        log.info(f" └ Cols = {self.plate['cols']}")
        
    def parse_wells(self):
        """Get all PE Wells"""
        log.info(f"[+] Reading Wells metadata")
        self.wells = []
        for well in self.xml.findall("./PE:Wells/PE:Well", self.NS):
            w = {
                "id": well.find("./PE:id", self.NS).text,
                "row": int(well.find("./PE:Row",self.NS).text),
                "col": int(well.find("./PE:Col",self.NS).text),
                "images": [self.images[wi.attrib["id"]] for wi in well.findall("./PE:Image",self.NS)]
            }
            self.wells.append(w)
        log.info(f" └ Wells: {len(self.wells)}")

    def save(self, output_path):

        if not os.path.exists(output_path):
            os.makedirs(output_path, exist_ok=True)
            
        out_json = os.path.join(output_path,f"{self.plate['name']}.json")
        log.info(f"[+] Saving JSON to : {out_json}")
        
        d = {
                "index": self.index_file,
                "plate": self.plate,
                "planes": self.planes,
                "wells": self.wells,
                "channels": self.channels
        }
        with open(out_json, "w") as f:
            json.dump(d, f, sort_keys=True, indent=1)
                     
    def write_manifest(self, output_path):
        if not os.path.exists(output_path):
            os.makedirs(output_path, exist_ok=True)
            
        manifest = open(f"{output_path}/manifest.tsv", 'w')
        
        for well in self.wells:
            iq = ImageQuery("", well["row"], well["col"], "")
            
            line = f"{iq.get_well_id()}\t{iq.get_well_id()[0]}\t{well['col']}\t{self.plate['name']}\t{self.index_file}\n"
            
            manifest.write(line)
        
        manifest.flush()
        manifest.close()
            
        
        

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse Perkin Elmer index XML file to JSON file')
    parser.add_argument('--input_file', type=str, required=False, help='Path to the PE Index file')
    parser.add_argument('--output_path', type=str, required=True, help='Path to the output directory where to storer the JSON file')
    parser.add_argument('--to_manifest', required=False, action='store_true', default=False, help='Store the output as a well,plate,index csv file')
    
    
    try:
       args = parser.parse_args()
    except:
       parser.print_help()
       exit(1)

    pe_data = PerkinElmerParser(args.input_file)
    
    if args.to_manifest:
        pe_data.write_manifest(args.output_path)
    else:
        pe_data.save(args.output_path)

