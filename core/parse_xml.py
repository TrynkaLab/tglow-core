#!/usr/bin/env python 
#
# Copied from Martin Prete's scripts
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

logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

class PerkinElmerParser(object):
    """Parse PerkinElmer Index XML into a easy to access dictionary"""
    def __init__(self, index_path):
        self.index_file = index_path
        log.debug(f"[+] Reading Index file: {index_path}")
        self.xml = ET.parse(index_path)
        self.NS = {"PE": re.findall(r'^{(.*)}',self.xml.getroot().tag)[0]}
        log.debug(f"[+] Using namespace '{self.NS}'")
        self.parse_plate()
        self.parse_images()
        self.parse_wells()
        self.parse_channels()
        self.parse_planes()
        
    def parse_channels(self):
        """Get PE channels"""
        log.debug(f"[+] Reading Channels metadata")
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
                channel["image_resolution"]: {
                    "x": {"unit":x.attrib["Unit"], "value":float(x.text),},
                    "y": {"unit":y.attrib["Unit"], "value":float(y.text),},
                }
                # more info inside :
                # ff_selector = f"./PE:Maps/PE:Map/PE:Entry[@ChannelID='{channel['id']}']/PE:FlatfieldProfile"
                # channel["flatfield"] = self.xml.find(ff_selector,self.NS).text
                self.channels.append(channel)
                log.debug(f" ├ {channel['id']} = {channel['name']}")
        log.debug(f" └ Channel count = {len(self.channels)}")


    def parse_planes(self):
        """Get PE planes as a sorted list"""
        log.debug(f"[+] Reading PlaneIDs from Images")
        all_planes = self.xml.findall("./PE:Images/PE:Image/PE:PlaneID", self.NS)
        self.planes = sorted(list(set([p.text for p in all_planes])))
        log.debug(f" └ Found planes: {self.planes}")

    def parse_images(self):
        """Get all PE images"""
        log.debug(f"[+] Reading images metadata")
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
                "sequence": image.find("./PE:SequenceID",self.NS).text,
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
        log.debug(f" └ Images: {len(self.images)}")
    
    def parse_plate(self):
        """Get all PE Plate"""
        log.debug(f"[+] Reading Plate metadata")
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
        log.debug(f" ├ ID   = '{self.plate['id']}'")
        log.debug(f" ├ Type = '{self.plate['type']}'")
        log.debug(f" ├ Time = {self.plate['time']}")
        log.debug(f" ├ Rows = {self.plate['rows']}")
        log.debug(f" └ Cols = {self.plate['cols']}")
        

    def parse_wells(self):
        """Get all PE Wells"""
        log.debug(f"[+] Reading Wells metadata")
        self.wells = []
        for well in self.xml.findall("./PE:Wells/PE:Well", self.NS):
            w = {
                "id": well.find("./PE:id", self.NS).text,
                "row": int(well.find("./PE:Row",self.NS).text),
                "col": int(well.find("./PE:Col",self.NS).text),
                "images": [self.images[wi.attrib["id"]] for wi in well.findall("./PE:Image",self.NS)]
            }
            self.wells.append(w)
        log.debug(f" └ Wells: {len(self.wells)}")

    def save(self, output_path):
        if not os.path.exists(output_path):
            os.makedirs(output_path, exist_ok=True)
        out_json = os.path.join(output_path,f"{self.plate['name']}.json")
        log.debug(f"[+] Saving JSON to : {out_json}")
        d = {
                "index": self.index_file,
                "plate": self.plate,
                "planes": self.planes,
                "wells": self.wells,
                "channels": self.channels
        }
        with open(out_json, "w") as f:
            json.dump(d, f, sort_keys=True, indent=1)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse Perkin Elmer index XML file to JSON file')
    parser.add_argument('--input_file', type=str, required=False, help='Path to the PE Index file')
    parser.add_argument('--output_path', type=str, required=True, help='Path to the output directory where to storer the JSON file')
    try:
       args = parser.parse_args()
    except:
       parser.print_help()
       exit(1)

    pe_data = PerkinElmerParser(args.input_file)
    pe_data.save(args.output_path)
