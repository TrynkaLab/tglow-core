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
from tglow.io.perkin_elmer_parser import PerkinElmerParser
#from tglow.io.image_query import ImageQuery

# Setup logging
logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

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

