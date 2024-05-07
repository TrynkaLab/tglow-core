import re
from tglow_utils import build_well_index

class ImageQuery:
    """Bean that stores data relevant for fetching and saving multiwell images"""

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

