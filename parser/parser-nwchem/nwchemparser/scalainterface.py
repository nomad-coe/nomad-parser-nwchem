"""
This is the access point to the parser for the scala layer in the
nomad project.
"""
from __future__ import absolute_import
import sys
import setup_paths
from nomadcore.parser_backend import JsonParseEventsWriterBackend
from nwchemparser import NWChemParser


if __name__ == "__main__":

    # Initialise the parser with the main filename and a JSON backend
    main_file = sys.argv[1]
    parser = NWChemParser(main_file, backend=JsonParseEventsWriterBackend)
    parser.parse()