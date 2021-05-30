#!/usr/bin/env python

import argparse, os, re, json

parser = argparse.ArgumentParser(description='Convert OBJ to legacy VTK series')
parser.add_argument('dir', help='directory where OBJ files are placed')
args = parser.parse_args()

times = [ float(f) for f in os.listdir('.') if re.match(r'[0-9]+[0-9.e+-]*', f) ]
times.sort()
del times[0]

files = [ f[:-3] + 'vtk' for f in os.listdir(args.dir) if re.match(r'.*\.obj', f) ]

data = {
    "file-series-version": "1.0",
    "files" : [ { 'name': f, 'time': t } for f, t in zip(files, times)]
}

print(json.dumps(data, sort_keys=True, indent=4))
