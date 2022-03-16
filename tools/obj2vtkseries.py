#!/usr/bin/env python

import argparse, os, sys, re, json

parser = argparse.ArgumentParser(description='Convert OBJ to legacy VTK series')
parser.add_argument('dir', help='directory where OBJ files are placed')
args = parser.parse_args()

for tdir in ['.', 'processor0']:
    times = [ float(f) for f in os.listdir(tdir) if re.match(r'[0-9]+[0-9.e+-]*', f) ]
    times.sort()
    del times[0]
    if len(times):
        break

if not len(times):
    print('Empty list of times!', file=sys.stderr)

files = [ f[:-3] + 'vtk' for f in os.listdir(args.dir) if re.match(r'.*\.obj', f) ]

data = {
    "file-series-version": "1.0",
    "files" : [ { 'name': f, 'time': t } for f, t in zip(files, times)]
}

print(json.dumps(data, sort_keys=True, indent=4))
