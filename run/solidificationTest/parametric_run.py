#!/usr/bin/env python3

import os, sys, argparse
import re, glob, shutil, subprocess, fileinput
from distutils import file_util

str2list = lambda s: [float(item) for item in s.split(',')]

parser = argparse.ArgumentParser(description='Script for running parametric simulations')
parser.add_argument('-G', type=str2list, default='5e5,1e6,2e6', help='delimited list of temperature gradients')
parser.add_argument('-V', type=str2list, default='0.012,0.03', help='delimited list of pulling speeds')
parser.add_argument('-f', '--hostfile', type=str, default=os.environ['WCOLL'], help='path to the hostfile')
parser.add_argument('-s', '--scale', type=float, default=2, help='rescale all dimensions')
parser.add_argument('-y', '--yscale', type=float, default=2, help='rescale along the temperature gradient')
parser.add_argument('-W', '--Wscale', type=float, default=0.5, help='rescale the interface width')
parser.add_argument('-n', '--dry-run', action='store_true', help='do not launch simulations')
parser.add_argument('-v', '--verbose', action='store_true', help='increase output verbosity')
args = parser.parse_args()

if args.verbose:
    from termcolor import colored

regex_float = re.compile('-?[0-9]+\.?[0-9]*[Ee]?[+-]?[0-9]*')

def read_property(filename, name):
    with open(filename, 'r') as f:
        for line in f.readlines():
            if name in line:
                return line.split()[1].replace(';', '')

def overwrite_line(line, func):
    if args.verbose:
        print(colored(f'-{line}', 'red'), file=sys.stderr, end='')
    line = func(line)
    if args.verbose:
        print(colored(f'+{line}', 'green'), file=sys.stderr, end='')
    return line

def change_property(filename, name, value):
    if args.verbose:
        print(f' -- Change file {filename}')
    regex = re.compile(f'^ *{name} ')
    with fileinput.input(filename, inplace=True) as f:
        for line in f:
            if regex.findall(line):
                line = overwrite_line(line, lambda line: re.sub(regex_float, value, line, count=1))
            print(line, end='')

def overwrite_property(filename, name, func):
    value = float(read_property(filename, name))
    change_property(filename, name, f'{func(value)}')

def change_blockMesh(filename):
    ivert = 0
    regex_3num = re.compile(f'\([-\d.]+ +[-\d.]+ +[-\d.]+\)')
    with fileinput.input(filename, inplace=True) as f:
        for line in f:
            if 'vertices' in line:
                ivert = 1
            if regex_3num.findall(line) and 0 < ivert and ivert < 9:
                ivert = ivert + 1
                func = lambda line: ' '.join(map(lambda x: str(float(x[1])*args.yscale) \
                    if x[0]==1 else x[1], enumerate(line.split()))) + '\n'
                line = overwrite_line(line, func)
            if 'simpleGrading' in line:
                cells = regex_3num.findall(line)[0].split()
                cellsX = int(int(cells[0].replace('(', ''))*args.scale/args.Wscale)
                cellsY = int(int(cells[1])*args.scale*args.yscale/args.Wscale)
                line2 = line.split()
                line2[9] = f'({cellsX}'; line2[10] = f'{cellsY}'
                line = overwrite_line(line, lambda _: ' '.join(line2) + '\n')
            print(line, end='')

app = read_property('system/controlDict', 'application')
cwd = os.getcwd()

if not args.dry_run:
    with open(args.hostfile) as f:
        hosts = f.read().splitlines()

for i, G in enumerate(args.G):
    for j, V in enumerate(args.V):
        k = j + len(args.V)*i
        dotT = G*V
        print(f"Case {k}: G = {G:.2g}, V = {V}, dotT = {dotT:.2g}")
        case = f"_case_G_{G:.2g}_dotT_{dotT:.2g}"

        if os.path.isdir(case):
            print(f" -- Directory {case} already exists")
        else:
            # 1. Copy all case files
            os.mkdir(case)
            paths = subprocess.check_output('git ls-files', shell=True).split();
            for path in map(lambda s: s.decode(), paths):
                path1 = os.path.dirname(path)
                path2 = os.path.join(case, path1)
                if path1:
                    os.makedirs(path2, exist_ok=True)
                shutil.copy(path, path2, follow_symlinks=False)

            # 2. Change parameters
            f = os.path.join(case, 'constant/problemProperties')
            overwrite_property(f, 'interfaceWidth', lambda x: x*args.Wscale)
            overwrite_property(f, 'nSeeds', lambda x: int(x*args.scale))
            overwrite_property(f, 'frontPosition', lambda x: x/args.scale/args.yscale)
            overwrite_property(f, 'tempGradient', lambda x: G)
            overwrite_property(f, 'coolingRate', lambda x: dotT)

            f = os.path.join(case, 'system/blockMeshDict')
            overwrite_property(f, 'scale', lambda x: x*args.scale)
            change_blockMesh(f)

        if os.path.exists(os.path.join(case, f'log.{app}')):
            print(f" -- Case {case} is already simulated")
            continue
        if args.dry_run:
            continue

        # 3. Run simulations
        path = os.path.join(cwd, case)
        cmd = f'ssh -f {hosts[k]} "bash -c \'. /etc/profile; ml openfoam paraview; \
            cd {path}; DISPLAY=2 ./Allrun -parallel > log\'"'
        print(cmd)
        os.system(cmd)

