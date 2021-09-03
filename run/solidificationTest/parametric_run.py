#!/usr/bin/env python3

import os, sys, argparse
import re, glob, shutil, subprocess, fileinput

str2list = lambda s: [float(item) for item in s.split(',')]

parser = argparse.ArgumentParser(description='Script for running parametric simulations')
parser.add_argument('-G', type=str2list, default='5e5,1e6,2e6', help='delimited list of temperature gradients')
parser.add_argument('-V', type=str2list, default='0.012,0.03', help='delimited list of pulling speeds')
parser.add_argument('-f', '--hostfile', type=str, default=None, help='path to the hostfile')
parser.add_argument('-s', '--scale', type=float, default=2, help='rescale all dimensions')
parser.add_argument('-y', '--yscale', type=float, default=2, help='rescale along the temperature gradient')
parser.add_argument('-W', '--Wscale', type=float, default=0.5, help='rescale the interface width')
parser.add_argument('-t', '--tscale', type=float, default=1.5, help='rescale the simulation time')
parser.add_argument('-N', '--Nwrites', type=int, default=100, help='number of writings during L/V')
parser.add_argument('--prefix', type=str, default='_case', help='prefix for directory names')
parser.add_argument('-p', '--progress', action='store_true', help='print progress instead of running')
parser.add_argument('-n', '--dry-run', action='store_true', help='do not launch simulations')
parser.add_argument('-v', '--verbose', action='store_true', help='increase output verbosity')
args = parser.parse_args()

if not args.hostfile:
    hostvar = 'WCOLL'
    print(f'Option --hostfile is missed. Trying to use {hostvar} variable...')
    args.hostfile = os.environ[hostvar]
if args.verbose or args.progress:
    from termcolor import colored
if args.progress:
    import time, humanize

class Regex:
    float1 = re.compile(r'-?[0-9]+\.?[0-9]*[Ee]?[+-]?[0-9]*')
    int3 = re.compile(r'\([-\d.]+ +[-\d.]+ +[-\d.]+\)')
    time = re.compile(r'^Time = ([0-9]+\.?[0-9]*[Ee]?[+-]?[0-9]*)')
    deltaT = re.compile(r'^deltaT = ([0-9]+\.?[0-9]*[Ee]?[+-]?[0-9]*)')

def reverse_readline(filename, buf_size=8192):
    """A generator that returns the lines of a file in reverse order

    Taken from https://stackoverflow.com/a/23646049/2531400
    """
    with open(filename) as fh:
        segment = None
        offset = 0
        fh.seek(0, os.SEEK_END)
        file_size = remaining_size = fh.tell()
        while remaining_size > 0:
            offset = min(file_size, offset + buf_size)
            fh.seek(file_size - offset)
            buffer = fh.read(min(remaining_size, buf_size))
            remaining_size -= buf_size
            lines = buffer.split('\n')
            if segment is not None:
                if buffer[-1] != '\n':
                    lines[-1] += segment
                else:
                    yield segment
            segment = lines[0]
            for index in range(len(lines) - 1, 0, -1):
                if lines[index]:
                    yield lines[index]
        if segment is not None:
            yield segment

def read_property(filename, name):
    regex_property = re.compile(f'^ *{name} +([a-zA-Z0-9+-.]+) *;')
    with open(filename, 'r') as f:
        for line in f.readlines():
            if (res := regex_property.findall(line)):
                return res[0]
    raise ValueError(f'Property `{name}` is not found in `{filename}`')

def read_last_regex(filename, regex):
    for line in reverse_readline(filename):
        if (res := regex.findall(line)):
            return res[0]
    raise ValueError(f'Regex `{regex}` is not found in `{filename}`')

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
    regex_property = re.compile(f'^ *{name} ')
    with fileinput.input(filename, inplace=True) as f:
        for line in f:
            if regex_property.findall(line):
                line = overwrite_line(line, lambda line: re.sub(Regex.float1, value, line, count=1))
            print(line, end='')

def overwrite_property(filename, name, func):
    value = float(read_property(filename, name))
    change_property(filename, name, f'{func(value)}')

def change_mesh(filename):
    Lmin = Lmax = ivert = 0
    with fileinput.input(filename, inplace=True) as f:
        for line in f:
            if 'vertices' in line:
                ivert = 1
            if Regex.int3.findall(line) and 0 < ivert and ivert < 9:
                ivert += 1
                func = lambda line: ' '.join(map(lambda x: str(float(x[1])*args.yscale) \
                    if x[0] == 1 else x[1], enumerate(line.split()))) + '\n'
                line = overwrite_line(line, func)
                y = float(line.split()[1])
                Lmin, Lmax = min(Lmin, y), max(Lmax, y)
            if 'simpleGrading' in line:
                cells = Regex.int3.findall(line)[0].split()
                cellsX = int(int(cells[0].replace('(', ''))*args.scale/args.Wscale)
                cellsY = int(int(cells[1])*args.scale*args.yscale/args.Wscale)
                line2 = line.split()
                line2[9] = f'({cellsX}'; line2[10] = f'{cellsY}'
                line = overwrite_line(line, lambda _: ' '.join(line2) + '\n')
            print(line, end='')
    return Lmax - Lmin

app = read_property('system/controlDict', 'application')
cwd = os.getcwd()

if args.progress:
    print(f'{"Case":30s}{"Vp":7s}{"deltaT":9s}{"Progress":14s}{"Elapsed time":20s}{"Finish time":20s}')
    now = time.time()
    for case in glob.glob(f'{args.prefix}*'):
        problemProperties = os.path.join(case, 'constant/problemProperties')
        controlDict = os.path.join(case, 'system/controlDict')
        videofile = os.path.join(case, 'video.avi')
        logfile = os.path.join(case, f'log.{app}')

        G = float(read_property(problemProperties, 'tempGradient'))
        dotT = float(read_property(problemProperties, 'coolingRate'))
        Vp = dotT/G
        end_time = float(read_property(controlDict, 'endTime'))
        if os.path.exists(logfile):
            curr_time = float(read_last_regex(logfile, Regex.time))
            deltaT = f'{float(read_last_regex(logfile, Regex.deltaT)):.1e}'
            progress = curr_time/end_time
            time_since_last_writing = now - os.path.getmtime(logfile)
            elapsed_time = os.path.getmtime(logfile) - os.path.getmtime(problemProperties)
            if progress < 1:
                remaining_time = elapsed_time/progress*(1 - progress)
                finish_time = humanize.naturaltime(remaining_time, future=True)
            else:
                finish_time = humanize.naturaltime(time_since_last_writing)
            progress = f'{progress:.0%}'
            color = 'yellow'
            if time_since_last_writing > 30:
                finish_time = 'stopped ' + humanize.naturaltime(time_since_last_writing)
                color = 'red'
        else:
            deltaT = ''
            elapsed_time = ''
            finish_time = ''
            progress = 'not started'
            color = 'red'
        if os.path.exists(videofile):
            finish_time = humanize.naturaltime(now - os.path.getmtime(videofile))
            progress += ' + video'
            color = 'green'
        elapsed_time = humanize.precisedelta(elapsed_time, minimum_unit='hours', format='%.0f')
        print(colored(f'{case:30s}{str(Vp):7s}{deltaT:9s}{progress:14s}{elapsed_time:20s}{finish_time:20s}', color))
    sys.exit()

if not args.dry_run:
    with open(args.hostfile) as f:
        hosts = f.read().splitlines()
        print(f'Hosts: {hosts}')
        if input(f'Are you sure to launch simulations? (y/n)') != 'y':
            sys.exit()

for i, G in enumerate(args.G):
    for j, V in enumerate(args.V):
        k = j + len(args.V)*i
        dotT = G*V
        print(f'Case {k}: G = {G:.2g}, V = {V}, dotT = {dotT:.2g}')
        case = f'{args.prefix}_G_{G:.2g}_dotT_{dotT:.2g}'

        if os.path.isdir(case):
            print(f' -- Directory {case} already exists')
        else:
            print(f' -- Directory {case} is created')
            # 1. Copy all case files
            os.mkdir(case)
            paths = subprocess.check_output('git ls-files', shell=True).split()
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
            scale = float(read_property(f, 'scale'))
            L = change_mesh(f)
            total_time = L*scale/V

            f = os.path.join(case, 'system/controlDict')
            overwrite_property(f, 'endTime', lambda x: total_time*args.tscale)
            overwrite_property(f, 'writeInterval', lambda x: total_time/args.Nwrites/args.yscale)

        if os.path.exists(os.path.join(case, f'log.{app}')):
            print(f' -- Directory {case} contains simulation results')
            continue
        if args.dry_run:
            continue

        # 3. Run simulations
        path = os.path.join(cwd, case)
        print(f' -- Simulation is running on {hosts[k]}')
        os.system(f'ssh -f {hosts[k]} "bash -c \'. /etc/profile; ml openfoam paraview; \
            cd {path}; DISPLAY=:2 ./Allrun -parallel > log\'"')
