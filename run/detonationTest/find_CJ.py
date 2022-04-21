#!/usr/bin/env python3

import sys, argparse
import numpy as np
from numpy import sqrt, exp
from scipy.integrate import solve_ivp

parser = argparse.ArgumentParser(description='Script for finding the Chapman--Jouguet solution')
parser.add_argument('-g', '--gamma', type=float, default=1.2, help='heat capacity ratio')
parser.add_argument('-E', type=float, default=25, help='activation energy')
parser.add_argument('-Q', type=float, default=50, help='heat release')
parser.add_argument('-k', type=float, default=35.956, help='rate constant')
parser.add_argument('-L', type=float, default=10, help='integration length')
parser.add_argument('-t', '--tol', type=float, default=1e-6, help='solver tolerance')
parser.add_argument('-v', '--verbose', action='store_true', help='increase output verbosity')
parser.add_argument('-o', '--output', default=sys.stdout, help='output filename')
args = parser.parse_args()

### Parameter-dependent constants
g = args.gamma
D = sqrt(g + (g**2-1)/2*args.Q) + sqrt((g**2-1)/2*args.Q)

if args.verbose:
    print(f'D = {D:.5g}')

_1l = lambda l: np.maximum(1-l, 1e-16)
_U = lambda l: 1/(g+1)*(D**2-g)/D*(1+sqrt(_1l(l)))
_p = lambda l: (1+D**2)/(g+1)*(1 + (D**2-g)/(1+D**2)*sqrt(_1l(l)))
_rho = lambda l: D**2*(g+1)/(g*(1+D**2) - (D**2-g)*sqrt(_1l(l)))
_dl = lambda t,l: args.k*_1l(l)/(_U(l)-D)*exp(-_rho(l)/_p(l)*args.E)

sol = solve_ivp(_dl, [0, -args.L], [0.], rtol=args.tol, dense_output=True)
X, Y = sol.t, sol.y[0]
if args.verbose:
    print(f'lambda(-1) = {sol.sol(-1)[0]:.5g}')

names = [ 'x', 'lambda', 'p', 'rho', 'U' ]
np.savetxt(args.output, np.transpose((X, Y, _p(Y), _rho(Y), _U(Y))), fmt='%.5g',
    header='%6s'*(len(names)) % tuple(names))

