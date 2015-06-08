#!/usr/bin/env python
import numpy as np
alpha = .5
beta  = 2./np.sqrt(np.pi)
N = 200
x = (np.arange(N)+.5)/N
T = 1-alpha*np.cos(2*np.pi*x)
dTdx = 2*np.pi*alpha*np.sin(2*np.pi*x)

K1 = -0.6463
u  = beta - K1*np.sqrt(T)*dTdx

def print_T(arr):
    print '\n'.join(['%f' % arr[i] for i in xrange(len(arr))])

def print_u(arr):
    print '\n'.join(['(%f 0 0)' % arr[i] for i in xrange(len(arr))])

#print_T(T)
#print '\n'
print_u(u)
