#!/usr/bin/env python
import numpy as np
tau =.5
N = 50
x = (np.arange(N)+.5)/N/2
T = 1-tau*np.cos(2*np.pi*x)
dTdx = 2*np.pi*tau*np.sin(2*np.pi*x)

K1 = -0.6463
u = -K1*np.sqrt(T)*dTdx

def print_T(arr):
	print '\n'.join(['%f' % arr[i] for i in xrange(len(arr))])

def print_u(arr):
	print '\n'.join(['(%f 0 0)' % arr[i] for i in xrange(len(arr))])

print_T(T)
# print '\n'
# print_u(u)
