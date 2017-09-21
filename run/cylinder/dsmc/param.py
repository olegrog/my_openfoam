#!/usr/bin/env python

import numpy as np
import os, glob

problem = glob.glob('*.geo')[0].split('.')[0]
solver = os.path.basename(os.getcwd())
k_B = 1.3806488e-23
mass, ndensity = 2*k_B, .5/k_B
R, r = 7.5, 0.5
volume = np.pi/2*(R**2-r**2)
ensemble = 50e5
total = 40
window = total*.25

H = [ 30, 50, 70, 90, 105 ]
Kn = [ 2.48e-5, 4.4460e-4, 5.5130e-3, 0.13338, 1.8903 ]
Ma = [ 2.278, 3.190, 3.032, 2.403, 1.33 ]
idx = 3
H, Kn, Ma = [H[idx]], [Kn[idx]], [Ma[idx]]

for h, kn, ma in zip(H, Kn, Ma):
    dirname = '../_{dir}/{h}km'.format(dir=solver, h=h)
    if os.path.exists('{dir}/{file}.geo'.format(dir=dirname, file=problem)):
        continue
    print 'H = {0}km, Kn = {1:.2e}, Ma = {2:.2e}, ensemble={3:.2e}'.format(h, kn, ma, ensemble)
    os.system('mkdir -p {dir}; cp -r * {dir}'.format(dir=dirname))
    U, Tb = ma*(5./6)**.5, 1+ma**2/3
    os.system('sed -i.sed s/_U_/{U}/g {dir}/{file}'.format(U=U, dir=dirname, file='0/boundaryU'))
    os.system('sed -i.sed s/_Tb_/{Tb}/g {dir}/{file}'.format(Tb=Tb, dir=dirname, file='0/boundaryT'))
    os.system('rm {dir}/{files}'.format(dir=dirname, files='0/*.sed'))
    os.system('sed -i.sed s/_nParticles_/{n}/g {dir}/{file}'.format(n=ndensity*volume/ensemble, dir=dirname, file='constant/dsmcProperties'))
    os.system('sed -i.sed s/_diameter_/{d}/g {dir}/{file}'.format(d=np.sqrt(mass/(np.sqrt(2)*np.pi*kn)), dir=dirname, file='constant/dsmcProperties'))
    os.system('sed -i.sed s/_mass_/{m}/g {dir}/{file}'.format(m=mass, dir=dirname, file='constant/dsmcProperties'))
    os.system('sed -i.sed s/_ndensity_/{n}/g {dir}/{file}'.format(n=ndensity, dir=dirname, file='constant/dsmcProperties'))
    os.system('sed -i.sed s/_ndensity_/{n}/g {dir}/{file}'.format(n=ndensity, dir=dirname, file='system/dsmcInitialiseDict'))
    os.system('sed -i.sed s/_U_/{U}/g {dir}/{file}'.format(U=U, dir=dirname, file='system/dsmcInitialiseDict'))
    os.system('sed -i.sed s/_w_/{w}/g {dir}/{file}'.format(w=window, dir=dirname, file='system/controlDict'))
    os.system('sed -i.sed s/_total_/{t}/g {dir}/{file}'.format(t=total, dir=dirname, file='system/controlDict'))

