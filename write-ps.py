import healpy as hp
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
from os import listdir
from os.path import isfile,join
import sys
import copy

NS = 2048
LMAX = 3*NS

def sort(value):
    n = float((value.split("zmax")[1]).split("_")[0])
    return(n)

path = '/scratch/r/rbond/remim/cal_sky_new/maps/lensed/GRF/'
lens = '/scratch/r/rbond/remim/remi_cib_lensing/'

maps = sorted([os.path.basename(f) for f in listdir(path) if isfile(join(path,f)) and 'fits' in f and 'alm' not in f and 'cib_field' in f],key=sort)
print(maps)
N = len(maps)

cl = np.zeros((LMAX,N+1))
cl[:,0] = np.arange(LMAX)

for i in range(N):
    m = hp.read_map(path+maps[i])
    t = hp.anafast(m)
    cl[:,i+1] = t

np.savetxt(path+'cib_field_GRF_lensed_'+str(NS)+'_powerspec_' + str(N) +'_l_bef.txt',cl)
