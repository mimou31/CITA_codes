import numpy as np
import healpy as hp
import os
from os import listdir
from os.path import isfile,join
import sys
sys.path.insert(0,'/scratch/r/rbond/remim/lenspower/py/')
import lenspowerWrapTools as lpw
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt


lenspowerDataDir ='../data/'
lenspowerExecDir ='/scratch/r/rbond/remim/lenspower/c/'
fileNameBase = 'testing'

LMAX = 6144
l = np.arange(LMAX)
l = np.insert(l,0,0)
print(len(l))

path = '/scratch/r/rbond/remim/cal_sky_new/maps/'
output = '/scratch/r/rbond/remim/cal_sky_new/maps/lensed_theo/'

mocks = '/scratch/r/rbond/remim/mocks/'

data = '/scratch/r/rbond/remim/lenspower/data/'

kappa = hp.read_map('/scratch/r/rbond/remim/mocks/cmb-maps/octant1_kappa.fits')

cl_k = hp.anafast(kappa)[:LMAX]

cl_dd = [4*cl_k[i]/(float(l[i+1])*float(l[i+1]+1)) for i in range(1,len(cl_k))]

cl_dd = np.insert(cl_dd,0,0)

unlensed = hp.read_map('/scratch/r/rbond/remim/mocks/cmb-maps/ffp1_unlensed.fits')
lensed = unlensed
    

cl_l = hp.anafast(lensed)[:LMAX]
cl_un = cl_l
   

    # UNLENSED

clTT = np.insert(cl_un,0,0)
clEE = np.insert(np.ones(LMAX),0,0)
clDD = np.insert(cl_dd,0,0)
clTE = np.insert(np.ones(LMAX),0,0)
clBB = np.insert(np.ones(LMAX),0,0)

    # LENSED

clTT_ = np.insert(cl_l,0,0)
clEE_ = np.insert(np.ones(LMAX),0,0)
clDD_ = np.insert(cl_dd,0,0)
clTE_ = np.insert(np.ones(LMAX),0,0)
clBB_ = np.insert(np.ones(LMAX),0,0)

    # NOISE

clTT_N = np.insert(np.zeros(LMAX),0,0)
clEE_N = np.insert(np.zeros(LMAX),0,0)
clDD_N = np.insert(np.zeros(LMAX),0,0)

np.savetxt(data+'unlensed_power_cmb_testing.txt',np.c_[l,clTT,clEE,clTE,clBB])
np.savetxt(data+'lensed_power_cmb_testing.txt',np.c_[l,clTT_,clEE_,clTE_,clBB_])
np.savetxt(data+'noise_power_T_testing.txt',np.c_[l,clTT_N])
np.savetxt(data+'noise_power_P_testing.txt',np.c_[l,clEE_N])
np.savetxt(data+'power_defl_testing.txt',np.c_[l,clDD])
np.savetxt(data+'noise_power_defl_testing.txt',np.c_[l,clDD_N])

execString = "cd " + lenspowerExecDir + '; ./lenspower ' + '"' + data + '" ' +  fileNameBase + ' power'
print '*** running: ', execString
os.system(execString)

cl= np.loadtxt(data+'lensed_power_cmb_exact_testing.txt',usecols=1)

np.savetxt(mocks+'cmb-maps/cmb_theo_august.txt',cl)
