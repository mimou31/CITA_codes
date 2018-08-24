import numpy as np
import healpy as hp
import os
from os import listdir
from os.path import isfile,join
import sys
sys.path.insert(0,'/scratch/r/rbond/remim/lenspower/py/')
import lenspowerWrapTools as lpw

lenspowerDataDir ='../data/'
lenspowerExecDir ='/scratch/r/rbond/remim/lenspower/c/'
fileNameBase = 'testing'

LMAX = 6144
l = np.arange(LMAX)
l = np.insert(l,0,0)
print(l)

path = '/scratch/r/rbond/remim/lenspower/data/'
output = '/scratch/r/rbond/remim/cal_sky_new/maps/lensed_theo/'
maps = '/scratch/r/rbond/remim/mocks/cib-maps/'

data = '/scratch/r/rbond/remim/lenspower/data/'

pathun = '/scratch/r/rbond/remim/mocks/test_lensing/maps/'

calsky = '/scratch/r/rbond/remim/cal_sky_new/maps/'

#maps_un = [os.path.basename(f) for f in listdir(pathun) if isfile(join(pathun,f)) and 'fits' in f and 'alm' not in f and 'pois' in f]
#maps_k = [os.path.basename(f) for f in listdir(kap) if isfile(join(kap,f)) and 'fits' in f and 'fwhm' in f]
#maps_k = [os.path.basename(f) for f in listdir(kap) if isfile(join(kap,f)) and 'fits' in f and 'fwhm' in f]
#name = [float((f.split("fwhm")[1]).split(".fits")[0]) for f in maps_k]
#print(name)
#kappa = hp.read_map('/scratch/r/rbond/remim/mocks/cmb-maps/octant1_kappa.fits')


#unlensed = hp.read_map('/scratch/r/rbond/remim/cal_sky_new/maps/cib/cib_field_nside2048_zmin0.4_zmax0.5_cib.fits')
#lensed = unlensed

    
cl_k1 = np.loadtxt(calsky+'kappa/total/kappa_total_2048_powerspec_23_l_bef.txt',usecols=np.arange(1,24))
cl_un1 = np.loadtxt(calsky+'cib/total/cib_total_2048_powerspec_20_l_bef.txt',usecols=np.arange(1,21))

z = np.arange(0.2,5.,0.2)
print(z)


for j in range(len(cl_un1)):
    cl_k = cl_k1[:,j]
    #cl_k = hp.anafast(hp.read_map(calsky+'kappa/kappa_field_nside2048_zmin0.0_zmax0.3_zsource0.25_kap_GRF.fits'))
    cl_dd = [4*cl_k[i]/(float(l[i+1]+1)*float(l[i+1]+1)) for i in range(1,len(cl_k))]
   # cl_dd = np.zeros(len(cl_dd))
    cl_dd = np.insert(cl_dd,0,0)
    #cl_k = cl_k1[:,j]
    cl_un = cl_un1[:,j]
    cl_l = cl_un
    #cl_dd = [4*cl_k[i]/(float(l[i+1]+1)*float(l[i+1]+1)) for i in range(1,len(cl_k))]
    #cl_psi = np.zeros(LMAX-1)
    

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
    print(len(cl))
    print(j)
    np.savetxt(calsky+'lensed_theo/total/cib_total_lens_zmin'+str(z[j])+'.txt',cl[1:])
    
sys.exit()
