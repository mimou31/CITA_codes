import healpy as hp
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

from matplotlib.gridspec import GridSpec 
plt.style.use('ggplot')
import utilities.savitzky_golay as sav

import os
from os import listdir
from os.path import isfile,join
import sys

def tri(s):
    return float((s.split("kap")[1]).split(".fits")[0])
def tritxt(s):
    return float((s.split("zmin")[1]).split(".txt")[0])


pathun = '/scratch/r/rbond/remim/cal-sky/test/cib_maps/slices/'
pathl = '/scratch/r/rbond/remim/mocks/cmb-maps/'

calsky = '/scratch/r/rbond/remim/cal_sky_new/maps/'
calsky2 = '/scratch/r/rbond/remim/cal_sky_new/maps_old_v2/'
maps_theo = sorted([os.path.basename(f) for f in listdir(calsky+'lensed_theo/GRF/') if isfile(join(calsky+'lensed_theo/GRF/',f))],key=tritxt)


zmax = [float((s.split("zmin")[1]).split(".txt")[0]) for s in maps_theo]
zmax = [1.2,1.24,1.28,1.32,1.36,1.4]
zmax = np.arange(2,50,2)
print(zmax)
maps_l=[]

cl_un=[]

print(maps_theo)

LMAX = 5000
LMIN = 10
n=len(maps_theo)
cl_theo1 = np.zeros((6144,n))
for i in range(len(maps_theo)):
   cl_theo1[:,i] = np.loadtxt(calsky+'lensed_theo/GRF/'+maps_theo[i],usecols=[0])
#np.savetxt(calsky+'lensed_theo/lensed_cib_no-l_4.txt',cl_theo1)


n=2
#cl_l1 = np.loadtxt(calsky+'lensed/GRF/cib_field_GRF_lensed_2048_powerspec_23_l_bef.txt',usecols=np.arange(1,n+1))
cl_l1 = np.loadtxt(calsky+'kappa/halos/2048/kappa_halos_2048_powerspec_23_l_bef.txt',usecols=np.arange(1,n+1))
cl_theo1 = np.loadtxt(calsky+'kappa/total/cl_kappa_remi_2_nonlinear.txt',usecols=np.arange(1,n+1))
cl_theo2 = np.loadtxt(calsky+'kappa/2048/cl_kappa_remi_2.txt',usecols=np.arange(1,n+1))

cl_un1 = np.loadtxt(calsky+'cib/GRF/cib_GRF_2048_powerspec_23_l_bef.txt',usecols=np.arange(1,n+1))
#cl_un1 = np.loadtxt(calsky+'kappa/total/kappa_total_2048_powerspec_23_l_bef.txt',usecols=np.arange(1,n+1))
cl_un2 = np.loadtxt(calsky+'kappa/2048/kappa_2048_powerspec_24_l_bef.txt',usecols=np.arange(1,n+1))
cl_l = np.zeros((LMAX-LMIN,n))
cl_theo = np.zeros((LMAX-LMIN,n))
cl_un = np.zeros((LMAX-LMIN,n))
cl_u2 = np.zeros((LMAX-LMIN,n))
cl_th2 = np.zeros((LMAX-LMIN,n))

for i in range(n):
    #cl_un[:,i] = cl_un1[LMIN:LMAX,i]
    print(i)
    cl_l[:,i] = cl_l1[LMIN:LMAX,i]
    cl_theo[:,i] = cl_theo1[LMIN:LMAX,i]
    cl_un[:,i] = cl_un1[LMIN:LMAX,i]
    cl_u2[:,i] = cl_un2[LMIN:LMAX,i]
    cl_th2[:,i] = cl_theo2[LMIN:LMAX,i]


cl_l = hp.anafast(hp.read_map(calsky+'lensed/2048/cib_totallensed_field_2048.fits'))[LMIN:LMAX]
cl_un = hp.anafast(hp.read_map(calsky+'cib/2048/cib_totalfield_2048.fits'))[LMIN:LMAX]
cl_theo = (np.loadtxt(pathl+'cmb_theo_august.txt')[1:])[LMIN:LMAX]
fig = plt.figure()#figsize=(22,14))

win = 201
deg = 5
width = 0.7
N = len(zmax)
l = np.arange(LMIN,LMAX)
mi,ma=10,0
ax_1 = []
ax_2 = []

row,col = 1,1
mx = 0
mn =10

lobs = np.asarray([53,114,187,320,502,684,890,1158,1505,1956,2649])
print(lobs)
cobs = np.asarray([2.36e5,1.24e5,6.63e4,3.34e4,1.91e4,1.25e4,9.17e3,6.83e3,5.34e3,4.24e3,3.42e3])
err = np.asarray([1.37e5,0.21e5,0.51e4,0.12e4,0.04e4,0.03e4,0.17e3,0.1e3,0.06e3,0.04e3,0.03e3])
ind = np.where((lobs >= LMIN) & (lobs <= LMAX))
maxp = cobs+err
minp = cobs-err

#row,col= 1,1

# FRACTIONAL DIFFERENCE 1 - NORMAL 2
plot = 2
n=1
if plot == 1:
    for i in range(n):
        ax2 = plt.subplot2grid((row,col),(i//col,i%col))
        ax_2.append(ax2)
        
        a,b,c = cl_theo,cl_un,cl_l
        
        cl = sav.savitzky_golay((a-b)/b,win,deg)
        #cl = (a-b)/b
    

        clt = sav.savitzky_golay((c-b)/b,win,deg)
       # clt = (c-b)/b
        ax2.loglog(l,np.abs(clt),label='simulation')
        ax2.loglog(lobs[ind],(maxp[ind]-minp[ind])/cobs[ind],label='Planck 2013 XXX')
        #ax2.plot(l,cl,label='theory')
        
        ax2.legend()
        ax2.set_xlabel('$l$')
        ax2.set_ylabel(r'$\Delta C_{l} / C_{l}$')
        #ax2.set_ylabel('$C_{l}$')
        #ax2.set_title(r'$\kappa$ : $ z_{max} = $'+str(float(zmax[i])/10))

if plot == 2:
    for i in range(n):
        ax1 = plt.subplot2grid((row,col),(i//col,i%col)) #power spec
        ax1 = plt.subplot2grid((row,col),(0,0))
        ax_1.append(ax1)
        
        a,b = cl_un,cl_l

        a,b = sav.savitzky_golay(a,win,deg),sav.savitzky_golay(b,win,deg)

        #ax1.loglog(l,cl_un[:,i],label='total')
        #ax1.loglog(l,cl_l[:,i],label='halos')
    
        #ax1.loglog(l,cl_u2[:,i],label='field')
        #ax1.loglog(l,cl_theo[:,i],label='linear theory')
        #ax1.loglog(l,cl_th2[:,i],label='theory + non-linear corr.')
        #ax1.plot(l,np.gradient(cl_un[:,i],np.log(l))-2)
        ax1.loglog(l,a*(10.**12),label='unlensed')
        ax1.loglog(l,b*(10.**12),label='lensed')
        ax1.errorbar(lobs[ind],cobs[ind],err[ind],label='Planck 13 XXX')
        ax1.set_xlabel('$l$')
        ax1.set_ylabel('$C_{l} [Jy^2 / Sr]$')
       # ax1.set_ylabel('$dln(C_l) / dln(l)$')
        ax1.legend()
        #ax1.set_title(r'$\kappa$ : $z_{min}=$'+str(zmin[i])+' / $z_{max}=$'+str(zmax[i]))
        #ax1.set_title(r'$z_{max} = $'+str(float(zmax[i])/10))

fig.tight_layout()
plt.savefig(calsky+'plot/cib_totalfield_powerspectrum_sm.png',format='png')
plt.show()
#for i in range(n):
    #i = 7
    #ax1 = plt.subplot2grid((row,col),(i//col,i%col)) #power spec
    #ax1 = plt.subplot2grid((row,col),(0,0))
    #ax_1.append(ax1)
    
    #ax1.loglog(l,cl_un[:,i],label='total')
    #ax1.loglog(l,cl_l[:,i],label='halos')
    
    #ax1.loglog(l,cl_u2[:,i],label='field')
    #ax1.loglog(l,cl_theo[:,i],label='linear theory')
    #ax1.loglog(l,cl_th2[:,i],label='theory + non-linear corr.')
    #ax1.set_xlabel('$l$')
    #ax1.set_ylabel('$C_{l}$')
    #ax1.legend()
    #ax1.set_title(r'$\kappa$ : $z_{min}=$'+str(zmin[i])+' / $z_{max}=$'+str(zmax[i]))
    #ax1.set_title(r'$z_{max} = $'+str(float(zmax[i])/10))
    #continue
    #ax2 = plt.subplot2grid((row,col),(i//col,i%col))
   # ax_2.append(ax2)
    #cl = sav.savitzky_golay(np.abs(cl_un[:,i]-cl_theo[:,i])/cl_theo[:,i],win,deg)
    #cl = cl_un[:,i]-cl_theo[:,i]/cl_theo[:,i]
    

    #clt = sav.savitzky_golay(np.abs(cl_th2[:,i]-cl_u2[:,i])/cl_th2[:,i],win,deg)
    #clt = cl_u2[:,i]-cl_th2[:,i]/cl_th2[:,i]
    #ax2.plot(l,clt,label='linear theory')
    #ax2.plot(l,cl,label='theory + non-linear corr.')
    #ax2.legend()
    #ax2.set_xlabel('$l$')
    #ax2.set_ylabel(r'$\Delta C_{l} / C_{l}$')
    #ax2.set_ylabel('$C_{l}$')
    #ax2.set_title(r'$\kappa$ : $ z_{max} = $'+str(float(zmax[i])/10))
   # ax2.set_title('Fractional difference')

    #mx_ = max(cl)
    #mn_ = min(cl)
    #if (mx_ > mx) :
    #    mx = mx_
    #if (mn_ < mn) :
    #    mn = mn_
    #mx_ = max(clt)
    #mn_ = min(clt)
    #if (mx_ > mx) :
    #    mx = mx_
    #if (mn_ < mn) :
    #    mn = mn_
#for ax in ax_2:
#    ax.set_ylim([mn,mx*10])


