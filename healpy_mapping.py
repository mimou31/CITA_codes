import healpy as hp
import numpy as np
from types import *
import matplotlib
#matplotlib.use('agg')
import matplotlib.pyplot as plt
import concurrent.futures
import multiprocessing
from  multiprocessing import Process,Pool,Manager,Value,Array
from itertools import product
import time
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys
def is_power2(num):
    return num != 0 and ((num & (num-1)) == 0)

def lowToHigh(m1,m2,pix):
    if(type(pix)!=int):
        print('TypeError: pix is not integer.')
        return
    n1 = hp.get_nside(m1)
    n2 = hp.get_nside(m2)

    if(is_power2(n1)==False):
        print("Nside 1 is not power of 2.")
        return
    if(is_power2(n2)==False):
        print("Nside 2 is not power of 2.")
        return

    r = n2/n1
    pixn = hp.ring2nest(n1,pix)
    ind = np.arange(r**2*pixn,r**2*(pixn+1))
    
    return(hp.nest2ring(n2,ind))

def highToLow(m1,m2,pix):
    if(type(pix)!=int):
        print('TypeError: pix is not integer.')
        return
    n1 = hp.get_nside(m1)
    n2 = hp.get_nside(m2)

    if(is_power2(n1)==False):
        print("Nside 1 is not power of 2.")
        return
    if(is_power2(n2)==False):
        print("Nside 2 is not power of 2.")
        return
    
    r = n1/n2
    pixn = hp.ring2nest(n1,pix)

    return hp.nest2ring(pixn//(r**2))



NK = 32
Np = multiprocessing.cpu_count()
#Np = 1
power = 3

print("Values of kappa : "+ str(hp.nside2npix(NK)) +" for "+str(Np)+" proc.")

kappaCMB = hp.read_map('/scratch/r/rbond/remim/mocks/cib-maps/kapcmb-half.fits')
kappa=hp.read_map('/scratch/r/rbond/remim/cal-sky/test/kappa_maps/kap_param.kappa_zmin0.0_zmax4.8_kap.fits')
ns = hp.get_nside(kappa)
kap_down = hp.ud_grade(kappa,NK)
cibG = hp.read_map('/scratch/r/rbond/remim/mocks/cib-maps/cib_sim_Pois_1024-lextend.fits')

cib = hp.read_map('/scratch/r/rbond/remim/cal-sky/test/cib_maps/sum/cib_ns2046_zmin_0.0_zmax_4.8.fits')
#cib = hp.read_map('/scratch/r/rbond/remim/mocks/cib-maps/cib_calsky_0-4.8_lensed_ip6.fits')
N = len(kap_down)


ti = time.time()

l = np.arange(6144)

def anafast_patch(i,ps,kap):
    t = time.time()
    #print "Starting " +str(i+1)+ " / " +str(N)
    n1 = hp.get_nside(kap_down)
    n2 = hp.npix2nside(len(cib))
    r = n2/n1
    pixn = hp.ring2nest(n1,i)
    ind = np.arange(r**2*pixn,r**2*(pixn+1))
    
    ind = hp.nest2ring(n2,ind)
    ps[i] = np.mean(cib[ind]**power) - (np.mean(cib[ind]))**power
    kap[i] = np.mean(kappa[ind])
    #print(ps[i])
    #print ps
    #print "Finished " +str(i+1)+ " / " +str(N) +" in " + str(time.time() - t) + " s"

man = Manager()
ps = man.list(range(N))
kap = man.list(range(N))

#p = Pool(Np)
#p.map(anafast_patch,range(N))
#p.close()
#p.join()
#p = Pool(Np)

jobs = []
for i in range(N):
    p = multiprocessing.Process(target=anafast_patch, args=(i,ps,kap))
    jobs.append(p)
    p.daemon = True
    p.start()

for proc in jobs:
    proc.join()


print("Finished all in " + str(time.time()-ti)+" s")
print(max(ps),min(ps))

fig,ax = plt.subplots()
ax.scatter(kap,ps,alpha=0.2)
#plt.yscale('log')

#ax.set_xlabel(r'$\langle \kappa \rangle$')
ax.set_ylabel(r'$\sigma(I^{CIB})_{patch}$')
ax.set_xlabel(r'$\sigma(\kappa^{'+str(power)+'})_{patch}$')
#ax.set_ylabel(r'$\langle I_{patch} \rangle$')
ax.ticklabel_format(axis='both',style='sci', scilimits=(0,0))
plt.savefig('/scratch/r/rbond/remim/mocks/i_patches_ns'+str(NK)+'_mean'+str(power)+'.png',format='png')
plt.show()
