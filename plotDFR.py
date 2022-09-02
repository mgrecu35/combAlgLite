fname="Chase_et_al_2021_NN-master/Unrimed_simulation_wholespecturm_train_V2.nc"
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib
import read_tables as cAlg
matplotlib.use('TkAgg')

cAlg.read_tables()
cAlg.initp2()
import numpy as np
nc=Dataset(fname)
temp=nc["T_env"][:]
a=np.nonzero(temp>-10)
zKa=nc["Z2"][:][a]
zKu=nc["Z"][:][a]
dm=nc["Dm"][:][a]*1e3
Nw=nc["Nw"][:][a]

plt.hist2d(dm,zKu-zKa,bins=(np.arange(40)*0.05,np.arange(40)*0.3-1),cmap='jet',norm=matplotlib.colors.LogNorm())
plt.plot(cAlg.tablep2.dms[0:253],cAlg.tablep2.dwr[0:253])
plt.xlabel("Dm (mm)")
plt.ylabel("DFR (dB)")
c=plt.colorbar()
c.ax.set_title("Counts")
plt.savefig("dfr_dm_distribution.png")

plt.figure()
plt.hist2d(dm,np.log10(Nw),bins=(np.arange(40)*0.05,np.arange(40)*0.2+4),cmap='jet',norm=matplotlib.colors.LogNorm())
plt.xlabel("Dm (mm)")
plt.ylabel("log10(Nw)")
c=plt.colorbar()
c.ax.set_title("Counts")
plt.savefig("nw_dm_distribution.png")

plt.figure()
zL=[]
dmL=[]
dm_stdL=[]
for i in range(40):
    a=np.nonzero((zKu-i)*(zKu-i-1)<0)
    if len(a[0])>0:
        zL.append(zKu[a].mean())
        dmL.append(dm[a].mean())
        dm_stdL.append(dm[a].std())



dmL=np.array(dmL)
dm_stdL=np.array(dm_stdL)
plt.figure(figsize=(8,8))
plt.subplot(211)
plt.plot(zL,dmL)
plt.plot(cAlg.tablep2.zkus[:253],cAlg.tablep2.dms[:253])
plt.ylabel('Dm mean (mm)')
plt.legend(['NN tables','CMB'])
plt.subplot(212)
plt.plot(zL,dm_stdL)
plt.ylabel('Dm std. (mm)')
plt.xlabel('Z(Ku) (dBZ)')
plt.savefig('zKu_dm.png')

for z1,d1,ds in zip(zL,dmL,dm_stdL):
    print(z1,d1,ds)
