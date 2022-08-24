import read_tables as cAlg
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from pyresample import kd_tree, geometry
cAlg.read_tables()
cAlg.initp2()
import glob
f1=sorted(glob.glob("Data2/2A*"))
f2=sorted(glob.glob("Data2/2B*"))
f3=sorted(glob.glob("Data2/1C*"))
f1=f1[1:3]
f2=f2[1:3]
f3=f3[1:3]

import numpy as np
from rteSimpleModule import rte,calcz
#from chaseRet import *


npix1,npix2=0,0
umu=np.cos(53/180*3.14)
def readData(f1_,f2_,f3_):
    fh=Dataset(f1_)
    zKu=fh["FS/PRE/zFactorMeasured"][:]
    lon=fh["FS/Longitude"][:]
    lat=fh["FS/Latitude"][:]
    sfcType=fh["FS/PRE/landSurfaceType"][:]
    pType=fh["FS/CSF/typePrecip"][:]
    pType=(pType/1e7).astype(int)
    fh_cmb=Dataset(f2_)
    qv=fh_cmb["KuKaGMI/vaporDensity"][:]
    press=fh_cmb["KuKaGMI/airPressure"][:]
    envNodes=fh_cmb["KuKaGMI/envParamNode"][:]
    airTemp=fh_cmb["KuKaGMI/airTemperature"][:]
    skTemp=fh_cmb["KuKaGMI/skinTemperature"][:]
    binNodes=fh_cmb["KuKaGMI/phaseBinNodes"][:]
    bcf=fh["FS/PRE/binClutterFreeBottom"][:]
    bsf=fh["FS/PRE/binRealSurface"][:]
    pwc=fh_cmb["KuKaGMI/precipTotWaterCont"][:]
    sfcEmiss=fh_cmb["KuKaGMI/surfEmissivity"][:]
    oe_sfcEmiss=fh_cmb["KuKaGMI/OptEst/OEsurfEmissivity"][:]
    dm=fh_cmb["KuKaGMI/precipTotDm"][:]
    cldw=fh_cmb["KuKaGMI/cloudLiqWaterCont"][:]
    sfcBin=fh_cmb["KuKaGMI/Input/surfaceRangeBin"][:,:,0]
    sfcPrecip=fh["FS/SLV/precipRateNearSurface"][:]
    bbPeak=fh["FS/CSF/binBBPeak"][:]
    fh.close()
    fh_cmb.close()
    a=np.nonzero(pType==1)
    fh_gmi=Dataset(f3_)
    gmi_lon=fh_gmi["S1/Longitude"][:]
    gmi_lat=fh_gmi["S1/Latitude"][:]
    tc=fh_gmi["S1/Tc"][:]
    gmi_lon2=fh_gmi["S2/Longitude"][:]
    gmi_lat2=fh_gmi["S2/Latitude"][:]
    tc2=fh_gmi["S2/Tc"][:]
    sat_lon=fh_gmi['S1/SCstatus/SClongitude'][:]
    dlon=(sat_lon-gmi_lon[:,110]).mean()
    swath_def = geometry.SwathDefinition(lons=gmi_lon, lats=gmi_lat)
    swath_def2 = geometry.SwathDefinition(lons=gmi_lon2, lats=gmi_lat2)
    target_def = geometry.SwathDefinition(lons=lon[:,:], lats=lat[:,:])
    tc_regrid = kd_tree.resample_gauss(swath_def, tc[:,:,:],
                                       target_def, radius_of_influence=25000, \
                                       sigmas=[12500 for k in range(9)])
    tc_regrid2 = kd_tree.resample_gauss(swath_def2, tc2[:,:,:],
                                        target_def, radius_of_influence=25000, \
                                    sigmas=[12500 for k in range(4)])

    umu=np.cos(53/180.0*np.pi)
    bsf1=bsf[:,:,0]
    n1=qv.shape[0]
    print(envNodes[132,14,:])
    print(airTemp.shape)
    print(envNodes[:,:,0].min(),envNodes[:,:,0].max())
    print(cldw.shape)
    print(umu)
    print(skTemp.shape)
    print(sfcEmiss.shape)
    a=np.nonzero(pType>0)
    print(a[0][0],a[1][0])
        
    return qv,press,sfcType,pType,airTemp,envNodes,binNodes,\
        bcf,bsf,pwc,sfcEmiss,dm,skTemp,lon,lat,tc_regrid,tc_regrid2,zKu,cldw,bbPeak


f1_=f1[0]
f2_=f2[0]
f3_=f3[0]
#stop

qv,press,sfcType,pType,airTemp,envNodes,binNodes,\
    bcf,bsf,pwc,sfcEmiss,dm,skTemp,lon,lat,tc_regrid,\
    tc_regrid2,zKu,cldw,bbPeak=\
        readData(f1_,f2_,f3_)
bbPeak=(bbPeak/2).astype(int)
zKuL=[]
tb13L=[]
tb13L_py=[]
bsf=(bsf[:,:,0]/2).astype(int)
import time
t3=time.time()
for i in range(200,400):
    for j in range(24,25):
        fint = cAlg.get_fint(bsf[i,j],binNodes[i,j,:])
        kexth1,salbh1,asymh1,zku1,zka1 = cAlg.calcz1(pwc[i,j,:],dm[i,j,:],fint,bsf[i,j],\
                                                     binNodes[i,j,:],bbPeak[i,j])
        #if bbPeak[i,j]>0:
        #    stop
        tb13=cAlg.rte1(binNodes[i,j,:],bsf[i,j],dm[i,j,:],\
                       kexth1,asymh1,salbh1,sfcEmiss[i,j,:],qv[i,j,:],\
                       airTemp[i,j,:],press[i,j,:],envNodes[i,j,:],skTemp[i,j],\
                       cldw[i,j,:],umu)
        #tb13_py=rte(binNodes,dm,pwc,kexth1,asymh1,salbh1,i,j,sfcEmiss,qv,airTemp,\
        #    press,envNodes,skTemp,cldw,umu,cAlg)
        #print(i,len(tb13L))
        if j==24:
            tb13L.append(tb13)
            tb13L_py.append(tb13)
            zKuL.append(zku1)

#stop
t4=time.time()
print('wall clock python t=%6.3f'%(t4-t3))
t1=time.time()
tbout2d = cAlg.calc_tb_f90(binNodes,pwc,dm,bsf,pType,envNodes,qv,airTemp,press,\
                           cldw,umu,skTemp,sfcEmiss,bbPeak)

tbout2d_f = cAlg.calc_tb_f90_flattened(binNodes[200:400,24,:],pwc[200:400,:],dm[200:400,:],\
                                       bsf[200:400,24],pType[200:400,24],envNodes[200:400,24,:],\
                                       qv[200:400,24,:],airTemp[200:400,24,:],press[200:400,24,:],\
                                       cldw[200:400,24,:],umu,skTemp[200:400,24],sfcEmiss[200:400,24,:],\
                                       bbPeak[200:400,24])
t2=time.time()
print('wall clock t=%6.3f'%(t2-t1))

zKuL=np.array(zKuL)
plt.subplot(211)
plt.pcolormesh(zKuL[:,::-1].T,vmin=10,vmax=45,cmap='jet')
plt.ylim(0,30)
plt.subplot(212)
plt.pcolormesh(zKu[200:400,24,::-1,0].T,vmin=10,vmax=45,cmap='jet')
plt.ylim(0,60)
plt.figure()
tb13L=np.array(tb13L)
tb13L_py=np.array(tb13L_py)
plt.plot(tc_regrid[200:400,24,8])
plt.plot(tb13L[:,8])
#plt.plot(tbout2d[200:400,24,8])
plt.show()
#plt.plot(tb13L_py[:,8])
