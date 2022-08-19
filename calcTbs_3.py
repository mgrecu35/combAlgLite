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
f1=f1[:3]
f2=f2[:3]
f3=f3[:3]

import numpy as np
from rteSimpleModule import rte,calcz
from chaseRet import *

dmL=[]
zKuL=[]
zKaL=[]
zL=[]
pRateL=[]
pwcL=[]
npix1,npix2=0,0
umu=np.cos(53/180*3.14)
tbsL_ocean=[]
tbsL_land=[]
tb_obsL_ocean=[]
tb_obsL_land=[]
tb2_obsL_ocean=[]
tb2_obsL_land=[]
pwcDist_NN=[[]for k in range(20)]
pwcDist=[[]for k in range(20)]
for i,f11 in enumerate(f1[:]):
    fh=Dataset(f11)
    zKu=fh["FS/PRE/zFactorMeasured"][:]
    lon=fh["FS/Longitude"][:]
    lat=fh["FS/Latitude"][:]
    sfcType=fh["FS/PRE/landSurfaceType"][:]
    #stop
    pType=fh["FS/CSF/typePrecip"][:]
    pType=(pType/1e7).astype(int)
    fh_cmb=Dataset(f2[i])
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
    fh.close()
    fh_cmb.close()
    a=np.nonzero(pType==1)
    fh_gmi=Dataset(f3[i])
    gmi_lon=fh_gmi["S1/Longitude"][:]
    gmi_lat=fh_gmi["S1/Latitude"][:]
    tc=fh_gmi["S1/Tc"][:]
    gmi_lon2=fh_gmi["S2/Longitude"][:]
    gmi_lat2=fh_gmi["S2/Latitude"][:]
    tc2=fh_gmi["S2/Tc"][:]
    sat_lon=fh_gmi['S1/SCstatus/SClongitude'][:]
    dlon=(sat_lon-gmi_lon[:,110]).mean()
    if dlon>0:
        idir=1
    else:
        idir=-1
    #a=[[],[]]
    #stop
    for i0,j0 in zip(a[0],a[1]):
        xn=[0,binNodes[i0,j0,1],binNodes[i0,j0,3],\
            binNodes[i0,j0,4]]
        fint=np.interp(range(88),xn,[0,0,1,1])
        if binNodes[i0,j0,1]-binNodes[i0,j0,0]>1 and binNodes[i0,j0,0]<binNodes[i0,j0,2]:
            zL=zKu[i0,j0,2*binNodes[i0,j0,0]:2*binNodes[i0,j0,1]+2:2,:]
            bins=np.arange(envNodes[i0,j0,0],envNodes[i0,j0,-1]+1)
            temp_int=np.interp(bins,envNodes[i0,j0,:],airTemp[i0,j0,:])
            tempC=temp_int[binNodes[i0,j0,0]-envNodes[i0,j0,0]:
                           binNodes[i0,j0,1]-envNodes[i0,j0,0]+1]-273.15
            a1=np.nonzero(zL[:,0]>12)
            
            if len(a1[0])>300:
                #print('*')
                #print(zL[:,0])
                #print(zL[:,1])
                Nw,Dm,IWC=retr(zL,tempC)
                for k in range(binNodes[i0,j0,1]-19,binNodes[i0,j0,1]+1):
                    if zKu[i0,j0,k,0]>10 and zKu[i0,j0,k,1]>10 and k>=binNodes[i0,j0,0]:
                        k0=k-binNodes[i0,j0,1]+19
                        pwcDist_NN[k0].append(IWC[k-binNodes[i0,j0,0]])
                        pwcDist[k0].append(pwc[i0,j0,k])
                pwc[i0,j0,binNodes[i0,j0,0]:binNodes[i0,j0,1]+1]=IWC
                dm[i0,j0,binNodes[i0,j0,0]:binNodes[i0,j0,1]+1]=Dm*1e3
                pwc[i0,j0,binNodes[i0,j0,3]:]*=1.07
                #pwc[i0,j0,binNodes[i0,j0,1]+\
                #    1:binNodes[i0,j0,3]]=np.interp(\
                #                                   range(binNodes[i0,j0,1]+1,\
                #                                         binNodes[i0,j0,3]),\
                #                                   [binNodes[i0,j0,1],binNodes[i0,j0,3]],\
                #                                   [pwc[i0,j0,binNodes[i0,j0,1]],\
                #                                    pwc[i0,j0,binNodes[i0,j0,3]]])
                #dm[i0,j0,binNodes[i0,j0,1]+\
                #    1:binNodes[i0,j0,3]]=np.interp(\
                #                                   range(binNodes[i0,j0,1]+1,\
                #                                         binNodes[i0,j0,3]),\
                #                                   [binNodes[i0,j0,1],binNodes[i0,j0,3]],\
                #                                   [dm[i0,j0,binNodes[i0,j0,1]],\
                #                                    dm[i0,j0,binNodes[i0,j0,3]]])
        zKu1D,zKa1D,pRate,kextH1,asymH1,salbH1=calcz(pwc,dm,cAlg,i0,j0,fint,sfcBin,binNodes)
        #if zKu[i0,j0,binNodes[i0,j0,-1],0]>10 and\
       #    pwc[i0,j0,binNodes[i0,j0,-1]]>1e-3:
        pwcL.append([pwc[i0,j0,binNodes[i0,j0,-1]],\
                     pwc[i0,j0,binNodes[i0,j0,-1]-1]])
        sfcEmiss[i0,j0,7:8]*=1.05
        sfcEmiss[i0,j0,7:8][sfcEmiss[i0,j0,7:8]>0.99]=0.99
        tb13=rte(binNodes,dm,pwc,kextH1,asymH1,salbH1,i0,j0,sfcEmiss,qv,airTemp,\
                 press,envNodes,skTemp,cldw,umu,cAlg)
      
        swp=((1-fint)*pwc[i0,j0,:]).sum()
        rwp=(fint*pwc[i0,j0,:]).sum()
        if sfcType[i0,j0]==0:
            npix1+=1
            tb13.append(swp)
            tb13.append(rwp)
            tbsL_ocean.append(tb13)
        else:
            npix2+=1
            tb13.append(swp)
            tb13.append(rwp)
            tbsL_land.append(tb13)

    nx=lon.shape[0]
    tb2d=np.zeros((nx,49,13),float)
    
        
    
    swath_def = geometry.SwathDefinition(lons=gmi_lon, lats=gmi_lat)
    swath_def2 = geometry.SwathDefinition(lons=gmi_lon2, lats=gmi_lat2)
    target_def = geometry.SwathDefinition(lons=lon[:,:], lats=lat[:,:])
    tc_regrid = kd_tree.resample_gauss(swath_def, tc[:,:,:],
                                       target_def, radius_of_influence=25000, \
                                       sigmas=[12500 for k in range(9)])
    tc_regrid2 = kd_tree.resample_gauss(swath_def2, tc2[:,:,:],
                                        target_def, radius_of_influence=25000, \
                                    sigmas=[12500 for k in range(4)])
    b=np.nonzero(sfcType[a]==0)
    tb_obsL_ocean.extend(tc_regrid[a][b])
    tb2_obsL_ocean.extend(tc_regrid2[a][b])
    b=np.nonzero(sfcType[a]!=0)
    tb_obsL_land.extend(tc_regrid[a][b])
    tb2_obsL_land.extend(tc_regrid2[a][b])
    #break
    #continue
    n1,n2=250,553 #45298
    n1,n2=200,400 #45313
    tb2d_NN=np.zeros((nx,49,13),float)
    for i0 in range(n1,n2):
        for j0 in range(49):
            xn=[0,binNodes[i0,j0,1],binNodes[i0,j0,3],\
                binNodes[i0,j0,4]]
            fint=np.interp(range(88),xn,[0,0,1,1])
            if pType[i0,j0]>0:
                zKu1D,zKa1D,pRate,kextH1,asymH1,salbH1=calcz(pwc,dm,cAlg,i0,j0,\
                                                             fint,sfcBin,binNodes)
            else:
                kextH1=np.zeros((88,8),float)
                salbH1=np.zeros((88,8),float)
                asymH1=np.zeros((88,8),float)
            tb13=rte(binNodes,dm,pwc,kextH1,asymH1,salbH1,i0,j0,sfcEmiss,qv,airTemp,\
                         press,envNodes,skTemp,cldw,umu,cAlg)
            swp=((1-fint)*pwc[i0,j0,:]).sum()
            tb2d[i0,j0,:]=tb13

    for i0 in range(n1,n2):
        for j0 in range(49):
            xn=[0,binNodes[i0,j0,1],binNodes[i0,j0,3],\
                binNodes[i0,j0,4]]
            fint=np.interp(range(88),xn,[0,0,1,1])
            if pType[i0,j0]>0:
                embed_NN(binNodes,envNodes,airTemp,zKu,pwc,dm,i0,j0)
                zKu1D,zKa1D,pRate,kextH1,asymH1,salbH1=calcz(pwc,dm,cAlg,i0,j0,\
                                                             fint,sfcBin,binNodes)
            else:
                kextH1=np.zeros((88,8),float)
                salbH1=np.zeros((88,8),float)
                asymH1=np.zeros((88,8),float)
            tb13=rte(binNodes,dm,pwc,kextH1,asymH1,salbH1,i0,j0,sfcEmiss,qv,airTemp,\
                         press,envNodes,skTemp,cldw,umu,cAlg)
            swp=((1-fint)*pwc[i0,j0,:]).sum()
            tb2d_NN[i0,j0,:]=tb13
            

    break
#



plt.figure()
plt.subplot(121)
plt.pcolormesh(lon[n1:n2,:],lat[n1:n2,:],tc_regrid[n1:n2,:,8],cmap='jet',vmin=200,vmax=280)
plt.subplot(122)
plt.pcolormesh(lon[n1:n2,:],lat[n1:n2,:],tb2d[n1:n2,:,8],cmap='jet',vmin=200,vmax=280)
plt.contour(lon[n1:n2,:],lat[n1:n2,:],sfcPrecip[n1:n2,:])
plt.figure()
#plt.subplot(121)
plt.pcolormesh(gmi_lon[:],gmi_lat[:],tc[:,:,5],cmap='jet')
plt.contour(lon,lat,sfcPrecip)

zKum=np.ma.array(zKu,mask=zKu<10)
fig=plt.figure(figsize=(8,10))
plt.subplot(311)
zim=plt.pcolormesh(range(n1,n2),np.arange(176)*0.125,\
               zKum[n1:n2,24,::-1,0].T,cmap='jet',vmin=10,vmax=45)
plt.ylabel("Height [km]")
plt.title("Observed Ku-band reflectivity")
plt.ylim(0,10.0)
plt.colorbar()
plt.subplot(312)
tp1=plt.plot(range(n1,n2),tc_regrid[n1:n2,24,8])
plt.plot(range(n1,n2),tb2d[n1:n2,24,8])
plt.plot(range(n1,n2),tb2d_NN[n1:n2,24,8])
plt.ylabel("Tb [K]")
plt.legend(['89-GHz H obs.','89-GHz H sim_1','89-GHz H sim_2'])
plt.xlim(n1,n2)
plt.subplot(313)
tp2=plt.plot(range(n1,n2),tc_regrid2[n1:n2,24,1])
plt.plot(range(n1,n2),tb2d[n1:n2,24,10])
plt.plot(range(n1,n2),tb2d_NN[n1:n2,24,10])
plt.legend(['166-GHz H obs.','166-GHz H sim_1','166-GHz H sim_2'])

plt.xlabel("Scan number")
plt.ylabel("Tb [K]")
plt.xlim(n1,n2)
plt.tight_layout()
l, b, w, h = zim.axes.get_position().bounds
l2, b2, w2, h2 = tp1[0].axes.get_position().bounds
tp1[0].axes.set_position([l2,b2,w,h2])
l3, b3, w3, h3 = tp2[0].axes.get_position().bounds
tp2[0].axes.set_position([l3,b3,w,h3])
fname_out='crossSect'+f11[37:-4]+'png'
plt.savefig(fname_out)
stop
tbsL_land=np.array(tbsL_land)
tb_obsL_land=np.array(tb_obsL_land)
print(np.corrcoef(tbsL_land[:,8],tb_obsL_land[:,-1]))
import matplotlib

from sklearn.neighbors import KNeighborsRegressor
X=np.array(tbsL_land)[:,5:-2]
y=np.array(tbsL_land)[:,-2:]
from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(X, y, \
                                                    test_size=0.33, random_state=42)
neigh = KNeighborsRegressor(n_neighbors=60,weights="distance")
neigh.fit(X_train, y_train)

tb2_obsL_land=np.array(tb2_obsL_land)
X_obs=[]

for i in range(len(tb_obsL_land)):
    x=[]
    x.extend(tb_obsL_land[i,5:])
    x.extend(tb2_obsL_land[i,:])
    X_obs.append(x)
    
X_obs=np.array(X_obs)
import matplotlib
