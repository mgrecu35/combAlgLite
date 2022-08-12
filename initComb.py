import read_tables as cAlg
import matplotlib.pyplot as plt
from netCDF4 import Dataset
cAlg.read_tables()
cAlg.initp2()
import glob
f1=sorted(glob.glob("Data/2A*"))
f2=sorted(glob.glob("Data/2B*"))
f3=sorted(glob.glob("Data/1C*"))
import numpy as np
from rteSimpleModule import rte,calcz

dmL=[]
zKuL=[]
zKaL=[]
zL=[]
pRateL=[]
for i,f11 in enumerate(f1):
    fh=Dataset(f11)
    zKu=fh["FS/PRE/zFactorMeasured"][:]
    lon=fh["FS/Longitude"][:]
    lat=fh["FS/Latitude"][:]
    sfcType=fh["FS/PRE/landSurfaceType"][:]
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
    dm=fh_cmb["KuKaGMI/precipTotDm"][:]
    cldw=fh_cmb["KuKaGMI/cloudLiqWaterCont"][:]
    fh.close()
    fh_cmb.close()
    a=np.nonzero(pType==1)
    for i0,j0 in zip(a[0],a[1]):
        if dm[i0,j0,binNodes[i0,j0,1]]>0 and dm[i0,j0,binNodes[i0,j0,3]]>0:
           
            #if dm[i0,j0,binNodes[i0,j0,0]]==dm[i0,j0,binNodes[i0,j0,:]].min():
            #    break
            xn=[0,binNodes[i0,j0,1],binNodes[i0,j0,3],\
                binNodes[i0,j0,4]]
            fint=np.interp(range(88),xn,[0,0,1,1])

            if fint[binNodes[i0,j0,0]]<0.01:
                zKu1D,zKa1D,pRate=calcz(pwc,dm,cAlg,i0,j0,fint)
                if zKu[i0,j0,2*binNodes[i0,j0,1],0]>0 and \
                   zKu[i0,j0,2*binNodes[i0,j0,1],1]>0 and\
                    zKu[i0,j0,2*binNodes[i0,j0,3],0]>0 and \
                    zKu[i0,j0,2*binNodes[i0,j0,3],1]>0:
                    dmL.append([dm[i0,j0,binNodes[i0,j0,1]],dm[i0,j0,binNodes[i0,j0,3]]])
                    zL.append([zKu[i0,j0,2*binNodes[i0,j0,1],0],\
                               zKu[i0,j0,2*binNodes[i0,j0,1],1],\
                               zKu1D[binNodes[i0,j0,1]],\
                               zKa1D[binNodes[i0,j0,1]],\
                               zKu[i0,j0,2*binNodes[i0,j0,3],0],\
                               zKu[i0,j0,2*binNodes[i0,j0,3],1],\
                               zKu1D[binNodes[i0,j0,3]],\
                               zKa1D[binNodes[i0,j0,3]]])
                    pRateL.append([pRate[binNodes[i0,j0,1]],\
                                   pRate[binNodes[i0,j0,3]],
                                   pwc[i0,j0,binNodes[i0,j0,1]]])
                #break
            #    print(fint[binNodes[i0,j0,:]],binNodes[i0,j0,:])
            #    print(xn)
    #break
    #print(ibin,swc1,swc_bin)

import pickle

import tensorflow as tf
d=pickle.load(open("scalerXy.pklz","rb"))
model = tf.keras.models.load_model('Chase_et_al_2021_NN-master/' + 'NN_6by8.h5',custom_objects=None,compile=True)
xm=d["scaler_X_mean_scale"][0]
xs=d["scaler_X_mean_scale"][1]
ym=d["scaler_y_mean_scale"][0]
ys=d["scaler_y_mean_scale"][1]
zL=np.array(zL)
X=zL[:,0:3]
X[:,1]=X[:,0]-X[:,1]
for k in range(2):
    X[:,k]=(X[:,k]-xm[k])/xs[k]
X[:,2]=(-3-xm[2])/xs[2]
yhat = model.predict(X[:,0:3])
#yhat = scaler_y.inverse_transform(yhat)
for k in range(3):
    yhat[:,k]=yhat[:,k]*ys[k]+ym[k]

yhat[:,1] = 10**yhat[:,1] 
Nw=yhat[:,0]
Dm=yhat[:,1]
Nw = 10**Nw #undo log, should be in m^-4
Dm = Dm/1000. # convert to m ^4
IWC = (Nw*(Dm)**4*1000*np.pi)/4**(4) # the 1000 is density of water (kg/m^3)
IWC = IWC*1000 #convert to g/m^3 
pRateL=np.array(pRateL)
umu=np.cos(53/180*3.14)
#tb13=rte(binNodes,dm,pwc,i0,j0,sfcEmiss,qv,airTemp,\
#        press,envNodes,skTemp,cldw,umu,cAlg)
import matplotlib
zL=np.array(zL)
fig=plt.figure()
ax=fig.add_subplot(111)
plt.hist2d(zL[:,1],zL[:,3],cmap='jet',norm=matplotlib.colors.LogNorm(),
           bins=np.arange(30)+10)
ax.set_aspect('equal')
fig=plt.figure()
ax=fig.add_subplot(111)
plt.hist2d(zL[:,5],zL[:,7],cmap='jet',norm=matplotlib.colors.LogNorm(),
           bins=np.arange(30)+10)
ax.set_aspect('equal')
