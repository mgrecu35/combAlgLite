import tensorflow as tf
import pickle
import numpy as np
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
d=pickle.load(open("Chase_et_al_2021_NN-master/scalerXy.pklz","rb"))
model = tf.keras.models.load_model('Chase_et_al_2021_NN-master/' + 'NN_6by8.h5',compile=True)
xm=d["scaler_X_mean_scale"][0]
xs=d["scaler_X_mean_scale"][1]
ym=d["scaler_y_mean_scale"][0]
ys=d["scaler_y_mean_scale"][1]
def retr(zL,temp):
    zL=np.array(zL)
    nz=zL.shape[0]
    X=np.zeros((nz,3),float)
    X[:,0]=zL[:,0]
    X[:,1]=zL[:,0]-zL[:,1]
    X[:,0][X[:,0]<0]=0
    X[:,1][X[:,1]<0]=0
    X[:,1][X[:,1]>10]=10
    for k in range(2):
        X[:,k]=(X[:,k]-xm[k])/xs[k]
        X[:,2]=(temp-xm[2])/xs[2]
    yhat = model.predict(X[:,0:3],verbose=0)
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
    
    return Nw,Dm,IWC

from bisectm import bisectm

def embed_NN(binNodes,envNodes,airTemp,zKu,pwc,dm,i0,j0,cAlg,\
             pwc1L,pwc2L,pwc3L,dm1L,dm2L,dm3L):
    if binNodes[i0,j0,1]-binNodes[i0,j0,0]>1 and \
       binNodes[i0,j0,0]<binNodes[i0,j0,2]:
        zL=zKu[i0,j0,2*binNodes[i0,j0,0]:2*binNodes[i0,j0,1]+2:2,:]
        bins=np.arange(envNodes[i0,j0,0],envNodes[i0,j0,-1]+1)
        temp_int=np.interp(bins,envNodes[i0,j0,:],airTemp[i0,j0,:])
        tempC=temp_int[binNodes[i0,j0,0]-envNodes[i0,j0,0]:
                       binNodes[i0,j0,1]-envNodes[i0,j0,0]+1]-273.15
        a1=np.nonzero(zL[:,0]>12)
        
        if len(a1[0])>3:
            Nw,Dm,IWC=retr(zL,tempC)
            dwr1=zL[-1,0]-zL[-1,1]
            zku1=zL[-1,0]
            if zL[-1,0]>15 and zL[-1,1]>15 and dwr1>1 and dwr1<9:
                dwr=(cAlg.tablep2.zkus-cAlg.tablep2.zkas)[:253]
                ibin=bisectm(dwr,253,dwr1)
                dn=(zku1-cAlg.tablep2.zkus[ibin])/10.
                pwc1=(cAlg.tablep2.swc[ibin])*10**dn
                pwc1L.append(pwc1)
                pwc2L.append(IWC[-1])
                pwc3L.append(pwc[i0,j0,binNodes[i0,j0,1]])
                dm1L.append(cAlg.tablep2.dms[ibin])
                dm2L.append(Dm[-1])
                dm3L.append(dm[i0,j0,binNodes[i0,j0,1]])
            pwc[i0,j0,binNodes[i0,j0,0]:binNodes[i0,j0,1]+1]=IWC
            dm[i0,j0,binNodes[i0,j0,0]:binNodes[i0,j0,1]+1]=Dm*1e3

            
