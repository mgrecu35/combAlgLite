from numpy import *
import numpy as np
from bisectm import *
from numba import jit
@jit(nopython=True)
def calc_pia(att,pia,dr,n):
    for k in range(n):
        if k==0:
            pia[k]=att[k]*dr
        else:
            pia[k]=pia[k-1]+(att[k]+att[k-1])*dr
            
def calcz(pwc,dm,cAlg,i0,j0,fint):
    swc=pwc[i0,j0,:]*(1-fint)
    rwc=pwc[i0,j0,:]*fint
    a1=np.nonzero(swc>0)
    zKus=np.zeros((88),float)-99
    zKas=np.zeros((88),float)-99
    pRate=np.zeros((88),float)
    att_ka=np.zeros((88),float)
    dr=0.25
    n=88
    pia_ka=np.zeros((88),float)
    for ik,swc1 in enumerate(swc[a1]):
        ibin=bisectm(cAlg.tablep2.dms[:253],253,dm[i0,j0,a1[0][ik]])
        swc_bin=cAlg.tablep2.swc[ibin]
        dn=np.log10(swc1/swc_bin)
        zKus[a1[0][ik]]=cAlg.tablep2.zkus[ibin]+10*dn
        zKas[a1[0][ik]]=cAlg.tablep2.zkas[ibin]+10*dn
        pRate[a1[0][ik]]=cAlg.tablep2.snowrate[ibin]*10**dn
        att_ka[a1[0][ik]]=cAlg.tablep2.attkas[ibin]*10**dn
    zKur=np.zeros((88),float)-99
    zKar=np.zeros((88),float)-99
    a1=np.nonzero(rwc>0)
    for ik,rwc1 in enumerate(rwc[a1]):
        ibin=bisectm(cAlg.tablep2.dmr[:289],289,dm[i0,j0,a1[0][ik]])
        rwc_bin=cAlg.tablep2.rwc[ibin]
        dn=np.log10(rwc1/rwc_bin)
        zKur[a1[0][ik]]=cAlg.tablep2.zkur[ibin]+10*dn
        zKar[a1[0][ik]]=cAlg.tablep2.zkar[ibin]+10*dn
        pRate[a1[0][ik]]+=cAlg.tablep2.rainrate[ibin]*10**dn
        att_ka[a1[0][ik]]+=cAlg.tablep2.attkar[ibin]*10**dn
    zKu1D=np.log10(10**(0.1*zKur)+10**(0.1*zKus))*10
    zKa1D=np.log10(10**(0.1*zKar)+10**(0.1*zKas))*10
    calc_pia(att_ka,pia_ka,dr,n)
    zKa1D-=pia_ka

    return zKu1D,zKa1D,pRate
def rte(binNodes,dm,pRate,i,j,emiss,qv,airTemp,\
        press,envNode,sfcTemp,cldw,umu,cAlg):
    freqs=[10.6,10.6,18.7,18.7,23.,37,37.,89,89.,166.,166.,186.3,190.3]
    r1L=[]
    r2L=[]
    kL=[]
    npol=[1,0,1,0,1,1,0,1,0,1,0,1,1]
    iFreq=[0,0,1,1,2,3,3,4,4,5,5,6,6,7,7]
    nfreq=8
    kextH=np.zeros((88,5,8),float)
    salbH=np.zeros((88,5,8),float)
    asymH=np.zeros((88,5,8),float)
    tbSim=[]
    iEnum=0
    for (pol,f) in zip(npol,freqs):
        kextL=[]
        for q,tk,pa,cld in zip(qv[i,j,:],airTemp[i,j,:],press[i,j,:],cldw[i,j,:]):
            if q<1e-3:
                q=1e-3
            absair,abswv=cAlg.gasabsr98(f,tk,q*1e-3,pa*1e2,1)
            if cld>0.001:
                z_clw=cAlg.gcloud(f,tk,cld)
            else:
                z_clw=0.
            kextL.append(absair+abswv+z_clw)
        bins=np.arange(envNode[i,j,0],envNode[i,j,-1]+1)
        kextInt=np.interp(bins,envNode[i,j,:],kextL)
        kextInt_copy=kextInt.copy()
        #kextInt=kextInt+kextH1[envNode[i,j,0]:envNode[i,j,-1]+1,iFreq[iEnum]]
        salb=kextInt.copy()*0.
        asym=kextInt.copy()*0.
        #salb=salbH1[envNode[i,j,0]:envNode[i,j,-1]+1,iFreq[iEnum]].copy()*\
        #    (1-kextInt_old/kextInt)
        #asym=asymH1[envNode[i,j,0]:envNode[i,j,-1]+1,iFreq[iEnum]].copy()
        tLayer=list(np.interp(bins,envNode[i,j,:],airTemp[i,j,:]))
        tLayer.append(sfcTemp[i,j])
        iEnum+=1
        emis=emiss[i,j,iEnum-1]
        ebar=emis
        nL=kextInt.shape[0]
        lambert=False
        tb1=cAlg.radtran(umu,sfcTemp[i,j],tLayer[::-1],\
                          np.arange(nL+1)*0.25,kextInt[::-1],\
                          salb[::-1],asym[::-1],2.7,emis,ebar,lambert)

        tbSim.append(tb1)
    return tbSim
        
