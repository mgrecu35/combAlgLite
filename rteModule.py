from numpy import *

def rte(psdNodes,binNodes,Nw,pRate,z13,i,j,emiss,qv,airTemp,\
        press,envNode,sfcTemp,umu,pyHB2,fract,idir,cldw):
    freqs=[10.6,10.6,18.7,18.7,23.,37,37.,89,89.,166.,166.,186.3,190.3]
    r1L=[]
    r2L=[]
    kL=[]
    npol=[1,0,1,0,1,1,0,1,0,1,0,1,1]
    iFreq=[0,0,1,1,2,3,3,4,4,5,5,6,6,7,7]
    nfreq=8
    kextH=zeros((88,5,8),float)
    salbH=zeros((88,5,8),float)
    asymH=zeros((88,5,8),float)
    tbSim=[]
    ns=psdNodes.shape[0]
    for dj in range(5):
        if idir==1:
            i1=min(ns-1,i+dj)
        else:
            i1=max(0,i+dj)
        #print(ns,i1,j,dj)
        if psdNodes[i1,j,0]>-999:
            Nwb=interp(arange(psdNodes[i1,j,0],psdNodes[i1,j,-1]+1),\
                       psdNodes[i1,j,:],Nw[i1,j,:,0])
        
            for k in range(max(0,binNodes[i1,j,0]),binNodes[i1,j,0]):
                kext,salb,asym=pyHB2.getsnowp(Nwb[0],\
                                              pRate[i1,j,binNodes[i1,j,0]],nfreq)
                kextH[k,dj,:]=kext
                salbH[k,dj,:]=salb
                asymH[k,dj,:]=asym
                
            for k in range(binNodes[i1,j,0],\
                           min(binNodes[i1,j,2],binNodes[i1,j,4])):
                kext,salb,asym=pyHB2.getsnowp(Nwb[k-psdNodes[i1,j,0]],\
                                              pRate[i1,j,k],nfreq)
                if z13[i1,j,k]>10 and binNodes[i1,j,3]>binNodes[i1,j,0]+2 and k<binNodes[i1,j,2]-1:
                    f= (k-binNodes[i1,j,0])/\
                       (0.0+binNodes[i1,j,3]-binNodes[i1,j,0])
                    kext,salb,asym,srate=pyHB2.getsnowp_z(Nwb[k-psdNodes[i1,j,0]]+\
                                                      fract,\
                                                      z13[i1,j,k],nfreq)
                    if k<binNodes[i,j,2]-1 and dj==0:
                        r1L.append(pRate[i,j,k])
                        r2L.append(srate)
                        kL.append(k)
                kextH[k,dj,:]=kext
                salbH[k,dj,:]=salb
                asymH[k,dj,:]=asym
                
            for k in range(binNodes[i1,j,2],binNodes[i1,j,4]):
                kext,salb,asym=pyHB2.getrainp(Nwb[k-psdNodes[i1,j,0]],\
                                              pRate[i1,j,k],nfreq)
                kextH[k,dj,:]=kext
                salbH[k,dj,:]=salb
                asymH[k,dj,:]=asym
    iEnum=0
    kextH1=kextH[:,0,:].copy()
    salbH1=salbH[:,0,:].copy()
    asymH1=asymH[:,0,:].copy()
    for k in range(88):
        dix=int((87-k)*0.25/5)
        if dix>=4:
            dix=4
            kextH1[k,:]=kextH[k,dix,:]
            salbH1[k,:]=salbH[k,dix,:]
            asymH1[k,:]=asymH[k,dix,:]
        else:
            f=(87-k)*0.25/5-dix
            dix1=dix+1
            #print(f,dix,dix1)
            if f<0.5:
                f=0
            else:
                f=1
            if dix1>4:
                dix1=4
            kextH1[k,:]=(1-f)*kextH[k,dix,:]+f*kextH[k,dix1,:]
            salbH1[k,:]=((1-f)*salbH[k,dix,:]*kextH[k,dix,:]+f*salbH[k,dix1,:]*kextH[k,dix1,:])/(kextH1[k,:]+1e-10)
            asymH1[k,:]=((1-f)*salbH[k,dix,:]*kextH[k,dix,:]*asymH[k,dix,:]+f*asymH[k,dix1,:]*salbH[k,dix1,:]*kextH[k,dix1,:])/\
                         ((1-f)*salbH[k,dix,:]*kextH[k,dix,:]+f*salbH[k,dix1,:]*kextH[k,dix1,:]+1e-10)
    for (pol,f) in zip(npol,freqs):
        kextL=[]
        for q,tk,pa,cld in zip(qv[i,j,:],airTemp[i,j,:],press[i,j,:],cldw[i,j,:]):
            absair,abswv=pyHB2.gasabsr98(f,tk,q*1e-3,pa*1e2,1)
            if cld>0.001:
                z_clw=pyHB2.gcloud(f,tk,cld)
            else:
                z_clw=0.
            kextL.append(absair+abswv+z_clw)
        bins=arange(envNode[i,j,0],envNode[i,j,-1]+1)
        kextInt=interp(bins,envNode[i,j,:],kextL)
        kextInt_old=kextInt.copy()
        kextInt=kextInt+kextH1[envNode[i,j,0]:envNode[i,j,-1]+1,iFreq[iEnum]]
        salb=kextInt.copy()*0.
        asym=kextInt.copy()*0.
        salb=salbH1[envNode[i,j,0]:envNode[i,j,-1]+1,iFreq[iEnum]].copy()*\
            (1-kextInt_old/kextInt)
        asym=asymH1[envNode[i,j,0]:envNode[i,j,-1]+1,iFreq[iEnum]].copy()
        tLayer=list(interp(bins,envNode[i,j,:],airTemp[i,j,:]))
        tLayer.append(sfcTemp[i,j])
        iEnum+=1
        emis=emiss[i,j,iEnum-1]
        ebar=emis
        nL=kextInt.shape[0]
        tb1=pyHB2.radtran(umu,nL,sfcTemp[i,j],tLayer[::-1],\
                          arange(nL+1)*0.25,kextInt[::-1],\
                          salb[::-1],asym[::-1],2.7,emis,ebar)
        tbSim.append(tb1)
    return tbSim, r1L,r2L,kL
                


def rte_noad(psdNodes,binNodes,Nw,pRate,z13,i,j,emiss,qv,airTemp,\
             press,envNode,sfcTemp,umu,pyHB2,idir,cldw):
    freqs=[10.6,10.6,18.7,18.7,23.,37,37.,89,89.,166.,166.,186.3,190.3]
    npol=[1,0,1,0,1,1,0,1,0,1,0,1,1]
    iFreq=[0,0,1,1,2,3,3,4,4,5,5,6,6,7,7]
    nfreq=8
    kextH=zeros((88,5,8),float)
    salbH=zeros((88,5,8),float)
    asymH=zeros((88,5,8),float)

    tbSim=[]
    ns=psdNodes.shape[0]
    for dj in range(5):
        if idir==1:
            i1=min(ns-1,i+dj)
        else:
            i1=max(0,i+dj)
        
        if psdNodes[i1,j,0]>-999:
            Nwb=interp(arange(psdNodes[i1,j,0],psdNodes[i1,j,-1]+1),\
                       psdNodes[i1,j,:],Nw[i1,j,:,0])
            for k in range(max(0,binNodes[i1,j,0]),binNodes[i1,j,0]):
                kext,salb,asym=pyHB2.getsnowp(Nwb[0],\
                                              pRate[i1,j,binNodes[i1,j,0]],nfreq)
                kextH[k,dj,:]=kext
                salbH[k,dj,:]=salb
                asymH[k,dj,:]=asym
            for k in range(binNodes[i1,j,0],\
                           min(binNodes[i1,j,2],binNodes[i1,j,4])):
                kext,salb,asym=pyHB2.getsnowp(Nwb[k-psdNodes[i1,j,0]],\
                                              pRate[i1,j,k],nfreq)
                
                kextH[k,dj,:]=kext
                salbH[k,dj,:]=salb
                asymH[k,dj,:]=asym
            for k in range(binNodes[i1,j,2],binNodes[i1,j,4]):
                kext,salb,asym=pyHB2.getrainp(Nwb[k-psdNodes[i1,j,0]],\
                                              pRate[i1,j,k],nfreq)
                kextH[k,dj,:]=kext
                salbH[k,dj,:]=salb
                asymH[k,dj,:]=asym
    iEnum=0
    kextH1=kextH[:,0,:].copy()
    salbH1=salbH[:,0,:].copy()
    asymH1=asymH[:,0,:].copy()
    for k in range(88):
        dix=int((87-k)*0.25/5)
        if dix>=4:
            dix=4
            kextH1[k,:]=kextH[k,dix,:]
            salbH1[k,:]=salbH[k,dix,:]
            asymH1[k,:]=asymH[k,dix,:]
        else:
            f=(87-k)*0.25/5-dix
            dix1=dix+1
            if dix1>4:
                dix1=4
            if f<0.5:
                f=0
            else:
                f=1.
            kextH1[k,:]=(1-f)*kextH[k,dix,:]+f*kextH[k,dix1,:]
            salbH1[k,:]=(1-f)*salbH[k,dix,:]+f*salbH[k,dix1,:]
            asymH1[k,:]=(1-f)*asymH[k,dix,:]+f*asymH[k,dix1,:]
            kextH1[k,:]=(1-f)*kextH[k,dix,:]+f*kextH[k,dix1,:]
            salbH1[k,:]=((1-f)*salbH[k,dix,:]*kextH[k,dix,:]+f*salbH[k,dix1,:]*kextH[k,dix1,:])/(kextH1[k,:]+1e-10)
            asymH1[k,:]=((1-f)*salbH[k,dix,:]*kextH[k,dix,:]*asymH[k,dix,:]+f*asymH[k,dix1,:]*salbH[k,dix1,:]*kextH[k,dix1,:])/\
                         ((1-f)*salbH[k,dix,:]*kextH[k,dix,:]+f*salbH[k,dix1,:]*kextH[k,dix1,:]+1e-10)
    for (pol,f) in zip(npol,freqs):
        kextL=[]
        for q,tk,pa,cld in zip(qv[i,j,:],airTemp[i,j,:],press[i,j,:],cldw[i,j,:]):
            absair,abswv=pyHB2.gasabsr98(f,tk,q*1e-3,pa*1e2,1)
            if cld>0.001:
                z_clw=pyHB2.gcloud(f,tk,cld)
            else:
                z_clw=0.
            kextL.append(absair+abswv+z_clw)
        bins=arange(envNode[i,j,0],envNode[i,j,-1]+1)
        kextInt=interp(bins,envNode[i,j,:],kextL)
        kextInt_old=kextInt.copy()
        kextInt=kextInt+kextH1[envNode[i,j,0]:envNode[i,j,-1]+1,iFreq[iEnum]]
        salb=kextInt.copy()*0.
        asym=kextInt.copy()*0.
        salb=salbH1[envNode[i,j,0]:envNode[i,j,-1]+1,iFreq[iEnum]].copy()*\
            (1-kextInt_old/kextInt)
        asym=asymH1[envNode[i,j,0]:envNode[i,j,-1]+1,iFreq[iEnum]].copy()
        tLayer=list(interp(bins,envNode[i,j,:],airTemp[i,j,:]))
        tLayer.append(sfcTemp[i,j])
        iEnum+=1
        emis=emiss[i,j,iEnum-1]
        ebar=emis
        nL=kextInt.shape[0]
        tb1=pyHB2.radtran(umu,nL,sfcTemp[i,j],tLayer[::-1],\
                          arange(nL+1)*0.25,kextInt[::-1],\
                          salb[::-1],asym[::-1],2.7,emis,ebar)
        tbSim.append(tb1)
    return tbSim
                
