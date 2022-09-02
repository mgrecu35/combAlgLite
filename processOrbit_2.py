import matplotlib.pyplot as plt
from netCDF4 import Dataset
import xarray as xr
from pyresample import kd_tree, geometry
#stop
import read_tables as cAlg
cAlg.read_tables()
cAlg.initp2()
cAlg.read_zku_dm()
import fkb
fkb.init_keras()


import glob
f1=[]
f2=[]
f3=[]
#stop
import pickle
if 1==0:
    for i in range(31):
        f11=sorted(glob.glob("/gpmdata/2018/12/%2.2i/radar/2A.GPM.DPR.V9-20211125.*"%(i+1)))
        f21=sorted(glob.glob("/gpmdata/2018/12/%2.2i/radar/2B.GPM.DPRGMI.CO*"%(i+1)))
        f31=sorted(glob.glob("/gpmdata/2018/12/%2.2i/1C/1C.GPM.GMI.*"%(i+1)))
        f1.extend(f11)
        f2.extend(f21)
        f3.extend(f31)
    for i in range(31):
        f11=sorted(glob.glob("/gpmdata/2019/01/%2.2i/radar/2A.GPM.DPR.V9-20211125.*"%(i+1)))
        f21=sorted(glob.glob("/gpmdata/2019/01/%2.2i/radar/2B.GPM.DPRGMI.CO*"%(i+1)))
        f31=sorted(glob.glob("/gpmdata/2019/01/%2.2i/1C/1C.GPM.GMI.*"%(i+1)))
        f1.extend(f11)
        f2.extend(f21)
        f3.extend(f31)
    for i in range(28):
        f11=sorted(glob.glob("/gpmdata/2019/02/%2.2i/radar/2A.GPM.DPR.V9-20211125.*"%(i+1)))
        f21=sorted(glob.glob("/gpmdata/2019/02/%2.2i/radar/2B.GPM.DPRGMI.CO*"%(i+1)))
        f31=sorted(glob.glob("/gpmdata/2019/02/%2.2i/1C/1C.GPM.GMI.*"%(i+1)))
        f1.extend(f11)
        f2.extend(f21)
        f3.extend(f31)

    import pickle

else:
    f1,f2,f3=pickle.load(open("djf_2018_fileList.pklz","rb"))


orbits=["2019/01/13 27701",
        "2019/01/08 27624",
        "2019/01/12 27686",
        "2019/01/11 27677",
        "2019/01/09 27645",
        "2019/01/08 27631"]

f1,f2,f3=[],[],[]
f4=[]
import os
for orbit in orbits:
    cal_date=orbit[:10]
    orb_numb=orbit[-5:]
    f11=sorted(glob.glob("/gpmdata/%s/radar/2A.GPM.DPR.V9-20211125.*%s*"%(cal_date,orb_numb)))
    f21=sorted(glob.glob("/gpmdata/%s/radar/2B.GPM.DPRGMI.CO*%s*"%(cal_date,orb_numb)))
    f31=sorted(glob.glob("/gpmdata/%s/1C/1C.GPM.GMI.*%s*"%(cal_date,orb_numb)))
    fpath="/PANFS/user/home/wsolson/data/MRMS_Matched/test2/*%s*%s*"%(cal_date.replace("/",""),\
                                                                      orb_numb)
    fmrms=sorted(glob.glob(fpath))
    #print(fmrms)
    for f in fmrms:
        #os.system("cp %s MRMS"%f)
        f4.append(f)
    #print(fmrms)
    f1.extend(f11)
    f2.extend(f21)
    f3.extend(f31)

#stop
import numpy as np
#from rteSimpleModule import rte,calcz
#from chaseRet import *


npix1,npix2=0,0

umu=np.cos(53/180*3.14)
def readData(f1_,f2_,f3_,n1,n2):
    fh=Dataset(f1_)
    lon=fh["FS/Longitude"][n1:n2,:]
    lat=fh["FS/Latitude"][n1:n2,:]
    #icount=0
    #n1=7300
    #n2=0
    #iscan=0
    #iscanL=[]
    #for lon1,lat1 in zip(lon[:,24],lat[:,24]):
    #    if lon1>-128 and lon1<-65 and lat1>30 and lat1<65:
    #        icount+=1
    #        iscanL.append(iscan)
    #        n1=min(iscan,n1)
    #        n2=max(iscan,n2)
    #    iscan+=1
    #print(iscanL)
    #print(icount,n1,n2)
    #plt.plot(lon[0:n1,24],lat[0:n1,24])
    #plt.plot(lon[n1:n2,24],lat[n1:n2,24])
    #return
    icount=1
    lon=fh["FS/Longitude"][n1:n2,:]
    lat=fh["FS/Latitude"][n1:n2,:]
    zKu=fh["FS/PRE/zFactorMeasured"][n1:n2,:,:,:]
    sfcType=fh["FS/PRE/landSurfaceType"][n1:n2,:]
    pType=fh["FS/CSF/typePrecip"][n1:n2,:]
    pType=(pType/1e7).astype(int)
    fh_cmb=Dataset(f2_)
    qv=fh_cmb["KuKaGMI/vaporDensity"][n1:n2,:,:]
    press=fh_cmb["KuKaGMI/airPressure"][n1:n2,:,:]
    envNodes=fh_cmb["KuKaGMI/envParamNode"][n1:n2,:,:]
    airTemp=fh_cmb["KuKaGMI/airTemperature"][n1:n2,:,:]
    skTemp=fh_cmb["KuKaGMI/skinTemperature"][n1:n2,:]
    binNodes=fh_cmb["KuKaGMI/phaseBinNodes"][n1:n2,:]
    bcf=fh["FS/PRE/binClutterFreeBottom"][n1:n2,:]
    bsf=fh["FS/PRE/binRealSurface"][n1:n2,:,:]
    pwc=fh_cmb["KuKaGMI/precipTotWaterCont"][n1:n2,:,:]
    sfcEmiss=fh_cmb["KuKaGMI/surfEmissivity"][n1:n2,:,:]
    dm=fh_cmb["KuKaGMI/precipTotDm"][n1:n2,:,:]
    cldw=fh_cmb["KuKaGMI/cloudLiqWaterCont"][n1:n2,:,:]
    sfcBin=fh_cmb["KuKaGMI/Input/surfaceRangeBin"][n1:n2,:,:]
    sfcPrecip=fh["FS/SLV/precipRateNearSurface"][n1:n2,:]
    bbPeak=fh["FS/CSF/binBBPeak"][n1:n2,:]
    #print(fh["FS/PRE"])
    #print(fh["FS/VER"])
    h0=fh["FS/VER/heightZeroDeg"][n1:n2,:]
    pRateCMB=fh_cmb["KuKaGMI/precipTotRate"][n1:n2,:,:]
    #print(fh_cmb["KuKaGMI"])
    #stop
    fh.close()
    fh_cmb.close()
    a=np.nonzero(pType==1)
    if icount>0:
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
    else:
        tc_regrid,tc_regrid2=[],[]
    umu=np.cos(53/180.0*np.pi)
    bsf1=bsf[:,:,0]
    n1=qv.shape[0]
    a=np.nonzero(pType>0)
    print(len(a[0]))
        
    return qv,press,sfcType,pType,airTemp,envNodes,binNodes,\
        bcf,bsf1,pwc,sfcEmiss,dm,skTemp,lon,lat,tc_regrid,\
        tc_regrid2,zKu,cldw,bbPeak,h0,pRateCMB,icount,n1,n2


f1_=f1[0]
f2_=f2[0]
f3_=f3[0]
#stop



def processOrbit(f1_,f2_,f3_):
    tboutL=[]
    tbout_dfrL=[]
    tbobsL=[]
    tbobs2L=[]
    sfcTypeL=[]
    h0L=[]
    swp=[]
    swp_dfr=[]
    swpL=[]
    ictot=0
    pRateDiagL=[]
    for i in range(10):
        qv,press,sfcType,pType,airTemp,envNodes,binNodes,\
            bcf,bsf,pwc,sfcEmiss,dm,skTemp,lon,lat,tc_regrid,\
            tc_regrid2,zObs,cldw,bbPeak,h0,pRateCMB,ic=\
                readData(f1_,f2_,f3_,i*700,(i+1)*700)
        if ic==0:
            continue
        else:
            ictot+=ic
            a=np.nonzero(pType>0)
            b=np.nonzero(h0.data[a]<0)
            print("#snow profiles",len(a[0][b]))
            continue
        zObs=zObs[:,:,::2,:]
        bbPeak=(bbPeak/2).astype(int)
        bsf=(bsf/2).astype(int)
        a=np.nonzero(pType>0)
        bL=[]
        if len(a[0])>0:
            for i1,i2 in zip(a[0],a[1]):
                t1=np.interp(range(88),envNodes[i1,i2],airTemp[i1,j1])
                if t1[binNodes[i1,j1,4]-1]<273.15:
                    print(pwc[i1,j1,binNodes[i1,j1,4]-1],pwc[i1,j1,binNodes[i1,j1,4]])
        if 1==0:
            tbout,tbout_dfr,\
                prate_out,\
                pwc_out,swp = \
                    cAlg.calc_tb_f90_flattened(zObs[a],binNodes[a],pwc[a],dm[a],\
                                               bsf[a],pType[a],envNodes[a],\
                                               qv[a],airTemp[a],press[a],\
                                               cldw[a],umu,skTemp[a],\
                                               sfcEmiss[a], \
                                               bbPeak[a])
            tboutL.extend(tbout)
            tbout_dfrL.extend(tbout_dfr)
            tbobsL.extend(tc_regrid[a])
            tbobs2L.extend(tc_regrid2[a])
            sfcTypeL.extend(sfcType[a])
            h0L.extend(h0[a])
            swpL.extend(swp)
            binNodes_=binNodes[a]
            bcf_=(bcf[a]/2).astype(int)
            for ip,pRate1 in enumerate(pRateCMB[a]):
                if binNodes_[ip,3]<=bcf_[ip]-1:
                    pRateDiagL.append([prate_out[ip,binNodes_[ip,1]],pRate1[binNodes_[ip,1]],\
                                       pRate1[binNodes_[ip,3]]])
                    
    
    return ictot
    tboutL=np.array(tboutL)
    tbout_dfrL=np.array(tbout_dfrL)
    tbobsL=np.array(tbobsL)
    tbobs2L=np.array(tbobs2L)
    
    
    swpX=xr.DataArray(swpL,dims=['nt','n2'])
    tbobsX=xr.DataArray(tbobsL,dims=['nt','n9'])
    tbobs2X=xr.DataArray(tbobs2L,dims=['nt','n4'])
    tboutX=xr.DataArray(tboutL,dims=['nt','n13'])
    tbout_dfrX=xr.DataArray(tbout_dfrL,dims=['nt','n13'])
    sfcTypeX=xr.DataArray(sfcTypeL,dims=['nt'])
    h0LX=xr.DataArray(h0L,dims=['nt'])
    
    ds=xr.Dataset({"tb_sim":tboutX,"tb_sim_dfr":tbout_dfrX,"tb_obs_S1":tbobsX,"tb_obs_S2":tbobs2X,\
                   "sfcType":sfcTypeX,"h0":h0LX,"swp":swpX})
    comp = dict(zlib=True, complevel=5)
    encoding = {var: comp for var in ds.data_vars}
    
    fnameout='simTb.'+f1_.split('.')[-4]+'.'+ f1_.split('.')[-3]+'.nc'
    ds.to_netcdf(fnameout, encoding=encoding)
    
for i in range(1,0):
    #print(f1_,f2_,f3_)
    if ('20190108' not in f1_) or ('27631' not in f1_):
        continue
    break
    ictot=processOrbit(f1_,f2_,f3_)
    print("ictot=",ictot)


import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')
pwcL=[]
zDBL=[]
sfcPrecipL=[]
tempL=[]
for f1_,f2_,f3_,f4_ in zip(f1[:],f2[:],f3[:],f4):
    i=0
    nc1=Dataset(f4_)
    lon_mrms=nc1["Longitude"][:]
    lat_mrms=nc1["Latitude"][:]
    scanN=nc1["scanNumber"][:]
    precipRateMRMS=nc1["precipRateMRMS"][:]
    precipRateSAT=nc1["nearSurfPrecipTotRateSAT"][:]
    pTypeMRMS=nc1["precipTypeMRMS"][:]
    pTypeMRMS_m=np.ma.array(pTypeMRMS,mask=pTypeMRMS<0)
    #plt.pcolormesh(lon_mrms,lat_mrms,pTypeMRMS_m)
    pTypeSAT=nc1["precipTypeSAT"][:]
    pTypeSAT_m=np.ma.array(pTypeSAT,mask=pTypeSAT<0)
    rqi=nc1["radarQualityIndex"][:]
    #plt.figure()
    #plt.pcolormesh(lon_mrms,lat_mrms,pTypeSAT_m)
    #plt.colorbar()
    #,stop
    qv,press,sfcType,pType,airTemp,envNodes,binNodes,\
        bcf,bsf,pwc,sfcEmiss,dm,skTemp,lon,lat,tc_regrid,\
        tc_regrid2,zObs,cldw,bbPeak,\
        h0,pRateCMB,ic,n1,n2=\
                              readData(f1_,f2_,f3_,scanN[0],scanN[-1]+1)
    zObs=zObs[:,:,::2,:]
    bbPeak=(bbPeak/2).astype(int)
    bsf=(bsf/2).astype(int)
    a=np.nonzero(pType>0)
    bL=[]
    
                
    tbout,tbout_dfr,\
        prate_out,\
        pwc_out,swp = \
                      cAlg.calc_tb_f90_flattened(zObs[a],binNodes[a],pwc[a],dm[a],\
                                                 bsf[a],pType[a],envNodes[a],\
                                                 qv[a],airTemp[a],press[a],\
                                                 cldw[a],umu,skTemp[a],\
                                                 sfcEmiss[a], \
                                                 bbPeak[a])
    if len(a[0])>0:
        for ic,i1 in enumerate(a[0]):
            i2=a[1][ic]
            t1=np.interp(range(88),envNodes[i1,i2],airTemp[i1,i2])
            if t1[binNodes[i1,i2,4]-1]<273.15 and pTypeMRMS[i1,i2]==3 and precipRateMRMS[i1,i2]>=-0.01 \
               and rqi[i1,i2]>0.8 and abs(i2-24)<24:
                #print(pwc[i1,i2,binNodes[i1,i2,4]-1],pwc[i1,i2,binNodes[i1,i2,4]])
                pwcL.append([pwc[i1,i2,binNodes[i1,i2,4]-1],pwc_out[ic,binNodes[i1,i2,4]-1]])
                sfcPrecipL.append([pRateCMB[i1,i2,binNodes[i1,i2,4]-1],pwc[i1,i2,binNodes[i1,i2,4]-1],\
                                   precipRateMRMS[i1,i2],precipRateSAT[i1,i2]])
                zDBL.append(zObs[i1,i2,binNodes[i1,i2,4]-1,:])
                tempL.append(t1[binNodes[i1,i2,4]-1])
    #stop
    #print(n1,n2)
    #if ic>0:
    #    plt.figure()
    #    plt.pcolormesh(lon,lat,h0,vmin=0)
    #    a=np.nonzero(pType>0)
    #    b=np.nonzero(h0.data[a]<0)
    #    print("# profiles",len(b[0]))
        

pwcL=np.array(pwcL)
sfcPrecipL=np.array(sfcPrecipL)
a1=np.nonzero(sfcPrecipL[:,1]>0)
sRateCoeff=np.polyfit(np.log10(sfcPrecipL[a1[0],1]),np.log10(sfcPrecipL[a1[0],2]),1)
zDBL=np.array(zDBL)
tempL=np.array(tempL)
zDBL[zDBL<0]=0

from chaseRet import *
Nw,dm,IWC= retr(zDBL,tempL-273.15)
Xf=np.zeros((zDBL.shape[0],3),float)
Xf[:,:2]=zDBL
Xf[:,2]=tempL-273.15
y_keras=fkb.call_keras(Xf)
iwc_dfr,dm_dfr=cAlg.iwc_from_dfr(zDBL,tempL-273)
zsim_out,iwc_dfr2,dm_dfr2,pRate_out=cAlg.iwc_from_dfr2(zDBL,tempL-273)
from sklearn.neighbors import KNeighborsRegressor

nc2=Dataset("collocatedZ_SSRGA_BF.nc")
z1=nc2["zKu_SSRGA"][:]
z2=nc2["zKa_SSRGA"][:]
tempC=nc2["tempC"][:]
iwc2=nc2["iwc_PSD"][:]
X=[]
y=[]
for x1,x2,x3,iwc21 in zip(z1,z2,tempC,iwc2):
    if x1>0 and x2>0 and iwc21>-0.1:
        X.append([x1,x2,x3])
        y.append(iwc21)
n_neigh=150
neigh_iwc=KNeighborsRegressor(n_neighbors=n_neigh,weights='distance')
neigh_iwc.fit(X,y)
y2=neigh_iwc.predict([[dbz1[0],dbz1[1],temp1-273.15] for dbz1,temp1 in zip(zDBL,tempL)])
