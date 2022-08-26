import read_tables as cAlg
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from pyresample import kd_tree, geometry
cAlg.read_tables()
cAlg.initp2()
import glob
f1=sorted(glob.glob("Data_2019/2A*"))
f2=sorted(glob.glob("Data_2019/2B*"))
f3=sorted(glob.glob("Data_2019/1C*"))
f1=f1[:]
f2=f2[:]
f3=f3[:]

import numpy as np
#from rteSimpleModule import rte,calcz
#from chaseRet import *


npix1,npix2=0,0

umu=np.cos(53/180*3.14)
def readData(f1_,f2_,f3_,n1,n2):
    fh=Dataset(f1_)
    zKu=fh["FS/PRE/zFactorMeasured"][n1:n2,:,:,:]
    lon=fh["FS/Longitude"][n1:n2,:]
    lat=fh["FS/Latitude"][n1:n2,:]
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
        bcf,bsf1,pwc,sfcEmiss,dm,skTemp,lon,lat,tc_regrid,tc_regrid2,zKu,cldw,bbPeak,h0,pRateCMB


f1_=f1[0]
f2_=f2[0]
f3_=f3[0]
#stop

import xarray as xr

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
    pRateDiagL=[]
    for i in range(10):
        qv,press,sfcType,pType,airTemp,envNodes,binNodes,\
            bcf,bsf,pwc,sfcEmiss,dm,skTemp,lon,lat,tc_regrid,\
            tc_regrid2,zObs,cldw,bbPeak,h0,pRateCMB=\
                readData(f1_,f2_,f3_,i*700,(i+1)*700)
        zObs=zObs[:,:,::2,:]
        bbPeak=(bbPeak/2).astype(int)
        bsf=(bsf/2).astype(int)
        a=np.nonzero(pType==1)
        if len(a[0])>0:
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
    

processOrbit(f1_,f2_,f3_)
