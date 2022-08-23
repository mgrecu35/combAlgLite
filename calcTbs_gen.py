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
from chaseRet import *


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
    return qv,press,sfcType,pType,airTemp,envNodes,binNodes,\
        bcf,bsf,pwc,sfcEmiss,dm,skTemp,lon,lat
