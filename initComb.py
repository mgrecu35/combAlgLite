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
dmL=[]
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
    temp=fh_cmb["KuKaGMI/airTemperature"][:]
    nodes=fh_cmb["KuKaGMI/phaseBinNodes"][:]
    bcf=fh["FS/PRE/binClutterFreeBottom"][:]
    bsf=fh["FS/PRE/binRealSurface"][:]
    pwc=fh_cmb["KuKaGMI/precipTotWaterCont"][:]
    sfcEmiss=fh_cmb["KuKaGMI/surfEmissivity"][:]
    dm=fh_cmb["KuKaGMI/precipTotDm"][:]
    a=np.nonzero(pType==1)
    for i0,j0 in zip(a[0],a[1]):
        if dm[i0,j0,nodes[i0,j0,1]]>0 and dm[i0,j0,nodes[i0,j0,3]]>0:
            dmL.append([dm[i0,j0,nodes[i0,j0,1]],dm[i0,j0,nodes[i0,j0,3]]])
    break