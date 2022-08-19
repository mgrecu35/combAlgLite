import numpy as np
from pyresample import kd_tree, geometry
swath_def = geometry.SwathDefinition(lons=lons, lats=lats)
result = kd_tree.resample_gauss(swath_def, data,
                                area_def, radius_of_influence=25000, sigmas=12500)
