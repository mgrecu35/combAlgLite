gfortran -c -fPIC -O3 -fopenmp -finit-real=zero -finit-integer=0 -fcheck=all src_f90/band_dble.f90
#gfortran -c -fPIC -O3  src_f90/band_dble.f90
#gfortran -c -g -fPIC -O src_f90/gcloud.f
#gfortran -c -g -fopenmp -finit-real=zero -finit-integer=0 -fcheck=all src_f90/radtran_tau_dble.f
gfortran -c -g -O3 -fPIC -fopenmp src_f90/radtran_tau_dble.f 
f2py -c --debug -m read_tables --f90flags="-fopenmp -finit-real=zero -finit-integer=0" src_f90/readTables_new.f90 src_f90/readTables_nonsph.f90\
    src_f90/bisection.f90 src_f90/hbprof_new.f90 src_f90/absorption3D.f90\
    src_f90/rosen.f  src_f90/emissivity-sp.f  radtran_tau_dble.o \
    band_dble.o src_f90/gcloud.f src_f90/eddington.f90 tbCalc.f90 src_f90/radtran_py_wrap.f90 -lgomp
