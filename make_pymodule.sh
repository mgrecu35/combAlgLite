#gfortran -c -fPIC -O3 src_f90/band_dble.f90
gfortran -c -fPIC -O src_f90/gcloud.f
f2py -c -m read_tables src_f90/readTables_new.f90 src_f90/readTables_nonsph.f90\
    src_f90/bisection.f90 src_f90/hbprof_new.f90 src_f90/absorption3D.f90\
    src_f90/rosen.f  src_f90/emissivity-sp.f src_f90/radtran_tau_dble.f \
    band_dble.o src_f90/gcloud.f src_f90/eddington.f90 tbCalc.f90
