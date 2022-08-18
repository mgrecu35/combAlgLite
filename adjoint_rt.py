import numpy as np

def tb_back(k,a,g,T,eps,dz,Ts,inc_ang,rt):		
    abig,b = rt.seteddington1d(k[:],a[:],g[:],T,eps,dz,Ts)
    #inc_ang=53.
    i01p= np.linalg.solve(abig,b)
    tb2 = rt.tbf90(i01p,T,inc_ang,k[:],a[:],g[:],eps,dz)
    tbb=1.0
    lam2=i01p
    i01pb,kb,ab,gb = rt.tbf90_b(tb2,tbb,i01p,T,inc_ang,k,a,g,eps,dz)
    y=0.
    yb=1.
    lam1=i01pb
    lam11=np.dot(i01pb,np.linalg.inv(abig))
    kb1,ab1,gb1 = rt.suml1al2_b(y,yb,lam11,lam2,k,a,g,T,eps,dz)
    y=0.
    yb=1
    kb2,ab2 = rt.suml1b_b(y,yb,lam11,k,a,g,T,eps,dz)
    y=0
    yb=1
    
    epsb = rt.tbf90_eps_b(tb2,tbb,i01p,T,inc_ang,k,a,g,eps,dz)
    epsb1,abigb = rt.suml1al2_eps_b(y,yb,lam11,lam2,k,a,g,T,eps,dz)
    dT_da=ab-ab1+ab2
    dT_dk=kb-kb1+kb2
    dT_dg=gb-gb1
    dT_deps=epsb-epsb1
    return tb2,dT_dk,dT_da,dT_dg, dT_deps

# A^-1*b=I
# L^T*A^-1*b=L^T*I=Tb
# dL^T*A^-1*b+L^T*d(A^-1)*b+dL^T*A^-1*db=dTb
# dL= from backwards tb90_b dL_dkb=kb
# d(A^-1)=-A^-1*dA*A^-1


#abig,b = rt.seteddington1d(k[:],a[:],g[:],T,eps,dz,Ts)
#incang=53.
#i01p=linalg.solve(abig,b)
##tb2 = rt.tbf90(i01p,T,incang,k[:],a[:],g[:],eps,dz)
#print(tb,tb2)
#tbb=1.0
#lam2=i01p
#i01pb,kb,ab,gb = rt.tbf90_b(tb2,tbb,i01p,T,incang,k,a,g,eps,dz)
#y=0.
#yb=1.
#lam1=i01pb
#lam11=dot(i01pb,linalg.inv(abig))
#kb1,ab1,gb1 = rt.suml1al2_b(y,yb,lam11,lam2,k,a,g,T,eps,dz)
#y=0.
#yb=1
#kb2,ab2 = rt.suml1b_b(y,yb,lam11,k,a,g,T,eps,dz)
#dT_da=ab-ab1+ab2
#dT_dk=kb-kb1+kb2
#dT_dg=gb-gb1
