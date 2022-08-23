using PyCall
push!(PyVector(pyimport("sys")."path"), "/home/grecu/combAlgLite")

py_read=pyimport("calcTbs_gen")

#qv,press,sfcType,pType,airTemp,envNodes,binNodes,\
#       bcf,bsf,pwc,sfcEmiss,dm=

f1_=py_read.f1[1]
f2_=py_read.f2[1]
f3_=py_read.f3[1]

satData=py_read.readData(f1_,f2_,f3_)

qv,press,sfcType,pType,airTemp,envNodes,binNodes,
bcf,bsf,pwc,sfcEmiss,dm,sfcTemp=satData

n1=size(qv)[1]

n1_=Array{Int32}(undef,1)
n1_[1]=Int32(n1);
tbout=Array{Float32}(undef,(n1,49,13));
pType=Int32.(pType);
binNodes=Int32.(binNodes);

bsf=Int32.(round.(bsf./2));
pType=Int32.(pType);
envNodes=Int32.(envNodes);
umu=Array{Float32}(undef,(1));
umu[1]=Float32(cos(53.0/180.0*3.1415));
cldw=zeros(Float32,n1,49,88);
#subroutine calc_tb_f90(n1,binNodes,pwc,dm,sfcBin,pType,envNodes,qv,airTemp,press,&
#                       cldw,umu,sfcTemp,emiss,tbout)

ccall((:calc_tb_f90_,"read_tables.cpython-39-x86_64-linux-gnu"),
      Cvoid,(Ref{Int32},Ref{Int32},
             Ref{Float32},Ref{Float32},
             Ref{Int32},Ref{Int32},Ref{Int32},Ref{Float32},
             Ref{Float32},Ref{Float32},Ref{Float32},Ref{Float32},Ref{Float32},Ref{Float32},Ref{Float32}),
      n1_,binNodes,pwc,dm,bsf,pType,envNodes,qv,airTemp,press,cldw,umu,sfcTemp,
      sfcEmiss,tbout)


