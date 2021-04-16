function  [D C R T] = crds(E )
% Function to determine the critical , redundant  and intermittent  nodes.
%=========================================================================
% first to claculate critical nodes using Minimun dominating set
   tic
   clc
% addpath 'c:\Program Files\mosek\8\toolbox\r2014a';
 mosekParams.MSK_DPAR_PRESOLVE_TOL_ABS_LINDEP= 1.0000e-006;
 mosekParams.MSK_DPAR_PRESOLVE_TOL_ABS_LINDEP=10^-6;
 mosekParams.MSK_IPAR_PRESOLVE_USE=1;
  NodSet=1:size(E,1);
  MDS= find(MDS_ILP(E));
 % MDS=find(sol1.xx>10^-6);
   D=MDS;
  
m=length(D);

n=size(E,1)
cds=ones(1,m);
clear mds;

mds.a=E;
mds.c= ones(n,1);
mds.blc=ones(n,1);
mds.buc=inf* ones(n,1);
mds.blx= zeros(n,1);
mds.bux = ones(n,1);

mds.a(n+1,:)=ones(1,n);
mds.blc(n+1)=m;
mds.buc(n+1)=m;

mds.ints.sub=1:n;

for i=1:m
       
    mds1=mds;  
    mds1.bux(MDS(i))=0;
   [Res sol2] = mosekopt('minimize echo(0)',mds1,mosekParams);
 
   if ~(strcmp(sol2.sol.int.prosta,'PRIMAL_FEASIBLE') )
       %|| sol2.sol.int.pobjval==m)
       cds(i)=0;
   end
  clear sol2;  
end
d=find(cds==0);
C=MDS(d);

%============================================================
 %% Calculate the of redundant nodes 
 
    REM=setdiff(NodSet,MDS);
    l=length(REM);
    rds=ones(1,l);

      for i=1:l
         mds1=mds;  
         mds1.blx(REM(i))=1;
        [Res sol2] = mosekopt('minimize echo(0)',mds1,mosekParams);
 
            if ~(strcmp(sol2.sol.int.prosta,'PRIMAL_FEASIBLE') || sol2.sol.int.pobjval==m)
                rds(i)=0; 
             end
         clear sol2;  
     end
t=find(rds==0);
t1=REM(t);
R=NodSet(t1);
%============================================================
% determine the intermittent  nodes.
T= setdiff(NodSet,union(C,R));
save info.mat
toc
end

