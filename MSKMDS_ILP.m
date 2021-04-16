function MDS= MSKMDS_ILP(E)
% This function to solve graph minimun dominating set by ILP model using 
% Mosek solver form matlab 
% addpath 'c:\Program Files\mosek\8\toolbox\r2014a';
 tic
   [r,res] = mosekopt('symbcon');
   sc = res.symbcon;
 
    n=size(E,1); 
    mds.a=sparse(E);
    mds.c= ones(n,1);
    mds.blc=ones(n,1);
    mds.buc=inf* ones(n,1);
    mds.blx= zeros(n,1);
    mds.bux = ones(n,1);

    mds.ints.sub=1:n;

 mosekParams.MSK_DPAR_PRESOLVE_TOL_ABS_LINDEP= 1.0000e-006;
 mosekParams.MSK_DPAR_PRESOLVE_TOL_ABS_LINDEP=10^-6;
 mosekParams.MSK_IPAR_PRESOLVE_USE=1;
 
 [ res sol1] = mosekopt('minimize echo(0) nokeepenv',mds,mosekParams);
 
 MDS=find(sol1.sol.int.xx>=10^-6)
 
  toc 
end