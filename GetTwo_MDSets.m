function [MDS1 MDS2] = GetTwo_MDSets( A )
% This Function developed to solve BigMatrix by ILP solver to generate
% Multiple minimum dominating set
% Remember use 'addpath c:\Program Files\mosek\8\toolbox\r2014a'; to
% connect matlab with mosek.

% Mosek intionalize
%*****************************************************************%
    mosekParams.MSK_DPAR_PRESOLVE_TOL_ABS_LINDEP= 1.0000e-006;
    mosekParams.MSK_DPAR_PRESOLVE_TOL_ABS_LINDEP=10^-6;
    mosekParams.MSK_IPAR_PRESOLVE_USE=1;
    
    mosekParams.MSK_IPAR_OPTIMIZER = 'MSK_OPTIMIZER_PRIMAL_SIMPLEX';
    
% generate the  Big Matrix  representing the ILP-Model
   
tic
n=size(A,1);
   MDS= find(MDS_ILP(A));
   nMDS=length(MDS);
   RS=zeros(1,nMDS);  % ***** RS is matrix for target multiple minimum dominating sets  ****
   
   %****  BigA is big matrix (4*n+1,3*n) to represent ILP model  
   [BigA BigC BigLC BigUC BigLX BigUX] =makeBigProblemR(sparse(A),nMDS);

    clear mds;
    mds.a=BigA;
    mds.c=BigC';
    mds.blc=BigLC';
    mds.buc=BigUC';
    mds.blx=BigLX;
    mds.bux=BigUX;

    mds.ints.sub=1:3*n;

   
   [Res R] = mosekopt('maximize echo(0) ',mds);

   
   MDS1=find(R.sol.int.xx(1:n)>=10^-6);
   MDS2=find(R.sol.int.xx(n+1:2*n)>=10^-6);
   %MDS1=R.sol.int.xx(1:n);
   %MDS2=R.sol.int.xx(n+1:2*n);
   
end
