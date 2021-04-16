function [  MMDS CT] = Get_MMDSets( A )
% This Function developed to solve BigMatrix by ILP solver to generate
% Multiple minimum dominating sets From the grpah adjancy Matrix.
% The function also find the critical nodes as the intersection between the
% first generated Two MDSets.
% Remember use 'addpath c:\Program Files\mosek\8\toolbox\r2014a'; to
% connect matlab with mosek.
% A is the adjancy matrix

% Mosek intionalize
%*****************************************************************%
    mosekParams.MSK_DPAR_PRESOLVE_TOL_ABS_LINDEP= 1.0000e-006;
    mosekParams.MSK_DPAR_PRESOLVE_TOL_ABS_LINDEP=10^-6;
    mosekParams.MSK_IPAR_PRESOLVE_USE=1;
% generate the  Big Matrix  representing the ILP-Model
   
tic
n=size(A,1);
   MDS= MDS_ILP(A);
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

   
    X=find(R.sol.int.xx(1:n)>=10^-6);
    Y=find(R.sol.int.xx(n+1:2*n)>=10^-6);
    Z=find(R.sol.int.xx(2*n+1:3*n)>=10^-6);
    
    co=1;
    RS(co,:)=X;    % first solution (fMDSet)
    RS(co+1,:)=Y;  % second solution (sMDSet)
    co=co+1;
    
    CT=intersect(X,Y);  % Critical Nodes
    XUY=unique(union(X,Y));
    f=1;
    clear BigA;
while f
    [Xnew XZ]= Gen_NewMDSet(XUY,nMDS,A);
    co=co+1;
    RS(co,:)=Xnew;
    
    if isempty(setdiff(Xnew,XUY))
          f=0;
     else
         XUY=unique(union(XUY,Xnew));
     end
end
   
toc
MMDS=RS';
[ crit_set red_set intr_set ] =nodesanalyze(MMDS,n );
toc
end
%============== Generate  New MDset ==============
%=================================================
function [Xnew XZ ]= Gen_NewMDSet(UXY,nMDS, A)

% This function to genearte new MDSet based with latge differnce from the
% vector of nodes(previous solutions)
n=size(A,1);
BigA = sparse(3*n+1,2*n);
for i=1:n
for j=1:n
BigA(i,j)=A(i,j);
end
end

i=1;
for k=1:n
    BigA(n+i,i)=1; %yi
    BigA(n+i,n+i)=1; %zi
    
    BigA(2*n+i,i)=1; %yi
    BigA(2*n+i,n+i)=-1; %zi
    i=i+1;
 end
BigA(3*n+1,1:2*n)=zeros(1,2*n);
BigA(3*n+1,1:n)=ones(1,n);

BigC=ones(1,2*n);

BigLC=zeros(1,3*n+1); %1 is to for sum X = nMDS 
BigUC=zeros(1,3*n+1); %1 is to for sum X = nMDS 
BigLC(1:n)=1; % 1<=AY<=Inf
BigUC(1:n)=inf; % 1<=AY<=Inf

BigLC(n+1:2*n)=-inf;
BigLC(3*n+1)=nMDS;
BigUC(3*n+1)=nMDS;

BigLX=zeros(1,2*n);
BigUX=ones(1,2*n);

mds.ints.sub=1:2*n;


for i =1:n
  if ismember(i,UXY)
      BigLC(n+i)=-inf;
      BigUC(n+i) =1;
  else
      BigLC(n+i)=-inf;
      BigUC(n+i)=2;
  end

  if ismember(i,UXY)
      BigLC(2*n+i)=-1;
      BigUC(2*n+i) =inf;
  else
      BigLC(2*n+i)=0;
      BigUC(2*n+i)=inf;
  end
end

mds.a=sparse(BigA);
mds.c=BigC';
mds.blc=BigLC';
mds.buc=BigUC';
mds.blx=BigLX;
mds.bux=BigUX;

   
[Res sol2] = mosekopt('maximize echo(0) ',mds);

X1=find(sol2.sol.int.xx(1:n)>=10^-6);
XZ=find(sol2.sol.int.xx(n+1:2*n)>=10^-6);
Xnew=X1;

end


