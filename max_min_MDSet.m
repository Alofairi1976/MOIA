function [Xnew XZ ]= max_min_MDSet(UXY,A)

% this function to genearte new MDSet based with latge differnce from the
% vector of nodes(previous solutions)
% UXY is the vector of valuses that using for determine new MDSet
% A is the network adjencey matrix
% S is selection variable( 0 min , 1 max )
tic
n=size(A,1);
MDS= MDS_ILP(A);
nMDS=length(MDS);
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

BigC=zeros(1,2*n);
BigC=ones(n+1,2*n);

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

[Res solv] = mosekopt('minimize echo(0) ',mds);

       
X1=find(solv.sol.int.xx(1:n)>=10^-6);
XZ=find(solv.sol.int.xx(n+1:2*n)>=10^-6);
Xnew=X1;
toc
end