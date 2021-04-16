function [BigA BigC BigLC BigUC BigLX BigUX] =makeBigProblemR(A,nMDS)

% This function tend to generate the Big Matrix that represent the problem
% to find the Most different MDSets 

n=size(A,1);
Z=zeros(n,n);
BigA=sparse([A Z Z; Z A Z]);
i=1;
for k=1:n
    BigA(2*n+i,i)=1; %xi
    BigA(2*n+i,n+i)=1; %yi
    BigA(2*n+i,2*n+i)=1; %zi
    i=i+1;
 end
 i=1;
 for k=1:n
    BigA(3*n+i,i)=1; %xi
    BigA(3*n+i,n+i)=1; %yi
    BigA(3*n+i,2*n+i)=-1; %zi
    i=i+1;
 end
 
 BigA(4*n+1,1:n)=ones(1,n);
 BigA(4*n+2,n+1:n+n)=ones(1,n);
 
BigC=ones(1,3*n); 
%BigC=[zeros(1,2*n), ones(1,n)];

BigLC=zeros(1,3*n+2); %2 is to for sum X = nMDS and sum Y = nmDS
BigLC(1:n)=1; % 1<=AX<=Inf
BigLC(n+1:2*n)=1;% 1<=AY<=Inf
BigLC(2*n+1:3*n)=-inf;  %-inf =< xi +yi+zi <=2
BigLC(3*n+1:4*n)=0;  %inf>= xi + yi-zi >=0
BigLC(4*n+1) =nMDS;
BigLC(4*n+2) =nMDS;

BigUC=zeros(1,3*n+2);
BigUC(1:n)=Inf; % 1<=AX<=Inf
BigUC(n+1:2*n)=Inf;% 1<=AY<=Inf
BigUC(2*n+1:3*n)=2;  %-inf =< xi +yi+zi <=2
BigUC(3*n+1:4*n)=inf;  %inf >= xi + yi-zi >=0
BigUC(4*n+1) = nMDS;
BigUC(4*n+2) = nMDS;

BigLX=zeros(1,n*3);
BigUX=ones(1,3*n);
%BigUX(2*n+1,n*3)=2;



%===================================================

end