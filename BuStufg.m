%% Init
clear
close all
clc

%% 
C=ff2n(18);
% C=C(sum(C,2)==3 | sum(C,2)==4,:); 
C(sum(C,2)~=4,:)=[];
D=pdist(C,'hamming')*size(C,2);
Dsqr=squareform(D); 
% tabulate(D)

A=squareform(D)<3;
for i=1:size(A,1)
    A(i,i)=0;
end
X=findMIS(A,randperm(size(A,1)));
smlDsqr=Dsqr(X>0,X>0); 
sum(X)
tabulate(smlDsqr(:))

%%
