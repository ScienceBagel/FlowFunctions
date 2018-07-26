function [r_GR_condB,Vdecomposition,Vtot,Etot,V_E_G_condB,V_E_G_condBR,E_V_G_condBR] = calcConditionalCorrelation(B,G,R,varargin)
% function output the percentile of each point in Query population in
% comparison to the Entire Population conditioned on a Cond variable 

arg.grid = 0:5:100;
arg.corrtype = 'spearman'; 
arg = parseVarargin(varargin,arg); 
arg.nbins = numel(arg.grid)-1; 
BlueGrid = prctile(B(:),arg.grid);
RedGrid = prctile(R(:),arg.grid); 
R=R(:);
B=B(:); 
G=G(:); 

Vdecomposition=nan(1,3); 
Variance_abs=nan(1,3); 
Gmean_cond_B = nan(arg.nbins,1); 
Gvar_cond_B = nan(arg.nbins,1); 
P_blue = nan(arg.nbins,1); 
cov_GB_perbin = nan(arg.nbins,1); 
Gmean_cond_BR=nan(arg.nbins,arg.nbins); 
Gvar_cond_BR = nan(arg.nbins,arg.nbins);
P_red = nan(arg.nbins,arg.nbins);

for i=2:numel(arg.grid)
    % find indexes of B in the bin
    ix_B_inBin = find(B>BlueGrid(i-1) & B<=BlueGrid(i)); 
    
    cv = cov(G(ix_B_inBin),R(ix_B_inBin)); 
    cov_GB_perbin(i-1)=cv(1,2); 
    P_blue(i-1)=numel(ix_B_inBin)/numel(B);
    Gmean_cond_B(i-1)=mean(G(ix_B_inBin));
    Gvar_cond_B(i-1) = var(G(ix_B_inBin));
    for j=2:numel(arg.grid)
        ix_R_inBin = R(ix_B_inBin)>RedGrid(j-1) & R(ix_B_inBin)<=RedGrid(j);
        P_red(i-1,j-1)=mean(ix_R_inBin); 
        Gmean_cond_BR(i-1,j-1)=nanmean(G(ix_B_inBin(ix_R_inBin))); 
        Gvar_cond_BR(i-1,j-1)=nanvar(G(ix_B_inBin(ix_R_inBin)));
    end
end

cov_GR_condB=sum(P_blue.*cov_GB_perbin); % a "fancy" numerical integral.
r_GR_condB=cov_GR_condB./std(G)./std(R); % estimate pearson corr coef

% 
E_V_G_condB = nansum(P_blue.*Gvar_cond_B); % just to test it sums...
V_E_G_condB = nansum(P_blue.*(Gmean_cond_B-sum(P_blue.*Gmean_cond_B)).^2); 

E_V_G_condBR = nansum(P_blue.*nansum(P_red.*Gvar_cond_BR,2)); 
V_E_G_condBR = nansum(P_blue.*nansum(P_red.*(Gmean_cond_BR-repmat(Gmean_cond_B(:),1,arg.nbins)).^2,2)); 

% 
% 
Vtot=var(G); 
Etot=mean(G); 


Vdecomposition(1)=V_E_G_condB/Vtot; % explained by blue
Vdecomposition(2)=V_E_G_condBR/Vtot; % explained by red (and not by blue)
Vdecomposition(3)=E_V_G_condBR/Vtot; % unexplained


