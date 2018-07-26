function [Percntile,Residuals] = percentileVsEntirePopUsingConditionalReadout(QueryPop_Readout,QueryPop_Cond,EntirePop_Readout,EntirePop_Cond,varargin)
% function output the percentile of each point in Query population in
% comparison to the Entire Poppulation conditioned on a Cond variable 


arg.grid = 0:5:100;
arg = parseVarargin(varargin,arg); 
QueryPop_Readout=QueryPop_Readout(:); 
QueryPop_Cond=QueryPop_Cond(:); 
EntirePop_Readout=EntirePop_Readout(:); 
EntirePop_Cond=EntirePop_Cond(:); 


Percntile=nan(size(QueryPop_Readout)); 
Residuals=nan(size(Percntile)); 
CondGrid = prctile(EntirePop_Cond(:),arg.grid);  

for i=2:numel(arg.grid)
    % find indexes of Entire based on Cond and their Readout vector
    ix_entire_cond = EntirePop_Cond>CondGrid(i-1) & EntirePop_Cond<=CondGrid(i); 
    
    % find indexes of cells with Query that have similar value
    % edges cases in both sides... 
    if i==2
        ix_query_cond = QueryPop_Cond<=CondGrid(i); 
    elseif i==numel(arg.grid)
        ix_query_cond = QueryPop_Cond > CondGrid(i-1); 
    else
        ix_query_cond = QueryPop_Cond > CondGrid(i-1) & QueryPop_Cond<=CondGrid(i); 
    end
    
    % replicate both vectors orthogonally
    RepMatQueryReadout = repmat(QueryPop_Readout(ix_query_cond),1,nnz(ix_entire_cond)); 
    RepMatEntireReadout = repmat(EntirePop_Readout(ix_entire_cond)',nnz(ix_query_cond),1); 
    Percntile(ix_query_cond) = mean(RepMatQueryReadout>RepMatEntireReadout,2); 
    Residuals(ix_query_cond) = QueryPop_Readout(ix_query_cond)-mean(QueryPop_Readout(ix_query_cond)); 
    
end
