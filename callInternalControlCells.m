function [QueryPop,RefPop,FRthrsh] = callInternalControlCells(ChaseData,varargin)

arg.channelcol = 3; 
arg.refpopcol = 8; 
arg.percentile = 8; 
arg.plot = false; 
arg = parseVarargin(varargin,arg); 

QueryPop=cell(size(ChaseData)); 
RefPop=cell(size(QueryPop)); 
FRthrsh=nan(size(QueryPop,1),1); 

for i = 1:size(ChaseData,1)
    FRthrsh(i) = prctile(ChaseData{i,arg.refpopcol}(:,arg.channelcol),arg.percentile );
    for j=1:size(ChaseData,2)
        if i == 1 % no internal ref in day 1
            RefPop{1,j} = ChaseData{1,arg.refpopcol};
            QueryPop{1,j} = ChaseData{1,j};
        else
            isref = ChaseData{i,j}(:,arg.channelcol)>FRthrsh(i); 
            RefPop{i,j} = ChaseData{i,j}(isref,:);
            QueryPop{i,j} =  ChaseData{i,j}(~isref,:);
        end
    end
end
          