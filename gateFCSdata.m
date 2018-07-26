function FCSmatrix = gateFCSdata(FCSmatrix,Gates,varargin)
% funciton to gate FCS matrix based on Gates. 
% FCSmatrix - Matrix of cell data
% Gates - a struct array with fields channels / gate
% examples
%    single cutoff
%         Gates(i).channel = 'SSC'; % must map to ColumHeader.
%         Gates(i).gate = @(ssc) ssc < 1.5; 
%         
%    polygon (a 3 verteces one)
%         Gates(i).channel = {'SSC','FCS'}
%         Gates(i).gate = @(sf) inpolygon(sf(:,1),sf(:,2),[0.1 1.5 0.7],[0.7 0.2 0.8]); 

% get inputs
arg.channels = {'BFP','GFP','FarRed','SSC','FSC'}; 
arg = parseVarargin(varargin,arg); 

includeInGate = true(size(FCSmatrix,1),1); 

for i=1:numel(Gates)
    % make into a cell array (
    if ~iscell(Gates(i).channel) && ischar(Gates(i).channel)
        Gates(i).channel={Gates(i).channel}; 
    end
    % gate data to gate on (could be multiple channels
    ToGateOn=nan(size(includeInGate,1),numel(Gates(i).channel));
    for j=1:size(ToGateOn,2)
        ToGateOn(:,j)=FCSmatrix(:,ismember(arg.channels,Gates(i).channel{j})); 
    end
    % add gate
    includeInGate = Gates(i).gate(ToGateOn); 
end
    
FCSmatrix=FCSmatrix(includeInGate,:);        