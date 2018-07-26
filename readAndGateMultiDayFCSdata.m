function [ChaseData,Populations] = readAndGateMultiDayFCSdata(ExperimentID,Gates,varargin)

arg.channelmapping = struct('BFP',{{'Pac Blue-A','DAPI-A'}},...
                            'GFP',{{'Alexa Fluor 488-A','GFP-A'}},...
                            'FarRed',{{'APC-A'}},...
                            'SSC',{{'SSC-A'}},...
                            'FSC',{{'FSC-A'}},...
                            'RFP',{{'PI'}});
arg.channels = {'BFP','GFP','FarRed','SSC','FSC'}; 
arg.basepth=''; 
arg = parseVarargin(varargin,arg); 

% get file names for this experiment ID
[FileNameCellArray,Populations] = getExperiment(ExperimentID); 

% load all data
ChaseData = cell(size(FileNameCellArray)); 
for i=1:numel(ChaseData)
    [~, fcshdr, ~, fcsdatcomp] = fca_readfcs(fullfile(arg.basepth,FileNameCellArray{i})); 
    for j=1:numel(arg.channels)
        chnl_ix = ismember({fcshdr.par.name},arg.channelmapping.(arg.channels{j}));
        ChaseData{i}(:,j)=fcsdatcomp(:,chnl_ix); 
    end
end

% gate
for i=1:numel(ChaseData)
    ChaseData{i}=gateFCSdata(ChaseData{i},Gates); 
end