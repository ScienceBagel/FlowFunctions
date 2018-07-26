function [FileNameCellArray,Populations] = getExperiment(ExperimentID,varargin)

arg.docid = '1RerIt9CSs5-gXcVf-qPLIdYZXFXqhJEiPO-a4ql6-a0'; 
arg.list = false; 
arg.filename = 'Data.csv'; 
arg = parseVarargin(varargin,arg); 

if ~arg.list 
    % find the GID for this experiment
    ExperimentList = getExperiment([],'list',true,'filename','List.csv');
    gid = ExperimentList.GID(ismember(ExperimentList.ExperimentID,ExperimentID));
else
    gid=0; % get the list sheet (first one)
end

% find which experiment we are going to use
url = sprintf('https://docs.google.com/spreadsheets/d/%s/export?format=csv&gid=%d',arg.docid,gid); 
tmpfilename = websave(arg.filename,url); 

SheetData = readtable(tmpfilename); 

% if getting list asigne and returm
if arg.list
    FileNameCellArray=SheetData; 
else % if not, process filename to add local paths and remove first col
    FileNameCellArray = table2cell(SheetData);
    Populations = SheetData.Properties.VariableNames;
    Populations(1)=[]; 
    
    for i=1:size(FileNameCellArray,1)
        for j=2:size(FileNameCellArray,2)
            FileNameCellArray(i,j)=fullfile(FileNameCellArray(i,1),FileNameCellArray(i,j));
        end
    end
    FileNameCellArray(:,1)=[]; 
end
