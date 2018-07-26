function [x_subplot,y_subplot] = makeMatFile(pth,savename,varargin)
% function output the percentile of each point in Query population in
% comparison to the Entire Population conditioned on a Cond variable 


%% Argument Inputs
arg.prefix = '';


%% Add Path
addpath('/home/maeve2/Desktop/Repositories/Flow_Analysis_Functions')

%% Convert files

list = dir(fullfile(pth,'*.fcs'));           
fcsname = struct2cell(list);                 
fcsfilenames=fcsname(1,1:size(fcsname,2));   
n=size(fcsfilenames,2);
wellname=cell(n,1);
% 
% for i=1:n
%     name = strsplit(fcsfilenames{i},'.');
%     if exist('var','arg.prefix') > 0;
%         name2 = strsplit(char(name(1)),arg.prefix);
%         wellname(i) = name2(2);
%     else   
%         wellname(i) = name(1);  
%     end
% end

for i=1:n
    name = strsplit(fcsfilenames{i},'.');
    if isempty(arg.prefix) == 0
        wellname(i) = name(1);
    end
    if isempty(arg.prefix) == 1
         name2 = strsplit(char(name(1)),arg.prefix);
         wellname(i) = name2(2);
    end
end    

x_subplot = round(sqrt(n))+1;
y_subplot = round(sqrt(n));


%% Load FCS Data

for i = 1:n
    
    % Read Data
    [fcsdat, fcshdr] = fca_readfcs(fullfile(pth,[fcsfilenames{i}]));
    channels = {fcshdr.par.name};
    
    % Determine number of channels
    num_channels = length({fcshdr.par.name})-3;
    
    % Pull forward and side scatter data
    fsc_ix = ismember({fcshdr.par.name},'FSC-A');
    ssc_ix = ismember({fcshdr.par.name},'SSC-A');
    fsc = fcsdat(:,fsc_ix);
    ssc = fcsdat(:,ssc_ix);
    
    if num_channels > 0
            fluor1_ix = ismember({fcshdr.par.name},channels(4));
            fluor1 = fcsdat(:,fluor1_ix);
    end    
    if num_channels > 1
        fluor2_ix = ismember({fcshdr.par.name},channels(5));
        fluor2 = fcsdat(:,fluor2_ix);
    end
    if num_channels > 2
        fluor3_ix = ismember({fcshdr.par.name},channels(6));
        fluor3 = fcsdat(:,fluor3_ix);
    end
    if num_channels > 3
        fluor4_ix = ismember({fcshdr.par.name},channels(7));
        fluor4 = fcsdat(:,fluor4_ix);
    end
    
% Create ARD to store single cell fluorescence data

    if num_channels == 0
        ARD{i} = [fsc,ssc];
    end
    if num_channels == 1
        ARD{i} = [fsc,ssc,fluor1];
    end
    if num_channels == 2
        ARD{i} = [fsc,ssc,fluor1,fluor2];
    end
    if num_channels == 3
        ARD{i} = [fsc,ssc,fluor1,fluor2,fluor3];
    end
    if num_channels == 4
        ARD{i} = [fsc,ssc,fluor1,fluor2,fluor3,fluor4];
    end
end


    
%% Rename indexes
fsc_ix = 1;
ssc_ix = 2;
fluor1_ix = 3;
fluor2_ix = 4;
fluor3_ix = 5;
fluor4_ix = 6;

%% Make gate for live cells and calculate percentage live cells
% Use the Ginput function to create gates around live cells

figure('Name','Comparing FSC and SSC')
clf
hold on
set(gcf, 'Position', get(0, 'Screensize'));
for i=1
    tmp = ARD{i};
    ScatterPlotFACS_v1([tmp(:,fsc_ix),tmp(:,ssc_ix)]);
    xlabel('FSC');
    ylabel('SSC');
    title(wellname{i})
end

output = ginput();
gate_live_x = output(:,1);
gate_live_y = output(:,2);

%% Calculate Percent Live and Percent Fluorescent

dead_or_alive=cell(n,1);
n_live=zeros(n,1);
p_live=zeros(n,1);

alive = cell(1,n);

for i=1:n
    [dead_or_alive{i},~]=inpolygon(ARD{i}(:,fsc_ix),ARD{i}(:,ssc_ix),gate_live_x,gate_live_y);
    n_live(i,1)=numel(ARD{i}(dead_or_alive{i}));
    p_live(i,1)=100*numel(ARD{i}(dead_or_alive{i}))./(size(ARD{i},1));
    alive{i}=ARD{1,i}(dead_or_alive{i},:); 
end

%% Calculate percentage with flourescence
% check that the number in the brackets corresponds to the blank sample
%% Calculate percentage with flourescence
% check that the number in the brackets corresponds to the blank sample

if exist('fluor1','var')
    fluor1gate=prctile(alive{1}(:,fluor1_ix),99.9);
    n_fluor1 = zeros(n,1);
    p_fluor1 = zeros(n,1);
    fluor1_alive = cell(n,1);
    for i = 1:n
        fluor1_alive{i} = dead_or_alive{i} & (ARD{i}(:,fluor1_ix)>fluor1gate);
        n_fluor1(i,1) = numel(ARD{i}(:,fluor1_ix)); % source of pain
        p_fluor1(i,1) = 100 * numel(ARD{i}(:,fluor1_ix))./numel(ARD{i}(dead_or_alive{i}));
    end
end


if exist('fluor2','var')
    fluor2gate=prctile(alive{1}(:,fluor2_ix),99.9);
    n_fluor2 = zeros(n,1);
    p_fluor2 = zeros(n,1);
    fluor2_alive = cell(n,1);
    for i = 1:n
        fluor2_alive{i} = dead_or_alive{i} & (ARD{i}(:,fluor2_ix)>fluor2gate);
        n_fluor2(i,1) = numel(ARD{i}(:,fluor2_ix));
        p_fluor2(i,1) = 100 * numel(ARD{i}(:,fluor2_ix))./numel(ARD{i}(dead_or_alive{i}));
    end
end

if exist('fluor3','var')
    fluor3gate=prctile(alive{1}(:,fluor3_ix),99.9);
    n_fluor3 = zeros(n,1);
    p_fluor3 = zeros(n,1);
    fluor3_alive = cell(n,1);
    for i = 1:n
        fluor3_alive{i} = dead_or_alive{i} & (ARD{i}(:,fluor3_ix)>fluor3gate);
        n_fluor3(i,1) = numel(ARD{i}(:,fluor3_ix));
        p_fluor3(i,1) = 100 * numel(ARD{i}(:,fluor3_ix))./numel(ARD{i}(dead_or_alive{i}));
    end
end

if exist('fluor4','var')
    fluor4gate=prctile(alive{1}(:,fluor4_ix),99.9);
    n_fluor4 = zeros(n,1);
    p_fluor4 = zeros(n,1);
    fluor4_alive = cell(n,1);
    for i = 1:n
        fluor4_alive{i} = dead_or_alive{i} & (ARD{i}(:,fluor4_ix)>fluor4gate);
        n_fluor4(i,1) = numel(ARD{i}(:,fluor4_ix));
        p_fluor4(i,1) = 100 * numel(ARD{i}(:,fluor4_ix))./numel(ARD{i}(dead_or_alive{i}));
    end
end

%% Save Mat File

cd('/home/maeve2/Desktop/Flow_Cytometry_Data/Mat Files')
save(savename)

