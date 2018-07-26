function [x_subplot,y_subplot] = makeAccuriMatFile(pth,savename,varargin)
% function output the percentile of each point in Query population in
% comparison to the Entire Population conditioned on a Cond variable 


%% Argument Inputs
arg.prefix = '';


%% Add Path
addpath('/home/maeve2/Desktop/Repositories/Flow_Analysis_Functions')

%% Load FCS Data
list = dir(fullfile(pth,'*.fcs'));           
fcsname = struct2cell(list);                 
fcsfilenames=fcsname(1,1:size(fcsname,2));   
n=size(fcsfilenames,2);
wellname=cell(n,1);

for i=1:n
    name = strsplit(fcsfilenames{i},'.');
    if exist('prefix','var') > 0;
        name2 =strsplit(char(name(1)),prefix);
        wellname(i) = name2(2);
    else   
        wellname(i) = name(1);  
    end
end

x_subplot = round(sqrt(n))+1;
y_subplot = round(sqrt(n));

%%

ARD = cell(1,n);
for i = 1:n
    
    [fcsdat, fcshdr] = fca_readfcs(fullfile(pth,[fcsfilenames{i}]));

    fscA_ix = ismember({fcshdr.par.name},'FSC-A');
    sscA_ix = ismember({fcshdr.par.name},'SSC-A');
    fl1A_ix = ismember({fcshdr.par.name},'FL1-A');
    fl2A_ix = ismember({fcshdr.par.name},'FL2-A');
    fl3A_ix = ismember({fcshdr.par.name},'FL3-A');
    fl4A_ix = ismember({fcshdr.par.name},'FL4-A');
    fscH_ix = ismember({fcshdr.par.name},'FSC-H');
    sscH_ix = ismember({fcshdr.par.name},'SSC-H');
    fl1H_ix = ismember({fcshdr.par.name},'FL1-H');
    fl2H_ix = ismember({fcshdr.par.name},'FL2-H');
    fl3H_ix = ismember({fcshdr.par.name},'FL3-H');
    fl4H_ix = ismember({fcshdr.par.name},'FL4-H');   
    
    fscA = fcsdat(:,fscA_ix);
    sscA = fcsdat(:,sscA_ix);
    fl1A = fcsdat(:,fl1A_ix);
    fl2A = fcsdat(:,fl2A_ix);
    fl3A = fcsdat(:,fl3A_ix);
    fl4A = fcsdat(:,fl4A_ix);  
    fscH = fcsdat(:,fscH_ix);
    sscH = fcsdat(:,sscH_ix);
    fl1H = fcsdat(:,fl1H_ix);
    fl2H = fcsdat(:,fl2H_ix);
    fl3H = fcsdat(:,fl3H_ix);
    fl4H = fcsdat(:,fl4H_ix);      
  
    ARD{i} = [fscA,sscA,fl1A,fl2A,fl3A,fl4A,fscH,sscH,fl1H,fl2H,fl3H,fl4H];
end
    
%% Make Live Cell Gate

figure('Name','Make a Live Cell Gate')
clf
hold on
set(gcf, 'Position', get(0, 'Screensize'));
for i=1
    tmp = ARD{i};
    ScatterPlotFACS_v1([tmp(:,fscA_ix),tmp(:,sscA_ix)]);
    xlabel('FSC');
    ylabel('SSC');
    title(wellname{i})
end

output = ginput();
gate_live_x = output(:,1);
gate_live_y = output(:,2);
    
%% Put Gate Into Action
dead_or_alive=cell(n,1);
n_live=zeros(n,1);
p_live=zeros(n,1);
alive = cell(1,n);

for i=1:n
    [dead_or_alive{i},~]=inpolygon(ARD{i}(:,fscA_ix),ARD{i}(:,sscA_ix),gate_live_x,gate_live_y);
    n_live(i,1)=numel(ARD{i}(dead_or_alive{i}));
    p_live(i,1)=100*numel(ARD{i}(dead_or_alive{i}))./(size(ARD{i},1));
    alive{i}=ARD{1,i}(dead_or_alive{i},:); 
end
    
%% Check that the live cell gate looks good
figure('Name','Check FSC and SSC Gate')
clf
hold on
set(gcf, 'Position', get(0, 'Screensize'));
for i=1
    tmp = ARD{i};
    ScatterPlotFACS_v1([tmp(:,fscA_ix),tmp(:,sscA_ix)]);
    xlabel('FSC');
    ylabel('SSC');
    title(wellname{i})
    plot(gate_live_x,gate_live_y,'k')
    legend(sprintf('%.1f % %',p_live(i,1)),'location','northeast');
end
   
%% Calculate Percent Live and Percent Fluorescent
gate_live=cat(2,gate_live_x,gate_live_y);

dead_or_alive=cell(n,1);
n_live=zeros(n,1);
p_live=zeros(n,1);

for i=1:n
    [dead_or_alive{i},~]=inpolygon(ARD{i}(:,fscA_ix),ARD{i}(:,sscA_ix),gate_live_x,gate_live_y);
    n_live(i,1)=numel(ARD{i}(dead_or_alive{i}));
    p_live(i,1)=100*numel(ARD{i}(dead_or_alive{i}))./(size(ARD{i},1));
    alive{i}=ARD{1,i}(dead_or_alive{i},:); 
end

%% Calculate percentage with flourescence
% check that the number in the brackets corresponds to the blank sample


fl1gate=prctile(alive{1}(:,fl1A_ix),99.9);
n_fl1 = zeros(n,1);
p_fl1 = zeros(n,1);
fl1_alive = cell(n,1);
for i = 1:n
    fl1_alive{i} = dead_or_alive{i} & (ARD{i}(:,fl1A_ix)>fl1gate);
    n_fl1(i,1) = numel(ARD{i}(:,fl1A_ix)); % source of pain
    p_fl1(i,1) = 100 * numel(ARD{i}(:,fl1A_ix))./numel(ARD{i}(dead_or_alive{i}));
end


fl2gate=prctile(alive{1}(:,fl2A_ix),99.9);
n_fl2 = zeros(n,1);
p_fl2 = zeros(n,1);
fl2_alive = cell(n,1);
for i = 1:n
    fl2_alive{i} = dead_or_alive{i} & (ARD{i}(:,fl2A_ix)>fl2gate);
    n_fl2(i,1) = numel(ARD{i}(:,fl2A_ix)); % source of pain
    p_fl2(i,1) = 100 * numel(ARD{i}(:,fl2A_ix))./numel(ARD{i}(dead_or_alive{i}));
end

fl3gate=prctile(alive{1}(:,fl3A_ix),99.9);
n_fl3 = zeros(n,1);
p_fl3 = zeros(n,1);
fl3_alive = cell(n,1);
for i = 1:n
    fl3_alive{i} = dead_or_alive{i} & (ARD{i}(:,fl3A_ix)>fl3gate);
    n_fl3(i,1) = numel(ARD{i}(:,fl3A_ix)); % source of pain
    p_fl3(i,1) = 100 * numel(ARD{i}(:,fl3A_ix))./numel(ARD{i}(dead_or_alive{i}));
end

fl4gate=prctile(alive{1}(:,fl4A_ix),99.9);
n_fl4 = zeros(n,1);
p_fl4 = zeros(n,1);
fl4_alive = cell(n,1);
for i = 1:n
    fl4_alive{i} = dead_or_alive{i} & (ARD{i}(:,fl4A_ix)>fl4gate);
    n_fl4(i,1) = numel(ARD{i}(:,fl4A_ix)); % source of pain
    p_fl4(i,1) = 100 * numel(ARD{i}(:,fl4A_ix))./numel(ARD{i}(dead_or_alive{i}));
end


%% Save Mat File

cd('/home/maeve2/Desktop/Flow_Cytometry_Data/Mat Files')
save(savename)


