clear all; close all;
scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);

% raw data location, modify if not the same as the location of this script
dataDir = scriptPath;
cd(dataDir);

load updatedworkspace201117.mat

Activin_labels = {'SB10';'x';'A100';...
    'SB10';'x';'A100';...
    'SB10';'x';'A100';...
    'SB10';'x';'A100';...
    'SB10';'x';'A100';...
    'SB10';'x';'A100';...
    'A10';'A10';'A30';...
    'A30';'SB30A100';'SB30A100';...
    'SB100';'SB100';'SB10A100';...
    'SB10A100';'SB10A100';'SB10A100';...
    'SB10A100';'mtesr'};
    
  CHIR_labels = {'IWP2';'IWP2';'IWP2';...
    'x';'x';'x';...
    'C1';'C1';'C1';...
    'C3';'C3';'C3';...
    'C6';'C6';'C6';...
    'C10';'C10';'C10';...
    'C3';'C6';'C3';...
    'C6';'C3';'C6';...
    'C3';'C6';'x';...
    'C1';'C3';'C6';...
    'C10';'Mtesr'}; 

% some more definitions
CHIR_labels_order = {'C10','C6','C3','C1','x','IWP2','mtesr'};
Activin_labels_order = {'mtesr','SB10','x','SB30A100','SB10A100','A10','A30','A100'};

channelLabels = {'DAPI', 'BRA', 'EOMES', 'SOX2'};

nPositions = numel(positions);
nConditions = numel(CHIR_labels);
nChannels = numel(channelLabels);

% to ignore some positions, could be list
excludePositions = 19; 

% calculate condition means and std
conditionMeans = zeros([nConditions nChannels]);
conditionStd = zeros([nConditions nChannels]);

for cii = 1:numel(CHIR_labels)
    
    startci = 4*(cii-1)+1;
    conditionIdx = startci:startci+3;
    conditionIdx = setdiff(conditionIdx, excludePositions);
    
    conditionNucLevels = [];
    
    % concatenate the nuclear levels tables for the different positions of
    % the condition
    for ci = conditionIdx
        conditionNucLevels = cat(1, conditionNucLevels, positions(ci).cellData.nucLevel);
    end
    
    conditionMeans(cii, :) = mean(conditionNucLevels, 1);
    conditionStd(cii, :) = std(conditionNucLevels, 1);
end

% normalize
conditionMeansNorm = (conditionMeans - min(conditionMeans))./(max(conditionMeans,[],1)-min(conditionMeans,[],1));
conditionStdNorm = conditionStd./(max(conditionMeans,[],1)-min(conditionMeans,[],1));

for ci = [2,3]
    figure,
    levels = conditionMeansNorm(:,ci);
    my_table = table(Activin_labels, CHIR_labels, levels);
    h = heatmap(my_table,'Activin_labels','CHIR_labels','ColorVariable','levels');
    h.Colormap = bone;
    h.Title = ['Normalized means of ' channelLabels{ci}];
    h.YDisplayData = CHIR_labels_order;
    h.XDisplayData = Activin_labels_order;
end


