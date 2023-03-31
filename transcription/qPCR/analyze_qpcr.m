clear all; close all;

scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);

% data location, modify if not the same as the location of this script
dataDir = scriptPath; 
cd(dataDir);

filenames = {'200829_FGFqPCR_FGF2_4_8.txt'};

% filenamesNorm = {'170115 A New TimeSeries ATP5 Oct4 Nanog Sox2_dataX.txt',...
%     '170115 A New TimeSeries ATP5O Oct4 Nanog Sox2 B_data.txt'};

% Conditions: (samples)
conditionLabels = {
	'Ctrl 0h',...
	'Ctrl 2h',...
	'FGF2',...
    'FGF2 + SB + LDN + iwp2',...
    'FGF2 + SB + LDN + iwp2 + MEKi',...
    'BMP',...
	'BMP + SB + iwp2',...
	'BMP + SB + iwp2 + MEKi',...
	'BMP + SB + iwp2 + FGFRi',...
	'WNT3A',...
	'WNT3A + SB + LDN',...
	'WNT3A + SB + LDN + MEKi',...
	'WNT3A + SB + LDN + FGFRi',...
	'Activin',...
	'Activin + LDN + iwp2',...
	'Activin + LDN + iwp2 + MEKi',...
	'Activin + LDN + iwp2 + FGFRi',...
	'Ctrl 0h dup',...
    'Ctrl 2h dup',...
	'SB',...
	'LDN',...
	'iwp2',...
	'MEKi',...
	'FGFRi',...
    };

normName = 'ATP5O';
normName = 'APT5O'; % Bohan misspelled

%% read data

tolerance = 1;
data = combineQPCR(dataDir, filenames, tolerance);

Ntargets = data.Ntargets;
Nsamples = data.Nsamples;
targets = data.targets;
CTmean = data.CTmean;
CTstd = data.CTstd;

% normalization
refIdx = find(strcmp(targets,normName));
if isempty(refIdx)
    disp('using normalization from other file');
	refData =  combineQPCR(dataDir, filenamesNorm);
    refIdx = find(strcmp(refData.targets,normName));
    CTref = refData.CTmean(:,refIdx);
    refIdx = 0; % to make sure this idx isn't excluded from plot later
else
    disp('using normalization from same file');
    CTref = CTmean(:,refIdx);
end

barefname = [targets{:}];

%% normalize

% index into data.samples for ctrl sample (e.g. untreated)
ctrlSampleIdx = 4; 

CTnorm = zeros([Nsamples Ntargets]);
% normalize to normalizing / housekeeping gene
for targeti = 1:Ntargets
    CTnorm(:,targeti) = -(CTmean(:,targeti) - CTref);
end
% normalize to ctrl sample
for targeti = 1:Ntargets
    CTnorm(:,targeti) = CTnorm(:,targeti) - CTnorm(ctrlSampleIdx,targeti);
end

%% special for Bohan, sort by sample number
% so we can use the ordered list of conditions above

[~,p] = sort(str2double(data.samples));
CTnorm = CTnorm(p,:);

%% make individual bar graphs for each target

labels = categorical(conditionLabels);

targetIndices = 2:data.Ntargets; % skipping 1 here because that is ATP5O

for targetIdx = targetIndices 

    %labels = categorical(data.samples);
	%labels = reordercats(labels,data.samples);
    
    bar(labels, CTnorm(:, targetIdx))
    ylabel('normalized CT');
    title(data.targets(targetIdx));
    fs = 16;
    set(gca, 'LineWidth', 2);
    set(gca,'FontSize', fs)
    
    saveas(gcf, ['normalizedCT_' data.targets{targetIdx} '.png']);
end
    
%% make a combined bar graph

bar(labels, CTnorm(:,targetIndices))
ylabel('normalized CT');
title('FGF 2,4,8 response at 2h');
legend(data.targets(targetIndices))
fs = 16;
set(gca, 'LineWidth', 2);
set(gca,'FontSize', fs)
    
saveas(gcf, ['normalizedCT_' data.targets{targetIndices} '.png']);



