clear all; close all;
scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);

% raw data location, modify if not the same as the location of this script
dataDir = scriptPath; 

%-----------metadata------------------------

% if actual stitching posPerCondition should be number of montages /
% colonies, noth the actual number of positions
manualMeta = struct();
manualMeta.nWells = 18;
manualMeta.posPerCondition = 11;

% specify if actual stitching:
manualMeta.montageGridSize = [3 3];      % [3 3] for 40x micropattern, [2 2] for 20x
manualMeta.montageOverlap = 20;   % percentage

meta = Metadata(dataDir, manualMeta);

save(fullfile(dataDir,'meta.mat'),'meta');

%-----------stitch options------------------------

stitchOptions = struct(     'saveStitchFullSize', true,...
                            'FusionMontage', false,...
                            'dichroicCorrectionChannels', [],...
                            'stitchStacks', true); 
                        
% empty stitch options if not stitching
% stitchOptions = struct();

% ------------MIP options--------------------------

% channels: channels to process, counting start from 0 here
% saveidx: for channels containing nuclear marker if multiple z-slices 
% tmax: cutoff on time points to process (e.g. when it went out of focus)
MIPoptions = struct(        'channels', 0:1,...    
                            'saveidx',  [true false],... 
                            'tmax',     []  );

%MIPoptions.videoOptions = struct('format','mp4', 'FrameRate', 5);
                        
makeStitchedPreviews = true;

%%
inputdir = dataDir;
outputdir = fullfile(dataDir, 'MIP');
if ~exist(outputdir,'dir'), mkdir(outputdir); end

batchMIP(inputdir, outputdir, MIPoptions);

if makeStitchedPreviews
    stitchedPreviews(dataDir, meta, stitchOptions, MIPoptions);
end
