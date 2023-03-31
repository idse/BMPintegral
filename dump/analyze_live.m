clear all; close all;
scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);

% raw data location, modify if not the same as the location of this script
dataDir = scriptPath;

meta = Metadata(dataDir);

% manual metadata
%-------------------

treatmentTime = 4;

meta.timeInterval = '10 min';
meta.nWells = 6;
meta.posPerCondition = 3;

% specify these when working with a grid:
% montageOverlap = [2 2]     % percent overlap of montage
% montageGridSize = 20;    % grid size n by m locations

if meta.nPositions ~= meta.nWells*meta.posPerCondition
    warning('total position doesnt match nWells x posPerCondition');
end

% in order of imaging (clockwise, which is not the layout of the dish)
meta.conditions = { 'ctrl', 'PD17=FGFR1i',  'PD03=MEK1/2i',...
                    'E6',   'E6+bFGF',      'SU54=RTKi'};

% set this to true if making an '8-well' loop through a 4-well
loop4well = false;

meta.nucChannel = 1;
ERKKTRChannel = 0;

%% save stitched previews of the MIPs

stitchedPreviews(dataDir, meta); 

%% extract nuclear and cytoplasmic levels

% Andor channel counting starts at 0
opts = struct(...
                    'dataChannels',     [ERKKTRChannel meta.nucChannel],...
                    'fgChannel',        ERKKTRChannel,...
                    'cytoplasmicLevels',true,... 
                    'segFG',            0,... % change if result below looks inverted
                    'segmentationDir',  fullfile(dataDir,'MIP'),...
                    'MIPidxDir',        fullfile(dataDir,'MIP'),...
                    'nucShrinkage',     2,...
                    'cytoSize',         3,...
                    'bgMargin',         10,...
                    'NCRcutoff',        [3 Inf]);

opts.cleanupOptions = struct('separateFused', true,...
    'clearBorder',true, 'minAreaStd', 1, 'minSolidity',0, 'minArea',250);


%% check that the options are set right

pi = 2;
P = Position(meta);
P.setID(pi);
time = 2;
opts.tMax = time;

% try out the nuclear cleanup settings on some frame:
% bla = nuclearCleanup(seg(:,:,time), opts.cleanupOptions);
% imshow(bla)

dbInfo = P.extractData(dataDir, opts);

bg = P.cellData(time).background
nucl = P.cellData(time).nucLevelAvg
cytl = P.cellData(time).cytLevelAvg
(nucl-bg)./(cytl - bg)

im = P.loadImage(dataDir, ERKKTRChannel, time);
%labels = P.loadTrackingLabels(MIPdir, nucChannel, time);
MIP = max(im,[],3);
A = imadjust(mat2gray(MIP));
s = 0.4;
figure,
imshow(cat(3, A + dbInfo.nucmaskraw*s, A + s*dbInfo.nucmask + s*dbInfo.bgmask, A + s*dbInfo.cytmask));

opts.tMax = meta.nTime;

%% run the analysis on all time points

tic
positions(meta.nPositions) = Position();

for pi = 1:meta.nPositions

    positions(pi) = Position(meta);
    positions(pi).setID(pi);
    positions(pi).extractData(dataDir, opts);
    positions(pi).makeTimeTraces();
    save(fullfile(dataDir,'positions'), 'positions');
end
toc

%% load results if above block was run previously

load(fullfile(dataDir,'positions'));

% % to concatenate two files (if microscope had to be stopped)
%
% positions_before = load(fullfile(dataDir,'positions_before'));
% positions_before = positions_before.positions;
% 
% for pi = 1:meta.nPositions
%     positions(pi).cellData = cat(2, positions_before(pi).cellData,...
%                                             positions(pi).cellData);
%     positions(pi).ncells = cat(2, positions_before(pi).ncells,...
%                                             positions(pi).ncells);
% 	positions(pi).nTime = positions_before(pi).nTime + positions(pi).nTime;
%     positions(pi).makeTimeTraces;
% end
% 
% save(fullfile(dataDir,'positions'), 'positions');

%% make figure of the time traces

options = struct('wellnrs', 1, 'channel', S4Channel,...
                'mode','N:C',...
                'yLabel', 'Smad4 (N:C)',...
                'loop4well', loop4well,...
                'minNCells', 10,...
                'dataChannels', opts.dataChannels,...
                'treatmentTime', treatmentTime,...
                'SNRcutoff', 1.5,...
                'signalingRange', [0.6 1.2]);

% individual conditions
wellnrs = options.wellnrs;
for wellnr = wellnrs
    options.wellnrs = wellnr;
    makeSignalingGraphs(positions, meta, options);
    export_fig(fullfile(dataDir,['timeTrace_well' num2str(wellnr) '.png']),'-native -m2');
end

% multiple conditions 
options.wellnrs = wellnrs;
makeSignalingGraphs(positions, meta, options)
legend(meta.conditions);
export_fig(fullfile(dataDir,'timeTrace_multipleConditions.png'),'-native -m2');

