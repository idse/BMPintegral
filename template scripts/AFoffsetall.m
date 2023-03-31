clear all; close all;
scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);

% raw data location, modify if not the same as the location of this script
protocolFile = '/Users/idse/data_tmp/190918_EF_DAPI_TripleReporter AW8.bkp';

% um/pixel
resolution = 0.301;      %0.603; %20x         %1.26; 10x?           0.301; % 40x
micropattern  = true;

gridSize = [3 3]; % for 40x
montageOverlap = 20; % percentage

% 3x3 grid will have length (2*(1024-overlap) + 1024-2*overlap)*resolution
% at 40x that is 868 micron so large enough for a 700 micron colony
overlapPixel = round(1024*montageOverlap/100);

%% make grid centered around single position (or average of two)

protocol = readFusionProtocol(protocolFile);

newProtocol = protocol;
newProtocol.AFoffset(2:end) = newProtocol.AFoffset(1);

%% generate new protocol file

writeFusionProtocol(protocolFile, newProtocol);

