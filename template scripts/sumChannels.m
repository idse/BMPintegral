clear all; close all;
scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);

% raw data location, modify if not the same as the location of this script
dataDir = scriptPath; 
MIPdir = fullfile(dataDir,'MIP');
cd(MIPdir);

channels = [0 1]; % which channels to sum

% make an array of the filenames for each channel
nChannels = numel(channels);
listing = {};
for ci = 1:nChannels
    listing{ci} = dir(fullfile(MIPdir,['*_MIP_*' sprintf('w%.4d',channels(ci)) '.tif']));
end

%%
profile on

% get numer of time points
fname = {};
for ci = 1:nChannels
    fname{ci} = fullfile(MIPdir,listing{ci}(1).name);
end
info = imfinfo(fname{1});

% start combining
for pi = 2:numel(listing{1})

    % filename to read
    for ci = 1:nChannels
        fname{ci} = fullfile(MIPdir,listing{ci}(pi).name);
    end
    % filename to write
    s = strsplit(listing{ci}(pi).name,'_w');
    sumfname = [s{1} '_SUM.tif'];
    
    % read data and sum
    for ti = 1:numel(info)
        im = {};
        for ci = 1:nChannels
            im{ci} = imread(fname{ci},ti);
        end
        imcat = cat(3,im{:});
        imsum = uint16(sum(imcat,3));
        
        if ti==1
            imwrite(imsum, fullfile(MIPdir,sumfname));
        else
            imwrite(imsum, fullfile(MIPdir,sumfname),'WriteMode','append');
        end
        % progress indicator
        fprintf('.');
        if mod(ti,80)==0
            fprintf('\n');
        end
    end
end

profile viewer
