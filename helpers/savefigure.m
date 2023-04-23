function savefigure(fname,f)

if ~exist('f','var')
    f = gcf;
end

[filepath,~,ext] = fileparts(fname);
if isempty(ext) || strcmp(ext,"")
    fname = strcat(fname,'.png');
end
if isempty(filepath) || length(filepath) == 1
    %if no path is specified, save to a default folder under today's date
    baseDir = 'DEFAULTPATHHERE';
    c = clock;
    saveDir = sprintf('%d%.2d%.2d',c(1),c(2),c(3));
    saveDir = fullfile(baseDir,saveDir(3:end));
    
    %if no such folder exists, make it
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    
    fname = fullfile(saveDir,fname);
else
    %if the specified folder does not exist, make it
    if ~exist(filepath,'dir'), mkdir(filepath); end
end

saveas(f,fname);

end