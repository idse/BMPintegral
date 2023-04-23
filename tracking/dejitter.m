function [xy, offsets] = dejitter(this, dataDir)
%stabilize calculated cell positions in shaky microscope time-lapse
%this = Position object with cell data
%dataDir = folder containing images used for alignment
%xy = cell array with new xy positions of each cell at each time
%offsets = array of cumulative frame-frame displacements

channel = this.nucChannel;
ntime = this.nTime;

imprev = loadImage(this, dataDir, channel, 1); %load first image
disps = zeros(ntime,2); %frame-frame displacements
for ti = 2:ntime
    fprintf('.')
    if mod(ti,45) == 0
        fprintf('\n')
    end
    %check if image is empty
    if this.ncells(ti) > 0
        im = loadImage(this, dataDir, channel, ti);
        [shiftx,shifty,~,~] = xcorr2fft(imprev, im);
        disps(ti,:) = [shiftx, shifty];
        imprev = im;
    end
end
fprintf('\n')

offsets = cumsum(disps,1); %cumulative displacements
xy = cell(ntime,1);
for ti = 1:ntime
    if this.ncells(ti) > 0
        xy{ti} = this.cellData(ti).XY - offsets(ti,2:-1:1);
    end
end


end