function positions_combined = concatenatePositions(positions_before, positions_after)
    % to concatenate two files (if microscope had to be stopped)

    for pi = 1:numel(positions_after)
        positions_after(pi).cellData = cat(2, positions_before(pi).cellData,...
                                                positions_after(pi).cellData);
        positions_after(pi).ncells = cat(2, positions_before(pi).ncells,...
                                                positions_after(pi).ncells);
        positions_after(pi).nTime = positions_before(pi).nTime + positions_after(pi).nTime;
        positions_after(pi).makeTimeTraces;
    end
    
    positions_combined = positions_after;
end