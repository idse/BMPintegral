function lines = writeFusionProtocol(protocolFile, protocol)
    
    newProtocolFile = [protocolFile(1:end-4) '_NEW.bkp'];
    fid = fopen(protocolFile);
    fido = fopen(newProtocolFile,'w');
    
    if fid == -1
        error('file not found');
    end

    % create new position lines
    newPositionLines = {};
    for i = 1:size(protocol.XYZ,1)
        
        newPositionLines{i} = {
            sprintf('{\n'),
            sprintf('"X": %f,\n', protocol.XYZ(i,1)),
            sprintf('"Y": %f,\n', protocol.XYZ(i,2)),
            sprintf('"ReferenceZ": %.2f,\n', protocol.XYZ(i,3)),
            sprintf('"FieldName": "%d",\n', i),
            sprintf('"DriftStabilisationOffset": %.2f,\n', protocol.AFoffset(i)),
            sprintf('"ScanZ": 0.0,\n'),
            sprintf('"IsValid": true,\n'),
            sprintf('"PropertyNames": [],\n'),
            sprintf('"IsDisposed": false,\n'),
            sprintf('},\n')
        };
    end
    
    % scan through file and write new positions 
    %--------------------------------------------
    lines = {};
    i = 1;
    lines{i} = fgets(fid);
    fwrite(fido, lines{i});
    
    while ischar(lines{i}) 
        % find section with positions
        if contains(lines{i},'MultiFieldDefinitions') && ~contains(lines{i},']')
            
            for j = 1:numel(newPositionLines)
                for k = 1:numel(newPositionLines{j})
                    fwrite(fido, newPositionLines{j}{k});
                end
            end
            
            % scan through section
            while ischar(lines{i}) && ~contains(lines{i},' ]')
                i = i + 1;
                lines{i} = fgets(fid);
            end
            fwrite(fido, lines{i});
        end
        i = i + 1;
        lines{i} = fgets(fid);
        fwrite(fido, lines{i});
    end
    fclose(fid);
    fclose(fido);

end