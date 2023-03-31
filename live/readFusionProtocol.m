function protocol = readFusionProtocol(protocolFile)
    % should eventually read all relevant features of protocol
    % for now just the positions

    protocol = struct('XYZ',[],'AFoffset',[]);

    fid = fopen(protocolFile);

    if fid == -1
        error('file not found');
    end

    % get positions
    %------------------
    tline = fgets(fid);
    while ischar(tline) 
        % find section with positions
        if contains(tline,'MultiFieldDefinitions') && ~contains(tline,']')
            % scan through section
            while ischar(tline) && ~contains(tline,' ]')
                % definition of individual positions
                if contains(tline,'{')
                    
                    xline = strsplit(fgets(fid),{'"X": ',','});
                    yline = strsplit(fgets(fid),{'"Y": ',','});
                    zline = strsplit(fgets(fid),{'"ReferenceZ": ',','});
                    XYZ = [str2num(xline{2}) str2num(yline{2}) str2num(zline{2})];
                    
                    fieldname = strsplit(fgets(fid),{'"FieldName": ',','});
                    AFline = strsplit(fgets(fid),{'"DriftStabilisationOffset": ',','});

                    protocol.XYZ = cat(1, protocol.XYZ, XYZ);
                    protocol.AFoffset = cat(1, protocol.AFoffset, str2num(AFline{2}));
                end
                tline = fgets(fid);
            end
        end
        tline = fgets(fid);
    end
    
    fclose(fid);
end

