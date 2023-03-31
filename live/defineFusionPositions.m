function lines = defineFusionPositions(protocol)
    
    lines = {};
    for i = 1:size(protocol.XYZ,1)
        
        lines{i} = {
            '{',
            sprintf('"X": %f,', protocol.XYZ(i,1)),
            sprintf('"Y": %f,', protocol.XYZ(i,2)),
            sprintf('"ReferenceZ": %f,', protocol.XYZ(i,3)),
            sprintf('"FieldName": "%d",',num2str(i)),
            sprintf('"DriftStabilisationOffset": %f,',num2str(protocol.AFoffset(i))),
            '"ScanZ": 0.0,',
            '"IsValid": true,',
            '"PropertyNames": [],',
            '"IsDisposed": false,',
            '},'
        };
    end
end